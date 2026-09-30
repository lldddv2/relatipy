"""Map independent geometric states through the public physical Orbit API.

Only units and ``dx^i/dt = u^i/u^t`` are handled here. Initial normalization,
physical transformations and evolution stay in production C. No private
binding or RHS is imported. This module is exclusively test infrastructure.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from importlib import metadata
import json
from pathlib import Path
import platform
import sys
from time import perf_counter

import numpy as np
from astropy import units as u
from astropy.constants import G, c

from relatipy.metrics.kerr import Kerr


REFERENCE_PATH = Path(__file__).resolve().parents[1] / "fixtures" / "orbit_peer_reference.json"
METHODS = ["radau", "dop853", "dp45"]
# Absolute regression budgets, not estimates of astrophysical uncertainty.
# t/r use T0/L0; angles use radians; u follows the normalized state convention.
ENDPOINT_BUDGET = np.array([1e-7, 1e-8, 1e-9, 1e-8, 1e-9, 1e-9, 1e-10, 1e-10])
TRAJECTORY_BUDGET = np.array([5e-6, 5e-8, 5e-8, 5e-8, 5e-9, 5e-9, 5e-10, 5e-10])


def load_references() -> dict:
    """Load independent, scientifically pending JSON references."""
    return json.loads(REFERENCE_PATH.read_text(encoding="utf-8"))


def scales(mass: u.Quantity) -> tuple[u.Quantity, u.Quantity]:
    """Return physical geometric scales using Astropy constants."""
    return (G * mass / c**2).to(u.m), (G * mass / c**3).to(u.s)


def public_orbit(case: dict, mass: u.Quantity = 4e6 * u.Msun):
    """Construct a public orbit sharing the reference's initial x and u."""
    length, time = scales(mass)
    state = np.asarray(case["initial_state_x_u"], dtype=float)
    metric = Kerr(mass=mass, spin=case["parameters"]["a_over_M"])
    orbit = metric.orbit(
        t=state[0] * time, R=state[1] * length,
        Theta=state[2] * u.rad, Phi=state[3] * u.rad,
        vR=state[5] / state[4] * c,
        vTheta=state[6] / state[4] * u.rad / time,
        vPhi=state[7] / state[4] * u.rad / time,
        tau=case["proper_time"][0] * time,
    )
    return orbit, length, time


def normalized_state(state, length: u.Quantity, time: u.Quantity) -> np.ndarray:
    """Read public BL x and contravariant u in the reference's unit order."""
    return np.stack((
        (state.t / time).to_value(u.one), (state.R / length).to_value(u.one),
        state.Theta.to_value(u.rad), state.Phi.to_value(u.rad),
        state.ut.to_value(u.one), (state.uR / c).to_value(u.one),
        (state.uTheta * time).to_value(u.rad),
        (state.uPhi * time).to_value(u.rad),
    ), axis=-1)


def residual(actual: np.ndarray, expected: np.ndarray) -> np.ndarray:
    """Subtract x/u after aligning only the equivalent azimuth branch."""
    result = np.asarray(actual) - np.asarray(expected)
    result[..., 3] -= 2 * np.pi * np.rint(result[..., 3] / (2 * np.pi))
    return result


def evaluate_case(
    case: dict, method: str, *, max_step: float | None = None,
    rtol: float = 1e-11, atol: float = 1e-13,
) -> dict:
    """Return component errors and elapsed public API time for one orbit.

    Endpoint integration uses stored samples. Trajectory comparison uses
    Solution.at and therefore includes separate interpolation error. The
    elapsed time includes object creation, integration and querying; this
    single run is not a performance benchmark or speed guarantee.
    """
    started = perf_counter()
    orbit, length, time = public_orbit(case)
    expected = np.asarray(case["state_x_u"])
    tau = np.asarray(case["proper_time"]) * time
    options = {"method": method, "rtol": rtol, "atol": atol}
    if max_step is not None:
        options["max_step"] = max_step * time
    initial_error = residual(normalized_state(orbit.initial, length, time), expected[0])
    solution = orbit.solve(tau_span=(tau[0], tau[-1]), **options)
    endpoint = normalized_state(solution[-1], length, time)
    endpoint_error = residual(endpoint, expected[-1])
    peer_error = residual(endpoint, np.asarray(case["pygro"]["endpoint_x_u"]))
    maxima = None
    if max_step is not None:
        sampled = solution.at(tau=tau)
        trajectory_error = residual(normalized_state(sampled, length, time), expected)
        maxima = np.max(np.abs(trajectory_error), axis=0)
    return {
        "case_id": case["id"], "method": method,
        "rtol": options["rtol"], "atol": options["atol"],
        "max_step_T0": max_step, "mass_Msun": 4e6,
        "status": solution.status, "steps": solution.integration.n_steps,
        "tau_final_T0": float((solution.tau[-1]/time).to_value(u.one)),
        "tau_final_error_T0": float(((solution.tau[-1]-tau[-1])/time).to_value(u.one)),
        "initial_roundtrip_error_x_u": initial_error.tolist(),
        "endpoint_error_x_u": endpoint_error.tolist(),
        "pygro_endpoint_error_x_u": peer_error.tolist(),
        "trajectory_max_abs_error_x_u": None if maxima is None else maxima.tolist(),
        "public_api_elapsed_seconds": perf_counter() - started,
        "endpoint_pass": bool(solution.status == 0 and np.all(np.abs(endpoint_error) <= ENDPOINT_BUDGET) and np.all(np.abs(peer_error) <= 2*ENDPOINT_BUDGET)),
        "trajectory_pass": None if maxima is None else bool(np.all(maxima <= TRAJECTORY_BUDGET)),
    }


def run() -> dict:
    """Run all methods/cases and return a JSON-friendly validation report."""
    cases = load_references()["cases"]
    endpoints = [evaluate_case(case, method) for case in cases for method in METHODS]
    trajectories = [evaluate_case(case, method, max_step=0.2)
                    for case in cases for method in METHODS]
    kerr = next(case for case in cases if case["id"] == "kerr_stable_bound")
    convergence = []
    for method in METHODS:
        coarse = evaluate_case(kerr, method, rtol=1e-8, atol=1e-10)
        fine = next(r for r in endpoints if r["case_id"] == kerr["id"] and r["method"] == method)
        half_step = evaluate_case(kerr, method, max_step=0.1)
        original_step = next(r for r in trajectories if r["case_id"] == kerr["id"] and r["method"] == method)
        endpoint_ratio = float(np.max(np.abs(fine["endpoint_error_x_u"])) /
                               np.max(np.abs(coarse["endpoint_error_x_u"])))
        time_interpolation_ratio = float(half_step["trajectory_max_abs_error_x_u"][0] /
                                         original_step["trajectory_max_abs_error_x_u"][0])
        convergence.append({
            "method": method, "coarse_endpoint": coarse,
            "fine_to_coarse_endpoint_error_ratio": endpoint_ratio,
            "half_step_trajectory": half_step,
            "half_to_original_time_interpolation_error_ratio": time_interpolation_ratio,
            "pass": endpoint_ratio < 0.1 and time_interpolation_ratio < 0.3,
        })
    return {
        "schema_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "versions": {name: metadata.version(name) for name in ("relatipy", "numpy", "scipy", "astropy")},
        "reference_sha256": hashlib.sha256(REFERENCE_PATH.read_bytes()).hexdigest(),
        "generation_command": "PYTHONPATH=relatipy/tests relatipy/.venv/bin/python -m reference.relatipy_peer_adapter --output relatipy/tests/fixtures/orbit_peer_reference_results.json",
        "manual_review": "pending", "reference_path": str(REFERENCE_PATH),
        "state_order": load_references()["state_order"],
        "endpoint_budget_x_u": ENDPOINT_BUDGET.tolist(),
        "trajectory_budget_x_u": TRAJECTORY_BUDGET.tolist(),
        "trajectory_includes_post_integration_interpolation": True,
        "interpolation_limitations": "Solver tolerance controls the integrated BL x/u, not post-integration Cartesian cubic splines or PCHIP t(tau). Circular BL motion may have very few accepted steps despite exact endpoints. Trajectory checks explicitly resolve stored spacing with max_step; no accuracy guarantee applies to an unconstrained sparse Solution.at query.",
        "budget_rationale": "Component regression budgets separate raw integration from post-integration interpolation. Endpoint accuracy is also checked by reducing rtol/atol 1000-fold; interpolation accuracy by halving max_step independently. Budgets are not certified numerical bounds, astrophysical errors, or a guarantee outside these four cases.",
        "status": "pass" if (all(r["endpoint_pass"] for r in endpoints)
                              and all(r["trajectory_pass"] for r in trajectories)
                              and all(r["pass"] for r in convergence)) else "fail",
        "endpoint_results": endpoints,
        "trajectory_results": trajectories,
        "kerr_convergence": convergence,
    }


def main() -> None:
    """Write a reproducible production-versus-peer validation report."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"{result['status']}: {args.output}")
    if result["status"] != "pass":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
