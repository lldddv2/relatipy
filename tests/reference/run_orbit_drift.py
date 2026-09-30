"""Measure public Orbit invariant drift on native accepted steps.

Run from the repository with ``python -m tests.reference.run_orbit_drift``.
The checked cases span five radial cycles of an inclined eccentric Kerr
orbit and five circular Schwarzschild revolutions.  Results describe these
cases and tolerances only; local solver tolerances are not global bounds.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sysconfig
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
from time import perf_counter

import numpy as np
from astropy import units as u
from astropy.constants import G, c

from relatipy.metrics.kerr import Kerr
from .timelike_invariants import public_canonical_state, public_orbit_from_canonical, timelike_invariants


RESULTS_PATH = Path(__file__).resolve().parents[1] / "fixtures" / "orbit_invariant_drift_results.json"
METHODS = ("radau", "dop853", "dp45")
TOLERANCES = {"coarse": (1e-6, 1e-8), "fine": (1e-10, 1e-12)}


def drift_cases() -> list[dict]:
    """Return fixed independently generated initial states and durations.

    KerrGeoPy 0.9.3 StableOrbit(a=.5,p=8,e=.2,x=.8), phases
    (0,.4,1.1,.2), supplies the Kerr initial state.  The duration maps five
    Mino radial periods using dτ/dλ=Σ, DOP853 rtol=2.3e-14, atol=1e-13,
    as implemented in peer_oracles.kerrgeopy_bound_reference.  Freeze those
    values so production validation does not require optional peer packages.
    Schwarzschild initial state and period are the exact circular solution.
    """
    radius = 10.0
    ut = 1 / np.sqrt(1 - 3 / radius)
    uphi = 1 / np.sqrt(radius**3 * (1 - 3 / radius))
    return [
        {
            "id": "kerr_inclined_eccentric_five_radial_periods",
            "spin": 0.5,
            "parameters": {"p": 8.0, "e": 0.2, "x": 0.8, "initial_phases": [0.0, 0.4, 1.1, 0.2]},
            "radial_periods": 5,
            "circular_periods": 0,
            "initial": [0.0, 6.747729739587268, 1.2951155714852487, 0.2,
                        1.331910837833271, 0.019334179047186116,
                        0.041179703326295995, 0.06784624037067309],
            "tau_final": 977.489445558449,
            "expected_initial": {"norm": -1.0, "energy": 0.9466042517710911,
                                 "Lz": 2.6976027753892358, "CarterQ": 4.102701297903183},
            "reference": "KerrGeoPy 0.9.3 analytic trajectory, independent initial x,u and proper-time mapping",
        },
        {
            "id": "schwarzschild_circular_r10_five_revolutions",
            "spin": 0.0,
            "parameters": {"radius": radius},
            "radial_periods": 0,
            "circular_periods": 5,
            "initial": [0.0, radius, np.pi / 2, 0.0, ut, 0.0, 0.0, uphi],
            "tau_final": float(5 * 2 * np.pi / uphi),
            "expected_initial": {"norm": -1.0, "energy": float((1-2/radius)*ut),
                                 "Lz": float(radius**2*uphi), "CarterQ": 0.0},
            "reference": "Exact equatorial circular Schwarzschild x,u and proper period 2*pi/uphi",
        },
    ]


def measure_case(case: dict, method: str, tolerance: str) -> dict:
    """Measure invariant errors on saved native steps without interpolation."""
    metric = Kerr(mass=1 * u.Msun, spin=case["spin"])
    orbit = public_orbit_from_canonical(metric, np.array(case["initial"]))
    time = (G * metric.mass / c**3).to(u.s)
    rtol, atol = TOLERANCES[tolerance]
    started = perf_counter()
    solution = orbit.solve(tau_span=(0 * u.s, case["tau_final"] * time),
                           method=method, rtol=rtol, atol=atol)
    runtime = perf_counter() - started
    states = public_canonical_state(solution, metric.mass)
    invariants = timelike_invariants(case["spin"], states)
    drift = {}
    for name, values in invariants.items():
        error = values - values[0]
        scale = max(1.0, abs(float(values[0])))
        drift[name] = {
            "initial": float(values[0]), "final": float(values[-1]),
            "max_abs_drift": float(np.max(np.abs(error))),
            "endpoint_signed_drift": float(error[-1]),
            "endpoint_abs_drift": float(abs(error[-1])),
            "max_scaled_drift": float(np.max(np.abs(error)) / scale),
            "scale": scale,
        }
    norm_error = invariants["norm"] + 1
    radial_velocity = states[:, 5]
    periapses = (int(np.count_nonzero((radial_velocity[:-1] < 0) & (radial_velocity[1:] >= 0)))
                 if case["radial_periods"] else None)
    return {
        "case_id": case["id"], "method": method, "tolerance": tolerance,
        "rtol": rtol, "atol": atol, "status": solution.status,
        "tau_requested": case["tau_final"],
        "tau_reached": float((solution.tau[-1] / time).to_value(u.one)),
        "runtime_seconds": runtime, "accepted_steps": solution.integration.n_steps,
        "rhs_evaluations": solution.integration.nfev, "stored_samples": len(solution),
        "sampling": "native accepted steps plus initial state; tau_eval omitted; no postprocessing interpolation",
        "invariants": drift,
        "max_abs_norm_plus_one": float(np.max(np.abs(norm_error))),
        "max_abs_initial_canonical_roundtrip_error": float(np.max(np.abs(states[0] - case["initial"]))),
        "observed_periapsis_crossings": periapses,
        "observed_azimuthal_revolutions": float((states[-1, 3] - states[0, 3]) / (2*np.pi)),
        "domain": {"radius_min": float(states[:,1].min()), "radius_max": float(states[:,1].max()),
                   "theta_min": float(states[:,2].min()), "theta_max": float(states[:,2].max()),
                   "minimum_abs_sin_theta": float(np.abs(np.sin(states[:,2])).min()),
                   "outer_horizon": float(1+np.sqrt(1-case["spin"]**2))},
    }


def measurement_pair_errors(case: dict, coarse: dict, fine: dict) -> list[str]:
    """Return failures of the declared fixed-case gates for one native method.

    The same gates are exposed to the standalone runner and pytest.  They
    assess measured invariant drift and domain coverage; they do not derive
    expected invariants or assert general global-error guarantees.
    """
    errors = []
    for row in (coarse, fine):
        label = row["tolerance"]
        if row["status"] != 0:
            errors.append(f"{label}: integration did not reach requested time")
        if not np.isclose(row["tau_reached"], case["tau_final"], rtol=2e-15, atol=0):
            errors.append(f"{label}: final proper time differs")
        if row["stored_samples"] != row["accepted_steps"] + 1:
            errors.append(f"{label}: samples do not match accepted native steps")
        if row["rhs_evaluations"] < row["accepted_steps"]:
            errors.append(f"{label}: invalid evaluation count")
        if row["max_abs_initial_canonical_roundtrip_error"] >= 5e-13:
            errors.append(f"{label}: initial canonical state changed")
        if row["domain"]["radius_min"] <= row["domain"]["outer_horizon"] + 1:
            errors.append(f"{label}: case left selected exterior domain")
        if row["domain"]["minimum_abs_sin_theta"] <= 0.79:
            errors.append(f"{label}: case left selected angular domain")
        if case["radial_periods"]:
            if row["observed_periapsis_crossings"] != case["radial_periods"]:
                errors.append(f"{label}: radial cycle count differs")
        elif not np.isclose(row["observed_azimuthal_revolutions"], case["circular_periods"],
                            rtol=2e-8, atol=2e-8):
            errors.append(f"{label}: circular revolution count differs")
    if fine["max_abs_norm_plus_one"] >= 2e-8:
        errors.append("fine: normalization exceeds gate")
    for name, measurement in fine["invariants"].items():
        if measurement["max_scaled_drift"] >= 2e-8:
            errors.append(f"fine: {name} drift exceeds gate")
        if measurement["endpoint_abs_drift"] > measurement["max_abs_drift"]:
            errors.append(f"fine: {name} endpoint exceeds maximum drift")
        # Exact circular solutions approach floating-point roundoff, so a
        # strict tolerance-convergence ratio is meaningful only for Kerr here.
        if case["radial_periods"] and measurement["max_abs_drift"] > max(
            3e-13, 0.05 * coarse["invariants"][name]["max_abs_drift"]
        ):
            errors.append(f"fine: {name} lacks required tolerance convergence")
    return errors


def run() -> dict:
    """Return JSON-friendly measurements and executable fixed-domain gates."""
    cases = drift_cases()
    rows = [measure_case(case, method, tolerance)
            for case in cases for method in METHODS for tolerance in TOLERANCES]
    from relatipy import _core
    extension = Path(_core.__file__)
    pairs = {(row["case_id"], row["method"], row["tolerance"]): row for row in rows}
    validations = [
        {"case_id": case["id"], "method": method,
         "errors": measurement_pair_errors(case, pairs[case["id"], method, "coarse"],
                                             pairs[case["id"], method, "fine"])}
        for case in cases for method in METHODS
    ]
    report = {
        "schema_version": "1.0", "manual_review": "pending",
        "status": "pass" if not any(item["errors"] for item in validations) else "fail",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "Selected-domain production Orbit invariant drift and tolerance convergence",
        "conventions": {"units": "G=c=M=mu=1", "state": "(t,r,theta,phi,ut,ur,utheta,uphi)",
                        "signature": "(-,+,+,+)", "CarterQ": "Fixed mu squared one, not K"},
        "formula_source": {"record_id": "INF-022", "claim_ids": ["C001", "C002", "C003", "C004"], "source_id": "S001",
                           "url": "https://arxiv.org/pdf/gr-qc/0202090",
                           "locator": "Schmidt 2002 Eq.(3),(4),(8),(17),(18); pp.2-5",
                           "operation": "Invert printed inverse metric; lower u; algebraically solve theta potential for Q"},
        "environment": {"python": platform.python_version(), "platform": platform.platform(),
                        "python_compiler": platform.python_compiler(),
                        "python_build_CC": sysconfig.get_config_var("CC"),
                        "python_build_CFLAGS": sysconfig.get_config_var("CFLAGS"),
                        "declared_extension_compile_args": ["-std=c11", "-O2"],
                        "compile_args_provenance": "setup.py declaration; actual build invocation is not introspectable",
                        "native_extension_sha256": hashlib.sha256(extension.read_bytes()).hexdigest(),
                        **{name: version(name) for name in ("numpy", "astropy", "scipy", "relatipy")}},
        "cases": cases, "results": rows, "validations": validations,
        "gates": {"fine_max_scaled_drift": 2e-8, "fine_max_abs_norm_plus_one": 2e-8,
                  "kerr_fine_to_coarse_drift_ratio_max": 0.05, "convergence_roundoff_floor": 3e-13,
                  "initial_canonical_roundtrip_max_abs": 5e-13,
                  "cycle_coverage": "Five periapsis crossings for eccentric Kerr; five azimuthal revolutions for circular Schwarzschild",
                  "scope": "Fixed cases only; no strict circular convergence gate at roundoff"},
        "limitations": ["Selected spins and regular bound BL orbits only; excludes horizon, axis and separatrix.",
                        "Invariant conservation alone does not bound phase or coordinate error.",
                        "Runtime is one execution including public output reconstruction, not a benchmark guarantee.",
                        "n_steps and nfev are available publicly; rejected-step and Jacobian counts are not exposed.",
                        "Scientific source and computed results remain pending human review."],
    }
    return report


def main() -> None:
    """Export the reproducible report; runtime is descriptive only."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=RESULTS_PATH)
    args = parser.parse_args()
    report = run()
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    for row in report["results"]:
        print(row["case_id"], row["method"], row["tolerance"],
              "max_scaled_drift=", max(item["max_scaled_drift"] for item in row["invariants"].values()),
              "steps=", row["accepted_steps"], "seconds=", round(row["runtime_seconds"], 4))
    if report["status"] != "pass":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
