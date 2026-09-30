"""Export independent x/u references without importing RelatiPy.

Run with the environment containing KerrGeoPy and PyGRO. JSON permits
production tests to use another Python ABI. States use BL order
``(t, r, theta, phi, ut, ur, utheta, uphi)`` with ``G = c = M = 1``.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import platform
from pathlib import Path
import sys

import numpy as np

from .peer_oracles import (
    FIXTURE_PATH, ReferenceTrajectory, exact_schwarzschild_circular_reference,
    integrate_with_pygro, kerrgeopy_bound_reference, load_cases,
    make_pygro_engine, software_versions,
)


def circular_reference(radius: float, samples: int = 257) -> ReferenceTrajectory:
    """Evaluate the exact circular Schwarzschild solution across one period."""
    ends = exact_schwarzschild_circular_reference(radius)
    tau = np.linspace(0.0, ends.proper_time[-1], samples)
    velocity = ends.four_velocity[:, 0]
    coordinates = ends.coordinates[:, :1] + velocity[:, None] * tau
    return ReferenceTrajectory(
        tau / radius**2, tau, coordinates,
        np.broadcast_to(velocity[:, None], coordinates.shape).copy(),
    )


def check_proper_time_mapping(case: dict, reference: ReferenceTrajectory) -> dict:
    """Cross-check the ODE clock with independent adaptive quadrature.

    Both calculations integrate KerrGeoPy's analytic Sigma, not a RelatiPy
    geodesic RHS. Seventeen declared mesh points include both endpoints.
    QUADPACK's estimate is recorded as an estimate, not a certified bound.
    """
    import kerrgeopy as kg
    from scipy.integrate import quad

    orbit = kg.StableOrbit(
        case["a_over_M"], case["p_over_M"], case["eccentricity"], case["x"],
        initial_phases=tuple(case["initial_phases"]),
    )
    trajectory = orbit.trajectory()
    indices = np.unique(np.linspace(0, reference.proper_time.size - 1, 17, dtype=int))
    differences, estimates = [], []
    for index in indices:
        value, estimate = quad(
            lambda lam: trajectory[1](lam)**2 + case["a_over_M"]**2 * np.cos(trajectory[2](lam))**2,
            0.0, reference.mino_time[index], epsabs=1e-10, epsrel=2.3e-14,
        )
        differences.append(value - reference.proper_time[index])
        estimates.append(estimate)
    return {
        "method": "scipy.integrate.quad / QUADPACK",
        "sample_indices": indices.tolist(),
        "proper_time_difference_at_checks": differences,
        "max_abs_difference_T0": float(np.max(np.abs(differences))),
        "max_quad_estimated_absolute_error_T0": float(np.max(estimates)),
        "epsabs": 1e-10, "epsrel": 2.3e-14,
    }


def export() -> dict:
    """Compute analytic references and independent PyGRO endpoint checks."""
    manifest = load_cases()
    metric, engine = make_pygro_engine()
    cases = []
    for case in manifest["cases"]:
        circular = case["kind"] == "circular"
        reference = (circular_reference(case["radius_over_M"])
                     if circular else kerrgeopy_bound_reference(case))
        goal = max(11, case["pygro_goal"])
        peer = integrate_with_pygro(metric, engine, reference,
                                    spin=case["a_over_M"], goal=goal)
        peer_endpoint = np.concatenate((
            reference.coordinates[:, -1] + peer.endpoint_error,
            peer.four_velocity[-1],
        ))
        final = np.concatenate((reference.coordinates[:, -1], reference.four_velocity[:, -1]))
        cases.append({
            "id": case["id"], "parameters": case,
            "primary_reference": "exact_schwarzschild_circular" if circular else "kerrgeopy",
            "proper_time": reference.proper_time.tolist(),
            "mino_time": reference.mino_time.tolist(),
            "state_x_u": np.concatenate((reference.coordinates, reference.four_velocity)).T.tolist(),
            "initial_state_x_u": reference.initial_state.tolist(),
            "pygro": {
                "integrator": "dp853", "accuracy_goal": goal, "precision_goal": goal,
                "initial_step": manifest["mapping"]["pygro_initial_step_M"],
                "endpoint_proper_time": float(peer.proper_time[-1]),
                "endpoint_x_u": peer_endpoint.tolist(),
                "endpoint_error_against_primary_x_u": (peer_endpoint-final).tolist(),
                "trajectory_coordinate_error_max_abs": np.max(np.abs(peer.trajectory_error), axis=1).tolist(),
                "normalization_error_max_abs": peer.max_abs_normalization_error,
                "steps": int(peer.proper_time.size - 1),
            },
            "reference_precision": {
                "representation": "IEEE-754 binary64 serialized with round-trip precision",
                "tau_mapping": None if circular else manifest["mapping"],
                "tau_mapping_quadrature_crosscheck": None if circular else check_proper_time_mapping(case, reference),
                "limits": "Analytic floating-point evaluation and numerical Mino-to-proper-time quadrature; no certified global error bound. PyGRO discrepancy is recorded independently, not treated as uncertainty.",
            },
        })
    return {
        "schema_version": "1.0", "manual_review": "pending",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "python": sys.version, "platform": platform.platform(),
        "versions": software_versions(),
        "generation_command": "PYTHONPATH=relatipy/tests MPLCONFIGDIR=/tmp/relatipy-matplotlib thesis_apply/.venv/bin/python -m reference.export_orbit_references --output relatipy/tests/fixtures/orbit_peer_reference.json",
        "manifest_sha256": hashlib.sha256(FIXTURE_PATH.read_bytes()).hexdigest(),
        "exporter_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "peer_oracles_sha256": hashlib.sha256(Path(__file__).with_name("peer_oracles.py").read_bytes()).hexdigest(),
        "conventions": manifest["conventions"],
        "state_order": ["t", "r", "theta", "phi", "ut", "ur", "utheta", "uphi"],
        "provenance": [
            {**manifest["provenance"][0], "source_ids": ["S001", "S002", "S003", "S004"],
             "locator": "INF-006 record.json: claims C001/C002/C004 evidence; KerrGeoPy Trajectory: Mino Time and Four Velocity"},
            {**manifest["provenance"][1], "source_ids": ["S001"],
             "locator": "INF-007 record.json: claims C002/C003/C007 evidence; PyGRO paper pp. 6-7 sec. 3.1-3.2 and pp. 12-13 sec. 3.6"},
            {
            "record_id": "INF-013", "manual_review": "pending",
            "source_ids": ["S002", "S003", "S004"],
            "claim_ids": ["C001", "C002", "C003", "C004"],
            "locator": "peer_oracles.py: exact_schwarzschild_circular_reference, kerrgeopy_bound_reference, integrate_with_pygro; record.json: claims C001-C004 and sources S002-S004",
            },
        ],
        "independence": "The exporter never imports RelatiPy. KerrGeoPy supplies analytic trajectories; PyGRO integrates its bundled metric. Circular trajectories use exact Schwarzschild formulas. No production RHS is reused.",
        "cases": cases,
    }


def main() -> None:
    """Write the reproducible cross-ABI reference artifact."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(export(), indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
