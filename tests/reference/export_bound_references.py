"""Freeze independent KerrGeoPy states for the bound-orbit input contract.

Run from the repository parent with ``thesis_apply/.venv/bin/python``.  This
exporter never imports RelatiPy, and KerrGeoPy is not a runtime dependency.
All saved states use ``G = c = M = 1`` and Boyer--Lindquist ``(x, u)`` order.
"""

from __future__ import annotations

import argparse
from importlib import metadata
import json
from pathlib import Path

import kerrgeopy as kg
import numpy as np


CASES = (
    ("inclined_prograde", 0.7, 12.0, 0.3, 0.7, 0.4, 1.2, 0.2),
    ("inclined_retrograde", 0.7, 12.0, 0.2, -0.6, 1.0, 0.7, -0.5),
    ("inward_southbound", 0.6, 13.0, 0.35, 0.65, 4.2, 4.0, 1.7),
    ("near_separatrix", 0.7, 4.8, 0.3, 0.7, 0.8, 0.9, 0.4),
    ("circular_inclined", 0.5, 10.0, 0.0, 0.8, 0.0, 0.6, 0.3),
    ("eccentric_equatorial", 0.5, 12.0, 0.3, 1.0, 0.6, 0.0, -0.2),
    ("schwarzschild_inclined", 0.0, 20.0, 0.1, 0.5, 0.3, 0.9, 0.4),
    ("weak_field", 0.5, 10000.0, 0.3, 0.8, 0.4, 0.5, 0.6),
)


def export() -> dict:
    """Evaluate each stable bound orbit at Mino time zero."""
    cases = []
    for name, spin, p, e, x, q_r0, q_theta0, q_phi0 in CASES:
        phases = (0.0, q_r0, q_theta0, q_phi0)
        orbit = kg.StableOrbit(spin, p, e, x, initial_phases=phases)
        coordinates = np.array([float(component(0.0)) for component in orbit.trajectory()])
        velocity = np.array([float(component(0.0)) for component in orbit.four_velocity()])
        state = np.concatenate((coordinates, velocity))
        if not np.all(np.isfinite(state)):
            raise ValueError(f"non-finite reference state: {name}")
        np.testing.assert_allclose(coordinates[0], 0.0, rtol=0, atol=5e-14)
        np.testing.assert_allclose(coordinates[3], q_phi0, rtol=0, atol=5e-14)
        cases.append({
            "id": name,
            "parameters": {
                "spin": spin, "p_over_M": p, "e": e, "x": x,
                "q_r0_rad": q_r0, "q_theta0_rad": q_theta0,
                "q_phi0_rad": q_phi0,
            },
            "initial_x_u": state.tolist(),
            "constants": {"energy": float(orbit.E), "Lz": float(orbit.L),
                          "CarterQ": float(orbit.Q)},
        })
    return {
        "schema_version": "1.0",
        "manual_review": "pending",
        "reference": "KerrGeoPy StableOrbit.trajectory and four_velocity at Mino lambda=0",
        "kerrgeopy_version": metadata.version("kerrgeopy"),
        "generation_command": "MPLCONFIGDIR=/tmp/relatipy-matplotlib thesis_apply/.venv/bin/python relatipy/tests/reference/export_bound_references.py --output relatipy/tests/fixtures/bound_orbit_reference.json",
        "conventions": {
            "units": "G = c = M = 1; p is measured in gravitational radii GM/c^2",
            "coordinates": "Boyer-Lindquist (t, r, theta, phi)",
            "state_order": ["t", "r", "theta", "phi", "ut", "ur", "utheta", "uphi"],
            "velocity": "contravariant four-velocity dx^mu/dtau",
            "phase_order": ["q_t0", "q_r0", "q_theta0", "q_phi0"],
            "phase_origin": "q_t0=0; KerrGeoPy subtracts radial and polar periodic offsets so t(0)=0 and phi(0)=q_phi0",
            "x": "cosine of inclination; sign selects prograde or retrograde branch",
        },
        "provenance": {
            "record_id": "INF-010", "claim_ids": ["C002", "C003", "C005"],
            "source_ids": ["S004", "S005"], "manual_review": "pending",
            "source_url": "https://kerrgeopy.readthedocs.io/en/latest/notebooks/Trajectory.html",
            "locator": "KerrGeoPy 0.9.3 stable.py stable_trajectory: initial phases and t/phi normalization; StableOrbit.trajectory/four_velocity; INF-010 record.json claims C002/C003/C005",
        },
        "independence": "Exporter imports KerrGeoPy and NumPy only; no RelatiPy production implementation is used.",
        "cases": cases,
    }


def main() -> None:
    """Write the checked reference artifact as JSON."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(export(), indent=2, allow_nan=False) + "\n",
                           encoding="utf-8")


if __name__ == "__main__":
    main()
