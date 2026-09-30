"""Public projected Radau checks using independent Kerr invariants."""

from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest
from astropy import units as u
from astropy.constants import G, c

from relatipy.metrics.kerr import Kerr

from .run_orbit_domain_extremes import (
    horizon_failure_case,
    near_horizon_case,
    near_polar_case,
)
from .run_orbit_drift import drift_cases
from .timelike_invariants import (
    public_canonical_state,
    public_orbit_from_canonical,
    timelike_invariants,
)


def solve_case(case: dict, method: str, rtol: float, atol: float):
    """Integrate one fixed Kerr state without interpolated output samples."""
    metric = Kerr(mass=1 * u.Msun, spin=case["spin"])
    orbit = public_orbit_from_canonical(metric, np.asarray(case["initial"]))
    time = (G * metric.mass / c**3).to(u.s)
    solution = orbit.solve(
        tau_span=(0 * u.s, case["tau_final"] * time),
        method=method, rtol=rtol, atol=atol,
    )
    return solution, public_canonical_state(solution, metric.mass)


@pytest.mark.parametrize("case", drift_cases(), ids=lambda case: case["id"])
def test_projected_steps_preserve_invariants_and_match_tight_reference(case):
    """Check saved steps and phase against tighter explicit integration."""
    solution, states = solve_case(case, "projection_radau", 1e-10, 1e-12)
    reference, reference_states = solve_case(case, "dop853", 1e-13, 1e-15)
    assert solution.success and reference.success
    assert solution.integration.method == "projection_radau"
    assert len(states) == solution.integration.n_steps + 1
    invariants = timelike_invariants(case["spin"], states)
    assert np.max(np.abs(invariants["norm"] + 1.0)) < 2e-12
    for name in ("energy", "Lz", "CarterQ"):
        scale = max(1.0, abs(float(invariants[name][0])))
        assert np.max(np.abs(invariants[name] - invariants[name][0])) < 2e-12 * scale
    endpoint_tolerance = 5e-7 if case["radial_periods"] else 1e-9
    assert np.max(np.abs(states[-1] - reference_states[-1])) < endpoint_tolerance
    if case["radial_periods"]:
        radial_velocity = states[:, 5]
        turns = np.count_nonzero(
            (radial_velocity[:-1] < 0) & (radial_velocity[1:] >= 0)
        )
        assert turns == case["radial_periods"]


def test_projected_method_preserves_domain_and_horizon_failure_contract():
    """Exercise near-axis and near-horizon cases plus partial failure output."""
    assert near_polar_case(1.0, "north", 1e-8, "projection_radau")["status"] == 0
    assert near_horizon_case(1.0, "projection_radau")["status"] == 0
    horizon_failure_case("projection_radau")


def test_projected_calls_keep_contexts_independent():
    """Parallel Kerr contexts reproduce their serial accepted trajectories."""
    cases = drift_cases()

    def solve(case):
        return solve_case(case, "projection_radau", 1e-9, 1e-12)[1]

    serial = [solve(case) for case in cases]
    with ThreadPoolExecutor(max_workers=2) as pool:
        concurrent = list(pool.map(solve, cases))
    for actual, expected in zip(concurrent, serial, strict=True):
        np.testing.assert_array_equal(actual, expected)
