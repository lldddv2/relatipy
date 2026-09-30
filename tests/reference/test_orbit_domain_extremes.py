"""Public BL domain regressions with independent metric and exact infall checks.

These finite cases cover valid moving near-axis states, extreme spin, and
near-horizon conditioning. They do not promise accuracy on a singular chart.
"""

import numpy as np
import pytest

from .run_orbit_domain_extremes import (
    METHODS,
    POLAR_DISTANCES,
    SPINS,
    horizon_failure_case,
    near_horizon_case,
    near_polar_case,
    radial_infall_case,
    rejected_initial_conditions,
    stage_failure_case,
)


@pytest.mark.parametrize("spin", SPINS)
@pytest.mark.parametrize("pole", ("north", "south"))
@pytest.mark.parametrize("distance", POLAR_DISTANCES)
def test_moving_near_polar_orbit_has_independent_invariants(spin, pole, distance):
    results = [near_polar_case(spin, pole, distance, method) for method in METHODS]
    endpoints = np.array([item["endpoint_normalized"] for item in results])
    # Shared RHS: this comparison is numerical, not an independent physics oracle.
    np.testing.assert_allclose(endpoints, np.broadcast_to(endpoints[0], endpoints.shape),
                               rtol=2e-9, atol=2e-8)


@pytest.mark.parametrize("spin", SPINS)
def test_short_near_horizon_orbits_remain_timelike_and_exterior(spin):
    results = [near_horizon_case(spin, method) for method in METHODS]
    endpoints = np.array([item["endpoint_normalized"] for item in results])
    np.testing.assert_allclose(endpoints, np.broadcast_to(endpoints[0], endpoints.shape),
                               rtol=2e-9, atol=2e-8)


@pytest.mark.parametrize("method", METHODS)
def test_near_horizon_schwarzschild_infall_matches_exact_radius(method):
    radial_infall_case(method)


@pytest.mark.parametrize("method", METHODS)
def test_true_polar_stage_crossing_retains_initial_state_as_failure(method):
    stage_failure_case(method)


@pytest.mark.parametrize("method", METHODS)
def test_horizon_chart_failure_returns_reconstructible_partial_solution(method):
    horizon_failure_case(method)


def test_public_axis_spin_and_horizon_rejections():
    assert len(rejected_initial_conditions()) == 45
