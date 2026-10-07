"""Null Schwarzschild trajectories against frozen INF-030:D001/D002.

Manual scientific review remains pending. Units are G = c = M = 1. The
radial comparison excludes the final time-cut interpolant, especially for
DOP853; its endpoint error is bounded separately. Numerical budgets record
actual calibration measurements in fixtures/null/schwarzschild.json.
"""

import numpy as np
import pytest
from astropy import units as u

from .null_reference import (
    canonical, constants_ray, load_fixture, make_metric, scales, solve_ray,
)


REFERENCE = load_fixture("schwarzschild.json")
METHODS = ("radau", "dop853", "dp45")
CONSTANTS = {
    row["quantity"]: float(row["value"])
    for row in REFERENCE["provenance"][0]["rows"]
}


def _radial_time(radius):
    """INF-030:C005 ingoing primitive, evaluated with extended precision."""
    radius = np.asarray(radius, dtype=np.longdouble)
    r0 = np.longdouble(REFERENCE["radial"]["r0"])
    return r0 - radius + 2 * np.log((r0 - 2) / (radius - 2))


def _tangent_ray(metric, radius):
    """Release an equatorial photon with zero radial coordinate velocity."""
    length, _ = scales(metric)
    return metric.null(
        R=radius * length, Theta=np.pi / 2 * u.rad, Phi=0 * u.rad,
        vPhi=1 * u.rad / u.s,
    )


@pytest.mark.parametrize("method", METHODS)
def test_null_radial_infall_matches_schwarzschild_time(method):
    """INF-030:C004/C005/C017, D002: ingoing t(r) from r0=50 to r=2.1."""
    case = REFERENCE["radial"]
    measured = case["methods"][method]
    ray = constants_ray(make_metric(0), case["r0"], 0, 0)
    solution = solve_ray(ray, case["t_final"], method=method, **REFERENCE["solver"])
    assert solution.success and solution.status == 0
    assert solution.termination is None  # Stop before the chart's horizon margin.
    states = canonical(solution)
    assert len(states) > 2
    assert np.all(np.diff(states[:, 1]) < 0)
    assert np.all(states[:, 1] > CONSTANTS["r_horizon"])

    # No t_eval. All rows except the last are accepted steps; the last row is
    # a quintic Hermite time-cut state, whose error is a distinct measurement.
    accepted = states[:-1]
    error = np.max(np.abs(accepted[:, 0].astype(np.longdouble)
                          - _radial_time(accepted[:, 1])))
    assert error <= measured["accepted_t_tolerance"]
    assert abs(states[-1, 1] - case["r_stop"]) <= measured["endpoint_r_tolerance"]

    # Freeze the D002 oracle as decimal source rows, rather than requiring the
    # brain in CI. Check their formula with the same physical error budget.
    for row in REFERENCE["provenance"][1]["rows"]:
        difference = abs(_radial_time(np.longdouble(row["r_over_M"]))
                         - np.longdouble(row["t_minus_t0_over_M"]))
        assert difference <= measured["accepted_t_tolerance"]


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("side, reason", [(-1, "horizon"), (1, "escape")])
def test_null_critical_impact_separates_capture_and_escape(method, side, reason):
    """INF-030:C002/C003, D001: inward b=bc(1±0.001) from r=100."""
    case = REFERENCE["capture"]
    b = CONSTANTS["b_critical"] * (1 + side * case["delta"])
    ray = constants_ray(make_metric(0), case["r0"], b, 0)
    solution = solve_ray(
        ray, case["t_final"], method=method, r_escape=case["r_escape"],
        **REFERENCE["solver"],
    )
    assert solution.success and solution.status == 1
    assert solution.termination.reason == reason
    states = canonical(solution)
    assert states[0, 5] < 0
    if reason == "horizon":
        # The horizon event retains the preceding exterior state. Its distance
        # budget is 10× the largest measured excess over r_h, not a crossing.
        excess = states[-1, 1] - CONSTANTS["r_horizon"]
        assert 0 < excess <= case["horizon_excess_tolerance"]
        assert states[-1, 5] < 0
    else:
        assert np.min(states[:, 1]) > CONSTANTS["r_photon"]
        assert states[-1, 1] >= case["r_escape"]
        assert states[-1, 5] > 0


@pytest.mark.parametrize("method", METHODS)
def test_null_tangential_photon_stays_on_schwarzschild_sphere(method):
    """INF-030:C001/C011, D001: r=3 tangent stays circular for t≤10."""
    case = REFERENCE["circle"]
    ray = _tangent_ray(make_metric(0), CONSTANTS["r_photon"])
    solution = solve_ray(ray, case["t_final"], method=method, **REFERENCE["solver"])
    assert solution.success and solution.status == 0
    error = np.max(np.abs(canonical(solution)[:, 1] - CONSTANTS["r_photon"]))
    # Measured initial effective perturbation and radial drift are both zero.
    # exp(lambda*10)=6.85 amplifies roundoff; the budget uses 10× the largest
    # perturbed-companion discrepancy (5.14e-12), above that roundoff floor.
    assert error <= case["methods"][method]["exact_r_tolerance"]


@pytest.mark.parametrize("method", METHODS)
def test_null_photon_sphere_perturbation_has_expected_short_time_growth(method):
    """INF-030:C001/C011: own linearization gives δr''=δr/27 in t."""
    case = REFERENCE["circle"]
    measured = case["methods"][method]
    dr = case["perturbed_dr"]
    ray = _tangent_ray(make_metric(0), CONSTANTS["r_photon"] + dr)
    solution = solve_ray(ray, case["t_final"], method=method, **REFERENCE["solver"])
    assert solution.success and solution.status == 0
    states = canonical(solution)
    # Tangential release has δr'(0)=0: both stable and unstable linear modes
    # contribute, giving cosh(lambda*t), rather than a pure exponential.
    # lambda=1/(3*sqrt(3)); the fixture records the C011-based derivation and
    # the measured effective perturbation (~2e-7 including radial forcing).
    expected = CONSTANTS["r_photon"] + dr * np.cosh(case["lyapunov"] * states[:, 0])
    assert np.max(np.abs(states[:, 1] - expected)) <= measured["linearization_tolerance"]
    assert states[-1, 1] - CONSTANTS["r_photon"] > 3 * dr
    assert np.max(np.abs(states[:, 1] - CONSTANTS["r_photon"])) <= measured["perturbed_r_tolerance"]
