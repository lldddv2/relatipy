"""Public production orbits against independent analytic and PyGRO states.

Four cases run for Radau, DOP853 and DP45 without optional dependency skips.
Peers run during export, potentially under another Python ABI. Endpoint and
post-integration interpolation errors are checked separately.
"""

from __future__ import annotations

import numpy as np
import pytest
from astropy import units as u

from .relatipy_peer_adapter import (
    ENDPOINT_BUDGET, TRAJECTORY_BUDGET, METHODS, evaluate_case,
    load_references, normalized_state, public_orbit,
)

REFERENCES = load_references()
CASES = REFERENCES["cases"]


def test_reference_provenance_and_shared_proper_time():
    """Keep cross-ABI provenance and pending scientific review explicit."""
    assert REFERENCES["manual_review"] == "pending"
    assert REFERENCES["versions"]["kerrgeopy"] == "0.9.3"
    assert REFERENCES["versions"]["pygro"] == "1.0.3"
    assert {item["record_id"] for item in REFERENCES["provenance"]} == {"INF-006", "INF-007", "INF-013"}
    for case in CASES:
        tau = np.asarray(case["proper_time"])
        assert np.all(np.diff(tau) > 0)
        assert np.asarray(case["state_x_u"]).shape == (len(tau), 8)
        assert abs(case["pygro"]["endpoint_proper_time"] - tau[-1]) < 1e-10
        errors = np.abs(case["pygro"]["endpoint_error_against_primary_x_u"])
        assert np.all(errors <= ENDPOINT_BUDGET), errors
        crosscheck = case["reference_precision"]["tau_mapping_quadrature_crosscheck"]
        if crosscheck is not None:
            assert crosscheck["max_abs_difference_T0"] < 1e-8


@pytest.mark.parametrize("case", CASES, ids=lambda case: case["id"])
@pytest.mark.parametrize("method", METHODS)
def test_public_orbit_endpoint_against_independent_references(case, method):
    """Check raw integrated endpoints independently of trajectory splines."""
    result = evaluate_case(case, method)
    assert result["status"] == 0
    assert result["tau_final_error_T0"] == 0
    assert np.max(np.abs(result["initial_roundtrip_error_x_u"])) <= 2e-12
    assert result["endpoint_pass"], result


@pytest.mark.parametrize("case", CASES, ids=lambda case: case["id"])
@pytest.mark.parametrize("method", METHODS)
def test_public_orbit_trajectory_including_interpolation(case, method):
    """Resolve stored samples before testing Solution.at against the oracle.

    PCHIP coordinate time has its own interpolation error. max_step=0.2 T0
    controls sample spacing independently of integration tolerance.
    """
    result = evaluate_case(case, method, max_step=0.2)
    assert result["status"] == 0
    assert result["trajectory_pass"], result
    assert np.all(np.asarray(result["trajectory_max_abs_error_x_u"]) <= TRAJECTORY_BUDGET)


def test_physical_mass_and_units_preserve_normalized_kerr_orbit():
    """Exercise Astropy conversion at two substantially different masses."""
    case = next(case for case in CASES if case["id"] == "kerr_stable_bound")
    endpoints = []
    for mass in (1 * u.Msun, (4e6 * u.Msun).to(u.kg)):
        orbit, length, time = public_orbit(case, mass)
        end = (case["proper_time"][-1] * time).to(u.day)
        orbit.integrate(end, method="dop853", rtol=1e-11, atol=1e-13)
        assert orbit.tau == end
        endpoints.append(normalized_state(orbit, length, time))
    np.testing.assert_allclose(endpoints[0], endpoints[1], rtol=1e-10, atol=1e-10)


@pytest.mark.parametrize("method", METHODS)
def test_endpoint_tolerance_convergence_against_kerrgeopy(method):
    """Tighter integration controls reduce independent endpoint error."""
    case = next(case for case in CASES if case["id"] == "kerr_stable_bound")
    coarse = evaluate_case(case, method, rtol=1e-8, atol=1e-10)
    fine = evaluate_case(case, method)
    coarse_error = np.max(np.abs(coarse["endpoint_error_x_u"]))
    fine_error = np.max(np.abs(fine["endpoint_error_x_u"]))
    assert fine_error < 0.1 * coarse_error, (coarse_error, fine_error)


@pytest.mark.parametrize("method", METHODS)
def test_trajectory_time_interpolation_converges_with_sample_spacing(method):
    """Halving stored step spacing resolves PCHIP time error separately."""
    case = next(case for case in CASES if case["id"] == "kerr_stable_bound")
    coarse = evaluate_case(case, method, max_step=0.2)
    fine = evaluate_case(case, method, max_step=0.1)
    assert fine["trajectory_max_abs_error_x_u"][0] < 0.3 * coarse["trajectory_max_abs_error_x_u"][0]
