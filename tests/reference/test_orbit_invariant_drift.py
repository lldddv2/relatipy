"""Independent invariant and multi-period public Orbit regression checks.

These gates cover the declared fixed cases, not arbitrary spins, durations,
polar or horizon limits.  They use accepted native steps; no interpolated
``tau_eval`` rows enter a conservation measurement.  Conserved quantities
cannot independently bound orbital phase error.
"""

from __future__ import annotations

import numpy as np
import pytest

from .run_orbit_drift import METHODS, drift_cases, measure_case, measurement_pair_errors
from .timelike_invariants import kerr_metric, timelike_invariants


@pytest.mark.parametrize("case", drift_cases(), ids=lambda case: case["id"])
def test_independent_invariants_match_reference_initial_state(case: dict) -> None:
    """Match exact circular constants or independently generated KerrGeoPy constants."""
    measured = timelike_invariants(case["spin"], np.array(case["initial"]))
    for name, expected in case["expected_initial"].items():
        np.testing.assert_allclose(measured[name], expected, rtol=3e-13, atol=3e-13)


def test_independent_metric_reduces_to_schwarzschild() -> None:
    """Check all metric components against the exact diagonal Schwarzschild limit."""
    radius, theta = 10.0, 1.1
    expected = np.diag([-(1-2/radius), 1/(1-2/radius), radius**2,
                        radius**2 * np.sin(theta)**2])
    np.testing.assert_allclose(kerr_metric(0.0, np.array([0, radius, theta, 0])),
                               expected, rtol=3e-15, atol=3e-15)


def test_carter_q_uses_fixed_timelike_rest_mass() -> None:
    """A nonunit velocity must not substitute measured norm for mu squared."""
    case = drift_cases()[0]
    canonical = np.array(case["initial"])
    original = timelike_invariants(case["spin"], canonical)
    factor = 1.01
    canonical[4:] *= factor
    changed = timelike_invariants(case["spin"], canonical)
    np.testing.assert_allclose(changed["norm"], factor**2 * original["norm"],
                               rtol=1e-14, atol=1e-14)
    # Terms quadratic in momenta scale by factor squared; the fixed mass term
    # does not.  Using -norm for mu squared would incorrectly scale all of Q.
    mass_term = case["spin"]**2 * np.cos(canonical[2])**2
    expected = factor**2 * original["CarterQ"] + (1-factor**2) * mass_term
    np.testing.assert_allclose(changed["CarterQ"], expected, rtol=1e-14, atol=1e-14)
    assert not np.isclose(changed["CarterQ"], factor**2 * original["CarterQ"],
                          rtol=1e-8, atol=1e-10)


@pytest.fixture(scope="module")
def measurements() -> dict:
    """Run both tolerance levels once per fixed orbit and native method."""
    return {(case["id"], method, tolerance): measure_case(case, method, tolerance)
            for case in drift_cases() for method in METHODS for tolerance in ("coarse", "fine")}


@pytest.mark.parametrize("case", drift_cases(), ids=lambda case: case["id"])
@pytest.mark.parametrize("method", METHODS)
def test_multiperiod_invariants_remain_bounded_and_converge(measurements, case, method) -> None:
    """Bound tight-tolerance drift and require Kerr convergence above roundoff."""
    coarse = measurements[case["id"], method, "coarse"]
    fine = measurements[case["id"], method, "fine"]
    assert measurement_pair_errors(case, coarse, fine) == []
