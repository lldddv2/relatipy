"""IntegrationWarning flags sparse tau_eval sampling and default-tolerance failures."""

import warnings

import numpy as np
import pytest
from astropy import units as u

import relatipy as rp
from relatipy import IntegrationError, IntegrationWarning, Kerr


def _eccentric_orbit():
    # Newtonian estimate: e ~ 0.92, Kepler period ~ 2.4 days.
    metric = Kerr(mass=1 * u.M_sun, spin=0.9)
    return metric.orbit(x=1e10 * u.m, vy=30 * u.km / u.s)


def _plunging_orbit():
    metric = Kerr(mass=1e6 * u.M_sun, spin=0.0)
    return metric.orbit(
        R=6 * metric.r_g, Theta=np.pi / 2 * u.rad, Phi=0 * u.rad,
        vR=-3e4 * u.km / u.s,
    )


def test_warning_is_public_user_warning() -> None:
    assert issubclass(IntegrationWarning, UserWarning)
    assert rp.geodesic.IntegrationWarning is IntegrationWarning
    assert "IntegrationWarning" in rp.__all__


def test_sparse_tau_eval_warns_about_aliasing() -> None:
    orbit = _eccentric_orbit()
    taus = np.linspace(0, 1, 201) * u.yr
    with pytest.warns(IntegrationWarning, match="may alias the orbit"):
        orbit.solve(tau_eval=taus, method="dop853", rtol=1e-10, atol=1e-12)


def test_dense_tau_eval_with_tight_tolerances_is_silent() -> None:
    orbit = _eccentric_orbit()
    taus = np.linspace(0, 10, 2001) * u.day
    with warnings.catch_warnings():
        warnings.simplefilter("error", IntegrationWarning)
        solution = orbit.solve(
            tau_eval=taus, method="dop853", rtol=1e-10, atol=1e-12
        )
    assert solution.status == 0


def test_single_tau_eval_sample_never_warns_about_aliasing() -> None:
    orbit = _eccentric_orbit()
    with warnings.catch_warnings():
        warnings.simplefilter("error", IntegrationWarning)
        orbit.solve(tau_eval=[0] * u.s, rtol=1e-10, atol=1e-12)


def test_default_tolerance_failure_warns_in_solve() -> None:
    orbit = _plunging_orbit()
    with pytest.warns(IntegrationWarning, match="default tolerances"):
        solution = orbit.solve(tau_span=(0 * u.s, 1e5 * u.s))
    assert solution.status == -1


def test_default_tolerance_failure_warns_before_integrate_raises() -> None:
    orbit = _plunging_orbit()
    with pytest.warns(IntegrationWarning, match="default tolerances"):
        with pytest.raises(IntegrationError):
            orbit.integrate(1e5 * u.s)


@pytest.mark.parametrize(
    "options",
    [{"rtol": 1e-4}, {"atol": 1e-7}, {"atol": np.full(8, 1e-6)}],
)
def test_failure_with_explicit_tolerances_does_not_warn(options) -> None:
    orbit = _plunging_orbit()
    with warnings.catch_warnings():
        warnings.simplefilter("error", IntegrationWarning)
        solution = orbit.solve(tau_span=(0 * u.s, 1e5 * u.s), **options)
    assert solution.status == -1
