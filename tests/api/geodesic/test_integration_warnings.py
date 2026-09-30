"""IntegrationWarning flags sparse tau_eval sampling and unfit tolerances."""

import warnings

import numpy as np
import pytest
from astropy import units as u

import relatipy as rp
from relatipy import IntegrationWarning, Kerr
from relatipy.geodesic.orbit import _state_scales


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


def _earth_like_orbit():
    metric = Kerr(mass=1 * u.M_sun, spin=0.01)
    return metric.orbit(x=1 * u.au, vy=30 * u.km / u.s)


def test_omitted_tolerances_are_chosen_from_the_state() -> None:
    orbit = _earth_like_orbit()
    with warnings.catch_warnings():
        warnings.simplefilter("error", IntegrationWarning)
        solution = orbit.solve(tau_eval=np.linspace(0, 1, 300) * u.yr)
    info = solution.integration
    assert info.rtol == 1e-10
    assert np.shape(info.atol) == (8,)
    expected = 1e-10 * _state_scales(np.asarray(orbit._initial_y))
    np.testing.assert_allclose(info.atol, expected)


@pytest.mark.parametrize("method", ["radau", "dop853", "dp45", "projection_radau"])
def test_automatic_tolerances_keep_weak_field_orbit_nearly_circular(method) -> None:
    # Regression: rtol=1e-3, atol=1e-6 turned this 1 au, 30 km/s orbit
    # (Newtonian e ~ 0.015) into a radius range of 0.015-2.4 au.
    solution = _earth_like_orbit().solve(
        tau_eval=np.linspace(0, 1, 300) * u.yr, method=method
    )
    radius = solution.r.to_value(u.au)
    assert solution.status == 0
    assert radius.min() > 0.999
    assert radius.max() < 1.031


def test_explicit_rtol_alone_scales_automatic_atol() -> None:
    orbit = _earth_like_orbit()
    solution = orbit.solve(tau_span=(0 * u.s, 1 * u.day), rtol=1e-8)
    expected = 1e-8 * _state_scales(np.asarray(orbit._initial_y))
    np.testing.assert_allclose(solution.integration.atol, expected)


@pytest.mark.parametrize(
    ("options", "match"),
    [
        ({"rtol": 1e-3}, "rtol=0.001 is looser"),
        ({"rtol": 1e-17}, "unattainable in double precision"),
        ({"rtol": 1e-10, "atol": 1e-3}, "scale of Theta, Phi, u\\^t"),
        ({"rtol": 1e-10, "atol": np.full(8, 1e-3)}, "scale of Theta, Phi, u\\^t"),
    ],
)
def test_unfit_explicit_tolerances_warn(options, match) -> None:
    orbit = _earth_like_orbit()
    with pytest.warns(IntegrationWarning, match=match):
        orbit.solve(tau_span=(0 * u.s, 1 * u.day), **options)
    with pytest.warns(IntegrationWarning, match=match):
        orbit.integrate(1 * u.day, **options)


@pytest.mark.parametrize(
    "options",
    [{"rtol": 1e-8, "atol": 1e-10}, {"rtol": 1e-10, "atol": 1e-12}, {"rtol": 0.0}],
)
def test_fit_explicit_tolerances_are_silent(options) -> None:
    orbit = _earth_like_orbit()
    with warnings.catch_warnings():
        warnings.simplefilter("error", IntegrationWarning)
        orbit.solve(tau_span=(0 * u.s, 1 * u.day), **options)


def test_default_failure_no_longer_suggests_tolerances() -> None:
    orbit = _plunging_orbit()
    with warnings.catch_warnings():
        warnings.simplefilter("error", IntegrationWarning)
        solution = orbit.solve(tau_span=(0 * u.s, 1e5 * u.s))
    assert solution.status == -1
