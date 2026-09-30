"""Coordinates are reconstructed once, only when a caller requests a view."""

import numpy as np
import pytest
from astropy import units as u
from astropy.constants import c

from relatipy import _core
from relatipy.geodesic._deferred import DeferredState
from relatipy.metrics.kerr import Kerr


def _metric() -> Kerr:
    return Kerr(mass=1 * u.Msun, spin=0.5)


def _orbit(metric: Kerr, family: str):
    radius = 12 * metric.r_g
    angular_velocity = (0.1 * c / radius) * u.rad
    if family == "cartesian":
        return metric.orbit(x=radius, vy=0.1 * c)
    if family == "spherical":
        return metric.orbit(
            r=radius, theta=90 * u.deg, phi=0 * u.rad,
            vphi=angular_velocity,
        )
    if family == "elements":
        return metric.orbit(a=100 * metric.r_g, e=0.2, inc=30 * u.deg)
    return metric.orbit(
        R=radius, Theta=90 * u.deg, Phi=0 * u.rad,
        vPhi=angular_velocity,
    )


@pytest.mark.parametrize("family", ["bl", "cartesian", "spherical", "elements"])
def test_solve_stores_bl_and_input_family(family: str) -> None:
    metric = _metric()
    orbit = _orbit(metric, family)
    solution = orbit.solve(
        tau_span=(0 * u.s, 0.01 * metric._time_scale), method="dp45"
    )
    assert isinstance(solution._state, DeferredState)
    assert set(solution._state._views) == {"bl", family}
    assert solution.tRQP.R.shape == solution.tau.shape
    assert set(solution._state._views) == {"bl", family}


def test_spherical_view_is_computed_once_and_shared_by_accessors(monkeypatch) -> None:
    metric = _metric()
    solution = _orbit(metric, "elements").solve(
        tau_span=(0 * u.s, 0.01 * metric._time_scale), method="dp45"
    )
    calls = []
    original = _core.reconstruct_canonical_family_batch

    def counted(spin, canonical, family):
        calls.append(family)
        return original(spin, canonical, family)

    monkeypatch.setattr(_core, "reconstruct_canonical_family_batch", counted)
    assert set(solution._state._views) == {"bl", "elements"}
    spherical = solution.trqp
    assert solution.trqp is spherical
    assert solution.vr.shape == solution.tau.shape
    assert calls == ["spherical"]
    assert solution.orbital_elements() is solution.orbital_elements()
    assert calls == ["spherical"]
    assert solution.txyz is solution.txyz
    assert calls == ["spherical", "cartesian"]


def test_requested_samples_and_interpolation_keep_other_views_deferred() -> None:
    metric = _metric()
    times = np.array([0.0, 0.005, 0.01]) * metric._time_scale
    solution = _orbit(metric, "elements").solve(tau_eval=times, method="dp45")
    assert set(solution._state._views) == {"bl", "elements"}
    exact = solution.at(tau=times[1])
    assert set(exact._views) == {"bl", "elements"}
    assert set(solution._state._views) == {"bl", "elements"}
    selected = solution[1]
    assert set(selected._views) == {"bl", "elements"}
    between = solution.at(tau=0.007 * metric._time_scale)
    assert set(between._views) == {"bl", "elements"}
    assert set(solution._state._views) == {"bl", "elements"}
    assert np.isfinite(between.trqp.r.to_value(u.m))
    assert "spherical" in between._views
    assert "spherical" not in solution._state._views


def test_lazy_views_match_full_native_reconstruction() -> None:
    from relatipy.geodesic.orbit import _state_from_native

    metric = _metric()
    solution = _orbit(metric, "bl").solve(
        tau_span=(0 * u.s, 0.01 * metric._time_scale), method="dp45"
    )
    lazy = solution._state
    eager = _state_from_native(
        metric, lazy._canonical,
        (lazy.tau / metric._time_scale).to_value(u.one),
        length_unit=lazy._length_unit, time_unit=lazy._time_unit,
    )
    for name in ("R", "Theta", "Phi", "vR", "vTheta", "vPhi", "ut",
                 "r", "theta", "phi", "vr", "vtheta", "vphi",
                 "x", "y", "z", "vx", "vy", "vz", "ux", "uy", "uz"):
        observed = getattr(lazy, name)
        expected = getattr(eager, name)
        np.testing.assert_allclose(observed.to_value(expected.unit), expected.value,
                                   rtol=2e-14, atol=0)
    for name in ("a", "e", "inc", "Omega", "omega", "f"):
        observed = getattr(lazy.orbital_elements(), name)
        expected = getattr(eager.orbital_elements(), name)
        if isinstance(observed, u.Quantity):
            np.testing.assert_allclose(observed.to_value(expected.unit),
                                       expected.value, rtol=2e-14, atol=0,
                                       equal_nan=True)
        else:
            np.testing.assert_allclose(observed, expected, rtol=2e-14, atol=0,
                                       equal_nan=True)
