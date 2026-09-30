"""Private Kerr reference-surface profile evaluated by the native core."""

import numpy as np
import pytest
from astropy import units as u

from relatipy import _core
from relatipy.metrics.kerr import Kerr


def test_schwarzschild_horizon_is_sphere_of_twice_r_g() -> None:
    bh = Kerr(mass=1 * u.Msun, spin=0.0)
    theta = np.linspace(0.0, np.pi, 17) * u.rad
    rho, z = bh._surface_profile("outer_horizon", theta)
    assert rho.shape == z.shape == (17,)
    assert rho.unit == u.m and z.unit == u.m
    radius = np.hypot(rho, z)
    np.testing.assert_allclose(
        (radius / bh.r_g).to_value(u.one), 2.0, rtol=1e-13
    )


@pytest.mark.parametrize("spin", [0.0, 0.5, 0.9, 1.0])
def test_ergosurface_meets_horizon_at_poles(spin: float) -> None:
    bh = Kerr(mass=4.3e6 * u.Msun, spin=spin)
    poles = [0.0, np.pi] * u.rad
    rho_h, z_h = bh._surface_profile("outer_horizon", poles)
    rho_e, z_e = bh._surface_profile("ergosurface", poles)
    scale = bh.r_g.to_value(u.m)
    np.testing.assert_allclose(rho_e.to_value(u.m), rho_h.to_value(u.m), atol=1e-13 * scale)
    np.testing.assert_allclose(z_e.to_value(u.m), z_h.to_value(u.m), rtol=1e-13)
    np.testing.assert_allclose(
        (z_h / bh.horizons.event).to_value(u.one), [1.0, -1.0], rtol=1e-13
    )


@pytest.mark.parametrize("spin", [0.0, 0.3, 0.9, 1.0])
def test_ergosurface_equator_is_hypot_two_spin(spin: float) -> None:
    bh = Kerr(mass=1 * u.Msun, spin=spin)
    rho, z = bh._surface_profile("ergosurface", np.pi / 2 * u.rad)
    assert rho.shape == z.shape == ()
    np.testing.assert_allclose(
        (rho / bh.r_g).to_value(u.one), np.hypot(2.0, spin), rtol=1e-13
    )
    assert abs((z / bh.r_g).to_value(u.one)) < 1e-14


def test_scalar_degrees_accepted_and_read_only() -> None:
    bh = Kerr(mass=1 * u.Msun, spin=0.7)
    rho, z = bh._surface_profile("ergosurface", [0.0, 45.0, 90.0] * u.deg)
    assert rho.shape == (3,)
    for value in (rho, z):
        with pytest.raises(ValueError):
            value[0] = 0 * u.m


@pytest.mark.parametrize(
    "theta",
    [[-0.1] * u.rad, [np.pi + 1e-6] * u.rad, [np.nan] * u.rad, [np.inf] * u.rad],
)
def test_invalid_theta_raises_value_error(theta: u.Quantity) -> None:
    bh = Kerr(mass=1 * u.Msun, spin=0.5)
    with pytest.raises(ValueError):
        bh._surface_profile("ergosurface", theta)


def test_non_angle_theta_rejected() -> None:
    bh = Kerr(mass=1 * u.Msun, spin=0.5)
    with pytest.raises(u.UnitConversionError):
        bh._surface_profile("outer_horizon", [1.0] * u.m)
    with pytest.raises(TypeError):
        bh._surface_profile("outer_horizon", [1.0])


@pytest.mark.parametrize("surface", ["inner_horizon", "", "ERGOSURFACE", 1, None])
def test_unknown_surface_raises_value_error(surface) -> None:
    bh = Kerr(mass=1 * u.Msun, spin=0.5)
    with pytest.raises(ValueError):
        bh._surface_profile(surface, [0.5] * u.rad)


def test_binding_validates_and_accepts_empty() -> None:
    rho, z = _core.kerr_surface_profile(0.5, 1, np.empty(0))
    assert rho.shape == z.shape == (0,)
    with pytest.raises(ValueError):
        _core.kerr_surface_profile(0.5, 2, np.array([0.5]))
    with pytest.raises(ValueError):
        _core.kerr_surface_profile(1.5, 0, np.array([0.5]))
    with pytest.raises(ValueError):
        _core.kerr_surface_profile(0.5, 0, np.array([4.0]))
