"""Orbit.get_keplerian_period: Newtonian period of the current osculating conic."""

from __future__ import annotations

import numpy as np
import pytest
from astropy import units as u
from astropy.constants import G, c

from relatipy import Kerr


def _expected(orbit, metric):
    a = orbit.orbital_elements().a
    return (2 * np.pi * np.sqrt(a**3 / (G * metric.mass))).to(orbit.tau.unit)


@pytest.fixture
def metric():
    return Kerr(mass=4e6 * u.Msun, spin=0.7)


@pytest.mark.parametrize("e", (0.0, 0.6))
def test_period_follows_kepler_third_law(metric, e):
    orbit = metric.orbit(a=200 * metric.r_g, e=e, inc=0.3 * u.rad, Omega=0 * u.rad,
                         omega=0 * u.rad, f=0.4 * u.rad)
    period = orbit.get_keplerian_period()
    assert period.isscalar and period.unit == orbit.tau.unit
    reference = 2 * np.pi * np.sqrt((200 * metric.r_g) ** 3 / (G * metric.mass))
    assert u.isclose(period, reference, rtol=1e-10)


def test_period_uses_the_current_state(metric):
    orbit = metric.orbit(a=20 * metric.r_g, e=0.5, inc=0.3 * u.rad, Omega=0 * u.rad,
                         omega=0 * u.rad, f=0 * u.rad)
    initial = orbit.get_keplerian_period()
    orbit.integrate(initial / 3, method="dp45")
    current = orbit.get_keplerian_period()
    assert u.isclose(current, _expected(orbit, metric), rtol=1e-12)
    assert current != initial  # relativistic precession changes the osculating a


def test_unbound_orbit_has_no_keplerian_period(metric):
    orbit = metric.orbit(x=100 * metric.r_g, vy=0.5 * c)
    with pytest.raises(ValueError, match="elliptic"):
        orbit.get_keplerian_period()
