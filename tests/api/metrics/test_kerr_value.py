"""Kerr is an immutable, hashable, picklable value object."""

import copy
import pickle

import pytest
from astropy import units as u

from relatipy import Kerr


def test_kerr_equality_hash_and_repr_follow_mass_and_spin() -> None:
    metric = Kerr(mass=1 * u.M_sun, spin=0.5)
    same = Kerr(mass=1 * u.M_sun, spin=0.5 * u.one)
    assert metric == same and hash(metric) == hash(same)
    assert metric != Kerr(mass=1 * u.M_sun, spin=0.4)
    assert repr(metric) == "Kerr(mass=<Quantity 1. solMass>, spin=0.5)"


def test_kerr_is_immutable_and_returns_read_only_radii() -> None:
    metric = Kerr(mass=1 * u.M_sun, spin=0.5)
    with pytest.raises(AttributeError, match="immutable"):
        metric._spin = 0.9
    with pytest.raises(AttributeError, match="immutable"):
        del metric._mass
    for radius in (metric.r_g, metric.r_isco_prograde, metric.r_isco_retrograde,
                   metric.r_photon_prograde, metric.r_photon_retrograde,
                   metric.horizons.event):
        assert not radius.flags.writeable


@pytest.mark.parametrize(
    "roundtrip",
    [copy.copy, copy.deepcopy, lambda value: pickle.loads(pickle.dumps(value))],
)
def test_kerr_orbit_and_solution_survive_copy_and_pickle(roundtrip) -> None:
    metric = Kerr(mass=1 * u.M_sun, spin=0.5)
    assert roundtrip(metric) == metric
    orbit = metric.orbit(a=50 * metric.r_g)
    solution = orbit.solve(tau_span=(0 * u.s, 100 * u.s))
    restored = roundtrip(solution)
    assert restored._metric == metric
    assert restored.tau.to_value(u.s).tolist() == solution.tau.to_value(u.s).tolist()
    assert roundtrip(orbit).initial.tau == orbit.initial.tau


def test_spin_rejects_dimensional_quantity_with_clear_message() -> None:
    with pytest.raises(TypeError, match="spin must be dimensionless; got unit m"):
        Kerr(mass=1 * u.M_sun, spin=0.5 * u.m)
    with pytest.raises(ValueError, match="mass must be a scalar; got shape"):
        Kerr(mass=[1] * u.M_sun, spin=0.5)
