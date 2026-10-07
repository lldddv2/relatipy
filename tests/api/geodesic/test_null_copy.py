"""Copy and serialization retain independent, immutable photon data."""

import copy
import pickle

import numpy as np
import pytest
from astropy import units as u

from relatipy import Kerr


def _pickle_roundtrip(value):
    return pickle.loads(pickle.dumps(value))


@pytest.fixture
def evolved_photon():
    metric = Kerr(mass=1 * u.Msun, spin=0.3)
    photon = metric.null(
        R=(10 * metric.r_g).to(u.km), Theta=1.1 * u.rad, Phi=0.2 * u.rad,
        b=3 * metric.r_g, eta=5 * metric.r_g**2,
        radial_sign=1, polar_sign=1, t=3 * u.ms,
    )
    photon.integrate(photon.t + 0.1 * metric._time_scale, method="dp45")
    return photon


@pytest.mark.parametrize("clone", [
    pytest.param(lambda photon: photon, id="original"),
    pytest.param(lambda photon: photon.copy(), id="copy-method"),
    pytest.param(copy.copy, id="shallow-copy"),
    pytest.param(copy.deepcopy, id="deepcopy"),
    pytest.param(_pickle_roundtrip, id="pickle"),
])
@pytest.mark.parametrize("name, replacement", [("b", 999 * u.m), ("eta", 999 * u.m**2)])
def test_null_invariants_reject_direct_quantity_assignment(evolved_photon, clone, name, replacement):
    quantity = getattr(clone(evolved_photon), name)
    with pytest.raises((ValueError, TypeError)):
        quantity[...] = replacement


@pytest.mark.parametrize("clone", [
    pytest.param(lambda photon: photon.copy(), id="copy-method"),
    pytest.param(copy.copy, id="shallow-copy"),
    pytest.param(copy.deepcopy, id="deepcopy"),
    pytest.param(_pickle_roundtrip, id="pickle"),
])
def test_null_copy_preserves_initial_current_units_and_exact_times(evolved_photon, clone):
    photon = clone(evolved_photon)
    assert type(photon) is type(evolved_photon)
    assert photon.t == evolved_photon.t
    assert photon.initial.t == evolved_photon.initial.t
    assert photon.t.unit == photon.initial.t.unit == u.ms
    assert photon.R.unit == photon.initial.R.unit == photon.b.unit == u.km
    assert photon.eta.unit == u.km**2
    np.testing.assert_array_equal(photon._initial_y, evolved_photon._initial_y)
    np.testing.assert_array_equal(photon._current_y, evolved_photon._current_y)
    assert photon.b == evolved_photon.b
    assert photon.eta == evolved_photon.eta

    current_time = evolved_photon.t.copy()
    photon.reset()
    assert photon.t == photon.initial.t
    assert evolved_photon.t == current_time
    photon.integrate(current_time, method="dp45")
    assert photon.t == current_time
    assert evolved_photon.t == current_time


@pytest.mark.parametrize("restore", [copy.deepcopy, _pickle_roundtrip])
def test_null_copy_of_restored_photon_has_independent_storage(evolved_photon, restore):
    restored = restore(evolved_photon)
    # Force lazy views before copying to detect aliasing of mutable caches.
    restored.xyz
    restored.initial.xyz
    photon = restored.copy()
    assert photon._current._views is not restored._current._views
    assert photon.initial._views is not restored.initial._views
    for name in ("_b", "_eta", "_initial_y", "_current_y"):
        assert not np.shares_memory(getattr(photon, name), getattr(restored, name))
    for name, replacement in (("b", 999 * u.m), ("eta", 999 * u.m**2)):
        with pytest.raises((ValueError, TypeError)):
            getattr(photon, name)[...] = replacement
        assert getattr(restored, name) == getattr(evolved_photon, name)
