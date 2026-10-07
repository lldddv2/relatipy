"""Preserve null termination contracts through copying and serialization."""

import copy
import pickle

import numpy as np
import pytest
from astropy import units as u
from astropy.constants import c

import relatipy as rp


def _pickle_roundtrip(value):
    return pickle.loads(pickle.dumps(value))


def _escaping_photon():
    metric = rp.Kerr(mass=1 * u.Msun, spin=0)
    photon = metric.null(
        R=10 * metric.r_g, Theta=np.pi / 2 * u.rad, Phi=0 * u.rad,
        vR=c,
    )
    return metric, photon, (metric.r_g / c).to(u.s)


@pytest.mark.parametrize("restore", [copy.copy, copy.deepcopy, _pickle_roundtrip])
def test_escape_solution_termination_remains_readonly_after_restore(restore):
    metric, photon, scale = _escaping_photon()
    solution = photon.solve(
        t_span=(photon.t, 100 * scale), method="dp45",
        r_escape=50 * metric.r_g,
    )
    restored = restore(solution)
    termination = restored.termination
    assert restored.status == 1
    assert termination.reason == "escape"
    assert termination.t == termination.state.t == solution.termination.t
    with pytest.raises((ValueError, TypeError)):
        termination.t[...] = 99 * u.s
    with pytest.raises((ValueError, TypeError)):
        termination.state.t[...] = 99 * u.s
    assert termination.t == termination.state.t


@pytest.mark.parametrize("restore", [copy.copy, copy.deepcopy, _pickle_roundtrip])
def test_null_terminal_exception_preserves_payload_after_restore(restore):
    metric, photon, scale = _escaping_photon()
    with pytest.raises(rp.IntegrationTerminated) as caught:
        photon.integrate(
            100 * scale, method="dp45", r_escape=50 * metric.r_g,
        )
    original = caught.value
    restored = restore(original)
    assert isinstance(restored, rp.IntegrationTerminated)
    assert type(restored) is type(original)
    assert restored.reason == original.reason == "escape"
    assert restored.t == restored.state.t == original.t
    assert str(restored) == str(original)
    np.testing.assert_array_equal(
        restored.state._canonical, original.state._canonical,
    )
    with pytest.raises((ValueError, TypeError)):
        restored.t[...] = 99 * u.s
    assert restored.t == restored.state.t
