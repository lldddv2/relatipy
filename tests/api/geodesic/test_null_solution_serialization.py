"""Regressions for continuous null azimuth and immutable solution copies."""

from __future__ import annotations

import copy
import pickle

import numpy as np
import pytest
from astropy import units as u
from astropy.constants import c

from relatipy import Kerr


def _pickle_roundtrip(value):
    return pickle.loads(pickle.dumps(value))


def test_null_solution_at_keeps_sparse_continuous_azimuth() -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0)
    time_scale = metric.r_g / c
    photon = metric.null(
        R=3 * metric.r_g, Theta=np.pi / 2 * u.rad, Phi=0 * u.rad,
        vPhi=1 * u.rad / time_scale,
    )
    solution = photon.solve(t_eval=np.array([0, 4, 8]) * np.sqrt(27) * time_scale)
    np.testing.assert_allclose(solution.Phi.to_value(u.rad), [0, 4, 8], atol=1e-6)

    # Close enough to the final sample that Cartesian interpolation agrees
    # with its position, while still exercising the interpolated branch.
    query = solution.t[-1] - 1e-8 * time_scale
    assert solution.at(query).Phi.to_value(u.rad) == pytest.approx(8, abs=1e-6)


@pytest.mark.parametrize("clone", [copy.copy, copy.deepcopy, _pickle_roundtrip])
def test_null_solution_copy_keeps_arrays_read_only(clone) -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0)
    time_scale = metric.r_g / c
    photon = metric.null(
        R=10 * metric.r_g, Theta=np.pi / 2 * u.rad, Phi=0 * u.rad, vR=c,
    )
    solution = photon.solve(t_eval=np.array([0, 0.1, 0.2]) * time_scale)
    # Prime lazy views too; restoration must not retain writable cached data.
    _ = solution.xyz, solution.vRQP
    restored = clone(solution)
    assert type(restored) is type(solution)
    assert restored.status == solution.status
    assert restored.message == solution.message
    assert restored.integration.method == solution.integration.method
    assert restored.integration.n_steps == solution.integration.n_steps
    assert restored.integration.nfev == solution.integration.nfev

    for name in ("t", "xyz", "vxyz", "R", "Phi", "vR"):
        quantity = getattr(restored, name)
        np.testing.assert_array_equal(quantity, getattr(solution, name))
        with pytest.raises((ValueError, TypeError)):
            quantity[...] = 999 * quantity.unit
    np.testing.assert_array_equal(restored.integration.atol, solution.integration.atol)
    with pytest.raises((ValueError, TypeError)):
        restored.integration.atol[...] = 999
