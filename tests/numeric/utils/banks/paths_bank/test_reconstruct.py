"""Tests for CoordinateBase reconstruction from NPZ."""

import json

import numpy as np
import pytest

from relatipy.numeric.coordinates import Cartesian
from relatipy.numeric.utils.banks.paths_bank import storage
from relatipy.numeric.utils.banks.paths_bank.reconstruct import (
    BankIntegrityError,
    rebuild,
)


def _make_cart_trajectory(N=5, a_kwarg=None):
    coord = Cartesian.__new__(Cartesian)
    xs = np.zeros((4, N)); xs[1] = np.linspace(1.0, 2.0, N)
    vs = np.zeros((3, N)); vs[0] = 0.1
    coord.xs = xs; coord.vs = vs
    coord.kwargs = {} if a_kwarg is None else {"a": float(a_kwarg)}
    coord.name_metric = "Cartesian"
    return coord


def test_rebuild_cartesian(tmp_path):
    N = 5
    coord = _make_cart_trajectory(N)
    taus = np.linspace(0.0, 1.0, N)

    npz = tmp_path / "p.npz"
    storage.save(npz, coord, taus, integrator="Yoshida6")
    rebuilt = rebuild(npz)

    assert type(rebuilt).__name__ == "Cartesian"
    assert np.array_equal(rebuilt.xs, coord.xs)
    assert np.array_equal(rebuilt.vs, coord.vs)
    assert np.array_equal(rebuilt.taus, taus)


def test_unknown_class_raises(tmp_path):
    coord = _make_cart_trajectory(3)
    taus = np.linspace(0.0, 1.0, 3)
    npz = tmp_path / "p.npz"
    storage.save(npz, coord, taus, integrator="Yoshida6")

    raw = storage.load_raw(npz)
    raw["coord_class"] = "TotallyBogusCoord"
    np.savez_compressed(
        npz,
        xs=raw["xs"],
        vs=raw["vs"],
        taus=raw["taus"],
        coord_class=np.asarray("TotallyBogusCoord"),
        coord_kwargs_json=np.asarray(json.dumps({})),
        meta=np.asarray(json.dumps(raw["meta"])),
    )
    with pytest.raises(BankIntegrityError):
        rebuild(npz)


def test_missing_npz_raises(tmp_path):
    with pytest.raises(BankIntegrityError):
        rebuild(tmp_path / "does_not_exist.npz")
