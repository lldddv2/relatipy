"""Tests for low-level NPZ storage."""

import numpy as np

from relatipy.numeric.coordinates import Cartesian
from relatipy.numeric.utils.banks.paths_bank import storage


def _make_trajectory(N=7):
    xs = np.zeros((4, N))
    xs[1] = np.linspace(1.0, 2.0, N)  # x
    vs = np.zeros((3, N))
    vs[0] = 0.1
    # Use Cartesian skeleton instance with batched arrays bolted on; storage
    # reads only .xs / .vs / .kwargs / name_metric / type, so we bypass the
    # CoordinateBase length-4 init by constructing a minimal stand-in.
    coord = Cartesian.__new__(Cartesian)
    coord.xs = xs
    coord.vs = vs
    coord.kwargs = {}
    coord.name_metric = "Cartesian"
    return coord


def test_roundtrip(tmp_path):
    N = 7
    coord = _make_trajectory(N)
    taus = np.linspace(0.0, 1.0, N)

    npz = tmp_path / "p.npz"
    meta = storage.save(npz, coord, taus, integrator="Yoshida6")
    raw = storage.load_raw(npz)

    assert np.array_equal(raw["xs"], coord.xs)
    assert np.array_equal(raw["vs"], coord.vs)
    assert np.array_equal(raw["taus"], taus)
    assert raw["coord_class"] == "Cartesian"
    assert raw["meta"]["integrator"] == "Yoshida6"
    assert meta["n_points"] == N
