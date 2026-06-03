"""End-to-end tests for PathsBank cache hit/miss behavior."""

import numpy as np
import pytest

from relatipy.numeric.geodesic.geodesic import Geodesic
from relatipy.numeric.utils.banks.paths_bank.reconstruct import BankIntegrityError


def _make_kw():
    return dict(
        integrator="Yoshida6",
        adaptative=False,
        steps_per_period=100,
        rtol=1e-9,
        atol=1e-12,
    )


def test_miss_then_hit_skips_integrator(bank, schw, ic_orbital, taus, monkeypatch):
    calls = {"n": 0}
    real = Geodesic.get_path

    def counted(self, *args, **kwargs):
        calls["n"] += 1
        return real(self, *args, **kwargs)

    monkeypatch.setattr(Geodesic, "get_path", counted)

    p1 = bank.get_path(schw, ic_orbital, taus, **_make_kw())
    assert calls["n"] == 1
    p2 = bank.get_path(schw, ic_orbital, taus, **_make_kw())
    assert calls["n"] == 1  # second call was a cache hit
    assert np.array_equal(p1.xs, p2.xs)
    assert np.array_equal(p1.vs, p2.vs)


def test_returned_object_carries_bank_hash(bank, schw, ic_orbital, taus):
    p = bank.get_path(schw, ic_orbital, taus, **_make_kw())
    assert isinstance(p._bank_hash, str)
    assert len(p._bank_hash) == 64


def test_has_and_hash_of(bank, schw, ic_orbital, taus):
    kw = _make_kw()
    assert bank.has(schw, ic_orbital, taus, **kw) is False
    bank.get_path(schw, ic_orbital, taus, **kw)
    assert bank.has(schw, ic_orbital, taus, **kw) is True
    h = bank.hash_of(schw, ic_orbital, taus, **kw)
    assert h in {e["hash"] for e in bank.list_paths()}


def test_get_path_by_hash(bank, schw, ic_orbital, taus):
    p = bank.get_path(schw, ic_orbital, taus, **_make_kw())
    p2 = bank.get_path_by_hash(p._bank_hash)
    assert np.array_equal(p.xs, p2.xs)


def test_delete_path_removes_file_and_entry(bank, schw, ic_orbital, taus):
    p = bank.get_path(schw, ic_orbital, taus, **_make_kw())
    h = p._bank_hash
    npz = bank.paths_dir / f"{h}.npz"
    assert npz.exists()
    bank.delete_path(h)
    assert not npz.exists()
    assert h not in bank.manifest.entries


def test_verify_detects_missing_npz(bank, schw, ic_orbital, taus):
    p = bank.get_path(schw, ic_orbital, taus, **_make_kw())
    (bank.paths_dir / f"{p._bank_hash}.npz").unlink()
    report = bank.verify()
    assert p._bank_hash in report["missing"]


def test_verify_prunes_missing(bank, schw, ic_orbital, taus):
    p = bank.get_path(schw, ic_orbital, taus, **_make_kw())
    (bank.paths_dir / f"{p._bank_hash}.npz").unlink()
    bank.verify(prune_missing=True)
    assert p._bank_hash not in bank.manifest.entries


def test_hit_with_missing_npz_raises(bank, schw, ic_orbital, taus):
    p = bank.get_path(schw, ic_orbital, taus, **_make_kw())
    (bank.paths_dir / f"{p._bank_hash}.npz").unlink()
    with pytest.raises(BankIntegrityError):
        bank.get_path(schw, ic_orbital, taus, **_make_kw())
