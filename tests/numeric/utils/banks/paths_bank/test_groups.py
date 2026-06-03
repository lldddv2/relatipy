"""Tests for the named-group API of PathsBank."""

import io
import numpy as np
import pytest

from relatipy.numeric.coordinates import OrbitalElements


def _kw():
    return dict(integrator="Yoshida6", adaptative=False, rtol=1e-9, atol=1e-12)


def _two_paths(bank, schw, taus):
    ic1 = OrbitalElements(a=50.0, e=0.2, mass=schw.mass)
    ic2 = OrbitalElements(a=60.0, e=0.1, mass=schw.mass)
    p1 = bank.get_path(schw, ic1, taus, **_kw())
    p2 = bank.get_path(schw, ic2, taus, **_kw())
    return p1, p2


def test_save_group_mixed_inputs(bank, schw, taus):
    p1, p2 = _two_paths(bank, schw, taus)
    g = bank.save_group("mix", p1, p2._bank_hash, notes="hello")
    assert g.path_hashes == [p1._bank_hash, p2._bank_hash]
    assert g.notes == "hello"
    assert "mix" in bank.list_groups()


def test_save_group_rejects_unknown_hash(bank, schw, ic_orbital, taus):
    bank.get_path(schw, ic_orbital, taus, **_kw())
    with pytest.raises(KeyError):
        bank.save_group("bad", "0" * 64)


def test_save_group_duplicate_name_raises(bank, schw, ic_orbital, taus):
    p = bank.get_path(schw, ic_orbital, taus, **_kw())
    bank.save_group("g", p)
    with pytest.raises(ValueError):
        bank.save_group("g", p)


def test_show_groups_prints(bank, schw, ic_orbital, taus, capsys):
    p = bank.get_path(schw, ic_orbital, taus, **_kw())
    bank.save_group("g", p)
    bank.show_groups()
    out = capsys.readouterr().out
    assert "g" in out
    assert p._bank_hash in out


def test_delete_group(bank, schw, ic_orbital, taus):
    p = bank.get_path(schw, ic_orbital, taus, **_kw())
    bank.save_group("g", p)
    bank.delete_group("g")
    assert "g" not in bank.list_groups()


def test_rename_group_collision_raises(bank, schw, taus):
    p1, p2 = _two_paths(bank, schw, taus)
    bank.save_group("a", p1)
    bank.save_group("b", p2)
    with pytest.raises(ValueError):
        bank.rename_group("a", "b")


def test_add_remove_from_group(bank, schw, taus):
    p1, p2 = _two_paths(bank, schw, taus)
    bank.save_group("g", p1)
    g = bank.add_to_group("g", p2)
    assert p2._bank_hash in g.path_hashes
    g2 = bank.remove_from_group("g", p1)
    assert p1._bank_hash not in g2.path_hashes


def test_delete_path_cascades_into_groups(bank, schw, taus):
    p1, p2 = _two_paths(bank, schw, taus)
    bank.save_group("g", p1, p2)
    bank.delete_path(p1._bank_hash)
    g = bank.get_group("g")
    assert p1._bank_hash not in g.path_hashes
    assert p2._bank_hash in g.path_hashes
