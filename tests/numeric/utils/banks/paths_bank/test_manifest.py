"""Tests for the Manifest dataclass and its atomic JSON I/O."""

import json

from relatipy.numeric.utils.banks.paths_bank.manifest import Manifest


def test_empty_round_trip(tmp_path):
    p = tmp_path / "manifest.json"
    Manifest.empty().save_atomic(p)
    m = Manifest.load(p)
    assert m.entries == {}
    assert m.groups == {}


def test_entries_persist(tmp_path):
    p = tmp_path / "manifest.json"
    m = Manifest.empty()
    m.add_entry("abc", {"npz_filename": "paths/abc.npz"})
    m.save_atomic(p)
    m2 = Manifest.load(p)
    assert "abc" in m2.entries
    assert m2.entries["abc"]["npz_filename"] == "paths/abc.npz"


def test_atomic_tmp_file_cleaned(tmp_path):
    p = tmp_path / "manifest.json"
    Manifest.empty().save_atomic(p)
    leftover = [f for f in tmp_path.iterdir() if f.name.endswith(".tmp")]
    assert leftover == []


def test_save_creates_parent_dirs(tmp_path):
    p = tmp_path / "nested" / "deeper" / "manifest.json"
    Manifest.empty().save_atomic(p)
    assert p.exists()
    assert json.loads(p.read_text())["bank_format_version"]
