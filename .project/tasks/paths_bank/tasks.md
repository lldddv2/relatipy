# PathsBank — Master Checklist

**Branch:** `feature/paths-bank`
**Started:** 2026-06-03
**Owner:** @luisdiaz
**Plan reference:** `/Users/ldiaz/.claude/plans/soft-rolling-trinket.md` (approved)

## Goal

Add a content-addressed cache (`PathsBank`) for `Geodesic.get_path` so that re-running the same configuration loads precomputed trajectories from disk instead of re-integrating. Bank also supports named **groups** for bundling related paths.

User-locked decisions:
- Storage: NPZ per path + central `manifest.json`
- Key: bit-for-bit exact (`repr()` floats → SHA256)
- `taus`: hash full array bytes + dtype + shape
- Default bank dir: `notebooks/paths_bank/` (any path accepted)

## Stages

- [x] Stage 0 — Branch + module skeleton + task files + `CLAUDE.md` Active Work block
- [x] Stage 1 — `key.py`, `manifest.py`, `storage.py`, `reconstruct.py` + unit tests
- [x] Stage 2 — `PathsBank` class core (`init`, `get_path`, `has`, `list_paths`, `get_path_by_hash`, `delete_path`, `verify`) + tests
- [x] Stage 3 — Groups API (`save_group`, `show_groups`, `list_groups`, `get_group`, `delete_group`, `rename_group`, `add_to_group`, `remove_from_group`) + tests
- [x] Stage 4 — Docs: notebook demo at `notebooks/paths_bank_demo.py`, module docstrings, `CLAUDE.md` Active Work block
- [x] Stage 5 — End-to-end integration test (`test_integration_end_to_end.py`) + full `pytest -q` green (all 33 bank tests pass; 1 unrelated pre-existing failure on main)

**Status:** functional. On merge: remove the Active Work block from `CLAUDE.md` and archive this directory.

## Quick resume for a fresh chat

1. `git checkout feature/paths-bank`
2. Read this file → find first unchecked stage.
3. Open the corresponding `stage_N_*.md` and follow it.
4. Commit per stage using `[FEAT]` / `[TEST]` / `[DOC]` per CLAUDE.md.

## Module target

```
src/relatipy/numeric/utils/banks/__init__.py            # re-export PathsBank, init_bank, Group
src/relatipy/numeric/utils/banks/paths_bank/
    __init__.py
    paths_bank.py    # orchestrator
    key.py           # canonical_key + hash_key (pure)
    manifest.py      # Manifest + atomic JSON I/O
    storage.py       # NPZ save/load
    reconstruct.py   # rebuild CoordinateSystem
    groups.py        # Group + group ops
    _version.py      # BANK_FORMAT_VERSION = "1.0"
```

Tests mirror under `tests/numeric/utils/banks/paths_bank/`.
