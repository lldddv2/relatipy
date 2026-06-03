# Stage 2 — PathsBank Core + Tests

## Goal
Wire the primitives together into the `PathsBank` class.

## Files

- `src/relatipy/numeric/utils/banks/paths_bank/paths_bank.py`
  - `class PathsBank` with:
    - `__init__(path)`
    - `init()` — create dir + empty manifest, or load
    - `get_path(metric, ic, taus, integrator="Radau", adaptative=True, steps_per_period=100, rtol=1e-12, atol=1e-13)` — miss → call `metric.geodesic.get_path(ic, taus, **kw)`, save, return; hit → load NPZ. Attach `coord._bank_hash = h`.
    - `has(metric, ic, taus, **kw)`
    - `hash_of(metric, ic, taus, **kw)`
    - `list_paths()` → list of `{hash, integrator, n_points, tau_span, coordinate_system, created_at}`
    - `get_path_by_hash(h)` — load NPZ
    - `delete_path(h)` — remove NPZ + manifest entry + cascade from any group (Stage 3 wires this)
    - `verify(prune_missing=False)` → `{"missing": [...], "orphan_files": [...]}`
  - Module-level `init_bank(bank)` helper
  - Re-exports: `from .paths_bank import PathsBank, init_bank`

- `src/relatipy/numeric/utils/banks/__init__.py` re-exports `PathsBank`, `init_bank`

## Tests

- `tests/numeric/utils/banks/paths_bank/test_bank_cache.py`
  - `test_miss_invokes_integrator_once` — monkeypatch `Geodesic.get_path` to count calls
  - `test_hit_skips_integrator`
  - `test_hit_returns_equal_arrays` — `np.array_equal(p1.xs, p2.xs)`
  - `test_returned_object_carries_bank_hash`
  - `test_delete_path_removes_file_and_entry`
  - `test_verify_detects_missing_npz`

Use Schwarzschild M=1.0 + tiny `taus` (≤200 points) + `Yoshida6` integrator (no C dep) for speed.

## Acceptance

`pytest tests/numeric/utils/banks/paths_bank/test_bank_cache.py -q` → green.

## Commit

```
[FEAT] paths-bank: add PathsBank class (init/get_path/has/list/delete/verify)
[TEST] paths-bank: cache hit/miss + delete + verify
```
