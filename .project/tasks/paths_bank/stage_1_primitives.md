# Stage 1 — Primitives + Tests

## Goal
Build the pure / I/O primitives that `PathsBank` depends on.

## Files

- `src/relatipy/numeric/utils/banks/paths_bank/key.py`
  - `_to_float(x)` — accepts plain float / `astropy.Quantity` (try `.to_value()`)
  - `_canonical_array(arr)` → `{"dtype", "shape", "sha256"}`
  - `_canonical_metric(metric)` → `{"type", "mass", ["a"]}` using `repr()`
  - `_canonical_ic(ic)` → `{"type", "system_name", "xs", "vs", "kwargs"}`
  - `_canonical_integrator(kw)` — merges defaults `("Radau", True, 100, 1e-12, 1e-13)`
  - `canonical_key(metric, ic, taus, integrator_kwargs)` → dict
  - `hash_key(key)` → SHA256 hex; excludes `relatipy_version` from hash

- `src/relatipy/numeric/utils/banks/paths_bank/manifest.py`
  - `Manifest` dataclass with `bank_format_version`, `relatipy_version`, `created_at`, `updated_at`, `entries: dict`, `groups: dict`
  - `Manifest.load(path)` / `Manifest.save_atomic(path)` (write `.tmp` → `os.replace`)
  - `Manifest.empty()` factory

- `src/relatipy/numeric/utils/banks/paths_bank/storage.py`
  - `save(npz_path, coord, taus, meta_dict)` — writes `xs`, `vs`, `taus`, `coord_class`, `coord_kwargs_json`, `meta`
  - Uses `np.savez_compressed`

- `src/relatipy/numeric/utils/banks/paths_bank/reconstruct.py`
  - `rebuild(npz_path) -> CoordinateBase` — lazy-imports `coordinate_systems` registry; raises `BankIntegrityError` for unknown class
  - Defines exception classes (`BankError`, `BankIntegrityError`)

## Tests

Under `tests/numeric/utils/banks/paths_bank/`:

- `test_key_canonicalization.py`
- `test_manifest.py`
- `test_storage.py`
- `test_reconstruct.py`

Cover: same-input collision, dtype-sensitivity, default-collapse, kwargs-order irrelevance, atomic-write recovery, NPZ round-trip, unknown-class raises.

## Acceptance

`pytest tests/numeric/utils/banks/paths_bank/ -q` → all pass.

## Commit

```
[FEAT] paths-bank: add key, manifest, storage, reconstruct primitives
[TEST] paths-bank: unit tests for primitives
```
