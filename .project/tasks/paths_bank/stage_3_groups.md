# Stage 3 — Groups API + Tests

## Goal
Add named grouping of cached paths.

## Files

- `src/relatipy/numeric/utils/banks/paths_bank/groups.py`
  - `@dataclass Group { name, created_at, updated_at, path_hashes: list[str], notes: str | None }`
  - Pure helpers: `add(group, hashes)`, `remove(group, hashes)`, etc. (mutate-in-place)

- Extend `PathsBank` (`paths_bank.py`) with:
  - `save_group(name, *paths_or_hashes, notes=None) -> Group` — accepts hash strings or objects with `_bank_hash`; validates each hash exists in manifest
  - `show_groups()` — pretty print to stdout
  - `list_groups() -> list[str]`
  - `get_group(name) -> Group`
  - `delete_group(name)`
  - `rename_group(old, new)` — raises if `new` exists
  - `add_to_group(name, *items) -> Group`
  - `remove_from_group(name, *items) -> Group`
- Wire `delete_path` (from Stage 2) to remove `h` from every group

- `banks/__init__.py` re-export `Group`

## Tests

- `tests/numeric/utils/banks/paths_bank/test_groups.py`
  - save with mix of CoordinateBase objects + hash strings
  - `show_groups` produces non-empty stdout
  - `delete_group`, `rename_group` collision raises
  - `add_to_group` / `remove_from_group`
  - `save_group` rejects unknown hash
  - `delete_path(h)` cascades into all groups containing `h`

## Acceptance

`pytest tests/numeric/utils/banks/paths_bank/test_groups.py -q` → green.

## Commit

```
[FEAT] paths-bank: add Groups (save/show/list/delete/rename/add/remove)
[TEST] paths-bank: groups + delete cascade
```
