# Stage 4 — Docs + Notebook

## Goal
Make the bank discoverable and easy to learn.

## Files

- Module docstrings on every public class/function in `paths_bank/*.py` (NumPy style)
- `notebooks/paths_bank/example_paths_bank.ipynb` — minimal demo: instantiate bank, run two identical `get_path` calls (timing diff), save a group, `show_groups`
- Top-level `src/relatipy/__init__.py` — consider re-exporting `PathsBank, init_bank` (verify it doesn't slow import). If circular import is a risk, skip and document explicit import path.

## Acceptance

- Notebook runs top-to-bottom without errors
- `help(PathsBank)` shows complete docstring with attributes + methods

## Commit

```
[DOC] paths-bank: notebook example + docstrings
```
