# Stage 0 — Branch + Skeleton

## Goal
Create the working branch, the empty module skeleton, the task files, and add the "Active Work" block to root `CLAUDE.md`.

## Files

- New: `src/relatipy/numeric/utils/banks/__init__.py`
- New: `src/relatipy/numeric/utils/banks/paths_bank/{__init__,_version}.py`
- New: `tests/numeric/utils/banks/__init__.py`, `tests/numeric/utils/banks/paths_bank/__init__.py`
- New: `.project/tasks/paths_bank/{tasks,stage_0..stage_5}.md`
- Modified: `CLAUDE.md` (Active Work block)

## Acceptance

- `git branch --show-current` → `feature/paths-bank`
- `python -c "import relatipy.numeric.utils.banks"` does not raise (empty module is fine)
- All task files exist

## Commit

```
[FEAT] paths-bank: scaffold module + stage tasks + CLAUDE.md note
```
