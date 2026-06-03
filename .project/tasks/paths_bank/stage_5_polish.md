# Stage 5 — End-to-End + Polish

## Goal
Confirm bank works under realistic conditions and is regression-free.

## Files

- `tests/numeric/utils/banks/paths_bank/test_integration_end_to_end.py`
  - Schwarzschild M=1.0 + `OrbitalElements(a=50, e=0.2, mass=1.0)` + `Yoshida6` + `taus = np.linspace(0, 1000, 200)`
  - Run twice; second run must NOT invoke integrator (monkeypatch counter)
  - `np.array_equal(p1.xs, p2.xs)` and `np.array_equal(p1.vs, p2.vs)`
  - `save_group("smoke", p1)` + `show_groups()`

- README / docs: short mention of `PathsBank` in user-facing summary

## Acceptance

- `pytest -q` (full suite) → green
- Manual smoke (from plan §Verification) works

## On merge

- Remove the "Active Work" block from root `CLAUDE.md`
- Delete or move `.project/tasks/paths_bank/` to a `done/` subfolder per project preference

## Commit

```
[TEST] paths-bank: end-to-end integration
[DOC] paths-bank: README mention
```
