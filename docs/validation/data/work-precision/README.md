# Work-precision: pericentre precession, Relatipy vs PyGRO

Generated 2026-10-28; rerun 2026-10-29 with all four Relatipy methods and two PyGRO integrators. Nothing in the relatipy repo was modified.

## Files
- `work_precision.py`   benchmark (helpers copied verbatim from
  `docs/validation/01-Theoretical-Validations.ipynb` cells 4, 6, 43, 45, 57)
- `work-precision.csv`  2800 rows = 100 orbits x 28 settings
- `run.log`       per-setting elapsed times, PyGRO engine compile time
- `validate_nb.py` / `validate_nb.log`  runs notebook cells 2, 4, 6, 43, 45 *verbatim*
  (`measure_precession`, rtol=1e-12 / goal=12) on orbits 0-4 for cross-check
- `summarize.py` / `summary.txt`  per-setting medians, status counts, cross-check

## Commands
```
PY=/home/lldddv/Documents/UdeA/thesis/relatipy/.venv/bin/python   # notebook kernel: Python 3.11.15
taskset -c 3 $PY work_precision.py --n 100 --pygro rkf78 dp45 --out work-precision.csv > run.log 2>&1
taskset -c 5 $PY validate_nb.py 5 > validate_nb.log
$PY summarize.py work-precision.csv > summary.txt
```

## Setup
- Machine: AMD Ryzen AI 7 350 w/ Radeon 860M (16 logical CPUs), Linux 7.2.6-arch2-1.
  Process pinned to one core (`taskset -c 3`); OMP/OPENBLAS/MKL/NUMEXPR/VECLIB_NUM_THREADS=1.
- Versions: relatipy 1.2.0 editable from `src/` (HEAD ebfe808 + uncommitted working tree, forces only),
  pygro 1.0.3, numpy 2.4.6, astropy 8.0.1.
- Sample: notebook cell 57 (seed 2026; draw 100 a log-uniform [30,1e4] M, then 100 e uniform
  [0.05,0.9]); orbit_id = index in that draw. Kerr chi=0.9, M=4.3e6 Msun, equatorial prograde
  (`bh.orbit(p=a(1-e^2), e=e, x=1.0)`), integrated over 100 Keplerian periods of proper time.
- eps = |dphi_int / dphi_orb - 1| (section 3.1): dphi_int from `hermite_pericentres` +
  `mean_advance`, dphi_orb from `exact_precession` (Bini et al. 2016 eqs. 25, 27).
- Relatipy: `orbit.solve(tau_span=..., method=m, rtol=tol)`, atol not passed (as in notebook):
  default atol = rtol x characteristic scale of each native state component at t=0.
  m in {projection_radau, radau, dop853, dp45}, tol in {1e-6, 1e-8, 1e-10, 1e-12, 1e-13}.
  n_steps = `sol.integration.n_steps` (accepted steps).
- PyGRO: rkf78 and dp45, `accuracy_goal = precision_goal = goal` (atol = rtol = 10^-goal), goal in {6,8,10,12},
  first step 1.0 M, horizon stop `Delta > 1e-4` (notebook cell 6). n_steps = len(geo.tau) - 1.
  Engine compile (autowrap) ~12 s, once per integrator, excluded.
- PyGRO integrator choice: pilot on the first 20 orbits with all five PyGRO 1.0.3 integrators
  (one process per integrator, separate cores; ranking only, not used in the figure). Medians at
  goal=12: dp45 eps 5.0e-6 / 1.08 s, rkf45 3.3e-5 / 1.00 s, dp853 3.5e-5 / 1.20 s,
  rkf78 2.2e-4 / 0.79 s; at goal=10 rkf78 is the cheapest for its eps. Kept: dp45 (lowest eps)
  and rkf78 (lowest cost at intermediate eps); dp853 and rkf45 cost more for the same eps.
  ck45 fails in PyGRO 1.0.3: `AttributeError: 'CashKarp45' object has no attribute 'verbose'`.
- Timing: `time.perf_counter` around the solve call only (Relatipy `orbit.solve`; PyGRO
  `engine.integrate`); orbit/Geodesic construction and post-processing excluded.
  One untimed warm-up solve per setting on orbit 0; then 3 repeats per orbit x setting, min kept.
- Timeout 120 s per orbit (SIGALRM for PyGRO; post-hoc check for Relatipy): never triggered.
- status: `ok`; `solver_status_-1` (Relatipy partial solution, numerical failure: "invariant
  projection failed; native status 10"); `exit_horizon` (PyGRO stopping criterion fired);
  `...;error:IndexError` = fewer than 2 pericentres, eps = NaN.
  For non-ok rows eps is still computed on the partial trajectory when possible; do not
  treat it as comparable.

## GYOTO (added 2026-10-29)
- `work_precision_gyoto.py`  GYOTO benchmark; runs with the Python where GYOTO is built
- `work-precision-gyoto.csv` 1000 rows = 100 orbits x 10 settings, same columns (library = gyoto)
- `run-gyoto.log`            per-setting elapsed times
- GYOTO: compiled from https://github.com/gyoto/Gyoto in a micromamba (conda-forge) environment,
  Python 3.12, numpy 2.5.4 (separate from the relatipy venv: Python 3.11).
- Commands:
```
$PY work_precision.py --n 100 --pygro --dump-states states.json   # initial states, tau_end, t_end, exact
taskset -c 3 $GYOTO_PY work_precision_gyoto.py states.json --n 100 \
    --integrators runge_kutta_fehlberg78 runge_kutta_dopri5 --out work-precision-gyoto.csv > run-gyoto.log
```
- Same sample, initial state (Boyer-Lindquist, M = 1), eps recipe (`hermite_pericentres` on stored
  nodes with tau and dx/dtau), timing (xFill only, best of 3, warm-up, one core, threads = 1).
- `Star` in `KerrBL`, RelTol = AbsTol = tol in {1e-6, 1e-8, 1e-10, 1e-12, 1e-13}, MaxIter = 1e8.
- `DeltaMaxOverR = 1e10`: the default 1 caps the step at r; with it the number of steps did not
  change between tol = 1e-6 and 1e-10 (smoke test, 3 orbits), so the cap, not the tolerance, set
  the cost. Lifted so the tolerance sets the step, as for Relatipy and PyGRO (no max step).
- GYOTO stops at a coordinate time: `t_end` = coordinate time of the Relatipy solution
  (dop853, rtol = 1e-13) at tau_end; status `short` if the last stored tau < 0.999 tau_end (none).
- Integrator choice: 20-orbit pilot of the Boost integrators (one process each, separate cores;
  ranking only). Medians at tol = 1e-13: dopri5 eps 4.2e-9 / 107 ms, cash_karp54 2.6e-8 / 83 ms
  (cash_karp54_classic identical), fehlberg78 2.0e-6 / 30 ms. Kept: dopri5 (lowest eps) and
  fehlberg78 (GYOTO default), the same pair as for PyGRO. `Legacy` is excluded: it warns that it
  does not compute proper time and its steps do not change with RelTol.
- Total wall time of the run: 136.0 s; all 1000 rows `ok`.

Total wall time of the run: 1894.1 s (first run, 14 settings: 940.8 s). The settings shared with
the first run reproduce it: identical eps, median wall times within 2%.
One `ok` row has eps = NaN: PyGRO dp45 goal=6, orbit 74 (fewer than 2 pericentres stored); the
notebook excludes it from the medians and draws that setting with an open marker.

## KerrGeoPy (added 2026-10-30)
- `work_precision_kerrgeopy.py`  time of the closed-form advance, 2 pi Omega_phi / Omega_r - 2 pi
  (`kg.StableOrbit(0.9, p, e, 1.0)` + `fundamental_frequencies()`, both timed); no integration, so no eps.
- `work-precision-kerrgeopy.csv` 100 rows, same columns (library = kerrgeopy, eps = NaN).
- KerrGeoPy 0.9.3; same protocol: one untimed warm-up, best of 3 repeats per orbit,
  `taskset -c 3 $PY work_precision_kerrgeopy.py` (2 s). Median 3.3 ms, range 3.0-4.9 ms per orbit.
- Notebook cell 61 draws it as a band (minimum to maximum) with a line at the median.
