Thin-disk MCMC roadmap
======================

This page records the scope of the thin-disk MCMC layer, what was
deliberately left out of it, where that work lives, and the suggested next
steps. It is a planning note for maintainers, not a user promise.

What this layer contains, and why
---------------------------------

The goal is one mergeable change that runs the five-star (S2, S29, S38,
S55, S301) Kerr MCMC with the thin disk of Piran et al. (2026), and nothing
else. Existing Kerr code (``_core``, ``_mcmc_core``, ``geodesic/``,
``metrics/kerr.py``, ``observables/kerr_mcmc.py``) is untouched; the layer
only adds files and additive exports.

- Native force layer (``native/src/force/``): the callback contract
  (``external_force.h``), the frame map and projected four-acceleration
  (``external.c``), the forced proper- and coordinate-time right-hand sides
  and their integrator adapters (``kerr_forced.c``), the generated Kerr
  kernel used when a force is active
  (``geodesic/physic/kerr_rhs_generated.c``) and the
  forced observable evaluator (``kerr_mcmc.c``).
- Thin-disk model (``native/src/force/models/``): ``thin_disk`` on the
  generic ``axisymmetric_table``, with ``elliptic`` and ``quadrature``.
- Binding: ``relatipy._force_core`` with the two forced observable
  evaluators, the thin-disk table and placed-model owners, and a callback
  evaluator used by the tests.
- Python: ``relatipy.forces`` (``ExternalForce``, ``ThinDisk``, the private
  bound handle), ``KerrMcmcMultiModel`` and ``KerrMcmcGaussianLikelihood``,
  ``PerturbedKerrMcmcModel`` and ``PerturbedKerrMcmcMultiModel``.
- Tests: native tests under ``native/tests/force/`` (see
  ``native/README.md``), with analytic test-only models in
  ``force_test_models.h``, and Python tests for the disk, the binding and
  the multi-star and perturbed models.
- Benchmark: ``tools/benchmarks/bench_perturbed_kerr_mcmc_multi.py``.

Compared with ``claude/kerr-external-forces``, which this layer was cut
from, the following were also stripped because the MCMC does not need them:

- ``rp_force_bl_to_cartesian`` and ``rp_force_four_acceleration``
  (validated direct evaluators) and ``rp_kerr_forced_integrator_jacobian``
  (analytic Jacobian of the inactive force);
- ``_force_core.integrate_kerr``, ``bl_to_cartesian``,
  ``four_acceleration``, ``CALLBACK_FAILURE`` and
  ``AxisymmetricModel.potential``;
- the flattened-grid and kink support of ``axisymmetric_table`` and
  ``rp_quadrature_merge_breakpoints``, used only by the spheroid
  (commits ``75749a1`` and ``13bac9a``); a disk uses a flattening of one and
  no kinks, so its tables are bitwise unchanged;
- the C benchmarks ``bench_forced.c`` and ``test_thin_disk_bench.c`` and the
  benchmarks ``bench_kerr_mcmc_multi.py`` and ``bench_mcmc_projection.py``.

The native tests that used the forced trajectory solve now integrate with
``rp_integrator_integrate`` and a step observer, with the same tolerances;
the thin-disk ``E'`` and ``L'`` check, a Python test of
``PerturbedKerr.orbit`` on the source branch, is a native test here.

What was left out, and where it lives
-------------------------------------

``claude/kerr-external-forces`` (head ``e36673d``)
   - ``relatipy._codegen``: the symbolic C99 optimizer, the Kerr RHS
     derivation and its verification; ``tools/codegen/`` with the
     specifications of ``kerr_rhs``, ``projection`` and ``ring_field``. The
     generator of ``native/src/geodesic/physic/kerr_rhs_generated.c`` lives
     there.
   - ``relatipy._jit``: the cached run-time C compiler.
   - ``rp.Force`` and ``rp.Potential``: forces defined at run time from
     SymPy expressions (``forces/compiled.py``), with the ``symbolic``
     optional dependency and its ``uv.lock`` entry.
   - ``rp.PerturbedKerr`` and ``PerturbedKerr.orbit``
     (``metrics/perturbed_kerr.py``, ``geodesic/_forced_orbit.py``) on the
     forced trajectory solve ``native/src/force/solution/solve.c``.
   - The generated ``projection`` and ``ring_field`` kernels (unused by the
     MCMC).
   - ``tools/benchmarks/bench_mcmc_projection.py`` (likelihood timing and
     emcee projection) and ``bench_kerr_mcmc_multi.py``.

``claude/kerr-generated-rhs``
   The homogeneous Kerr path (``_core`` and ``_mcmc_core``) switched to the
   generated kernel, about 2.7-2.9 times faster per MCMC call (five-star
   ``KerrMcmcMultiModel`` 20.6 to 7.5 ms), and the legacy-compatible Kerr
   connection removed. With it, an inactive and a zero active force share
   the homogeneous kernel, so ``test_active_zero_force`` checks equality.

``claude/force-extra-models``
   The built-in Plummer sphere and oblate power-law spheroid force models,
   which need the flattened-grid and kink support removed here.

``claude/funny-galileo-7lm47r``
   Research data under ``research/s301_friends``: S2 and S301 inputs, the
   element table of S301 and its neighbours, observation schedules, mock
   scenarios and the Kerr plus disk prototype report.

Suggested merge order
---------------------

1. This layer (``feature/thin-disk``) into ``main``.
2. The homogeneous switch of ``claude/kerr-generated-rhs``: it changes
   existing Kerr code, so it deserves its own review. Rebase it on this
   layer; ``evaluate_homogeneous`` in ``kerr_forced.c`` then always calls
   the shared kernel, and the documentation that the forced path is faster
   than the homogeneous one no longer holds.
3. The code generator, the run-time compiler, ``rp.Force``,
   ``rp.Potential`` and ``PerturbedKerr.orbit`` from
   ``claude/kerr-external-forces``, with the forced trajectory solve, the
   ``symbolic`` extra and the stripped direct evaluators.
4. The Plummer and spheroid models of ``claude/force-extra-models``,
   restoring the table flattening and kinks.

Next steps
----------

- Light propagation: the observable models use a straight-line Rømer
  delay. Shapiro delay and lensing (photon propagation in Kerr, and the
  disk's contribution) are not modelled.
- A 2PN check of the forced equation: ``test_newtonian_limit`` bounds the
  1PN residual of the force effect and its ``1/r`` scaling; a comparison
  with a 2PN integration of a point mass plus the same force would test the
  coupling one order further.
- An analytic Jacobian for Radau on the forced path: with an active force
  Radau uses finite differences (eight extra right-hand sides per
  Jacobian); the homogeneous analytic Jacobian plus the derivative of the
  projected force would make Radau competitive.
- ``relatipy._jit`` resolves the compiler, and puts its identity in the
  cache key, before it looks in the cache, so a process without a compiler
  cannot reuse a cached library. Look up the cache first (or key it without
  the compiler when loading) before run-time forces are relied on in pool
  workers.
- MCMC cost: ``bench_mcmc_projection.py`` measured 13.5 ms per five-star
  likelihood call with the disk (dop853, ``rtol = 1e-10``), about 73
  evaluations/s in one process and 177/s with four processes. For 43
  parameters, 172 walkers (four per parameter) and 2e4 steps, 3.4e6
  evaluations, that projects to about 13 h in one process and 5.4 h with
  four. The homogeneous switch of step 2 does not speed up the forced path,
  which already uses the generated kernel; Radau with an analytic Jacobian
  or a cheaper disk table lookup would.
