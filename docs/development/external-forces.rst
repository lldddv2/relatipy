External forces and the thin-disk MCMC
======================================

This page fixes the contract of the external-force layer used by the
perturbed MCMC models: a Kerr geodesic plus a Newtonian force. The forced
observable evaluator (``native/src/force/kerr_mcmc.c``) is also the
homogeneous one: ``rp_kerr_mcmc_evaluate`` in
``native/src/geodesic/kerr_mcmc.c``, behind ``_mcmc_core``, calls it with no
force, so both paths share the initial states, the Rømer arrival solve, the
observation terms and one generated Kerr kernel. The orbit integration of
``_core`` does not use the force layer. The user-facing description is
:doc:`/user-guide/forces`.

Equation of motion
------------------

In Boyer--Lindquist coordinates ``(t, r, theta, phi)``, signature
``(-,+,+,+)`` and geometric units with the black-hole mass set to one, the
perturbed geodesic is

.. math::

   \frac{du^\mu}{d\tau} + \Gamma^\mu{}_{\alpha\beta} u^\alpha u^\beta = f^\mu,
   \qquad f^\mu = \left(g^{\mu\nu} + u^\mu u^\nu\right) F_\nu .

The homogeneous term is the corrected Levi--Civita contraction of the Kerr
metric, the one used by every homogeneous production path:

the forced right-hand side always calls the generated kernel
``rp_kerr_geodesic_physic_evaluate_generated_rhs``
(``native/src/geodesic/physic/kerr_rhs_generated.c``), the kernel behind
``rp_kerr_null_geodesic_rhs``. An inactive force (a NULL force or a NULL
callback) therefore gives observables bitwise those of ``_mcmc_core``, and a
zero active force (a thin disk of zero mass) adds a zero four-acceleration,
so its right-hand side and its ``dop853`` observables are identical too.

The generated kernel is committed C, derived symbolically from the
corrected Kerr connection. Do not edit the file by hand.

The projector keeps ``u_mu f^mu = 0``, so ``g_(mu nu) u^mu u^nu = -1`` is
preserved. The covector ``F_nu`` is built from a Newtonian acceleration
``a`` evaluated in the black-hole Cartesian frame

.. math::

   x = \sqrt{r^2 + a_*^2}\,\sin\theta\cos\phi,\quad
   y = \sqrt{r^2 + a_*^2}\,\sin\theta\sin\phi,\quad
   z = r\cos\theta ,

as ``F_t = 0`` and ``F_i = a . dx/dx^i`` for ``i = r, theta, phi``, where
``a_*`` is the spin length. For a potential, ``a = -grad Phi`` and
``F_nu = -partial_nu Phi``. At ``M = 0`` these coordinates are flat oblate
spheroidal coordinates, so the coupling reduces exactly to Newtonian motion.

For a static potential that is axisymmetric about the spin axis,
``E' = -u_t exp(Phi)`` and ``L' = u_phi exp(Phi)`` are exact invariants of
this equation; the native tests use them to validate the projection.
Carter's constant is not conserved, so the forced evaluator rejects the
``projection_*`` methods when the force is active; it accepts them only with
an inactive force, which is how the homogeneous evaluator offers them
(``rp_kerr_forced_check_method`` lists the methods without projection). The
perturbed Python models reject the ``projection_*`` names for every force.
Radau uses finite differences for its Jacobian in this evaluator.

The observable evaluator integrates in coordinate time for arrival and
coordinate samples, every proper-time derivative (force term included) being
divided by ``u^t``, and in proper time for proper-time samples. The initial
states, Rømer arrival solve and observation terms are those of
``geodesic/kerr_mcmc.c``.

Native callback
---------------

Every force is one native function::

    typedef int (*rp_external_acceleration_fn)(
        const double position[3],  /* black-hole Cartesian frame */
        const double velocity[3],  /* dx/dt in the same frame */
        const void *model,         /* read-only model data */
        double acceleration[3]     /* Newtonian acceleration */
    );

It returns zero on success. A nonzero value is reported as
``RP_FORCE_STATUS_CALLBACK_FAILURE`` by the forced right-hand side, as an RHS
failure by the integrators and as an integration failure by the observable
evaluators (:class:`~relatipy.geodesic.IntegrationError`, ``-inf`` in the
likelihood). The tabulated thin disk fails at the centre ``r = 0`` and for a
NULL model.

Model data are immutable during an integration, own no global state and are
safe to share between threads. Units are geometric with the black-hole mass
set to one: lengths in ``G M / c^2``, times in ``G M / c^3``. The thin disk
reads an ``rp_axisymmetric_model``: a borrowed table, the mass, the length
scale and a rotation built from the symmetry axis. No Python callback runs
inside an integrator.

Layout
------

- ``native/src/force/``: callback type (``external_force.h``), the frame map
  and projection (``external.c``), the forced right-hand sides in proper and
  coordinate time (``kerr_forced.c``) and the forced observable evaluator
  (``kerr_mcmc.c``), parallel to ``geodesic/``.
- ``native/src/geodesic/physic/kerr_rhs_generated.{c,h}``: the generated
  Kerr kernel described above, beside the Kerr physics.
- ``native/src/force/models/``: the built-in ``thin_disk`` on the generic
  ``axisymmetric_table``, with its ``elliptic`` and ``quadrature`` helpers.
- ``native/tests/force/``: C tests, with analytic test-only models in
  ``force_test_models.h``; build commands are in ``native/README.md``.
- ``bindings/cython/_force_core.pyx``: the private ``relatipy._force_core``
  extension (see :doc:`/development/architecture`): the two forced
  observable evaluators, mirroring ``_mcmc_core``, the owners of thin-disk
  tables and their placed models, and a callback evaluator for tests.
- ``relatipy.forces``: the public ``ExternalForce`` base, ``ThinDisk``,
  ``SurfaceDensityDisk``, ``Force``, ``Potential``, ``get_force`` and the
  definition helpers ``elliptic_k``, ``elliptic_e`` and ``integral``, and
  the private bound-force handle.
- ``relatipy._codegen`` (private): the SymPy code generator that turns a
  force, a potential or a surface density into C.
- ``relatipy._jit`` (private): compiles that C with ``$CC`` or ``cc`` into a
  shared library cached in ``$RELATIPY_CACHE_DIR`` or
  ``~/.cache/relatipy/jit`` and loads it with ``ctypes``.
- ``relatipy.observables``: ``KerrMcmcMultiModel`` and
  ``KerrMcmcGaussianLikelihood``, and the perturbed models
  ``PerturbedKerrMcmcModel`` and ``PerturbedKerrMcmcMultiModel``, parallel to
  the homogeneous ones.

Python binding of a force
-------------------------

A force is an immutable value with named parameters. Binding it to a
black-hole mass gives a private handle with three things: the callback
address, the model-data address and a reference that keeps the table alive.
Parameters are converted to geometric units by their dimensions (a quantity
of length dimension ``L``, time ``T`` and mass ``M`` is divided by
``r_g**L t_g**T M**M``) and written into the model data in place, so an MCMC
sampler updates the mass, the disk mass or the normal between proposals
without rebuilding a table. Directions are given in observer axes
``(Dec, RA, away)`` and mapped through the proposal's spin frame. Handles
are never pickled with their addresses: tabulated forces pickle by their
constructor arguments, and tables are recreated on first use in each
process; ``Force`` and ``Potential`` pickle their generated C source and
parameters, and a process loads the compiled library from the cache or
compiles it on first use.

Thin disk
---------

The razor-thin power-law disk of Piran et al. (2026),
``Sigma(R) proportional to R^alpha`` with ``alpha = -1`` between ``R_in``
and ``R_out``, each ring softened with ``H = eps R'`` (``eps = 0.01``). Its
acceleration is a sum of ring fields with complete elliptic integrals,
tabulated once for unit mass and ``R_out = 1`` as the smooth combinations
``A = -r^3 a_d / d``, ``C = -r^3 a_z / z`` and ``P = -r Phi`` and
interpolated with tensor-product C2 cubic splines stored as bicubic Hermite
data, with a self-similar inner continuation and an exact multipole far
field. A table depends only on ``alpha``, ``R_in / R_out`` and ``eps``
(about 1 s and 11 MB to build); the disk mass, ``R_out`` (any length) and
the disk normal are model parameters. The potential and the acceleration are
separate interpolants that agree to about ``5e-7`` (median) and ``9e-6``
(worst) inside the softened plane, which bounds how well ``E'`` and ``L'``
computed from the tabulated potential are conserved.

The Python layer caches built tables per shape in each process (a few
shapes), and every placed model keeps its table referenced.

Known limitations
-----------------

- The observation model uses the straight-line Rømer delay. Shapiro delay
  and gravitational lensing are not modelled.
- The forced equation of motion has not been compared with a second
  post-Newtonian (2PN) integration. ``test_newtonian_limit``
  (``native/tests/force/test_newtonian_limit.c``) only bounds the residual
  against Newtonian motion by its first post-Newtonian (1PN) size.
- The observable evaluators, with or without a force, give Radau a
  finite-difference Jacobian, one extra right-hand-side evaluation per state
  component; the analytic Kerr Jacobian is not used there.
- ``relatipy._jit`` resolves the compiler and puts its identity in the cache
  key before it looks in the cache. A process without a C compiler therefore
  cannot reuse a library that is already cached.
- Proper-time samples requested with ``tau_eval`` can differ between
  integration methods by up to about ``1e-4``, with or without projection;
  ``tests/api/observables/test_kerr_mcmc_model.py`` notes this and compares
  the methods only for ``t_obs`` and ``t_eval`` samples.
