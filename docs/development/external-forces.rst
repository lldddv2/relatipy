External forces and potentials
==============================

This page fixes the contract of the external-force layer: a perturbed Kerr
geodesic integrated beside, not inside, the homogeneous Kerr path. Existing
Kerr modules, integrators and the ``_core`` and ``_mcmc_core`` extensions are
not modified; the new layer only calls them.

Equation of motion
------------------

In Boyer--Lindquist coordinates ``(t, r, theta, phi)``, signature
``(-,+,+,+)`` and geometric units with the black-hole mass set to one, the
perturbed geodesic is

.. math::

   \frac{du^\mu}{d\tau} + \Gamma^\mu{}_{\alpha\beta} u^\alpha u^\beta = f^\mu,
   \qquad f^\mu = \left(g^{\mu\nu} + u^\mu u^\nu\right) F_\nu .

The homogeneous term is the existing optimized Kerr right-hand side. The
projector keeps ``u_mu f^mu = 0``, so ``g_(mu nu) u^mu u^nu = -1`` is
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
this equation; tests use them to validate the projection. Carter's constant
is not conserved, so ``projection_radau`` is rejected when a force is set.

Native callback
---------------

Every force, built in or compiled at run time, is one native function::

    typedef int (*rp_external_acceleration_fn)(
        const double position[3],  /* black-hole Cartesian frame */
        const double velocity[3],  /* dx/dt in the same frame */
        const void *model,         /* read-only model data */
        double acceleration[3]     /* Newtonian acceleration */
    );

It returns zero on success; a nonzero value is reported as a right-hand-side
failure. Model data are immutable during an integration, own no global
state, and are safe to share between threads. Units are geometric with the
black-hole mass set to one: lengths in ``G M / c^2``, times in ``G M / c^3``.
No Python callback runs inside an integrator.

Layout
------

- ``native/src/force/``: callback types, the frame map and projection, the
  forced right-hand sides (proper and coordinate time), the forced trajectory
  solve and the forced observable model, parallel to ``geodesic/``.
- ``native/src/force/models/``: built-in force models. Generated sources
  carry a header naming the ``tools/codegen`` specification that regenerates
  them.
- ``native/tests/force/``: C tests.
- ``bindings/cython/_force_core.pyx``: the private ``relatipy._force_core``
  extension, added to ``setup.py`` beside the existing two.
- ``relatipy._codegen``: symbolic optimizer that turns SymPy expressions,
  matrices and tensors into optimized C99. SymPy is an optional dependency
  (``relatipy[symbolic]``) imported only when needed.
- ``relatipy._jit``: compiles generated C with the system compiler into a
  cached shared library and returns the callback address.
- ``relatipy.forces``: public ``Force``, ``Potential`` and built-in models.
- ``relatipy.metrics.PerturbedKerr``: Kerr plus one force, with the same
  orbit API as :class:`~relatipy.metrics.Kerr`.
- ``relatipy.observables``: multi-star observable models, homogeneous and
  perturbed, for MCMC.

Built-in model
--------------

The default model is the razor-thin power-law disk of Piran et al. (2026),
``Sigma(R) proportional to R^alpha`` with ``alpha = -1`` between ``R_in`` and
``R_out``, softened with ``H = 0.01 R'``. Its acceleration is a sum of ring
fields with complete elliptic integrals, tabulated once in the disk frame and
interpolated with a bicubic scheme that is at least C1; the disk mass scales the table and the
disk normal enters through a rotation, so both can change between MCMC
proposals without rebuilding the table.
