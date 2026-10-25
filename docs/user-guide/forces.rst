.. _forces-guide:

External forces and potentials
==============================

An extended mass around the black hole, such as a stellar disk or a cluster,
adds a Newtonian pull to each stellar orbit. :class:`relatipy.metrics.PerturbedKerr`
integrates the Kerr geodesic equation plus such a force, projected so that
the four-velocity stays normalized:

.. math::

   \frac{du^\mu}{d\tau} + \Gamma^\mu{}_{\alpha\beta} u^\alpha u^\beta
   = \left(g^{\mu\nu} + u^\mu u^\nu\right) F_\nu .

``F_nu`` is built from the Newtonian acceleration of the force in the
black-hole Cartesian frame (``z`` along the spin). The background metric is
unchanged: the force is a test-particle perturbation. The contract and the
native layout are described in :doc:`/development/external-forces`.

Forces
------

Every force lives in :mod:`relatipy.forces` (the common ones are also
``rp.<Name>``) and is an immutable value with named run-time parameters.

The built-in model, the thin disk of Piran et al. (2026), runs without any
extra dependency:

.. code-block:: python

   from astropy import units as u
   import relatipy as rp

   # Piran et al. (2026) fiducial disk: Sigma ~ 1/R out to the S301 apocentre
   disk = rp.ThinDisk(mass=1e3 * u.M_sun, r_out=1362.5 * u.au)
   disk.parameter_names        # ('disk_mass', 'normal_theta', 'normal_phi')

:class:`~relatipy.forces.ThinDisk` is a razor-thin disk with
``Sigma(R) ~ R^slope`` between ``r_in`` (default 0) and ``r_out``, each ring
softened by ``H = softening * R'`` (default 0.01). Other mass models, such as
a Plummer sphere, are written as a run-time
:class:`~relatipy.forces.Potential` (see below).

The disk reads its field from a table of the unit shape, built in native
code on first use (about 1 s and 11 MB) and cached per shape in each
process. The mass, the outer radius and the normal are written into a small
native model, so they change without a rebuild. The ``normal`` argument is a
vector; it is stored as the angles ``normal_theta`` (from the third axis) and
``normal_phi`` (from the first axis towards the second).

Orbits
------

:class:`~relatipy.metrics.PerturbedKerr` takes the Kerr mass and spin plus
one force, and its :meth:`~relatipy.metrics.PerturbedKerr.orbit` accepts the
initial conditions of :meth:`Kerr.orbit <relatipy.metrics.Kerr.orbit>`:

.. code-block:: python

   bh = rp.PerturbedKerr(mass=4.3e6 * u.M_sun, spin=0.5, force=disk)
   orbit = bh.orbit(a=1000 * u.au, e=0.88, inc=2.3 * u.rad)
   solution = orbit.solve(tau_span=(0 * u.s, 30 * u.yr), method="dop853")

   bh.kerr                     # the unperturbed background metric

The result is an ordinary :class:`~relatipy.geodesic.Solution` with the same
interpolation, coordinate views and plots. Here the disk normal is in the
black-hole frame, so the default ``(0, 0, 1)`` puts the disk in the
equatorial plane of the spin. Points to keep in mind:

- ``method`` is ``"radau"``, ``"dop853"`` or ``"dp45"``; ``projection_radau``
  raises :class:`ValueError` because its projector restores Carter's
  constant, which a force does not conserve. Radau uses a finite-difference
  Jacobian.
- ``get_E``, ``get_Lz`` and ``get_Q`` report the Kerr quantities, which a
  force changes. For a static potential axisymmetric about the spin,
  ``E' = -u_t exp(Phi)`` and ``L' = u_phi exp(Phi)`` are conserved instead.
- A failing force callback (for example the disk at its exact centre) is a
  numerical failure: ``Solution.status == -1``.
- A zero force reproduces :class:`~relatipy.metrics.Kerr` orbits to rounding
  error, not bit for bit, because the forced right-hand side uses its own
  generated Kerr kernel.

Forces defined at run time
--------------------------

:class:`~relatipy.forces.Potential` and :class:`~relatipy.forces.Force`
turn a definition into C, compile it with the system compiler and use it
like a built-in model. The definition is evaluated in geometric units of the
black hole: lengths in ``r_g = G M / c^2``, velocities in ``c``, ``Phi`` in
``c^2``. Parameters given as quantities are converted when the force is
bound; plain numbers are used as given.

.. code-block:: python

   # Miyamoto-Nagai disk in cylindrical coordinates of the black-hole frame
   mn = rp.Potential(
       lambda R, phi, z, m, a, b: -m / (R**2 + (a + (z**2 + b**2) ** 0.5) ** 2) ** 0.5,
       coords="cylindrical",
       parameters={"m": 1e3 * u.M_sun, "a": 500 * u.au, "b": 50 * u.au},
   )

   # The same kind of definition as a string, and a velocity-dependent force
   halo = rp.Potential("-k / sqrt(r**2 + b**2)", coords="spherical",
                       parameters={"k": 1e-4, "b": 2e4})
   drag = rp.Force(lambda R, phi, z, vR, vphi, vz, k: (-k * vR, -k * vphi, -k * vz),
                   coords="cylindrical", parameters={"k": 1e-12})

A callable receives the coordinates, then (for a
:class:`~relatipy.forces.Force`) the
velocity components in the same orthonormal basis, then the parameters, as
SymPy symbols. Strings use the names ``x, y, z, vx, vy, vz`` (Cartesian),
``R, phi, z, v_R, v_phi, v_z`` (cylindrical) or
``r, theta, phi, v_r, v_theta, v_phi`` (spherical). Cylindrical and
spherical coordinates are the flat coordinates of the black-hole Cartesian
frame, not Boyer--Lindquist ones.

Compilation happens on first use and needs:

- SymPy, installed with ``pip install relatipy[symbolic]``; without it
  :class:`ImportError` names the extra;
- a C11 compiler, ``$CC`` or ``cc``; without it :class:`RuntimeError` says
  so. The compiled libraries are shared objects, so run-time forces work on
  POSIX systems.

Libraries are cached in ``$RELATIPY_CACHE_DIR`` or ``~/.cache/relatipy/jit``
under a hash of the source, the compiler and the flags, so later processes
load them without compiling. A generated callback reports failure when an
acceleration component is not finite. String definitions are parsed with
``sympy.parse_expr``, which evaluates Python code: pass only trusted
strings.

Pickling stores the generated C source and the parameters, not the Python
definition, so a force built from a lambda can be sent to a multiprocessing
pool; a worker loads the cached library without SymPy.

MCMC with a force
-----------------

:class:`~relatipy.observables.PerturbedKerrMcmcMultiModel` extends the layout
of :class:`~relatipy.observables.KerrMcmcMultiModel` with a final force
block, ``force.<name>`` for each run-time parameter of the force in its
unit. Directions are read in observer axes ``(Dec, RA, away)``; with the
default normal a disk is seen face-on.

.. code-block:: python

   model = rp.PerturbedKerrMcmcMultiModel(
       astrometry_times, spectroscopy_times, reference_epoch=2000.0, force=disk
   )
   model.set_solver(method="dop853", rtol=1e-10, atol=1e-13)
   model.parameter_names[-3:]
   # ('force.disk_mass', 'force.normal_theta', 'force.normal_phi')

   likelihood = rp.KerrMcmcGaussianLikelihood(model, data)

One force handle is bound for the model and its parameters are rewritten in
place for every proposal with the proposal mass and orientation; nothing is
rebuilt or recompiled. Proposals outside a force's domain (a negative mass)
raise :class:`ValueError`, so the prior must exclude them. The model and the
likelihood are picklable; each pool worker builds its own disk table once.

For five S-stars with 100-150 epochs each, one likelihood call costs about
13 ms with the thin disk and 22 ms without a force (dop853, ``rtol=1e-10``).
At that tolerance the numerical error of the disk model is below 0.1 uas.
``tools/benchmarks/bench_mcmc_projection.py`` measures this and projects the
cost of a full emcee run.

:class:`~relatipy.observables.PerturbedKerrMcmcModel` is the single-star
model: :meth:`~relatipy.observables.PerturbedKerrMcmcModel.get_ra_dec_vr`
accepts ``force=...`` with physical inputs, and its legacy call takes the 13
:class:`~relatipy.observables.KerrMcmcModel` values followed by the force
block.

Accuracy notes
--------------

- The thin disk's potential and acceleration are separate interpolants that
  agree to about ``1e-5`` in the softened plane, so ``E'`` and ``L'`` computed
  from the tabulated potential hold to that level; for an analytic
  :class:`~relatipy.forces.Potential` they hold to the integration
  tolerance.
- Light travel uses the straight-line Rømer model of the Kerr observable
  models; no photon is traced through the perturbed field.
