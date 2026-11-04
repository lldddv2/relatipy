.. _forces-guide:

MCMC with a thin disk
=====================

An extended mass around the black hole, such as a stellar disk, adds a
Newtonian pull to each stellar orbit. The perturbed observable models
integrate the Kerr geodesic equation plus such a force, projected so that the
four-velocity stays normalized:

.. math::

   \frac{du^\mu}{d\tau} + \Gamma^\mu{}_{\alpha\beta} u^\alpha u^\beta
   = \left(g^{\mu\nu} + u^\mu u^\nu\right) F_\nu .

``F_nu`` is built from the Newtonian acceleration of the force in the
black-hole Cartesian frame (``z`` along the spin). The background metric is
unchanged: the force is a test-particle perturbation. The contract and the
native layout are described in :doc:`/development/external-forces`.

The thin disk
-------------

The built-in force is the razor-thin power-law disk of Piran et al. (2026):

.. code-block:: python

   from astropy import units as u
   import relatipy as rp

   # Piran et al. (2026) fiducial disk: Sigma ~ 1/R out to the S301 apocentre
   disk = rp.ThinDisk(mass=1e3 * u.M_sun, r_out=1362.5 * u.au)
   disk.parameter_names        # ('disk_mass', 'normal_theta', 'normal_phi')

:class:`~relatipy.forces.ThinDisk` has ``Sigma(R) ~ R^slope`` between
``r_in`` (default 0) and ``r_out``, each ring softened by
``H = softening * R'`` (default 0.01). Its field is read from a table of the
unit shape, built in native code on first use (about 1 s and 11 MB) and
cached per shape in each process. The mass, the outer radius and the normal
are written into a small native model, so they change without a rebuild.
The ``normal`` argument is a vector; it is stored as the angles
``normal_theta`` (from the third axis) and ``normal_phi`` (from the first
axis towards the second).

Forces written as text
----------------------

:func:`~relatipy.forces.get_force` builds a force from exactly one of a
potential, a force or a surface density written as text. Coordinates use
fixed names in the black-hole frame: ``x, y, z`` (Cartesian), ``R, phi, z``
(cylindrical) or ``r, theta, phi`` (spherical); a force may also use the
velocity components ``vx, vy, vz``, ``v_R, v_phi, v_z`` or
``v_r, v_theta, v_phi``. Every other name is a parameter and must be a key
of ``params``. The coordinate system is inferred from the names when
``coords`` is omitted.

.. code-block:: python

   from astropy import units as u
   import relatipy as rp

   # Plummer sphere; Phi is in units of c^2, so k = G m / c^2 is a length.
   plummer = rp.get_force(
       potential="-k / sqrt(r**2 + b**2)",
       params={"k": 1e-6 * u.au, "b": 1e3 * u.au},
   )
   plummer.parameter_names     # ('k', 'b')

   # Linear drag, one tuple of three components.
   drag = rp.get_force(force="(-k*vx, -k*vy, -k*vz)", params={"k": 1e-3})

   # Exponential thin disk; R is in units of r_out.
   disk = rp.get_force(
       surface_density="exp(-R/h)", params={"h": 0.25},
       mass=1e3 * u.M_sun, r_out=1362.5 * u.au,
   )
   disk.parameter_names        # ('disk_mass', 'normal_theta', 'normal_phi')

A potential or a force becomes a :class:`~relatipy.forces.Potential` or a
:class:`~relatipy.forces.Force`, whose ``params`` are run-time parameters:
they form the force block of the perturbed MCMC models and change without
recompiling. Expressions are evaluated in geometric units of the black hole
(``G = c = M = 1``); each parameter is converted by its unit when the force
is bound, and plain numbers are used as given. A surface density becomes a
:class:`~relatipy.forces.SurfaceDensityDisk`: ``Sigma(R)`` must be
axisymmetric, its ``params`` fix the shape of its table, and ``mass``,
``r_out`` and the optional ``r_in``, ``softening`` and ``normal`` are those
of the thin disk above, with the same run-time parameters.

Expressions may use numbers, ``pi`` and the functions ``sqrt exp log sin cos
tan asin acos atan atan2 sinh cosh tanh abs erf``. The classes also accept
Python callables that build SymPy expressions; there
:func:`~relatipy.forces.elliptic_k`, :func:`~relatipy.forces.elliptic_e` and
:func:`~relatipy.forces.integral` add complete elliptic integrals and
run-time quadrature.

These forces need two optional tools:

- SymPy (``pip install sympy``), to turn the text into C. A surface density
  is converted when it is constructed, a potential or a force when it is
  first used.
- A C11 compiler, ``$CC`` or ``cc``, to compile that C into a shared
  library the first time the force is bound, that is at the first prediction
  of a perturbed model that uses it.

A missing one raises :class:`ImportError` or :class:`RuntimeError`.
Libraries are cached in ``$RELATIPY_CACHE_DIR`` or, when it is unset,
``~/.cache/relatipy/jit``, so a later process with the same compiler loads
them without compiling. :class:`~relatipy.forces.Force` and
:class:`~relatipy.forces.Potential` pickle their generated C source, so
multiprocessing workers need no SymPy; a surface-density disk is rebuilt
from its text in each worker, which needs SymPy there.

Multi-star MCMC
---------------

:class:`~relatipy.observables.PerturbedKerrMcmcMultiModel` extends the layout
of :class:`~relatipy.observables.KerrMcmcMultiModel` (shared block, spin
block, one block per star) with a final force block, ``force.<name>`` for
each run-time parameter of the force in its unit. The disk normal is read in
observer axes ``(Dec, RA, away)``; with the default normal the disk is seen
face-on.

.. code-block:: python

   model = rp.PerturbedKerrMcmcMultiModel(
       astrometry_times, spectroscopy_times, reference_epoch=2000.0, force=disk
   )
   model.set_solver(method="dop853", rtol=1e-10, atol=1e-13)
   model.parameter_names[-3:]
   # ('force.disk_mass', 'force.normal_theta', 'force.normal_phi')

   likelihood = rp.KerrMcmcGaussianLikelihood(model, data)

``astrometry_times`` and ``spectroscopy_times`` map each star name to its
epochs; ``data`` maps each star to its observed ``alpha``, ``delta`` and
``v_los`` with their errors (see :doc:`/reference/observables`).

One force handle is bound for the model and its parameters are rewritten in
place for every proposal with the proposal mass and orientation; no table is
rebuilt. Proposals outside the disk's domain (a negative mass) raise
:class:`ValueError`, so the prior must exclude them; native failures give
``-inf``. The model and the likelihood are picklable, so they can be passed
to a ``multiprocessing`` pool of an emcee sampler; each worker builds its own
disk table once.

For five S-stars (S2, S29, S38, S55, S301) with 100-150 epochs each, one
likelihood call costs about 16 ms with the thin disk and 8 ms without a
force (dop853, ``rtol=1e-10``); both paths use the same generated Kerr
kernel, and the disk adds its table lookup and projection. At that
tolerance the numerical error of the disk model is below 0.1 uas.
``tools/benchmarks/bench_perturbed_kerr_mcmc_multi.py`` measures it.

:class:`~relatipy.observables.PerturbedKerrMcmcModel` is the single-star
model: :meth:`~relatipy.observables.PerturbedKerrMcmcModel.get_ra_dec_vr`
accepts ``force=...`` with physical inputs, and its legacy call takes the 13
:class:`~relatipy.observables.KerrMcmcModel` values followed by the force
block.

Accuracy notes
--------------

- A disk of zero mass evaluates the same Kerr right-hand side as the Kerr
  models and adds a zero four-acceleration, so ``dop853`` predictions are
  identical.
- The perturbed models accept ``dop853``, ``dp45`` and ``radau``. The
  ``projection_*`` methods of :class:`~relatipy.observables.KerrMcmcModel`
  are not available: their projector restores Carter's constant, which a
  force does not conserve. Radau uses a finite-difference Jacobian.
- The disk's potential and acceleration are separate interpolants that agree
  to about ``1e-5`` in the softened plane; the acceleration matches the
  direct ring sum to about ``1e-7``.
- Light travel uses the straight-line Rømer model of the Kerr observable
  models; no photon is traced through the perturbed field.
