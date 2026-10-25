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

The library ships one force, the razor-thin power-law disk of Piran et al.
(2026):

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
- ``projection_radau`` is not available: its projector restores Carter's
  constant, which a force does not conserve. Radau uses a finite-difference
  Jacobian.
- The disk's potential and acceleration are separate interpolants that agree
  to about ``1e-5`` in the softened plane; the acceleration matches the
  direct ring sum to about ``1e-7``.
- Light travel uses the straight-line Rømer model of the Kerr observable
  models; no photon is traced through the perturbed field.
