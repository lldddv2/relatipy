.. _observables-reference:

Observables
===========

.. py:module:: relatipy.observables

:mod:`relatipy.observables` provides :class:`KerrMcmcModel`, a direct model for
Kerr astrometry and redshift predictions, :class:`KerrMcmcMultiModel` for
several stars around one black hole, and :class:`KerrMcmcGaussianLikelihood`
for their Gaussian log-likelihood. It evolves a timelike Kerr geodesic
and applies a straight-line Rømer observation model. It does not trace photons
or include Shapiro delay or gravitational lensing.

General interface
-----------------

Use :meth:`KerrMcmcModel.get_ra_dec_vr` for physical inputs. ``kerr``
supplies:

- a positive mass quantity;
- dimensionless spin in ``[0, 1]``;
- a three-value observer-frame direction ``(Dec, RA, away)``.

``orbit`` accepts a supported
``Kerr.orbit`` initial-condition family. Cartesian, spherical, and classical
element inputs use observer axes; Boyer--Lindquist and bound-element inputs
use the spin-aligned frame.

``sol`` contains exactly one nonempty one-dimensional quantity named
``t_eval``, ``tau_eval``, or ``t_obs``. The input sequence may be unordered;
results retain its order. ``distance`` is a positive physical length. The
method returns read-only ``alpha`` and ``delta`` arrays in arcsec plus a
read-only ``v_los`` array in km/s. Each output has shape ``(n,)`` for the
``n`` supplied epochs.

The method raises :class:`TypeError`, :class:`ValueError`, or
:class:`astropy.units.UnitConversionError` for invalid mappings, values,
shapes, or units. It raises :class:`relatipy.geodesic.IntegrationError` when
native integration or observed-time matching fails.

Legacy interface
----------------

The callable interface takes a 13-value numeric proposal. It requires that:

- legacy astrometry and spectroscopy epochs have been supplied;
- :meth:`KerrMcmcModel.set_solver` has been supplied;
- the native extension is built; otherwise evaluation raises
  :class:`ImportError`.

Its parameter units and order are specified in the class documentation and
in :doc:`/user-guide/kerr-mcmc`. It returns read-only arcsec offsets and km/s
line-of-sight velocities.

Multi-star interface
--------------------

:class:`KerrMcmcMultiModel` evaluates one :class:`KerrMcmcModel` per star from
one numeric vector. It takes mappings from star name to astrometric and
spectroscopic epochs and one shared ``reference_epoch``;
:meth:`KerrMcmcMultiModel.set_solver` applies the same solver to every star.
The vector holds, in order:

.. list-table::
   :header-rows: 1
   :widths: 20 45 35

   * - Block
     - Names
     - Units
   * - ``shared``
     - ``D``, ``M``, ``xS0``, ``yS0``, ``vxS0``, ``vyS0``, ``v_LSR``
     - kpc, million solar masses, arcsec (two), arcsec/year (two), km/s
   * - ``spin``
     - ``spin_dec``, ``spin_ra``, ``spin_away``
     - dimensionless spin vector in observer axes ``(Dec, RA, away)``
   * - one per star ``S``
     - ``S.t_p``, ``S.a``, ``S.e``, ``S.inc``, ``S.Omega``, ``S.omega``
     - Julian years, arcsec, dimensionless, radians (three)

``parameter_names``, ``parameter_slices`` and ``parameter_units`` describe
this layout, and :meth:`~KerrMcmcMultiModel.pack` and
:meth:`~KerrMcmcMultiModel.unpack` convert between named values and vectors;
``pack`` also accepts convertible scalar quantities. Calling the model returns
a mapping from each star to its read-only ``(alpha, delta, v_los)`` arrays.
For one star the arrays equal those of :class:`KerrMcmcModel` called with
``(D, M, t_p, a, e, inc, Omega, omega, xS0, yS0, vxS0, vyS0, v_LSR)`` and the
spin vector.

:class:`KerrMcmcGaussianLikelihood` takes the model and, for each star, the
observed ``alpha``, ``delta`` and ``v_los`` with their errors ``alpha_err``,
``delta_err`` and ``v_los_err`` in arcsec and km/s. Calling it with a vector
returns the independent Gaussian log-likelihood, or ``-inf`` when native
evaluation fails. Proposals outside the model domain raise
:class:`ValueError`, so a prior must exclude them. Both classes are picklable,
for example for a ``multiprocessing`` pool passed to an emcee
sampler.

Perturbed models
----------------

:class:`PerturbedKerrMcmcModel` and :class:`PerturbedKerrMcmcMultiModel` are
the counterparts of :class:`KerrMcmcModel` and :class:`KerrMcmcMultiModel`
with one :class:`~relatipy.forces.ExternalForce`. They share initial
states, light model, units, solver settings and errors; the orbit adds the
force. :meth:`PerturbedKerrMcmcModel.get_ra_dec_vr` takes an optional
``force``; the legacy call appends the force block to the 13 values. The
multi-star vector appends one final block ``force.<name>`` for each of the
force's :attr:`~relatipy.forces.ExternalForce.parameter_names`, in its
:attr:`~relatipy.forces.ExternalForce.parameter_units`; for a
:class:`~relatipy.forces.ThinDisk` these are ``force.disk_mass`` in solar
masses (or the unit given) and ``force.normal_theta`` and
``force.normal_phi`` in radians, the normal in observer axes
``(Dec, RA, away)``. One bound force is shared by all stars and rewritten in
place for each proposal. :class:`KerrMcmcGaussianLikelihood` accepts the
perturbed multi-star model; both are picklable.

API
---

.. autoclass:: relatipy.observables.KerrMcmcModel
   :members: set_solver, get_ra_dec_vr, __call__
   :show-inheritance:

.. autoclass:: relatipy.observables.KerrMcmcMultiModel
   :members: stars, reference_epoch, astrometry_times, spectroscopy_times,
      n_parameters, parameter_names, parameter_units, parameter_slices,
      pack, unpack, set_solver, __call__
   :inherited-members:
   :special-members: __call__
   :no-show-inheritance:

.. autoclass:: relatipy.observables.KerrMcmcGaussianLikelihood
   :members: model
   :special-members: __call__

.. autoclass:: relatipy.observables.PerturbedKerrMcmcModel
   :members: force, n_force_parameters, set_solver, get_ra_dec_vr, __call__

.. autoclass:: relatipy.observables.PerturbedKerrMcmcMultiModel
   :members: force, stars, reference_epoch, astrometry_times, spectroscopy_times,
      n_parameters, parameter_names, parameter_units, parameter_slices,
      pack, unpack, set_solver, __call__
   :inherited-members:
   :special-members: __call__
   :no-show-inheritance:
