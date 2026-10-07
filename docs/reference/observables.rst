.. _observables-reference:

Observables
===========

.. py:module:: relatipy.observables

:mod:`relatipy.observables` provides :class:`KerrMcmcModel`, a direct model for
Kerr astrometry and redshift predictions. It evolves a timelike Kerr geodesic
and applies a straight-line Rømer observation model. It does not trace photons
or include Shapiro delay or gravitational lensing.

General interface
-----------------

Use :meth:`KerrMcmcModel.get_ra_dec_vr` for physical inputs. ``kerr`` supplies
a positive mass quantity, dimensionless spin in ``[0, 1]``, and a three-value
observer-frame direction ``(Dec, RA, away)``. ``orbit`` accepts a supported
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

The callable interface takes a 13-value numeric proposal only after legacy
astrometry and spectroscopy epochs and :meth:`KerrMcmcModel.set_solver` have
been supplied. Its parameter units and order are specified in the class
documentation and in :doc:`kerr-mcmc`. It returns read-only arcsec offsets and
km/s line-of-sight velocities. The native extension must be built; otherwise
evaluation raises :class:`ImportError`.

API
---

.. autoclass:: relatipy.observables.KerrMcmcModel
   :members: set_solver, get_ra_dec_vr, __call__
   :show-inheritance:
