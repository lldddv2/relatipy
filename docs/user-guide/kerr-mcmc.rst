.. _kerr-mcmc:

Kerr MCMC observables
=====================

:class:`relatipy.observables.KerrMcmcModel` predicts sky-plane offsets and a
redshift-equivalent velocity from a timelike Kerr geodesic. It does not create
a public :class:`~relatipy.metrics.Kerr`, :class:`~relatipy.geodesic.Orbit`, or
:class:`~relatipy.geodesic.Solution` object. The native evaluator receives a
complete proposal, integrates it, and evaluates all requested observables.

General evaluator
-----------------

Construct an empty model and call ``get_ra_dec_vr`` with physical inputs. The
``kerr`` mapping requires:

- ``mass``, a positive mass quantity;
- ``spin``, a dimensionless scalar in ``[0, 1]``;
- ``vec``, three finite observer-frame components. ``vec`` gives the spin
  direction in ``(Dec, RA, away)`` axes. It is normalized for nonzero spin; a
  zero vector is valid only when ``spin`` is zero.

``distance`` is a positive observer distance quantity.

The ``orbit`` mapping accepts the same one-family initial conditions as
:meth:`~relatipy.metrics.Kerr.orbit`:

- classical elements;
- Cartesian inputs;
- spherical inputs;
- Boyer--Lindquist inputs;
- bound Kerr inputs.

Cartesian, spherical, and classical-element inputs use observer axes.
Boyer--Lindquist and bound inputs use the spin-aligned body frame. Optional
``t`` and ``tau`` set the physical coordinate-time and proper-time origins.

The ``sol`` mapping requires exactly one nonempty one-dimensional quantity:

- ``t_eval``, coordinate time;
- ``tau_eval``, proper time;
- ``t_obs``, arrival time at the observer.

Input times may be unordered. The returned arrays retain that order. The
mapping can also set:

- ``method`` (``"dop853"`` or ``"radau"``);
- positive scalar ``rtol`` and ``atol``;
- positive integer ``max_steps``.

Without an override, the general evaluator uses ``dop853``, ``rtol=1e-9``,
``atol=1e-12``, and ``max_steps=100000``.

.. code-block:: python

   import numpy as np
   from astropy import units as u
   from relatipy import KerrMcmcModel

   model = KerrMcmcModel()
   alpha, delta, v_los = model.get_ra_dec_vr(
       kerr={
           "mass": 4e6 * u.M_sun,
           "spin": 0.4,
           "vec": (0.0, 0.0, 1.0),
       },
       orbit={
           "a": 1000 * u.au,
           "e": 0.5,
           "inc": 0.6 * u.rad,
           "Omega": 0.4 * u.rad,
           "omega": 0.3 * u.rad,
       },
       sol={"t_eval": np.array([0.0, 1e5]) * u.s},
       distance=8 * u.kpc,
   )

``alpha`` and ``delta`` are read-only NumPy arrays in arcsec. ``v_los`` is a
read-only NumPy array in km/s. All three have shape ``(len(requested_times),)``.
They are numeric arrays, not Astropy quantities. ``alpha`` and ``delta`` are
offsets, not absolute right ascension and declination: observer-frame ``Y``
maps to ``alpha`` and ``X`` maps to ``delta``. Positive observer-frame ``Z``
points away from the observer.

Legacy fixed-epoch evaluator
----------------------------

The constructor also supports the fixed-epoch API. Supply astrometric and
spectroscopic Julian-year epochs plus ``reference_epoch``, then configure a
solver with :meth:`~relatipy.observables.KerrMcmcModel.set_solver` before calling the
model. A constructor ``spin_vector`` supplies the default observer-frame
spin; a call-time ``spin_vector`` overrides it for one proposal.

.. code-block:: python

   import numpy as np
   from relatipy import KerrMcmcModel

   model = KerrMcmcModel(
       [2002.33, 2003.0], [2002.33], reference_epoch=2000.0
   )
   model.set_solver(method="dop853", rtol=1e-8, atol=1e-10)
   params = np.array([
       8.33, 4.35, 2002.33, 0.1255, 0.8839,
       np.deg2rad(134.18), np.deg2rad(226.94), np.deg2rad(65.51),
       0.001, -0.002, 0.0001, -0.0002, 12.0,
   ])
   alpha, delta, v_los = model(params)

This API accepts exactly 13 finite values:

.. list-table:: Fixed-epoch parameter vector
   :header-rows: 1
   :widths: 22 20 58

   * - Name
     - Unit
     - Meaning
   * - ``D``
     - kpc
     - Positive distance to the system.
   * - ``M``
     - million solar masses
     - Positive central mass.
   * - ``t_p``
     - Julian year
     - Coordinate-time reference at the input periapsis state.
   * - ``a``
     - arcsec
     - Positive angular semimajor-axis input.
   * - ``e``
     - dimensionless
     - Eccentricity in ``[0, 1)``.
   * - ``inc``, ``Omega``, ``omega``
     - rad
     - Classical orbital angles.
   * - ``xS0``, ``yS0``
     - arcsec
     - Constant offsets added to ``alpha`` and ``delta``.
   * - ``vxS0``, ``vyS0``
     - arcsec/year
     - Linear drifts added to ``alpha`` and ``delta``.
   * - ``v_LSR``
     - km/s
     - Offset subtracted from ``v_los``.

The returned ``alpha`` and ``delta`` have length ``len(astrometry_times)``;
``v_los`` has length ``len(spectroscopy_times)``. Their respective input
orders are preserved. The input periapsis state uses an osculating Kepler
convention, so its ``a`` and ``e`` are not exact Kerr turning-point parameters
when spin is nonzero.

Numerical and physical limits
-----------------------------

The observation model uses small-angle sky projection and straight-line light
travel with Rømer delay. It does not trace photons or include Shapiro delay or
gravitational lensing. ``v_los`` is redshift-equivalent velocity, not the
star's coordinate velocity along the line of sight. Solver tolerances apply
to the internal geometric state. Check observable convergence for each fit.

Errors are reported as follows:

- invalid mappings, units, shapes, values, spins, or solver settings raise
  ``TypeError``, ``ValueError``, or an Astropy unit-conversion error;
- native integration or observed-time matching failures raise
  :class:`~relatipy.geodesic.IntegrationError`;
- a missing compiled extension raises :class:`ImportError`.
