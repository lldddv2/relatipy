Usage
=====

Kerr geometry and one orbit
---------------------------

A :class:`~relatipy.metrics.Kerr` object fixes the black-hole mass and
its dimensionless spin. The supported spin interval is ``0 <= spin <= 1``.
The mass sets the physical length scale :attr:`~relatipy.metrics.Kerr.r_g`
and time scale used inside the native solver. Public coordinates and times are
Astropy quantities.

This example constructs a timelike orbit from a Cartesian position and
coordinate velocity, then integrates from the saved initial state:

.. code-block:: python

   from astropy import units as u
   from astropy.constants import c
   from relatipy import Kerr

   bh = Kerr(mass=1 * u.Msun, spin=0.5)
   orb = bh.orbit(x=12 * bh.r_g, vy=0.1 * c)
   sol = orb.solve(tau_span=(0 * u.s, 2e-8 * u.s), method="dp45")

   assert sol.status == 0
   assert sol.success
   assert orb.tau == orb.initial.tau  # solve leaves the orbit unchanged
   print(sol.tau, sol.xyz.shape)

To request particular output times, pass a nonempty, one-dimensional, finite,
strictly increasing ``tau_eval`` sequence. A single requested time is valid.

.. code-block:: python

   import numpy as np

   times = np.array([0, 1e-8, 2e-8]) * u.s
   sampled = orb.solve(tau_eval=times, method="dp45")
   assert len(sampled) == len(times)

Only one initial-condition family may be supplied. They are classical elements
``(a, e, inc, Omega, omega, f)``, Cartesian position and velocity, spherical
position and coordinate velocity, Boyer--Lindquist position and coordinate
velocity, or bound Kerr parameters ``(p, e, x, q_r0, q_theta0, q_phi0)``.
The bound family also accepts a
:class:`~relatipy.coordinates.KerrOrbitalElements` instance. In that family,
``x`` is the dimensionless Kerr inclination parameter, so it is not a
Cartesian coordinate. Missing velocity components default to zero. The
classical-element angles use a right-handed Cartesian frame centered on the
black hole, with ``z`` aligned to the spin axis. The initial Boyer--Lindquist
radius must exceed the outer horizon.

To advance the mutable current point, call :meth:`~relatipy.geodesic.Orbit.integrate`
with an absolute proper time. :meth:`~relatipy.geodesic.Orbit.reset`
restores the initial point, while :meth:`~relatipy.geodesic.Orbit.copy`
creates an independently evolving orbit. A target before the current point
restarts from the saved initial conditions; a target before the initial time
is invalid.

Numerical controls and outcomes
-------------------------------

``integrate`` and ``solve`` share five controls: ``method``, ``rtol``,
``atol``, ``first_step``, and ``max_step``. The default method is
``"radau"``; ``"dop853"``, ``"dp45"``, and ``"projection_radau"``
are explicit alternatives. Method names are case-sensitive.
The defaults are ``rtol=None`` and ``atol=None``, which select tolerances
automatically (see :ref:`integration-tolerances`), and ``None`` for both
step controls. Supplied step controls are positive proper-time quantities.
``atol`` is a nonnegative number or a numeric array of shape ``(8,)``,
ordered as the native state
``(t/T0, R/r_g, Theta, Phi, u^t, u^R, u^Theta, u^Phi)``.
These tolerances weight local error in the normalized state; they do not
bound the global trajectory error.

.. _integration-tolerances:

Automatic and explicit tolerances
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Omitted tolerances are chosen automatically for every method.
``rtol=None`` selects ``1e-10``. ``atol=None`` selects the shape ``(8,)``
array ``rtol * s``, where ``rtol`` is the effective relative tolerance
(automatic or supplied) and ``s`` holds one characteristic scale per native
component of the starting state: ``s_i = max(abs(y0_i), floor_i)``. With
``R0`` the starting radius in units of ``r_g``, the floors are Newtonian
orders of magnitude:

.. list-table::
   :header-rows: 1

   * - Component
     - ``t/T0``
     - ``R/r_g``
     - ``Theta``
     - ``Phi``
     - ``u^t``
     - ``u^R``
     - ``u^Theta``
     - ``u^Phi``
   * - Floor
     - ``R0**1.5``
     - ``R0``
     - ``1``
     - ``1``
     - ``1``
     - ``R0**-0.5``
     - ``R0**-1.5``
     - ``R0**-1.5``

``rtol`` and ``atol`` are resolved independently, so supplying only
``rtol`` still scales the automatic ``atol``. A solution records the
effective values in the ``rtol`` and ``atol`` fields of its
:class:`~relatipy.geodesic.IntegrationInfo` (``Solution.integration``); an
automatic ``atol`` is a read-only array of shape ``(8,)``.

Explicit tolerances are used as given. Before integrating, both
:meth:`~relatipy.geodesic.Orbit.integrate` and
:meth:`~relatipy.geodesic.Orbit.solve` issue
:class:`~relatipy.geodesic.IntegrationWarning` when an explicit tolerance is
unfit for the starting state:

* ``rtol > 1e-6``;
* ``0 < rtol < 100 * eps``, which double precision cannot attain;
* ``atol_i > 1e-6 * s_i`` for any of ``t/T0``, ``R/r_g``, ``Theta``,
  ``Phi``, or ``u^t``.

The message names each cause and suggests omitting ``rtol`` and ``atol``.
The warning does not change the tolerances used, the returned status, or
the exceptions raised.

.. code-block:: python

   import numpy as np
   from astropy import units as u
   from relatipy import Kerr

   orbit = Kerr(mass=1 * u.Msun, spin=0.01).orbit(x=1 * u.au, vy=30 * u.km / u.s)
   taus = np.linspace(0, 1, 300) * u.yr

   solution = orbit.solve(tau_eval=taus)    # automatic tolerances
   solution.integration.rtol                # 1e-10
   orbit.solve(tau_eval=taus, rtol=1e-3)    # IntegrationWarning: unfit tolerances

The automatic choice is a heuristic based on the starting state. It does
not guarantee the global trajectory error, for example on unbound
trajectories that travel far from their starting radius or near the
horizon.

``"projection_radau"`` uses the Radau integrator and projects each accepted
native step outside the horizon toward the initial timelike norm
``-1``, energy, axial angular momentum, and Carter constant. It keeps the
step position fixed while correcting the four-velocity. A failed projection
is a numerical failure; ``solve`` reports a partial result when possible,
and ``integrate`` raises ``IntegrationError``. Both retain the last valid
state. Projection acts on accepted native steps, not on intermediate stages
or :meth:`~relatipy.geodesic.Solution.at` interpolation. Requested
``tau_eval`` samples between accepted steps are interpolated and do not
inherit the projection guarantee. Conserved quantities at accepted steps
remain subject to floating-point error. Projection does not guarantee
global trajectory accuracy.

The returned :class:`~relatipy.geodesic.Solution` has read-only samples,
``status``, ``success``, ``message``, and numerical settings and counts in
``integration``. ``status == 0`` means the requested endpoint was reached;
``status == 1`` means an internal terminal event stopped the solver;
``status == -1`` means numerical failure with a partial result. For a
confirmed outer-horizon crossing, ``termination.reason`` is
``"outer_horizon"``. Its state and proper time are the last stored valid
point outside the horizon, not the exact crossing point. An attempted step
that fails internally before a crossing is confirmed is a numerical failure.
The native event check runs after each accepted step and tests
``R <= r_+``; no crossing point is localized.
With ``tau_eval``, a partial result keeps only the requested samples reached
before the failure or event; if none was reached, ``solve`` raises
:class:`~relatipy.geodesic.IntegrationError` or
:class:`~relatipy.geodesic.IntegrationTerminated` instead of returning
unrequested samples.
``integrate`` instead raises :class:`~relatipy.geodesic.IntegrationTerminated`
for a terminal event or :class:`~relatipy.geodesic.IntegrationError`
for numerical failure, retaining its last valid point.

Querying a solution
-------------------

:meth:`relatipy.geodesic.Solution.at` accepts exactly one Astropy time
quantity, ``tau`` (proper time) or ``t`` (coordinate time). A scalar query
returns a scalar :class:`~relatipy.geodesic.State` with ``xyz.shape == (3,)``;
a one-dimensional query returns a state series with
``xyz.shape == (m, 3)``. Exact stored rows are selected without
reconstruction. Between-sample queries interpolate Cartesian position and
coordinate velocity separately, then reconstruct the remaining state views
in one native batch call. They do not reintegrate the orbit.

.. code-block:: python

   point = sol.at(tau=sol.tau[0])
   assert point.xyz.shape == (3,)

The intermediate azimuth follows the stored angular phase, but samples whose
angles differ by at least ``pi`` radians cannot uniquely determine a full
turn count. The Boyer--Lindquist chart policy requires a domain error outside
``0 < Theta < pi`` or when ``abs(sin(Theta)) <= 64 * DBL_EPSILON``. The guard
indicates numerical proximity to the axis, not a confirmed physical crossing.
Its full application across integration paths remains under validation; an
integration stage may fail before the threshold is reached.

Osculating elements
--------------------

``orb.orbital_elements()``, ``point.orbital_elements()``, and
``sol.orbital_elements()`` expose the instantaneous classical Kepler conic
computed from Cartesian position and coordinate velocity. Its radial or
parabolic character does not classify the complete Kerr trajectory.
An exactly zero computed Kepler energy produces ``a = +inf`` as a length
quantity. Undefined angles for radial or numerically near-radial states are
``NaN`` angular quantities. Other valid physical state fields remain finite;
these auxiliary values do not prevent integration or solution lookup.

Use the read-only ``defined`` mask to identify available element fields:

.. code-block:: python

   elements = point.orbital_elements()
   assert elements.defined.shape == (6,)
   # Columns: a, e, inc, Omega, omega, f.
   available_angles = elements.defined[2:]
   assert sol.orbital_elements().defined.shape == (len(sol), 6)

The mask marks ``+inf`` semi-major axis as ``True`` and undefined ``NaN``
fields as ``False``. Scalar elements use shape ``(6,)``; series use
``(n, 6)``. Copying and selection retain the representation, including
queries between solution samples. No new parabolic tolerance is used, so
roundoff can produce finite semi-major axis near zero energy. Ordinary
circular and equatorial conventions remain unchanged. Interpolation
reconstructs each query's own elements and can change conic classification.
Use a supported position/velocity initial-condition family for a parabolic
state. The bound Kerr ``p`` input family represents stable bound geodesics;
it is distinct from this instantaneous Kepler classification.

Plotting
--------

Matplotlib is optional and is imported only when a plotting method is called.
Install it with ``pip install "relatipy[plot]"``.
:meth:`~relatipy.geodesic.Orbit.preview` draws the current
osculating Kepler conic and marks the current position. It is a preview of
the instantaneous classical orbit, **not** the integrated Kerr path. The
optional outer-horizon and prograde/retrograde ISCO circles are equatorial
references in the oblate Cartesian ``xy`` plane; they are not 3D surfaces.
``preview`` raises ``ValueError`` for undefined required angles or a
parabolic osculating conic.
:meth:`~relatipy.geodesic.Solution.plot` draws the stored numerical samples
of an actual integrated trajectory. Both methods accept ``projection="3d"``,
``"xy"``, ``"xz"``, ``"yz"``, or ``"views"``. ``views`` is the default and
produces three orthogonal planes. With ``interactive=False`` they return
Matplotlib figures; with ``interactive=True`` they return Plotly figures. If
``interactive`` is omitted, ``"3d"`` selects Plotly and other projections
select Matplotlib.

.. code-block:: python

   preview_fig, preview_ax = orb.preview(projection="xy")
   path_fig, path_ax = sol.plot(projection="xy")

See :doc:`api` for signatures and :doc:`integration-roadmap` for the
implementation limits.
