.. _geodesic-reference:

Geodesics and integration results
=================================

.. py:module:: relatipy.geodesic

The geodesic frontend separates a mutable :class:`~relatipy.geodesic.Orbit`
from immutable state and result records.  It evolves scalar timelike Kerr
orbits.  Native calculations use a mass-normalized, eight-component
Boyer--Lindquist state; public positions, velocities, and times are Astropy
quantities.

State records
-------------

:class:`~relatipy.geodesic.State` represents one state or a one-dimensional
series of states.  For a scalar state, scalar components have shape ``()``
and Cartesian vector properties such as ``xyz``, ``vxyz``, and ``uxyz`` have
shape ``(3,)``.  For ``n`` samples, scalar components have shape ``(n,)`` and
these vectors have shape ``(n, 3)``.  All stored public arrays and quantities
are read-only.

``tau`` and ``t`` are proper and coordinate time.  Cartesian positions and
velocities use length and length-per-time units; spherical and
Boyer--Lindquist angular coordinates are angles, and their corresponding
angular velocity and four-velocity components have angle-per-time units.
``ut`` is dimensionless.
:meth:`~relatipy.geodesic.State.orbital_elements` returns instantaneous
Kepler osculating elements computed from Cartesian position and coordinate
velocity.  Those elements describe a local conic, not the full Kerr
trajectory.  Their ``defined`` mask identifies degenerate fields.

Integration controls and outcomes
---------------------------------

:meth:`~relatipy.geodesic.Orbit.integrate` advances the mutable current
state to an absolute proper time and returns ``None``.
:meth:`~relatipy.geodesic.Orbit.solve` leaves the orbit unchanged and
returns a read-only :class:`~relatipy.geodesic.Solution`.
Both accept the following controls:

- ``method``: the default is ``"radau"``; ``"dop853"``, ``"dp45"``, and
  ``"projection_radau"`` are supported alternatives.
- ``rtol``: the default is ``None``, which selects ``1e-10``.
- ``atol``: the default is ``None``, which selects ``rtol`` times a
  characteristic scale of each native component at the starting state
  (see :ref:`integration-tolerances`).  A supplied ``atol`` is a finite,
  non-negative scalar or a read-only numeric vector of shape ``(8,)``; the
  native order is ``(t/T0, R/r_g, Theta, Phi, u^t, u^R, u^Theta, u^Phi)``.
- ``first_step`` and ``max_step``: the default is ``None``.  A supplied step
  control is a finite, positive proper-time quantity.

Explicit tolerances are used as given, but an unfit value issues
:class:`~relatipy.geodesic.IntegrationWarning` before integration.  These
controls scale local error estimates and do not guarantee global error.

Integration can end early in two ways:

- The outer horizon is an internal terminal event.  ``solve`` records it as
  ``status == 1`` and a :class:`~relatipy.geodesic.Termination`;
  ``integrate`` raises :class:`~relatipy.geodesic.IntegrationTerminated`.  In
  either case, the recorded terminal state is the last valid exterior state,
  not a localized crossing.
- A numerical failure uses ``status == -1`` in a solution or raises
  :class:`~relatipy.geodesic.IntegrationError` from ``integrate``.

Stored-result access
--------------------

``Solution[index]`` selects stored samples without interpolation.  An integer
returns a scalar ``State``; a slice, one-dimensional integer array, or
one-dimensional Boolean mask returns a series.  ``Solution.at(tau=...)`` or
``Solution.at(t=...)`` accepts exactly one scalar or one-dimensional time
quantity.  Exact stored times are copied exactly.  Between samples, Cartesian
position and coordinate velocity are interpolated independently; a native
batch reconstructs the Kerr-coordinate state without reintegration.  Queries
are subject to these conditions:

- Queries must lie within the stored domain.
- Coordinate-time queries also require that the stored ``t(tau)`` relation
  be uniquely monotone.
- Angular winding between adjacent samples that differ by at least ``pi``
  radians cannot be recovered uniquely.

API
---

.. autoclass:: relatipy.geodesic.Orbit
   :members:
   :show-inheritance:

.. autoclass:: relatipy.geodesic.InitialState
   :members:
   :exclude-members: state_vector
   :show-inheritance:

.. autoclass:: relatipy.geodesic.State
   :members:
   :exclude-members: state_vector
   :show-inheritance:

.. autoclass:: relatipy.geodesic.Solution
   :members:
   :exclude-members: success
   :special-members: __len__, __getitem__
   :show-inheritance:

.. autoclass:: relatipy.geodesic.IntegrationInfo
   :members:
   :show-inheritance:

.. autoclass:: relatipy.geodesic.Termination
   :members:
   :show-inheritance:

.. autoexception:: relatipy.geodesic.IntegrationError

.. autoexception:: relatipy.geodesic.IntegrationTerminated
   :exclude-members: termination, reason, tau, state

.. autoexception:: relatipy.geodesic.IntegrationWarning
