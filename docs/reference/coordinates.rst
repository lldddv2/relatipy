.. _coordinates-reference:

Coordinates and orbital elements
================================

.. py:module:: relatipy.coordinates

The :mod:`relatipy.coordinates` package provides immutable value objects for
public coordinate data, coordinate velocities, and orbital-element inputs.
All dimensional inputs are :class:`astropy.units.Quantity` objects. Scalar
values and one-dimensional series are accepted where indicated by each class;
series components must have the same shape. Stored arrays and quantities are
read-only.

Coordinate conventions
----------------------

Cartesian positions use ``(x, y, z)`` in a right-handed frame centered on the
black hole. Spherical coordinates use ``(r, theta, phi)``. Boyer--Lindquist
coordinates use ``(R, Theta, Phi)``. The containers store supplied values
without converting between these families.

Values produced by :class:`~relatipy.geodesic.Orbit` and
:class:`~relatipy.geodesic.Solution` use spin-aligned *oblate* Cartesian axes:

.. math::

   x = \sqrt{R^2 + a^2}\,\sin\Theta\cos\Phi, \qquad
   y = \sqrt{R^2 + a^2}\,\sin\Theta\sin\Phi, \qquad
   z = R\cos\Theta,

with ``a = spin * r_g``. Their spherical view is the Euclidean one,
``r = sqrt(x**2 + y**2 + z**2)`` and ``theta = atan2(sqrt(x**2 + y**2), z)``,
with ``phi = Phi``. For nonzero spin, ``r`` and ``theta`` therefore differ from
the Boyer--Lindquist ``R`` and ``Theta``. Velocities are derivatives with
respect to coordinate time ``t``; four-velocities are derivatives with respect
to proper time. Units must be compatible as follows:

- lengths with metres;
- times with seconds;
- angles with radians;
- radial velocities with length per time;
- angular velocities with angle per time.

``CartesianCoordinates.xyz`` has shape ``(3,)`` for one position and
``(n, 3)`` for a series. ``CartesianStateVector.xyz`` and ``vxyz`` have the
same two alternatives. Validation errors are raised as follows:

- invalid quantity types raise :class:`TypeError`;
- incompatible dimensions raise :class:`astropy.units.UnitConversionError`;
- invalid or mismatched shapes raise :class:`ValueError`.

Classical elements
------------------

:class:`OrbitalElements` stores an already computed classical Keplerian conic.
It does not derive elements from a state vector or define a dynamical
convention. Its fields are ``(a, e, inc, Omega, omega, f)``. ``a`` is a length,
``e`` is dimensionless, and the remaining fields are angles. The read-only
``defined`` mask has field order ``(a, e, inc, Omega, omega, f)`` and shape
``(6,)`` or ``(n, 6)``. A positive infinite semi-major axis is defined for an
exactly parabolic conic; ``NaN`` and other infinities are not.

:class:`KerrOrbitalElements` is a separate scalar input for a stable bound
timelike Kerr geodesic. Its fields satisfy:

- ``p`` is a positive physical semilatus rectum;
- ``e`` is in ``[0, 1)``;
- ``x`` is in ``[-1, 1]``;
- the three Mino phases are finite angles.

It validates only public units, scalar shapes, and basic parameter
ranges. Native construction checks stability and the Boyer--Lindquist chart
domain.

API
---

Cartesian
~~~~~~~~~

.. autoclass:: relatipy.coordinates.CartesianCoordinates
   :exclude-members: xyz
   :show-inheritance:

.. autoclass:: relatipy.coordinates.CartesianStateVector
   :show-inheritance:

Spherical and Boyer--Lindquist
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: relatipy.coordinates.SphericalCoordinates
   :show-inheritance:

.. autoclass:: relatipy.coordinates.SphericalVelocity
   :show-inheritance:

.. autoclass:: relatipy.coordinates.SphericalFourVelocity
   :show-inheritance:

.. autoclass:: relatipy.coordinates.BoyerLindquistCoordinates
   :show-inheritance:

.. autoclass:: relatipy.coordinates.BoyerLindquistVelocity
   :show-inheritance:

.. autoclass:: relatipy.coordinates.BoyerLindquistFourVelocity
   :show-inheritance:

Orbital elements
~~~~~~~~~~~~~~~~

.. autoclass:: relatipy.coordinates.OrbitalElements
   :exclude-members: defined
   :show-inheritance:

.. autoclass:: relatipy.coordinates.KerrOrbitalElements
   :show-inheritance:
