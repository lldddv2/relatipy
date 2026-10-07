.. _metrics-reference:

Kerr metric
===========

.. py:module:: relatipy.metrics

RelatiPy currently provides one public metric frontend,
:class:`~relatipy.metrics.Kerr`.  A metric is an immutable value: its mass
sets the conversion between the native geometric calculation and public
physical quantities, while ``spin`` is the dimensionless Kerr spin parameter.
The public frontend accepts ``0 <= spin <= 1``.  It does not presently expose
negative-spin parameterizations.

Units and reference radii
-------------------------

``mass`` must be a finite, positive scalar :class:`astropy.units.Quantity`
convertible to mass.  ``spin`` must be finite, dimensionless, and in the
closed interval ``[0, 1]``.  The length and time scales held by a metric are

.. math::

   r_g = \frac{G M}{c^2}, \qquad t_g = \frac{G M}{c^3}.

All returned radii are immutable scalar quantities.  ``r_g`` and ``r_s`` are
the gravitational and Schwarzschild length scales.  The ISCO, photon-orbit,
and horizon properties are Boyer--Lindquist coordinate radii; they are not
Cartesian distances.  :meth:`~relatipy.metrics.Kerr.r_ergosurface` accepts a scalar or
one-dimensional polar-angle quantity in ``[0, pi]`` and returns a radius with
the same shape.

Initial orbit construction
--------------------------

:meth:`~relatipy.metrics.Kerr.orbit` constructs one scalar, timelike
:class:`~relatipy.geodesic.Orbit`.
It accepts exactly one initial-condition family:

- stable bound Kerr elements;
- osculating Kepler elements;
- Cartesian coordinates and coordinate velocity;
- spherical coordinates and coordinate velocity;
- Boyer--Lindquist coordinates and coordinate velocity.

See :doc:`/user-guide/usage` for supported field
combinations and units.  The first public frontend is scalar; batch orbit
construction is not implemented.

API
---

.. autoclass:: relatipy.metrics.Kerr
   :members:
   :show-inheritance:

.. autoclass:: relatipy.metrics.Horizons
   :members:
   :show-inheritance:
