"""
Symbolic spherical coordinates for Cartesian conversion rules.

This module provides :class:`Spherical`, a subclass of
:class:`~relatipy.symbolic.coordinates.base.BaseCoordinate` that expresses
Cartesian coordinates :math:`(x, y, z)` as functions of the generalized
coordinates :math:`x^1 = r`, :math:`x^2 = \\theta`, and :math:`x^3 = \\phi`,
where :math:`r \\ge 0` is the radial distance, :math:`\\theta` is the polar
angle measured from the :math:`z` axis, and :math:`\\phi` is the azimuth about
:math:`z`. The map implemented here is

.. math::

    x &= r \\sin\\theta \\cos\\phi,\\\\
    y &= r \\sin\\theta \\sin\\phi,\\\\
    z &= r \\cos\\theta.

The symbols :math:`x^i(t)` and the array ``xs_p`` are defined in
``relatipy.symbolic.coordinates.base`` and are used consistently with the
numeric spherical chart in ``relatipy.numeric.coordinates.spherical``.

See Also
--------
relatipy.numeric.coordinates.spherical :
    Numeric spherical convention :math:`(t, r, \\theta, \\phi)` and velocity
    maps.

Examples
--------
>>> from relatipy.symbolic.coordinates.spherical import Spherical
>>> s = Spherical()
>>> M = s.cartesian_conversion_rule()
>>> M.shape
(3, 1)
"""

from __future__ import annotations

import sympy as sp

from .base import BaseCoordinate, xs_p


class Spherical(BaseCoordinate):
    """
    Symbolic spherical chart with standard :math:`(r, \\theta, \\phi)` map.

    Generalized coordinates are :math:`x^1 = r`, :math:`x^2 = \\theta`,
    :math:`x^3 = \\phi` as functions of an affine parameter :math:`t` (see
    ``xs_p`` in the parent module). This class only supplies the static
    conversion to Cartesian components; velocity relations are obtained via
    base-class helpers.

    Notes
    -----
    The constructor takes no arguments. No attributes are defined beyond those
    of :class:`~relatipy.symbolic.coordinates.base.BaseCoordinate`.

    Examples
    --------
    >>> from relatipy.symbolic.coordinates.spherical import Spherical
    >>> type(Spherical()).__name__
    'Spherical'
    """

    def __init__(self) -> None:
        """
        Construct a symbolic spherical coordinate helper.

        Examples
        --------
        >>> from relatipy.symbolic.coordinates.spherical import Spherical
        >>> isinstance(Spherical(), Spherical)
        True
        """
        super().__init__()

    def cartesian_conversion_rule(self) -> sp.Matrix:
        """
        Return Cartesian coordinates as functions of :math:`(r, \\theta, \\phi)`.

        Components are built from ``xs_p[0]``, ``xs_p[1]``, and ``xs_p[2]``
        following

        .. math::

            x &= x^1 \\sin(x^2) \\cos(x^3),\\\\
            y &= x^1 \\sin(x^2) \\sin(x^3),\\\\
            z &= x^1 \\cos(x^2).

        Returns
        -------
        sympy.Matrix
            Column vector ``Matrix([[x], [y], [z]])`` of shape ``(3, 1)``.

        Examples
        --------
        >>> from relatipy.symbolic.coordinates.spherical import Spherical
        >>> s = Spherical()
        >>> out = s.cartesian_conversion_rule()
        >>> out.shape
        (3, 1)
        """
        x = xs_p[0] * sp.sin(xs_p[1]) * sp.cos(xs_p[2])
        y = xs_p[0] * sp.sin(xs_p[1]) * sp.sin(xs_p[2])
        z = xs_p[0] * sp.cos(xs_p[1])

        return sp.Matrix([x, y, z])
