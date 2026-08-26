"""
Symbolic cylindrical coordinates for Cartesian conversion rules.

This module provides :class:`Cylindrical`, a subclass of
:class:`~relatipy.symbolic.coordinates.base.BaseCoordinate` that expresses
Cartesian coordinates :math:`(x, y, z)` as functions of the generalized
coordinates :math:`x^1 = \\rho`, :math:`x^2 = \\phi`, and :math:`x^3 = z`,
where :math:`\\rho \\ge 0` is the radial distance from the :math:`z` axis,
:math:`\\phi` is the azimuth about that axis, and :math:`z` is the usual
Cartesian height. The map implemented here is

.. math::

    x &= \\rho \\cos\\phi,\\\\
    y &= \\rho \\sin\\phi,\\\\
    z &= z.

The symbols :math:`x^i(t)` and the array ``xs_p`` are defined in
``relatipy.symbolic.coordinates.base`` and are used consistently with the
numeric cylindrical chart in ``relatipy.numeric.coordinates.cilindrical``.

See Also
--------
relatipy.numeric.coordinates.cilindrical :
    Numeric cylindrical convention :math:`(t, \\rho, \\phi, z)` and velocity
    maps.

Examples
--------
>>> from relatipy.symbolic.coordinates.cylindrical import Cylindrical
>>> c = Cylindrical()
>>> M = c.cartesian_conversion_rule()
>>> M.shape
(3, 1)
"""

from __future__ import annotations

import sympy as sp

from .base import BaseCoordinate, xs_p


class Cylindrical(BaseCoordinate):
    """
    Symbolic cylindrical chart with standard :math:`(\\rho, \\phi, z)` map.

    Generalized coordinates are :math:`x^1 = \\rho`, :math:`x^2 = \\phi`,
    :math:`x^3 = z` as functions of an affine parameter :math:`t` (see
    ``xs_p`` in the parent module). This class only supplies the static
    conversion to Cartesian components; velocity relations are obtained via
    base-class helpers.

    Notes
    -----
    The constructor takes no arguments. No attributes are defined beyond those
    of :class:`~relatipy.symbolic.coordinates.base.BaseCoordinate`.

    Examples
    --------
    >>> from relatipy.symbolic.coordinates.cylindrical import Cylindrical
    >>> type(Cylindrical()).__name__
    'Cylindrical'
    """

    def __init__(self) -> None:
        """
        Construct a symbolic cylindrical coordinate helper.

        Examples
        --------
        >>> from relatipy.symbolic.coordinates.cylindrical import Cylindrical
        >>> isinstance(Cylindrical(), Cylindrical)
        True
        """
        super().__init__()

    def cartesian_conversion_rule(self) -> sp.Matrix:
        """
        Return Cartesian coordinates as functions of :math:`(\\rho, \\phi, z)`.

        Components are built from ``xs_p[0]``, ``xs_p[1]``, and ``xs_p[2]``
        following

        .. math::

            x &= x^1 \\cos(x^2),\\\\
            y &= x^1 \\sin(x^2),\\\\
            z &= x^3.

        Returns
        -------
        sympy.Matrix
            Column vector ``Matrix([[x], [y], [z]])`` of shape ``(3, 1)``.

        Examples
        --------
        >>> from relatipy.symbolic.coordinates.cylindrical import Cylindrical
        >>> c = Cylindrical()
        >>> out = c.cartesian_conversion_rule()
        >>> out.shape
        (3, 1)
        """
        x = xs_p[0] * sp.cos(xs_p[1])
        y = xs_p[0] * sp.sin(xs_p[1])
        z = xs_p[2]

        return sp.Matrix([x, y, z])
