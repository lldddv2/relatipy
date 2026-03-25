"""
Symbolic base chart for Cartesian maps and velocity relations.

This module defines the affine parameter :data:`t`, the generalized
coordinates :math:`x^i(t)` as SymPy functions, and the array :data:`xs_p`
that collects them. Subclasses of :class:`BaseCoordinate` implement
:meth:`~BaseCoordinate.cartesian_conversion_rule` to express Cartesian
components :math:`(x, y, z)` as functions of :math:`(x^1, x^2, x^3)`.
The base class then builds velocity relations used consistently across
symbolic coordinate charts (e.g. spherical, cylindrical).

Notes
-----
Generalized velocities are :math:`\\dot{x}^i = \\mathrm{d}x^i/\\mathrm{d}t`.
The method :meth:`BaseCoordinate.calculate_cartesian_conversion_rule_velocity`
forms the Jacobian :math:`J_{ij} = \\partial x^i / \\partial x^j` and returns
:math:`\\dot{\\mathbf{x}}_{\\mathrm{cart}} = J\\,\\dot{\\mathbf{x}}`.

The method :meth:`BaseCoordinate.calculate_relation_vi_dxi_dt` constructs
an orthonormal frame from normalized coordinate basis vectors
:math:`\\partial \\mathbf{x}/\\partial x^i` and projects the Cartesian
velocity onto that frame.

Examples
--------
>>> import sympy as sp
>>> from relatipy.symbolic.coordinates.base import t, xs_p
>>> sp.simplify(xs_p[0].diff(t))  # doctest: +ELLIPSIS
Derivative(x^1(t), t)
"""

from __future__ import annotations

import numpy
import sympy as sp
from itertools import product

t = sp.symbols("t")
x1 = sp.Function("x^1", positive=True)(t)
x2 = sp.Function("x^2", positive=True)(t)
x3 = sp.Function("x^3", positive=True)(t)

xs_p = numpy.array([x1, x2, x3])


class BaseCoordinate:
    """
    Abstract symbolic coordinate chart with Cartesian conversion hooks.

    Subclasses must implement :meth:`cartesian_conversion_rule` returning a
    ``sympy.Matrix`` column vector ``(x, y, z)`` in terms of ``xs_p``.
    Velocity helpers on this class use that map and the global symbols
    :data:`t` and :data:`xs_p`.

    Notes
    -----
    The default constructor takes no arguments and defines no instance
    attributes; concrete charts supply :meth:`cartesian_conversion_rule`.

    Examples
    --------
    >>> from relatipy.symbolic.coordinates.base import BaseCoordinate
    >>> from relatipy.symbolic.coordinates.spherical import Spherical
    >>> issubclass(Spherical, BaseCoordinate)
    True
    """

    def __init__(self) -> None:
        """
        Initialize the coordinate helper.

        Examples
        --------
        >>> from relatipy.symbolic.coordinates.base import BaseCoordinate
        >>> BaseCoordinate()  # doctest: +ELLIPSIS
        <relatipy.symbolic.coordinates.base.BaseCoordinate object at ...>
        """
        pass

    def calculate_relation_vi_dxi_dt(self) -> sp.Matrix:
        """
        Components of Cartesian velocity in a normalized coordinate frame.

        Let :math:`\\mathbf{x}(x^1,x^2,x^3)` be the Cartesian position from
        :meth:`cartesian_conversion_rule`. For each :math:`i`, define the
        (unnormalized) coordinate basis vector
        :math:`\\mathbf{e}_i = \\partial \\mathbf{x} / \\partial x^i`,
        normalize it, and form the Cartesian velocity via the chain rule

        .. math::

            \\dot{\\mathbf{x}} = \\sum_{j}
            \\frac{\\partial \\mathbf{x}}{\\partial x^j}\\,\\dot{x}^j.

        This method returns the three dot products
        :math:`\\dot{\\mathbf{x}}\\cdot\\hat{\\mathbf{e}}_i` as a column
        matrix, with ``Abs`` replaced by the identity for simplification.

        Returns
        -------
        sympy.Matrix
            Column vector of shape ``(3, 1)`` with the projected velocity
            components.

        Raises
        ------
        AttributeError
            If :meth:`cartesian_conversion_rule` is not implemented on the
            subclass.

        Examples
        --------
        >>> from relatipy.symbolic.coordinates.spherical import Spherical
        >>> vs = Spherical().calculate_relation_vi_dxi_dt()
        >>> vs.shape
        (3, 1)
        """
        xs = self.cartesian_conversion_rule()
        # Derivatives of generalized coordinates with respect to t
        xs_prime_dot = sp.Matrix([xs_p[0], xs_p[1], xs_p[2]]).diff(t)

        es = [0, 0, 0]

        for i in range(3):
            es_temp = xs.diff(xs_p[i])
            es_norm = sp.sqrt((es_temp.T @ es_temp)[0, 0]).simplify()
            es[i] = sp.simplify(es_temp / es_norm)

        # Cartesian velocity (chain rule)
        v_cart = sp.zeros(1, 3)
        for i, j in product(range(3), repeat=2):
            v_cart[i] += sp.diff(xs[i], xs_p[j]) * xs_prime_dot[j]

        # Projections v·e_i
        v0 = sp.simplify(v_cart.dot(es[0]))
        v1 = sp.simplify(v_cart.dot(es[1]))
        v2 = sp.simplify(v_cart.dot(es[2]))

        vs = sp.Matrix([v0, v1, v2]).replace(sp.Abs, lambda x: x)
        vs.simplify()

        return vs

    def calculate_cartesian_conversion_rule_velocity(self) -> sp.Matrix:
        """
        Cartesian velocity from generalized coordinates and velocities.

        Computes :math:`\\dot{x}`, :math:`\\dot{y}`, :math:`\\dot{z}` using
        the Jacobian of the map
        :math:`(x^1,x^2,x^3)\\mapsto(x,y,z)` given by
        :meth:`cartesian_conversion_rule` and the derivatives
        :math:`\\dot{x}^i = \\mathrm{d}x^i/\\mathrm{d}t`.

        .. math::

            \\begin{pmatrix}\\dot{x}\\\\ \\dot{y}\\\\ \\dot{z}\\end{pmatrix}
            = J\\,
            \\begin{pmatrix}\\dot{x}^1\\\\ \\dot{x}^2\\\\ \\dot{x}^3
            \\end{pmatrix},
            \\quad
            J_{ij} = \\frac{\\partial x^i}{\\partial x^j}.

        Returns
        -------
        sympy.Matrix
            Simplified column vector ``Matrix([[dx/dt], [dy/dt], [dz/dt]])``
            of shape ``(3, 1)``.

        Raises
        ------
        AttributeError
            If :meth:`cartesian_conversion_rule` is not implemented on the
            subclass.

        Examples
        --------
        >>> from relatipy.symbolic.coordinates.spherical import Spherical
        >>> v = Spherical().calculate_cartesian_conversion_rule_velocity()
        >>> v.shape
        (3, 1)
        """
        # Cartesian coordinates as functions of x^i(t)
        xs = self.cartesian_conversion_rule()

        # Partial derivatives with respect to each generalized coordinate
        J = sp.Matrix(
            [[sp.diff(xs[i], xs_p[j]) for j in range(3)] for i in range(3)]
        )

        # Time derivatives of x^i (generalized velocities)
        xs_prime_dot = sp.Matrix([sp.diff(q, t) for q in xs_p])

        # Cartesian velocity: v = J * q_dot
        v_cart = J * xs_prime_dot

        return sp.simplify(v_cart)
