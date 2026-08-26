"""
Symbolic Boyer–Lindquist coordinates for Kerr-related charts.

This module defines :class:`BoyerLindquist`, a symbolic coordinate system whose
generalized coordinates :math:`(x^1, x^2, x^3)` are identified with the radial,
polar, and azimuthal Boyer–Lindquist coordinates :math:`(r, \\theta, \\phi)`.
The map to Cartesian spatial coordinates uses the spin parameter :math:`a` in
the usual way, with auxiliary radius :math:`\\varrho = \\sqrt{r^2 + a^2}`.

.. math::

    x &= \\varrho \\sin\\theta \\cos\\phi,\\\\
    y &= \\varrho \\sin\\theta \\sin\\phi,\\\\
    z &= r \\cos\\theta.

Notes
-----
The symbols :math:`x^i(t)` and the array ``xs_p`` are shared with
:mod:`relatipy.symbolic.coordinates.base`; see that module for the convention
that :math:`x^1`, :math:`x^2`, and :math:`x^3` are positive functions of proper
or coordinate time.

Examples
--------
>>> import sympy as sp
>>> from relatipy.symbolic.coordinates.boyer_lindquist import BoyerLindquist
>>> bl = BoyerLindquist()
>>> M = bl.cartesian_conversion_rule()
>>> M.shape
(3, 1)
"""

from .base import BaseCoordinate, xs_p
import sympy as sp


class BoyerLindquist(BaseCoordinate):
    """
    Boyer–Lindquist chart with symbolic map to Cartesian coordinates.

    Generalized coordinates are :math:`(x^1, x^2, x^3) \\equiv (r, \\theta,
    \\phi)` in the sense of the parent class, together with a constant spin
    parameter :math:`a > 0` used in the oblate spheroidal radius
    :math:`\\sqrt{r^2 + a^2}`.

    Notes
    -----
    The constructor takes no arguments. The spin parameter ``a`` is created as
    a positive SymPy symbol.

    Attributes
    ----------
    a : sympy.Symbol
        Kerr spin parameter :math:`a`, declared positive.

    Examples
    --------
    >>> bl = BoyerLindquist()
    >>> bl.a.is_positive
    True
    """

    def __init__(self):
        super().__init__()
        self.a = sp.symbols("a", positive=True)

    def cartesian_conversion_rule(self):
        """
        Return the Cartesian components :math:`(x, y, z)` from Boyer–Lindquist
        coordinates.

        Uses :math:`x^i` from ``xs_p`` as :math:`(r, \\theta, \\phi)` and the
        instance attribute ``a`` for the spin parameter.

        Returns
        -------
        sympy.Matrix
            Column vector :math:`[x, y, z]^\\mathsf{T}` with entries

            .. math::

                x &= \\sqrt{(x^1)^2 + a^2}\\,\\sin x^2\\,\\cos x^3,\\\\
                y &= \\sqrt{(x^1)^2 + a^2}\\,\\sin x^2\\,\\sin x^3,\\\\
                z &= x^1 \\cos x^2.

        Examples
        --------
        >>> import sympy as sp
        >>> from relatipy.symbolic.coordinates.boyer_lindquist import BoyerLindquist
        >>> bl = BoyerLindquist()
        >>> out = bl.cartesian_conversion_rule()
        >>> out[2]  # z = r cos(theta)
        x^1(t)*cos(x^2(t))
        """
        xa = sp.sqrt(xs_p[0] ** 2 + self.a**2)
        sin_norm = xa * sp.sin(xs_p[1])
        x = sin_norm * sp.cos(xs_p[2])
        y = sin_norm * sp.sin(xs_p[2])
        z = xs_p[0] * sp.cos(xs_p[1])

        return sp.Matrix([x, y, z])
