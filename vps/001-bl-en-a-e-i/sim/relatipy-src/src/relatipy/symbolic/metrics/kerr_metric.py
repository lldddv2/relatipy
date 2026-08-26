"""
Symbolic Kerr metric in Boyer–Lindquist–type coordinates.

This module defines :class:`Kerr`, which builds an ``einsteinpy`` symbolic
metric tensor and Christoffel symbols for the Kerr line element using
coordinates :math:`(x^0, x^1, x^2, x^3)` mapped to :math:`(t, r, \\theta, \\phi)`.
The mass and angular momentum appear only indirectly through the parameters
:math:`R_s` (Schwarzschild radius) and :math:`a` (spin length), which are
independent SymPy symbols in :meth:`Kerr._compute_metric`.

The metric components follow the block structure

.. math::

    g_{\\mu\\nu} = \\begin{pmatrix}
        A & 0 & 0 & E/2 \\\\
        0 & B & 0 & 0 \\\\
        0 & 0 & C & 0 \\\\
        E/2 & 0 & 0 & D
    \\end{pmatrix},

with :math:`\\Sigma`, :math:`\\Delta`, and :math:`A,\\ldots,E` constructed as in
:meth:`Kerr._compute_metric`.

Notes
-----
Christoffel symbols are obtained via
``einsteinpy.symbolic.ChristoffelSymbols.from_metric`` and each component is
expanded and simplified with :func:`_simplify_expr`.

References
----------
Bardeen, J. M., Press, W. H., & Teukolsky, S. A. (1972). Rotating black holes:
locally nonrotating frames, energy extraction, and scalar synchrotron radiation.
*The Astrophysical Journal*, 178, 347–370.

Examples
--------
>>> from relatipy.symbolic.metrics import Kerr
>>> k = Kerr()
>>> g = k.metric()
>>> g.tensor().shape
(4, 4)
"""

from itertools import product

import einsteinpy.symbolic as es
import numpy as np
import sympy as sp


def _simplify_expr(expr):
    """
    Expand and simplify a SymPy expression.

    Parameters
    ----------
    expr : sympy.Expr
        Expression to simplify.

    Returns
    -------
    sympy.Expr
        Expanded and simplified expression.

    Examples
    --------
    >>> import sympy as sp
    >>> x = sp.Symbol("x")
    >>> _simplify_expr((x + 1)**2 - 1)
    x*(x + 2)
    """
    return expr.expand().simplify()


# Coordinates x^mu ~ (t, r, theta, phi)
x0, x1, x2, x3 = sp.symbols("x^0 x^1 x^2 x^3")
xs = [x0, x1, x2, x3]


class Kerr:
    """
    Lazy symbolic Kerr geometry (metric and Christoffel symbols).

    The metric and Christoffel symbols are computed once on first access and
    cached on the instance.

    Parameters
    ----------
    None
        This class takes no constructor arguments.

    Attributes
    ----------
    metic_has_been_computed : bool
        Whether the metric tensor has been built and cached.
    christoffel_has_been_computed : bool
        Whether the Christoffel symbols have been built and cached.
    _metric : einsteinpy.symbolic.metric.MetricTensor or None
        Cached metric tensor, or ``None`` before first call to :meth:`metric`.
    _christoffel_symbols : einsteinpy.symbolic.christoffel.ChristoffelSymbols or None
        Cached Christoffel symbols, or ``None`` before first call to
        :meth:`christoffel_symbols`.

    Examples
    --------
    >>> from relatipy.symbolic.metrics import Kerr
    >>> k = Kerr()
    >>> k.metic_has_been_computed
    False
    >>> _ = k.metric()
    >>> k.metic_has_been_computed
    True
    """

    def __init__(self):
        self.metic_has_been_computed = False
        self.christoffel_has_been_computed = False
        self._metric = None
        self._christoffel_symbols = None

    def metric(self):
        """
        Return the symbolic Kerr metric tensor.

        The first call computes and caches the tensor; later calls return the
        same object.

        Returns
        -------
        einsteinpy.symbolic.metric.MetricTensor
            Metric tensor in coordinates ``xs``.

        Examples
        --------
        >>> from relatipy.symbolic.metrics import Kerr
        >>> k = Kerr()
        >>> g = k.metric()
        >>> g.syms
        [x^0, x^1, x^2, x^3]
        """
        if not self.metic_has_been_computed:
            self._metric = self._compute_metric()
            self.metic_has_been_computed = True

        return self._metric

    def christoffel_symbols(self):
        """
        Return the symbolic Christoffel symbols for this metric.

        The first call computes and caches the symbols; later calls return the
        same object.

        Returns
        -------
        einsteinpy.symbolic.christoffel.ChristoffelSymbols
            Christoffel symbols :math:`\\Gamma^\\mu_{\\nu\\sigma}` associated with
            :meth:`metric`.

        Examples
        --------
        >>> from relatipy.symbolic.metrics import Kerr
        >>> k = Kerr()
        >>> ch = k.christoffel_symbols()
        >>> ch.tensor().shape
        (4, 4, 4)
        """
        if not self.christoffel_has_been_computed:
            self._christoffel_symbols = self._compute_christoffel_symbols()
            self.christoffel_has_been_computed = True
        return self._christoffel_symbols

    def _compute_metric(self):
        """
        Build the Kerr metric as a ``MetricTensor`` from symbolic components.

        Uses independent SymPy symbols ``R_s`` and ``a`` for the Schwarzschild
        radius and spin parameter, and coordinates ``xs[1]`` (:math:`r`),
        ``xs[2]`` (:math:`\\theta`) in the Boyer–Lindquist-type definitions of
        :math:`\\Sigma`, :math:`\\Delta`, and the diagonal and off-diagonal
        entries ``A``–``E``.

        Returns
        -------
        einsteinpy.symbolic.metric.MetricTensor
            Kerr metric in coordinates ``xs``.

        Examples
        --------
        >>> from relatipy.symbolic.metrics import Kerr
        >>> k = Kerr()
        >>> g = k._compute_metric()
        >>> g.tensor().shape
        (4, 4)
        """
        # constant sp.symbols
        G, c = sp.symbols("G c")

        # parameters
        M = sp.symbols("M")
        J = sp.symbols("J")

        # Derived parameters
        R_s = 2 * G * M / c**2
        a = J / (M * c)
        R_s, a = sp.symbols("R_s a")

        # Metric
        Sigma = xs[1] ** 2 * (1 + a**2 / xs[1] ** 2 * sp.cos(xs[2]) ** 2)
        Delta = xs[1] ** 2 * (1 - (R_s * xs[1] + a**2) / xs[1] ** 2)

        A = 1 - R_s * xs[1] / Sigma
        B = -Sigma / Delta
        C = -Sigma
        D = -(xs[1] ** 2 + a**2 + R_s * xs[1] * a**2 / Sigma * sp.sin(xs[2]) ** 2) * sp.sin(
            xs[2]
        ) ** 2
        E = 2 * R_s * xs[1] * a / Sigma * sp.sin(xs[2]) ** 2

        metric = np.diag([A, B, C, D])
        metric[0, 3] = metric[3, 0] = E / 2

        g = es.MetricTensor(metric, xs)

        return g

    def _compute_christoffel_symbols(self):
        """
        Compute Christoffel symbols from the metric and simplify each component.

        Returns
        -------
        einsteinpy.symbolic.christoffel.ChristoffelSymbols
            Christoffel symbols with each component expanded and simplified.

        Examples
        --------
        >>> from relatipy.symbolic.metrics import Kerr
        >>> k = Kerr()
        >>> ch = k._compute_christoffel_symbols()
        >>> ch.syms
        [x^0, x^1, x^2, x^3]
        """
        ch = es.ChristoffelSymbols.from_metric(self.metric())

        christoffel = np.zeros((4, 4, 4), dtype=object)

        for mu, nu, sigma in product(range(4), repeat=3):
            expr = _simplify_expr(ch.tensor()[mu, nu, sigma])
            christoffel[mu, nu, sigma] = expr

        christoffel = es.ChristoffelSymbols(christoffel, xs)

        return christoffel
