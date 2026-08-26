"""
Schwarzschild vacuum metric in spherical coordinates.

This module provides :class:`Schwarzschild`, a concrete implementation of
:class:`~relatipy.numeric.metrics.base.BaseMetric` for the static, spherically
symmetric black hole solution. Coordinates are spherical
:math:`(t, r, \\theta, \\phi)` in the same dimensionless geometric convention as
the rest of the numeric stack (lengths scaled by :math:`GM/c^2`, times by
:math:`GM/c^3` when interfacing through the base class).

The line element implemented by :meth:`Schwarzschild._metric_dimensionless` is

.. math::

    ds^2 = \\left(1 - \\frac{R_s}{r}\\right)\\, dt^2
        - \\left(1 - \\frac{R_s}{r}\\right)^{-1} dr^2
        - r^2\\, d\\theta^2
        - r^2 \\sin^2\\theta\\, d\\phi^2,

where :math:`R_s = 2GM/c^2` is stored as :attr:`~relatipy.numeric.metrics.base.BaseMetric.R_s`.

Notes
-----
The innermost stable circular orbit (ISCO) for massive test particles in the
equatorial plane is at :math:`r = 3 R_s` in these units; see :attr:`Schwarzschild.isco`.

References
----------
Schwarzschild, K. (1916). Über das Gravitationsfeld eines Massenpunktes nach der
Einsteinschen Theorie. Sitzungsberichte der Königlich Preussischen Akademie der
Wissenschaften, 189–196.

Examples
--------
>>> import numpy
>>> from relatipy.numeric.metrics import Schwarzschild
>>> bh = Schwarzschild(1.989e30)
>>> xs = numpy.array([0.0, 1.0e10, numpy.pi / 2, 0.0])
>>> g = bh.metric(xs, dimensionless=True)
>>> g.shape
(4, 4)
"""
import numpy

from .base import BaseMetric


class Schwarzschild(BaseMetric):
    """
    Schwarzschild metric in spherical coordinates.

    The central mass sets the Schwarzschild radius
    :math:`R_s = 2GM/c^2` on the base class. The metric is diagonal in
    :math:`(t, r, \\theta, \\phi)`.

    Parameters
    ----------
    mass : float
        Mass of the gravitating body in kilograms (validated by the base class).

    Attributes
    ----------
    isco : float
        Radius of the innermost stable circular orbit, :math:`3 R_s`, in the same
        length units as ``R_s``.

    Examples
    --------
    >>> from relatipy.numeric.metrics import Schwarzschild
    >>> bh = Schwarzschild(5.972e24)
    >>> bh.isco == 3 * bh.R_s
    True
    """

    def __init__(self, mass):
        """
        Construct a Schwarzschild metric for the given mass.

        Parameters
        ----------
        mass : float
            Mass in kilograms.

        Examples
        --------
        >>> from relatipy.numeric.metrics import Schwarzschild
        >>> Schwarzschild(1.989e30).R_s > 0
        True
        """
        super().__init__(mass, valid_coordinate="Spherical")
        self.isco = self._get_isco()

    def _get_isco(self):
        """
        Innermost stable circular orbit radius for equatorial geodesics.

        For the Schwarzschild metric, massive circular orbits in the equatorial
        plane are stable for :math:`r > 6GM/c^2` and unstable below; the ISCO
        where the radial epicyclic frequency vanishes lies at :math:`r = 3 R_s`.

        Returns
        -------
        float
            :math:`3 R_s` in geometric length units.

        Examples
        --------
        >>> from relatipy.numeric.metrics import Schwarzschild
        >>> bh = Schwarzschild(1.989e30)
        >>> bh._get_isco() == 3 * bh.R_s
        True
        """
        return 3 * self.R_s

    def _metric_dimensionless(self, xs):
        """
        Schwarzschild metric tensor in dimensionless spherical coordinates.

        Returns the diagonal components
        :math:`g_{tt} = 1 - R_s/r`,
        :math:`g_{rr} = -(1 - R_s/r)^{-1}`,
        :math:`g_{\\theta\\theta} = -r^2`,
        :math:`g_{\\phi\\phi} = -r^2 \\sin^2\\theta`.

        Parameters
        ----------
        xs : array_like
            Spherical coordinates. A single event has shape ``(4,)`` with
            ``[t, r, theta, phi]``. Multiple events have shape ``(N, 4)`` with
            one row per event.

        Returns
        -------
        ndarray
            If ``xs`` is one-dimensional, a ``(4, 4)`` diagonal matrix. If ``xs``
            has shape ``(N, 4)``, an ``(N, 4, 4)`` array whose ``[i]`` slice is
            the metric at the ``i``-th row of ``xs``.

        Examples
        --------
        >>> import numpy
        >>> from relatipy.numeric.metrics import Schwarzschild
        >>> bh = Schwarzschild(1.989e30)
        >>> g = bh._metric_dimensionless([0.0, 1.0e10, numpy.pi / 2, 0.0])
        >>> g.shape
        (4, 4)
        """
        xs = numpy.asarray(xs, dtype=float)

        # Single point: shape (4,)
        if xs.ndim == 1:
            A = 1 - self.R_s / xs[1]
            B = -1 / A
            C = -xs[1] ** 2
            D = -xs[1] ** 2 * numpy.sin(xs[2]) ** 2
            return numpy.diag([A, B, C, D])

        # Multiple points: shape (N, 4)
        r = xs[:, 1]
        theta = xs[:, 2]

        A = 1 - self.R_s / r
        B = -1 / A
        C = -r ** 2
        D = -r ** 2 * numpy.sin(theta) ** 2

        N = len(xs)
        metrics = numpy.zeros((N, 4, 4))
        metrics[:, 0, 0] = A
        metrics[:, 1, 1] = B
        metrics[:, 2, 2] = C
        metrics[:, 3, 3] = D

        return metrics

    def _get_christoffel_symbols(self, xs):
        """
        Christoffel symbols :math:`\\Gamma^k_{\\mu\\nu}` for the Schwarzschild metric.

        Indices follow ``(t, r, \\theta, \\phi)`` :math:`\\equiv (0,1,2,3)`.
        The array layout is ``Gamma[k, mu, nu]`` with symmetry in the lower pair
        (only lower-triangle values are guaranteed to be filled consistently with
        the analytic expressions used here).

        Parameters
        ----------
        xs : ndarray
            Coordinates with shape ``(4,)`` as ``[t, r, theta, phi]``, or
            ``(N, 4)`` for a batch of points.

        Returns
        -------
        ndarray
            Shape ``(4, 4, 4)`` for a single point, or ``(N, 4, 4, 4)`` when
            ``xs`` is two-dimensional.

        Notes
        -----
        Components involve :math:`r`, :math:`\\theta`, and :math:`R_s`. Polar
        terms use :math:`\\sin\\theta` and :math:`\\cos\\theta`; callers should
        avoid :math:`\\theta` values where the spherical chart is singular unless
        the corresponding Christoffel components are not needed.

        Examples
        --------
        >>> import numpy
        >>> from relatipy.numeric.metrics import Schwarzschild
        >>> bh = Schwarzschild(1.989e30)
        >>> x = numpy.array([0.0, 1.0e10, numpy.pi / 2, 0.0])
        >>> G = bh._get_christoffel_symbols(x)
        >>> G.shape
        (4, 4, 4)
        """
        xs = numpy.asarray(xs, dtype=float)

        if xs.ndim == 2:
            N = len(xs)
            Gammas = numpy.zeros((N, 4, 4, 4))
            for i, x in enumerate(xs):
                Gammas[i] = self._get_christoffel_symbols(x)
            return Gammas

        r_s = self.R_s
        r = xs[1]
        theta = xs[2]
        cos = numpy.cos
        sin = numpy.sin
        Gamma = numpy.array(
            [
                [
                    [0, r_s / (2 * r * (r - r_s)), 0, 0],
                    [r_s / (2 * r * (r - r_s)), 0, 0, 0],
                    [0, 0, 0, 0],
                    [0, 0, 0, 0],
                ],
                [
                    [-r_s * (-r + r_s) / (2 * r**3), 0, 0, 0],
                    [0, r_s * (-r + r_s) / (2 * r**3 * (1 - r_s / r) ** 2), 0, 0],
                    [0, 0, -r + r_s, 0],
                    [0, 0, 0, (-r + r_s) * sin(theta) ** 2],
                ],
                [
                    [0, 0, 0, 0],
                    [0, 0, 1 / r, 0],
                    [0, 1 / r, 0, 0],
                    [0, 0, 0, -sin(theta) * cos(theta)],
                ],
                [
                    [0, 0, 0, 0],
                    [0, 0, 0, 1 / r],
                    [0, 0, 0, cos(theta) / sin(theta)],
                    [0, 1 / r, cos(theta) / sin(theta), 0],
                ],
            ]
        )

        return Gamma
