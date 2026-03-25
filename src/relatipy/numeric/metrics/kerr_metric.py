"""
Kerr black-hole spacetime in Boyer–Lindquist coordinates (numeric).

This module provides :class:`Kerr`, a concrete metric used by the numeric stack
for geodesic integration and related calculations. Coordinates are
:math:`(t, r, \\theta, \\phi)` with the same dimensionless conventions as
:class:`~relatipy.numeric.metrics.base.BaseMetric` (geometric units with
:math:`G=c=1` in internal formulas, and masses handled as in the validator).

The line element in Boyer–Lindquist form is

.. math::

    ds^2 = -\\left(1 - \\frac{2Mr}{\\Sigma}\\right)\\, dt^2
    - \\frac{4Mra\\sin^2\\theta}{\\Sigma}\\, dt\\, d\\phi
    + \\frac{\\Sigma}{\\Delta}\\, dr^2 + \\Sigma\\, d\\theta^2
    + \\left(r^2 + a^2 + \\frac{2M r a^2 \\sin^2\\theta}{\\Sigma}\\right)
      \\sin^2\\theta\\, d\\phi^2,

with :math:`\\Sigma = r^2 + a^2 \\cos^2\\theta` and
:math:`\\Delta = r^2 - 2Mr + a^2`.

Notes
-----
Christoffel components in :meth:`Kerr._get_christoffel_symbols` were generated
from symbolic expressions (see ``utils/symbolic_to_numeric.py``) and then
manually optimized.

References
----------
Bardeen, J. M., Press, W. H., & Teukolsky, S. A. (1972). Rotating black holes:
locally nonrotating frames, energy extraction, and scalar synchrotron radiation.
*The Astrophysical Journal*, 178, 347–370.

Examples
--------
>>> from relatipy.numeric.metrics import Kerr
>>> bh = Kerr(1.0, 0.5)
>>> bool(bh.isco_prograde < bh.isco_retrograde)
True
"""

import numpy
from numpy import sin, cos, tan

from .base import BaseMetric


class Kerr(BaseMetric):
    """
    Kerr metric in Boyer–Lindquist coordinates.

    The constructor stores the dimensionless spin parameter and precomputes
    innermost stable circular orbit (ISCO) radii for prograde and retrograde
    equatorial motion.

    Parameters
    ----------
    mass : float or astropy.units.Quantity
        Black-hole mass, in the same convention as :class:`BaseMetric`
        (validated by ``validator.validate_scalar``; bare floats are treated as
        geometric mass in units of the reference solar mass).
    a : float
        Dimensionless Kerr spin parameter :math:`a/M` (often denoted
        :math:`a_*`), with :math:`|a| \\leq 1` for sub-extremal black holes.
        It is converted internally to the length-like spin :math:`a` used in
        Boyer–Lindquist formulas via ``a * R_s / 2``.

    Attributes
    ----------
    a : float
        Spin parameter in length units consistent with :attr:`BaseMetric.R_s`
        and coordinate :math:`r` (internal Boyer–Lindquist :math:`a`).
    isco_prograde : float
        ISCO radius for co-rotating equatorial orbits (geometric units).
    isco_retrograde : float
        ISCO radius for counter-rotating equatorial orbits (geometric units).

    Examples
    --------
    >>> from relatipy.numeric.metrics import Kerr
    >>> k = Kerr(1.0, 0.9)
    >>> bool(k.a > 0)
    True
    """

    def __init__(self, mass, a):
        """
        Set internal spin, register Boyer–Lindquist coordinates, and compute ISCO radii.

        Parameters
        ----------
        mass : float or astropy.units.Quantity
            Black-hole mass (same convention as :class:`BaseMetric`).
        a : float
            Dimensionless Kerr spin :math:`a/M` used to set :attr:`a`.

        Examples
        --------
        >>> from relatipy.numeric.metrics import Kerr
        >>> float(Kerr(1.0, 0.0).isco_prograde)
        6.0
        """
        super().__init__(mass, valid_coordinate="BoyerLindquist", kwargs={"a": a})
        self.a = a * self.R_s / 2
        self.isco_prograde = self._get_isco(prograde=True)
        self.isco_retrograde = self._get_isco(prograde=False)

    def _get_isco(self, prograde: bool = True) -> float:
        """
        Innermost stable circular orbit (ISCO) radius for Kerr in the equatorial plane.

        Uses the closed-form expression from Bardeen, Press & Teukolsky (1972)
        in terms of the dimensionless spin :math:`\\hat{a} = a_{\\mathrm{int}}/M`,
        where :math:`a_{\\mathrm{int}}` is the internal spin :attr:`a` and
        :math:`M` is :attr:`~BaseMetric.mass`.

        Parameters
        ----------
        prograde : bool, optional
            If True, return the ISCO for co-rotating (prograde) orbits; if
            False, for counter-rotating (retrograde) orbits.

        Returns
        -------
        float
            ISCO radius in the same geometric length units as coordinates
            (consistent with :attr:`BaseMetric.R_s` and :attr:`a`).

        References
        ----------
        Bardeen, J. M., Press, W. H., & Teukolsky, S. A. (1972). Rotating black
        holes: locally nonrotating frames, energy extraction, and scalar
        synchrotron radiation. *The Astrophysical Journal*, 178, 347–370.

        Examples
        --------
        >>> from relatipy.numeric.metrics import Kerr
        >>> k = Kerr(1.0, 0.0)
        >>> bool(abs(k._get_isco(True) - 6.0) < 1e-10)
        True
        """
        a_hat = self.a / self.mass  # spin adimensional ∈ [0, 1]

        Z1 = 1 + (1 - a_hat**2)**(1/3) * (
            (1 + a_hat)**(1/3) + (1 - a_hat)**(1/3)
        )
        Z2 = numpy.sqrt(3 * a_hat**2 + Z1**2)

        sign = -1 if prograde else +1
        return self.mass * (3 + Z2 + sign * numpy.sqrt((3 - Z1) * (3 + Z1 + 2 * Z2)))

    def _metric_dimensionless(self, xs):
        """
        Kerr metric tensor :math:`g_{\\mu\\nu}` in Boyer–Lindquist coordinates.

        Coordinates are :math:`x^\\mu = (t, r, \\theta, \\phi)`. The implementation
        matches the standard Kerr form with :math:`\\Sigma` and :math:`\\Delta` as
        above, using :attr:`BaseMetric.R_s` for :math:`2M`.

        Parameters
        ----------
        xs : array_like
            Coordinates ``[t, r, theta, phi]``. A single point has shape ``(4,)``;
            several points have shape ``(N, 4)``.

        Returns
        -------
        numpy.ndarray
            If ``xs`` is one-dimensional, an array of shape ``(4, 4)``. If
            ``xs`` has shape ``(N, 4)``, an array of shape ``(N, 4, 4)``.

        Notes
        -----
        Off-diagonal components implement the :math:`g_{t\\phi} = g_{\\phi t}`
        coupling. The matrix is symmetric; entries are set explicitly only where
        non-zero.

        Examples
        --------
        >>> import numpy as np
        >>> from relatipy.numeric.metrics import Kerr
        >>> k = Kerr(1.0, 0.5)
        >>> g = k._metric_dimensionless(np.array([0.0, 10.0, 1.0, 0.0]))
        >>> g.shape
        (4, 4)
        """
        xs = numpy.asarray(xs, dtype=float)

        if xs.ndim == 1:
            x1   = xs[1]
            x1_2 = x1 * x1
            s_x2 = numpy.sin(xs[2]) ** 2
            c_x2 = numpy.cos(xs[2]) ** 2
            a    = self.a
            a_2  = a * a
            R_s  = self.R_s

            Sigma = x1_2 + a_2 * c_x2
            Delta = x1_2 - R_s * x1 + a_2

            A = 1 - R_s * x1 / Sigma
            B = -Sigma / Delta
            C = -Sigma
            D = -(x1_2 + a_2 + R_s * x1 * a_2 / Sigma * s_x2) * s_x2
            E = 2 * R_s * x1 * a / Sigma * s_x2

            metric = numpy.diag([A, B, C, D])
            metric[0, 3] = metric[3, 0] = E / 2
            return metric

        # Multiple points: shape (N, 4)
        r    = xs[:, 1]
        r2   = r * r
        s_x2 = numpy.sin(xs[:, 2]) ** 2
        c_x2 = numpy.cos(xs[:, 2]) ** 2
        a    = self.a
        a_2  = a * a
        R_s  = self.R_s

        Sigma = r2 + a_2 * c_x2
        Delta = r2 - R_s * r + a_2

        A = 1 - R_s * r / Sigma
        B = -Sigma / Delta
        C = -Sigma
        D = -(r2 + a_2 + R_s * r * a_2 / Sigma * s_x2) * s_x2
        E = 2 * R_s * r * a / Sigma * s_x2

        N = len(xs)
        metrics = numpy.zeros((N, 4, 4))
        metrics[:, 0, 0] = A
        metrics[:, 1, 1] = B
        metrics[:, 2, 2] = C
        metrics[:, 3, 3] = D
        metrics[:, 0, 3] = metrics[:, 3, 0] = E / 2

        return metrics

    # Generated by utils/symbolic_to_numeric.py
    def _get_christoffel_symbols(self, xs):
        """
        Christoffel symbols :math:`\\Gamma^\\lambda_{\\mu\\nu}` of the Kerr metric.

        Indices are ordered as ``Gamma[lambda, mu, nu]``. Only non-zero components
        are filled; symmetry :math:`\\Gamma^\\lambda_{\\mu\\nu} =
        \\Gamma^\\lambda_{\\nu\\mu}` is enforced by copying mirrored entries.

        Parameters
        ----------
        xs : array_like of shape (4,) or (N, 4)
            Boyer–Lindquist coordinates ``(t, r, theta, phi)``. For ``(N, 4)``,
            Christoffels are computed pointwise in a loop.

        Returns
        -------
        numpy.ndarray
            Array of shape ``(4, 4, 4)`` for a single point, or ``(N, 4, 4, 4)``
            for ``N`` points.

        Notes
        -----
        Expressions were generated symbolically and then partially folded for
        performance. The polar angle :math:`\\theta` must avoid the coordinate
        singularities of :math:`\\tan\\theta` (e.g. :math:`\\theta = 0, \\pi`).

        Examples
        --------
        >>> import numpy as np
        >>> from relatipy.numeric.metrics import Kerr
        >>> k = Kerr(1.0, 0.5)
        >>> x = np.array([0.0, 10.0, np.pi / 4, 0.0])
        >>> G = k._get_christoffel_symbols(x)
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

        x0, x1, x2, x3 = xs

        ############### Auxiliary variables ##############
        a = self.a
        R_s = self.R_s

        cos_I2_x2I_ = cos(2 * x2)
        cos_I4_x2I_ = cos(4 * x2)
        cos_Ix2I_ = cos(x2)
        sin_I2_x2I_ = sin(2 * x2)
        sin_Ix2I_ = sin(x2)
        tan_Ix2I_ = tan(x2)

        R_s_pow2 = R_s**2
        _1_over_tan_Ix2I_ = 1 / tan_Ix2I_
        a_pow2 = a**2
        a_pow4 = a**4
        a_pow6 = a**6
        cos_Ix2I__pow2 = cos_Ix2I_**2
        cos_Ix2I__pow4 = cos_Ix2I_**4
        sin_Ix2I__pow2 = sin_Ix2I_**2
        sin_Ix2I__pow4 = sin_Ix2I_**4
        x1_pow2 = x1**2
        x1_pow3 = x1**3
        x1_pow4 = x1**4
        x1_pow5 = x1**5
        x1_pow6 = x1**6

        # Modified manually for optimization
        a_pow2_cos_Ix2I__pow2 = a_pow2 * cos_Ix2I__pow2
        R_s_a_pow2_cos_Ix2I__pow2_x1 = R_s * a_pow2_cos_Ix2I__pow2 * x1
        R_s_x1_pow3 = R_s * x1_pow3
        _2_x1_pow2 = 2 * x1_pow2
        _R_s_x1 = -R_s * x1
        _a_pow2_cos_Ix2I__pow2_x1_pow2 = -a_pow2_cos_Ix2I__pow2 * x1_pow2
        a_pow2_cos_I2_x2I_ = a_pow2 * cos_I2_x2I_
        _a_pow2_sin_I2_x2I__over__Ia_pow2_cos_I2_x2I___plus__a_pow2__plus__2_x1_pow2I_ = (
            -a_pow2 * sin_I2_x2I_ / (a_pow2_cos_I2_x2I_ + a_pow2 + _2_x1_pow2)
        )
        a_pow2_x1_pow2 = a_pow2 * x1_pow2
        _a_pow2_x1_pow2 = -a_pow2_x1_pow2
        _cos_I4_x2I_ = -cos_I4_x2I_
        _x1_pow4 = -x1_pow4
        a_pow4_cos_Ix2I__pow2 = a_pow4 * cos_Ix2I__pow2

        ##################################################

        ############### Christoffel symbols ##############
        Gamma = numpy.zeros((4, 4, 4))

        # Only non-zero components and symmetry properties

        # mu = 0
        Gamma[0, 0, 1] = (
            R_s
            * (
                a_pow2_x1_pow2 * sin_Ix2I__pow2
                + a_pow4 * sin_Ix2I__pow2
                - a_pow4
                + x1_pow4
            )
            / (
                2
                * (a_pow2_cos_Ix2I__pow2 + x1_pow2) ** 2
                * (_R_s_x1 + a_pow2 + x1_pow2)
            )
        )
        Gamma[0, 0, 2] = (
            2
            * R_s
            * _a_pow2_sin_I2_x2I__over__Ia_pow2_cos_I2_x2I___plus__a_pow2__plus__2_x1_pow2I_
            * x1
            / (_2_x1_pow2 + a_pow2 + a_pow2_cos_I2_x2I_)
        )
        Gamma[0, 1, 3] = (
            R_s
            * a
            * sin_Ix2I__pow2
            * (
                _a_pow2_cos_Ix2I__pow2_x1_pow2
                + _a_pow2_x1_pow2
                + 3 * _x1_pow4
                + a_pow4_cos_Ix2I__pow2
            )
            / (
                2
                * (a_pow2_cos_Ix2I__pow2 + x1_pow2) ** 2
                * (_R_s_x1 + a_pow2 + x1_pow2)
            )
        )
        Gamma[0, 2, 3] = (
            R_s
            * a**3
            * cos_Ix2I_
            * sin_Ix2I_**3
            * x1
            / (a_pow2_cos_Ix2I__pow2 + x1_pow2) ** 2
        )

        Gamma[0, 1, 0] = Gamma[0, 0, 1]
        Gamma[0, 2, 0] = Gamma[0, 0, 2]
        Gamma[0, 3, 1] = Gamma[0, 1, 3]
        Gamma[0, 3, 2] = Gamma[0, 2, 3]

        # mu = 1
        Gamma[1, 0, 0] = (
            R_s
            * (
                R_s_a_pow2_cos_Ix2I__pow2_x1
                - R_s_x1_pow3
                + _a_pow2_cos_Ix2I__pow2_x1_pow2
                + _a_pow2_x1_pow2
                + a_pow4_cos_Ix2I__pow2
                + x1_pow4
            )
            / (2 * (a_pow2_cos_Ix2I__pow2 + x1_pow2) ** 3)
        )
        Gamma[1, 0, 3] = (
            R_s
            * a
            * sin_Ix2I__pow2
            * (
                -R_s_a_pow2_cos_Ix2I__pow2_x1
                + R_s_x1_pow3
                + _x1_pow4
                + a_pow2_cos_Ix2I__pow2 * x1_pow2
                + a_pow2_x1_pow2
                - a_pow4_cos_Ix2I__pow2
            )
            / (2 * (a_pow2_cos_Ix2I__pow2 + x1_pow2) ** 3)
        )
        Gamma[1, 1, 1] = (
            -R_s * a_pow2_cos_Ix2I__pow2 / 2
            + R_s * x1_pow2 / 2
            + a_pow2 * x1
            + a_pow2_cos_Ix2I__pow2 * x1
        ) / (
            R_s_a_pow2_cos_Ix2I__pow2_x1
            + R_s_x1_pow3
            + _a_pow2_cos_Ix2I__pow2_x1_pow2
            + _x1_pow4
            + a_pow2_x1_pow2
            + a_pow4_cos_Ix2I__pow2
        )
        Gamma[1, 1, 2] = (
            _a_pow2_sin_I2_x2I__over__Ia_pow2_cos_I2_x2I___plus__a_pow2__plus__2_x1_pow2I_
        )
        Gamma[1, 2, 2] = (
            x1 * (R_s * x1 + a_pow2 - x1_pow2) / (a_pow2_cos_Ix2I__pow2 + x1_pow2)
        )
        Gamma[1, 3, 3] = (
            sin_Ix2I__pow2
            * (
                R_s * a_pow2 * sin_Ix2I__pow2 * x1_pow4 / 2
                + 2 * R_s * a_pow2_cos_Ix2I__pow2 * x1_pow4
                + R_s * a_pow4 * cos_Ix2I__pow4 * x1_pow2
                - R_s * a_pow4 * sin_Ix2I__pow2 * x1_pow2 / 2
                - R_s * a_pow4 * x1_pow2 * (_cos_I4_x2I_ + 1) / 16
                + R_s * a_pow6 * (_cos_I4_x2I_ + 1) / 16
                + R_s * x1_pow6
                - R_s_pow2 * a_pow2 * sin_Ix2I__pow2 * x1_pow3 / 2
                + R_s_pow2 * a_pow4 * x1 * (_cos_I4_x2I_ + 1) / 16
                + a_pow2 * x1_pow5
                - 2 * a_pow2_cos_Ix2I__pow2 * x1_pow5
                - a_pow4 * cos_Ix2I__pow4 * x1_pow3
                + 2 * a_pow4_cos_Ix2I__pow2 * x1_pow3
                + a_pow6 * cos_Ix2I__pow4 * x1
                - x1**7
            )
            / (a_pow2_cos_Ix2I__pow2 + x1_pow2) ** 3
        )

        Gamma[1, 2, 1] = Gamma[1, 1, 2]
        Gamma[1, 3, 0] = Gamma[1, 0, 3]

        # mu = 2
        Gamma[2, 0, 0] = (
            4
            * R_s
            * _a_pow2_sin_I2_x2I__over__Ia_pow2_cos_I2_x2I___plus__a_pow2__plus__2_x1_pow2I_
            * x1
            / (_2_x1_pow2 + a_pow2 + a_pow2_cos_I2_x2I_) ** 2
        )
        Gamma[2, 0, 3] = (
            4
            * R_s
            * a
            * sin_I2_x2I_
            * x1
            * (a_pow2 + x1_pow2)
            / (_2_x1_pow2 + a_pow2 + a_pow2_cos_I2_x2I_) ** 3
        )
        Gamma[2, 1, 1] = (
            -a_pow2
            * cos_Ix2I_
            * sin_Ix2I_
            / (
                R_s_a_pow2_cos_Ix2I__pow2_x1
                + R_s_x1_pow3
                + _a_pow2_cos_Ix2I__pow2_x1_pow2
                - _a_pow2_x1_pow2
                + _x1_pow4
                + a_pow4_cos_Ix2I__pow2
            )
        )
        Gamma[2, 1, 2] = x1 / (a_pow2_cos_Ix2I__pow2 + x1_pow2)
        Gamma[2, 2, 2] = (
            _a_pow2_sin_I2_x2I__over__Ia_pow2_cos_I2_x2I___plus__a_pow2__plus__2_x1_pow2I_
        )
        Gamma[2, 3, 3] = (
            cos_Ix2I_
            * sin_Ix2I_
            * (
                -2 * R_s_x1_pow3 * a_pow2 * sin_Ix2I__pow2
                - _2_x1_pow2 * a_pow4_cos_Ix2I__pow2
                + _R_s_x1 * a_pow4 * sin_Ix2I__pow4
                + _R_s_x1 * a_pow4 * (_cos_I4_x2I_ + 1) / 4
                + _x1_pow4 * a_pow2
                + 2 * _x1_pow4 * a_pow2_cos_Ix2I__pow2
                - a_pow4 * cos_Ix2I__pow4 * x1_pow2
                - a_pow6 * cos_Ix2I__pow4
                - x1_pow6
            )
            / (a_pow2_cos_Ix2I__pow2 + x1_pow2) ** 3
        )

        Gamma[2, 2, 1] = Gamma[2, 1, 2]
        Gamma[2, 3, 0] = Gamma[2, 0, 3]

        # mu = 3
        Gamma[3, 0, 1] = (
            R_s
            * a
            * (-a_pow2_cos_Ix2I__pow2 + x1_pow2)
            / (
                2
                * (a_pow2_cos_Ix2I__pow2 + x1_pow2) ** 2
                * (_R_s_x1 + a_pow2 + x1_pow2)
            )
        )
        Gamma[3, 0, 2] = (
            _1_over_tan_Ix2I_ * _R_s_x1 * a / (a_pow2_cos_Ix2I__pow2 + x1_pow2) ** 2
        )
        Gamma[3, 1, 3] = (
            R_s * _a_pow2_cos_Ix2I__pow2_x1_pow2
            + R_s * _a_pow2_x1_pow2 * sin_Ix2I__pow2 / 2
            + R_s * _x1_pow4
            + R_s * a_pow4 * (_cos_I4_x2I_ + 1) / 16
            + 2 * a_pow2_cos_Ix2I__pow2 * x1_pow3
            + a_pow4 * cos_Ix2I__pow4 * x1
            + x1_pow5
        ) / ((a_pow2_cos_Ix2I__pow2 + x1_pow2) ** 2 * (_R_s_x1 + a_pow2 + x1_pow2))
        Gamma[3, 2, 3] = (
            _1_over_tan_Ix2I_
            * (
                R_s * a_pow2 * sin_Ix2I__pow2 * x1
                + _2_x1_pow2 * a_pow2
                + 2 * _a_pow2_x1_pow2 * sin_Ix2I__pow2
                - 2 * a_pow4 * sin_Ix2I__pow2
                + a_pow4 * sin_Ix2I__pow4
                + a_pow4
                + x1_pow4
            )
            / (_2_x1_pow2 * a_pow2_cos_Ix2I__pow2 + a_pow4 * cos_Ix2I__pow4 + x1_pow4)
        )

        Gamma[3, 1, 0] = Gamma[3, 0, 1]
        Gamma[3, 2, 0] = Gamma[3, 0, 2]
        Gamma[3, 3, 1] = Gamma[3, 1, 3]
        Gamma[3, 3, 2] = Gamma[3, 2, 3]

        ##################################################

        return Gamma
