"""
Boyer–Lindquist coordinates for the Kerr geometry (numeric).

This module provides :class:`BoyerLindquist`, a chart with coordinates
:math:`(t, r, \\theta, \\phi)` and spatial velocities
:math:`(v_r, v_\\theta, v_\\phi)` related to coordinate derivatives by the
usual Kerr line-element conventions (geometric units). The class maps between
these velocities and :math:`\\mathrm{d}x^i/\\mathrm{d}\\tau`, converts to
Cartesian Kerr–Schild–type spatial coordinates, and evaluates Kerr conserved
quantities :math:`L_z` and Carter’s :math:`Q` when a compatible metric and
state layout are supplied.

The auxiliary radius appearing in the spatial map is
:math:`\\sqrt{r^2 + a^2}` and

.. math::

    \\Sigma = r^2 + a^2 \\cos^2\\theta.

Notes
-----
Velocities :math:`v_i` used here are the “physical” Boyer–Lindquist components
linked to :math:`\\dot r`, :math:`\\dot\\theta`, :math:`\\dot\\phi` by the
factors :math:`\\sqrt{\\Sigma}`, :math:`\\sqrt{r^2+a^2}`, and
:math:`\\sin\\theta` as implemented in :meth:`BoyerLindquist._get_dxs_dt_from_vs`.

See Also
--------
relatipy.numeric.metrics.kerr_metric.Kerr :
    Kerr metric in the same Boyer–Lindquist chart.
relatipy.numeric.coordinates.base.CoordinateBase :
    Base class for position, velocities, and state vectors.

Examples
--------
>>> import numpy as np
>>> from relatipy.numeric.coordinates import BoyerLindquist
>>> bl = BoyerLindquist(
...     np.array([0.0, 10.0, np.pi / 2, 0.0]),
...     vels=np.zeros(3),
...     a=0.5,
...     from_dxs_dt=False,
... )
>>> bl.name_metric
'BoyerLindquist'
"""

import numpy
from numpy import sin, cos, sqrt, arctan2, arccos

from .base import CoordinateBase


class BoyerLindquist(CoordinateBase):
    """
    Boyer–Lindquist chart for Kerr spacetime with spin parameter ``a``.

    Coordinates are :math:`x^\\mu = (t, r, \\theta, \\phi)` with
    :math:`\\theta` the polar angle and :math:`\\phi` the azimuth. Spatial
    velocities follow the same normalization as in the Kerr line element
    (see module notes).

    Parameters
    ----------
    xs : array_like of shape (4,)
        Event :math:`(t, r, \\theta, \\phi)` in geometric units.
    vels : array_like of shape (3,) or None, optional
        Either :math:`(v_r, v_\\theta, v_\\phi)` if ``from_dxs_dt`` is False,
        or :math:`(\\dot r, \\dot\\theta, \\dot\\phi)` if True. If None,
        velocities and derivatives default to zeros.
    a : float or None, optional
        Dimensionless spin parameter :math:`a` (required; must not be None).
    from_dxs_dt : bool, optional
        If True, interpret ``vels`` as coordinate time derivatives of
        :math:`(r, \\theta, \\phi)`; otherwise interpret them as
        :math:`v_r, v_\\theta, v_\\phi`.

    Attributes
    ----------
    a : float
        Spin parameter stored from ``kwargs`` (same as constructor argument).

    Raises
    ------
    ValueError
        If ``a`` is None.

    Examples
    --------
    >>> import numpy as np
    >>> from relatipy.numeric.coordinates import BoyerLindquist
    >>> BoyerLindquist(
    ...     np.array([0.0, 5.0, np.pi / 2, 0.1]),
    ...     vels=np.array([0.0, 0.0, 0.0]),
    ...     a=0.9,
    ... ).a
    0.9
    """

    def __init__(self, xs, vels=None, a=None, from_dxs_dt=False):
        super().__init__(
            xs, vels=vels, from_dxs_dt=from_dxs_dt, system_name="BoyerLindquist", a=a
        )
        if self.a is None:
            raise ValueError(
                "The spin parameter 'a' must be provided for Boyer-Lindquist coordinates."
            )

    def _get_dxs_dt_from_vs(self):
        """
        Map physical Boyer–Lindquist velocities to coordinate derivatives.

        With :math:`x^1 = r`, :math:`x^2 = \\theta`, :math:`x^3 = \\phi`,

        .. math::

            \\dot r &= v_r \\frac{\\sqrt{r^2 + a^2}}{\\sqrt{\\Sigma}},\\\\
            \\dot\\theta &= \\frac{v_\\theta}{\\sqrt{\\Sigma}},\\\\
            \\dot\\phi &= \\frac{v_\\phi}{\\sqrt{r^2 + a^2}\\,\\sin\\theta},

        where :math:`\\Sigma = r^2 + a^2 \\cos^2\\theta`.

        Returns
        -------
        numpy.ndarray
            Array :math:`(\\dot r, \\dot\\theta, \\dot\\phi)` with the same
            shape convention as ``self.vs``.

        Examples
        --------
        >>> import numpy as np
        >>> from relatipy.numeric.coordinates import BoyerLindquist
        >>> bl = BoyerLindquist(
        ...     np.array([0.0, 2.0, np.pi / 2, 0.0]),
        ...     vels=np.array([0.0, 0.0, 0.0]),
        ...     a=0.5,
        ...     from_dxs_dt=False,
        ... )
        >>> np.allclose(bl.dxs_dt, 0.0)
        True
        """
        a = self.a
        xs = self.xs
        vs = self.vs

        sqrt_cos = sqrt(a**2 * cos(xs[2]) ** 2 + xs[1] ** 2)
        sqrt_a = sqrt(a**2 + xs[1] ** 2)
        sin_x2 = sin(xs[2])

        dx1_dt = vs[0] / sqrt_cos * sqrt_a
        dx2_dt = vs[1] / sqrt_cos
        dx3_dt = vs[2] / (sqrt_a * sin_x2)

        return numpy.array([dx1_dt, dx2_dt, dx3_dt])

    def _get_vs_from_dxs_dt(self):
        """
        Map coordinate derivatives to physical Boyer–Lindquist velocities.

        Inverse of :meth:`_get_dxs_dt_from_vs` with the same definitions of
        :math:`\\Sigma` and :math:`\\sqrt{r^2+a^2}`.

        Returns
        -------
        numpy.ndarray
            Array :math:`(v_r, v_\\theta, v_\\phi)` with the same shape as
            ``self.dxs_dt``.

        Examples
        --------
        >>> import numpy as np
        >>> from relatipy.numeric.coordinates import BoyerLindquist
        >>> bl = BoyerLindquist(
        ...     np.array([0.0, 2.0, np.pi / 2, 0.0]),
        ...     vels=np.array([0.0, 0.0, 0.0]),
        ...     a=0.5,
        ...     from_dxs_dt=True,
        ... )
        >>> np.allclose(bl.vs, 0.0)
        True
        """
        a = self.a
        xs = self.xs
        dxs_dt = self.dxs_dt

        sqrt_cos = sqrt(a**2 * cos(xs[2]) ** 2 + xs[1] ** 2)
        sqrt_a = sqrt(a**2 + xs[1] ** 2)
        sin_x2 = sin(xs[2])

        v1 = dxs_dt[0] * sqrt_cos / sqrt_a
        v2 = dxs_dt[1] * sqrt_cos
        v3 = dxs_dt[2] * sqrt_a * sin_x2

        return numpy.array([v1, v2, v3])

    @staticmethod
    def _convert_to_cartesian(xs, vs, a):
        """
        Map Boyer–Lindquist position and velocities to Cartesian coordinates.

        Positions use :math:`\\sqrt{r^2+a^2}\\sin\\theta` in the equatorial
        plane factor and :math:`z = r\\cos\\theta`. Velocities are obtained by
        differentiating the map using :math:`\\dot r`, :math:`\\dot\\theta`,
        :math:`\\dot\\phi` implied by ``vs`` via the same relations as in
        :meth:`_get_dxs_dt_from_vs`.

        Parameters
        ----------
        xs : array_like of shape (4,)
            :math:`(t, r, \\theta, \\phi)`.
        vs : array_like of shape (3,)
            :math:`(v_r, v_\\theta, v_\\phi)`.
        a : float
            Spin parameter :math:`a`.

        Returns
        -------
        xs_p : numpy.ndarray
            Cartesian :math:`(t, x, y, z)`.
        vs_p : numpy.ndarray
            Cartesian velocity components :math:`(v_x, v_y, v_z)` (spatial
            part aligned with the position map).

        Examples
        --------
        >>> import numpy as np
        >>> from relatipy.numeric.coordinates import BoyerLindquist
        >>> t = np.array([0.0, 3.0, np.pi / 2, 0.0])
        >>> v = np.zeros(3)
        >>> xp, vp = BoyerLindquist._convert_to_cartesian(t, v, a=0.0)
        >>> float(np.round(xp[1], 6))
        3.0
        """
        xs_p = numpy.zeros_like(xs)
        vs_p = numpy.zeros_like(vs)

        xs_p[0] = xs[0]

        r, theta, phi = xs[1], xs[2], xs[3]
        xa = sqrt(r**2 + a**2)
        sqrt_cos = sqrt(a**2 * cos(theta) ** 2 + r**2)
        sin_t = sin(theta)
        cos_t = cos(theta)
        sin_p = sin(phi)
        cos_p = cos(phi)

        xs_p[1] = xa * sin_t * cos_p
        xs_p[2] = xa * sin_t * sin_p
        xs_p[3] = r * cos_t

        dr_dt = vs[0] * xa / sqrt_cos
        dtheta_dt = vs[1] / sqrt_cos
        dphi_dt = vs[2] / (xa * sin_t)

        vs_p[0] = (
            (r * dr_dt * sin_t * cos_p / xa)
            + (xa * cos_t * cos_p * dtheta_dt)
            - (xa * sin_t * sin_p * dphi_dt)
        )
        vs_p[1] = (
            (r * dr_dt * sin_t * sin_p / xa)
            + (xa * cos_t * sin_p * dtheta_dt)
            + (xa * sin_t * cos_p * dphi_dt)
        )
        vs_p[2] = (dr_dt * cos_t) - (r * sin_t * dtheta_dt)

        return xs_p, vs_p

    @staticmethod
    def _convert_from_cartesian(xs_p, vs_p, a):
        """
        Build a :class:`BoyerLindquist` instance from Cartesian state.

        Recovers :math:`r`, :math:`\\theta`, :math:`\\phi` from
        :math:`(x,y,z)` and maps Cartesian velocities to Boyer–Lindquist
        :math:`v_r`, :math:`v_\\theta`, :math:`v_\\phi` using the inverse
        construction consistent with :meth:`_convert_to_cartesian`.

        Parameters
        ----------
        xs_p : array_like of shape (4,)
            Cartesian :math:`(t, x, y, z)`.
        vs_p : array_like of shape (3,)
            Cartesian velocities :math:`(v_x, v_y, v_z)`.
        a : float
            Spin parameter :math:`a`.

        Returns
        -------
        BoyerLindquist
            Coordinate object with ``from_dxs_dt=True`` (``vels`` are
            coordinate derivatives of :math:`r,\\theta,\\phi`).

        Examples
        --------
        >>> import numpy as np
        >>> from relatipy.numeric.coordinates import BoyerLindquist
        >>> cart_t = np.array([0.0, 3.0, 0.0, 0.0])
        >>> cart_v = np.zeros(3)
        >>> bl = BoyerLindquist._convert_from_cartesian(cart_t, cart_v, a=0.0)
        >>> float(bl.xs[1])
        3.0
        """
        xs = numpy.zeros_like(xs_p)
        vs = numpy.zeros_like(vs_p)

        xs[0] = xs_p[0]

        w = (xs_p[1] ** 2 + xs_p[2] ** 2 + xs_p[3] ** 2) - (a**2)
        xs[1] = sqrt(0.5 * (w + sqrt((w**2) + (4 * (a**2) * (xs_p[3] ** 2)))))
        xs[2] = arccos(xs_p[3] / xs[1])
        xs[3] = arctan2(xs_p[2], xs_p[1])

        w = (xs_p[1] ** 2 + xs_p[2] ** 2 + xs_p[3] ** 2) - (a**2)
        dw_dt = 2 * (xs_p[1] * vs_p[0] + xs_p[2] * vs_p[1] + xs_p[3] * vs_p[2])

        vs[0] = (1 / (2 * xs[1])) * (
            (dw_dt / 2)
            + (
                (w * dw_dt + 4 * (a**2) * xs_p[3] * vs_p[2])
                / (2 * sqrt((w**2) + (4 * (a**2) * (xs_p[3] ** 2))))
            )
        )
        vs[1] = (-1 / sqrt(1 - (xs_p[3] / xs[1]) ** 2)) * (
            (vs_p[2] * xs[1] - vs[0] * xs_p[3]) / (xs[1] ** 2)
        )
        vs[2] = (1 / (1 + (xs_p[2] / xs_p[1]) ** 2)) * (
            (vs_p[1] * xs_p[1] - vs_p[0] * xs_p[2]) / (xs_p[1] ** 2)
        )

        coordinate = BoyerLindquist(xs, vels=vs, from_dxs_dt=True, a=a)
        return coordinate

    def _get_Lz(self, metric, mass_particle=1.0):
        """
        Azimuthal conserved momentum :math:`L_z = p_\\phi` for Kerr.

        For geodesic motion, the covariant momentum component

        .. math::

            L_z = p_\\phi = g_{\\phi\\mu} u^\\mu
                = g_{\\phi\\phi} u^\\phi + g_{t\\phi} u^t

        is conserved. The :math:`g_{t\\phi} u^t` term encodes frame dragging;
        omitting it (as in a naive Newtonian :math:`\\rho v_\\phi`) does not
        yield a conserved quantity on Kerr geodesics.

        This implementation uses the ansatz :math:`u^0 = 1`,
        :math:`u^i = \\mathrm{d}x^i/\\mathrm{d}\\tau` and sums
        :math:`g_{3\\mu} u^\\mu` with index :math:`3` corresponding to
        :math:`\\phi`.

        Parameters
        ----------
        metric : BaseMetric
            Metric providing :meth:`metric` compatible with ``self.xs``
            (typically :class:`~relatipy.numeric.metrics.kerr_metric.Kerr`).
        mass_particle : float, optional
            Multiplicative mass factor applied to the contraction (default 1).

        Returns
        -------
        numpy.ndarray
            Values of :math:`L_z` per column; shape ``(N,)`` when ``metric``
            returns ``(4, 4, N)`` and each ``self.dxs_dt[i]`` has length ``N``.

        Notes
        -----
        The implementation indexes ``g[3, :, :]`` and uses
        ``len(self.dxs_dt[0])`` to size the four-velocity; callers must use a
        batched layout consistent with :meth:`BaseMetric.metric` returning
        shape ``(4, 4, N)`` (see that method’s batching convention).

        Examples
        --------
        >>> import numpy as np  # doctest: +SKIP
        >>> from relatipy.numeric.coordinates import BoyerLindquist  # doctest: +SKIP
        >>> from relatipy.numeric.metrics import Kerr  # doctest: +SKIP
        >>> xs = np.array([[0.0], [10.0], [np.pi / 2], [0.0]])  # doctest: +SKIP
        >>> vs = np.zeros((3, 1))  # doctest: +SKIP
        >>> bl = BoyerLindquist(xs, vels=vs, a=0.5, from_dxs_dt=False)  # doctest: +SKIP
        >>> lz = bl._get_Lz(Kerr(1.0, 0.5))  # doctest: +SKIP
        """
        from numpy import einsum, ones

        g = metric.metric(self.xs)

        u = ones((4, len(self.dxs_dt[0])))
        u[1:, :] = self.dxs_dt

        return mass_particle * einsum("jn,jn->n", g[3, :, :], u)

    def _get_Q(self, metric):
        """
        Carter constant :math:`Q` from the Kerr metric and state.

        Uses the conserved energy :math:`E` from :meth:`CoordinateBase._get_E`,
        :math:`L_z` from :meth:`_get_Lz`, the metric component
        :math:`g_{\\theta\\theta}`, and :math:`\\dot\\theta` to form

        .. math::

            Q = p_\\theta^2 + \\cos^2\\theta\\left(
                a^2(1 - E^2) + \\frac{L_z^2}{\\sin^2\\theta}
                \\right),

        with :math:`p_\\theta = g_{\\theta\\theta}\\, \\dot\\theta`.

        Parameters
        ----------
        metric : BaseMetric
            Kerr (or compatible) metric instance.

        Returns
        -------
        numpy.ndarray
            Carter :math:`Q` per batch column, same length convention as
            :meth:`_get_Lz`.

        Notes
        -----
        Batching requirements are the same as for :meth:`_get_Lz`.

        Examples
        --------
        >>> import numpy as np  # doctest: +SKIP
        >>> from relatipy.numeric.coordinates import BoyerLindquist  # doctest: +SKIP
        >>> from relatipy.numeric.metrics import Kerr  # doctest: +SKIP
        >>> xs = np.array([[0.0], [10.0], [np.pi / 2], [0.0]])  # doctest: +SKIP
        >>> vs = np.zeros((3, 1))  # doctest: +SKIP
        >>> bl = BoyerLindquist(xs, vels=vs, a=0.5, from_dxs_dt=False)  # doctest: +SKIP
        >>> q = bl._get_Q(Kerr(1.0, 0.5))  # doctest: +SKIP
        """
        a = self.a
        r, theta = self.xs[1], self.xs[2]

        E = self._get_E(metric)
        Lz = self._get_Lz(metric)

        g = metric.metric(self.xs)

        g_thth = g[2, 2, :]
        p_theta = g_thth * self.dxs_dt[1]

        return p_theta**2 + numpy.cos(theta) ** 2 * (
            a**2 * (1 - E**2) + Lz**2 / numpy.sin(theta) ** 2
        )

