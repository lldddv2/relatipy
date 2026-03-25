"""
Abstract base class for numerical spacetime metrics.

This module defines :class:`BaseMetric`, which ties together a mass parameter,
the Schwarzschild radius :math:`R_s = 2GM/c^2`, optional coordinate-system
metadata, and a :class:`~relatipy.numeric.geodesic.geodesic.Geodesic` helper.
Concrete metrics implement the dimensionless line element in subclasses.

Quantities are expressed in **geometric units** (:math:`G=c=1`) for internal
components unless a method explicitly converts to SI via
:meth:`BaseMetric._metric_geom_to_si` or
:meth:`BaseMetric._christoffel_dimensionless_to_si`.

Notes
-----
Subclasses must implement :meth:`BaseMetric._metric_dimensionless` and
:meth:`BaseMetric._get_christoffel_symbols` to supply :math:`g_{\\mu\\nu}` and
Christoffel symbols :math:`\\Gamma^\\rho_{\\mu\\nu}` in the chosen coordinates.

References
----------
See subclass modules (e.g. Schwarzschild, Kerr) for explicit metric forms.

Examples
--------
>>> import numpy
>>> from relatipy.numeric.metrics.schwarzschild_metric import Schwarzschild
>>> m = Schwarzschild(mass=2e30)
>>> x = numpy.array([0.0, 1e11, numpy.pi / 2, 0.0])
>>> g = m.metric(x)
>>> g.shape
(4, 4)
"""

import numpy
from itertools import product

from ..constants import _c, _c_SI, _G
from ..utils.dimensions import validator
from ..geodesic.geodesic import Geodesic


class BaseMetric:
    """
    Base class for numerical metrics and related geometric quantities.

    Stores the central mass, derived Schwarzschild radius, coordinate-system
    label, auxiliary keyword arguments, and a bound :class:`~relatipy.numeric.geodesic.geodesic.Geodesic`
    instance. Public methods delegate the actual tensor components to
    subclasses.

    Parameters
    ----------
    mass : float
        Gravitational source mass :math:`M` (kg). Used to compute
        :math:`R_s = 2GM/c^2`.
    valid_coordinate : str, optional
        Name of the coordinate system implemented by the subclass (e.g.
        ``"Cartesian"``, ``"Spherical"``). Default is ``"Cartesian"``.
    kwargs : dict, optional
        Extra keyword arguments stored on the instance for coordinate or
        integrator code.

    Attributes
    ----------
    mass : float
        Validated scalar mass (kg).
    valid_coordinate : str
        Coordinate system label.
    R_s : float
        Schwarzschild radius :math:`2GM/c^2` (m).
    kwargs : dict
        Additional keyword arguments.
    geodesic : Geodesic
        Geodesic integrator bound to this metric.

    Examples
    --------
    >>> from relatipy.numeric.metrics.schwarzschild_metric import Schwarzschild
    >>> m = Schwarzschild(mass=1e31)
    >>> m.R_s > 0
    True
    """

    def __init__(self, mass, valid_coordinate="Cartesian", kwargs={}):
        """
        Initialize the metric with mass and coordinate metadata.

        Parameters
        ----------
        mass : float
            Source mass :math:`M` (kg).
        valid_coordinate : str, optional
            Coordinate system name. Default ``"Cartesian"``.
        kwargs : dict, optional
            Optional keyword arguments stored in :attr:`BaseMetric.kwargs`.

        Examples
        --------
        >>> from relatipy.numeric.metrics.schwarzschild_metric import Schwarzschild
        >>> m = Schwarzschild(mass=2e30)
        >>> m.valid_coordinate
        'Spherical'
        """
        self.mass = validator.validate_scalar(mass)
        self.valid_coordinate = valid_coordinate
        self.R_s = 2 * _G * self.mass / _c**2  # Schwarzschild radius
        self.kwargs = kwargs
        self.geodesic = Geodesic(self)

    def metric(self, xs, dimensionless=True):
        """
        Metric tensor components :math:`g_{\\mu\\nu}` at one or many events.

        For a single point, returns a :math:`4 \\times 4` matrix. For a batch,
        the input is interpreted as **four rows and** ``N`` **columns**
        (shape ``(4, N)``); the returned array has shape ``(4, 4, N)`` after
        restoring the leading dimensions expected by calling code.

        Parameters
        ----------
        xs : array_like
            Coordinates. One point: shape ``(4,)``. Multiple points: shape
            ``(4, N)`` with columns :math:`(x^0, x^1, x^2, x^3)^\\top` per event.
        dimensionless : bool, optional
            If ``True`` (default), return :math:`g_{\\mu\\nu}` in geometric units.
            If ``False``, convert to SI using :meth:`BaseMetric._metric_geom_to_si`.

        Returns
        -------
        numpy.ndarray
            If ``xs`` is 1D, shape ``(4, 4)``. If ``xs`` is 2D with shape
            ``(4, N)``, shape ``(4, 4, N)``.

        Raises
        ------
        ValueError
            If ``xs`` is not 1D or 2D.

        Examples
        --------
        >>> import numpy
        >>> from relatipy.numeric.metrics.schwarzschild_metric import Schwarzschild
        >>> m = Schwarzschild(mass=2e30)
        >>> x = numpy.array([0.0, 1e11, numpy.pi / 2, 0.0])
        >>> g = m.metric(x)
        >>> g[0, 0] < 1
        True
        """
        xs = numpy.asarray(xs, dtype=object)

        if xs.ndim == 1:
            xs = validator.validate_vector(xs)
            metric = self._metric_dimensionless(xs)
            return metric if dimensionless else self._metric_geom_to_si(metric)

        if xs.ndim == 2:
            xs = xs.T  # (4, N) -> (N, 4)
            metrics = self._metric_dimensionless(xs)
            metrics = metrics if dimensionless else numpy.array([self._metric_geom_to_si(g) for g in metrics])
            return metrics.T  # (N, 4, 4) -> (4, 4, N)

        raise ValueError(f"xs must be 1D (single point) or 2D (N points), got shape {xs.shape}")

    def _metric_dimensionless(self, xs):
        """
        Metric :math:`g_{\\mu\\nu}` in geometric units at given coordinates.

        Parameters
        ----------
        xs : numpy.ndarray
            Coordinates in the subclass's chart. Single event, shape ``(4,)``;
            or batch, shape ``(N, 4)``.

        Returns
        -------
        numpy.ndarray
            Shape ``(4, 4)`` or ``(N, 4, 4)`` consistent with ``xs``.

        Raises
        ------
        NotImplementedError
            Always, unless overridden by a subclass.

        Examples
        --------
        >>> import numpy
        >>> from relatipy.numeric.metrics.schwarzschild_metric import Schwarzschild
        >>> m = Schwarzschild(mass=2e30)
        >>> x = numpy.array([0.0, 1e10, numpy.pi / 2, 0.0])
        >>> g = m._metric_dimensionless(x)
        >>> g.shape
        (4, 4)
        """
        raise NotImplementedError("This method should be implemented by subclasses.")

    def get_4velocity(self, coordinate):
        """
        Construct a four-velocity :math:`u^\\mu` from coordinate-time derivatives.

        Assumes parametrization by coordinate time :math:`t = x^0` with
        :math:`u^i = (u^0/c)\\,\\mathrm{d}x^i/\\mathrm{d}t` for :math:`i \\in \\{1,2,3\\}`.
        The time component :math:`u^0` is fixed so that
        :math:`g_{\\mu\\nu} u^\\mu u^\\nu = -c^2` (Lorentzian signature with the
        conventions used in this codebase).

        Parameters
        ----------
        coordinate : object
            Must provide ``xs`` (length-4 position) and ``dxs_dt`` (length-3
            array of :math:`\\mathrm{d}x^i/\\mathrm{d}t`), as in
            :class:`~relatipy.numeric.coordinates.base.CoordinateBase`.

        Returns
        -------
        numpy.ndarray
            Four-vector :math:`u^\\mu`, shape ``(4,)``.

        Examples
        --------
        >>> import numpy
        >>> from relatipy.numeric.coordinates.cartesian import Cartesian
        >>> from relatipy.numeric.metrics.Minkowski_metric import Minkowski
        >>> m = Minkowski()
        >>> coord = Cartesian(
        ...     numpy.array([0.0, 0.0, 0.0, 0.0]),
        ...     vels=numpy.array([1e5, 0.0, 0.0]),
        ...     from_dxs_dt=False,
        ... )
        >>> u = m.get_4velocity(coord)
        >>> u.shape
        (4,)
        """
        us = numpy.zeros(4)
        g = self.metric(coordinate.xs)

        u_t2 = g[0, 0]
        for i in range(1, 4):
            u_t2 += 2 / _c * g[0, i] * coordinate.dxs_dt[i - 1]
            for j in range(1, 4):
                u_t2 += 1 / _c**2 * g[i, j] * coordinate.dxs_dt[i - 1] * coordinate.dxs_dt[j - 1]

        # u_t2 = 1/(1-(vs[1]*vs[1] + vs[2]*vs[2] + vs[3]*vs[3])/_c**2)
        us[0] = _c * numpy.sqrt(1 / u_t2)

        for i in range(1, 4):
            us[i] = coordinate.dxs_dt[i - 1] * us[0] / _c
        return us

    def get_4state_vector(self, coordinate):
        """
        Concatenate position and four-velocity into an 8-component state vector.

        Parameters
        ----------
        coordinate : object
            Object with ``xs`` and compatible with :meth:`BaseMetric.get_4velocity`.

        Returns
        -------
        numpy.ndarray
            Eight-component state, shape ``(8,)``: ``xs`` followed by
            :math:`u^\\mu` from :meth:`BaseMetric.get_4velocity`.

        Examples
        --------
        >>> import numpy
        >>> from relatipy.numeric.coordinates.cartesian import Cartesian
        >>> from relatipy.numeric.metrics.Minkowski_metric import Minkowski
        >>> m = Minkowski()
        >>> coord = Cartesian(
        ...     numpy.array([0.0, 0.0, 0.0, 0.0]),
        ...     vels=numpy.array([0.0, 0.0, 0.0]),
        ...     from_dxs_dt=False,
        ... )
        >>> y = m.get_4state_vector(coord)
        >>> y.shape
        (8,)
        """
        return numpy.concatenate((coordinate.xs, self.get_4velocity(coordinate)))

    def get_dxs_dt_from_4velocity(self, us):
        """
        Spatial coordinate derivatives :math:`\\mathrm{d}x^i/\\mathrm{d}t` from :math:`u^\\mu`.

        Uses :math:`\\mathrm{d}x^i/\\mathrm{d}t = c\\,u^i/u^0`.

        Parameters
        ----------
        us : array_like
            Four-velocity components :math:`u^\\mu` with nonzero :math:`u^0`.

        Returns
        -------
        numpy.ndarray
            Same shape as ``us``; each spatial row scaled by :math:`c/u^0`.

        Raises
        ------
        ZeroDivisionError
            If ``us[0] == 0`` (division by zero).

        Examples
        --------
        >>> import numpy
        >>> from relatipy.numeric.metrics.Minkowski_metric import Minkowski
        >>> m = Minkowski()
        >>> u = numpy.array([2.0, 0.1, 0.0, 0.0])
        >>> dx = m.get_dxs_dt_from_4velocity(u)
        >>> numpy.isclose(dx[1], 0.05)
        True
        """
        return _c * us / us[0]

    def get_christoffel_symbols(self, xs, dimensionless=True):
        """
        Christoffel symbols :math:`\\Gamma^\\rho_{\\mu\\nu}` of the Levi-Civita connection.

        For a batch of points, ``xs`` must have shape ``(N, 4)`` (unlike
        :meth:`BaseMetric.metric`, which uses ``(4, N)`` for the 2D case).

        Parameters
        ----------
        xs : array_like
            Single point, shape ``(4,)``; or ``N`` points, shape ``(N, 4)``.
        dimensionless : bool, optional
            If ``True`` (default), return symbols in geometric units. If ``False``,
            map to SI via :meth:`BaseMetric._christoffel_dimensionless_to_si`.

        Returns
        -------
        numpy.ndarray
            Shape ``(4, 4, 4)`` for one point, or ``(N, 4, 4, 4)`` for a batch.

        Raises
        ------
        ValueError
            If ``xs`` is not 1D or 2D.

        Examples
        --------
        >>> import numpy
        >>> from relatipy.numeric.metrics.schwarzschild_metric import Schwarzschild
        >>> m = Schwarzschild(mass=2e30)
        >>> x = numpy.array([0.0, 1e11, numpy.pi / 2, 0.0])
        >>> Gamma = m.get_christoffel_symbols(x)
        >>> Gamma.shape
        (4, 4, 4)
        """
        xs = numpy.asarray(xs, dtype=object)

        if xs.ndim == 1:
            xs = validator.validate_vector(xs)
            christoffel = self._get_christoffel_symbols(xs)
            return christoffel if dimensionless else self._christoffel_dimensionless_to_si(christoffel)

        if xs.ndim == 2:
            christoffels = self._get_christoffel_symbols(xs)
            return christoffels if dimensionless else numpy.array([self._christoffel_dimensionless_to_si(G) for G in christoffels])

        raise ValueError(f"xs must be 1D (single point) or 2D (N points), got shape {xs.shape}")

    @staticmethod
    def _christoffel_dimensionless_to_si(Gamma_geom):
        """
        Map Christoffel symbols from geometric to SI coordinates.

        Uses the reference length :math:`L_\\mathrm{ref}` and speed of light
        :math:`c_\\mathrm{SI}` from project constants. Under
        :math:`t_\\mathrm{geom} = t_\\mathrm{SI}/T_\\mathrm{ref}` and
        :math:`r_\\mathrm{geom} = r_\\mathrm{SI}/L_\\mathrm{ref}`,

        .. math::

            \\Gamma^{\\rho,\\mathrm{SI}}_{\\mu\\nu}
            = c_\\mathrm{SI}^{\\,(\\delta_\\mu^0 + \\delta_\\nu^0 - \\delta_\\rho^0)}
              L_\\mathrm{ref}^{\\,f(\\rho,\\mu,\\nu)}
              \\Gamma^{\\rho,\\mathrm{geom}}_{\\mu\\nu},

        where the exponent :math:`f` distinguishes time/radial versus angular
        indices as implemented in code.

        Parameters
        ----------
        Gamma_geom : array_like
            Christoffel symbols in geometric units, shape ``(4, 4, 4)`` with
            indices ``(rho, mu, nu)``.

        Returns
        -------
        numpy.ndarray
            Christoffel symbols in SI convention, same shape as ``Gamma_geom``.

        Examples
        --------
        >>> import numpy
        >>> Gamma = numpy.zeros((4, 4, 4))
        >>> G_si = BaseMetric._christoffel_dimensionless_to_si(Gamma)
        >>> numpy.allclose(G_si, 0)
        True
        """
        from ..constants import _L_ref

        Gamma_geom = numpy.asarray(Gamma_geom, dtype=float)
        Gamma_si = numpy.zeros_like(Gamma_geom)

        for rho in range(4):
            for mu, nu in product(range(4), repeat=2):
                exp_c = (mu == 0) + (nu == 0) - (rho == 0)
                in_tr = lambda a: (a == 0) or (a == 1)
                exp_L = in_tr(rho) - in_tr(mu) - in_tr(nu)
                factor = (_c_SI**exp_c) * (_L_ref**exp_L)
                Gamma_si[rho, mu, nu] = factor * Gamma_geom[rho, mu, nu]

        return Gamma_si

    @staticmethod
    def _metric_geom_to_si(g_geom):
        """
        Map metric components from geometric to SI coordinates.

        Temporal components pick up factors of :math:`c_\\mathrm{SI}` per index;
        angular blocks (indices 2 and 3) pick up powers of :math:`L_\\mathrm{ref}`
        because quantities such as :math:`g_{\\phi\\phi}` involve squared radii
        expressed in geometric length units.

        Parameters
        ----------
        g_geom : array_like
            Metric :math:`g_{\\mu\\nu}` in geometric units, shape ``(4, 4)``.

        Returns
        -------
        numpy.ndarray
            Metric in SI convention, shape ``(4, 4)``.

        Examples
        --------
        >>> import numpy
        >>> g = numpy.eye(4)
        >>> g[1, 1] = g[2, 2] = g[3, 3] = -1
        >>> g_si = BaseMetric._metric_geom_to_si(g)
        >>> g_si.shape
        (4, 4)
        """
        from ..constants import _L_ref

        g_geom = numpy.asarray(g_geom, dtype=float)
        g_si = numpy.zeros_like(g_geom)

        for mu, nu in product(range(4), repeat=2):
            n_zero = (mu == 0) + (nu == 0)
            n_r = (mu in [2, 3]) + (nu in [2, 3])
            g_si[mu, nu] = (_c_SI**n_zero) * (_L_ref**n_r) * g_geom[mu, nu]

        return g_si
