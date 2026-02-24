import numpy
from itertools import product

from ..constants import _c, _c_SI, _G
from ..utils.dimensions import validator
from ..geodesic.geodesic import Geodesic

class BaseMetric:
    def __init__(self, mass, valid_coordinate="Cartesian", kwargs={}):
        self.mass = validator.validate_scalar(mass)
        self.valid_coordinate = valid_coordinate
        self.R_s = 2 * _G * self.mass / _c**2  # Schwarzschild radius
        self.kwargs = kwargs
        self.geodesic = Geodesic(self)

    def metric(self, xs, dimensionless=True):
        """
        Returns the metric tensor components g_{mu nu} as a 2D array.

        Parameters
        ----------
        xs : list
            List of coordinates [t, x, y, z] or array of shape (4,) or (4, N).
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
        Returns the metric tensor components g_{mu nu} as a 2D array.

        Parameters
        ----------
        xs : list
            List of coordinates dimensionless [t, x, y, z].
        """
        raise NotImplementedError("This method should be implemented by subclasses.")

    @staticmethod
    def _metric_geom_to_si(g_geom):
        """
        Convert metric from geometric convention (x0 = ct)
        to SI convention (x0 = t).

        Parameters
        ----------
        g_geom : array-like (4x4)
            Metric tensor in geometric units.

        Returns
        -------
        numpy.ndarray
            Metric tensor in SI convention.
        """
        g_geom = numpy.asarray(g_geom, dtype=float)
        g_si = numpy.zeros_like(g_geom)

        for mu, nu in product(range(4), repeat=2):
            n_zero = (mu == 0) + (nu == 0)
            g_si[mu, nu] = (_c_SI ** n_zero) * g_geom[mu, nu]

        return g_si

    def get_4velocity(self, coordinate):
        """
        Returns the four-velocity of a test particle in the given metric.

        """
        us = numpy.zeros(4)
        g = self.metric(coordinate.xs)

        u_t2 = g[0, 0]
        for i in range(1, 4):
            u_t2 += 2 / _c * g[0, i] * coordinate.dxs_dt[i - 1] 
            for j in range(1, 4):
                u_t2 += 1 / _c**2 * g[i, j] * coordinate.dxs_dt[i - 1] * coordinate.dxs_dt[j - 1]

        # u_t2 = 1/(1-(vs[1]*vs[1] + vs[2]*vs[2] + vs[3]*vs[3])/_c**2)
        us[0] = _c * numpy.sqrt(1/u_t2)
        
        for i in range(1, 4):
            us[i] = coordinate.dxs_dt[i - 1] * us[0] / _c
        return us

    def get_4state_vector(self, coordinate):
        return numpy.concatenate((coordinate.xs, self.get_4velocity(coordinate)))
    
    def get_dxs_dt_from_4velocity(self, us):
        return _c * us/us[0]

    def get_christoffel_symbols(self, xs, dimensionless=True):
        """
        Returns the Christoffel symbols of the metric.

        Parameters
        ----------
        xs : list or array of shape (4,) or (N, 4)
            List of coordinates [x0, x1, x2, x3].
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
        Gamma_geom = numpy.asarray(Gamma_geom, dtype=float)
        Gamma_si = numpy.zeros_like(Gamma_geom)

        for rho in range(4):
            for mu, nu in product(range(4), repeat=2):
                exponent = (mu == 0) + (nu == 0) - (rho == 0)
                Gamma_si[rho, mu, nu] = (_c_SI ** exponent) * Gamma_geom[rho, mu, nu]

        return Gamma_si
