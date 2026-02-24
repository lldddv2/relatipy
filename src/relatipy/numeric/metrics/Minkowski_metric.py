import numpy
import itertools

from ..constants import _G, _c
from .base import BaseMetric
from ..geodesic.geodesic import Geodesic
from ..utils.dimensions import validator

class Minkowski(BaseMetric):
    def __init__(self):
        pass

    def _metric_dimensionless(self, xs):
        """
        Returns the Minkowski metric tensor.

        Parameters
        ----------
        xs : list
            List of coordinates [t, r, theta, phi] of spherical coordinates.
        mass : float
            Mass of the black hole in kg.
        """
        xs = numpy.asarray(xs, dtype=float)
        
        if xs.ndim == 1:
            return numpy.diag([1, -1, -1, -1])
        
        # Multiple points: shape (N, 4) — métrica constante, se repite N veces
        N = len(xs)
        metrics = numpy.zeros((N, 4, 4))
        metrics[:, 0, 0] =  1
        metrics[:, 1, 1] = -1
        metrics[:, 2, 2] = -1
        metrics[:, 3, 3] = -1
        return metrics

    def _get_christoffel_symbols(self, xs):
        """
        Returns the Christoffel symbols of the Minkowski metric (all zero).

        Parameters
        ----------
        xs : array of shape (4,) or (N, 4)
            Coordinates [x0, x1, x2, x3].
        """
        xs = numpy.asarray(xs, dtype=float)
        if xs.ndim == 2:
            N = len(xs)
            return numpy.zeros((N, 4, 4, 4))
        return numpy.zeros((4, 4, 4))
