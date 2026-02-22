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
        metric = numpy.diag([1, -1, -1, -1])
        return metric

    def get_christoffel_symbols(self, xs):
        """
        Returns the Christoffel symbols of the Kerr metric.

        Parameters
        ----------
        xs : list
            List of coordinates [t, r, theta, phi] of boyer-lindquist coordinates.
        """
        Gamma = numpy.zeros((4, 4, 4))
        return Gamma
