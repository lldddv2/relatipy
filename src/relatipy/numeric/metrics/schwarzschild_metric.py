import numpy
import itertools

from .base import BaseMetric

class Schwarzschild(BaseMetric):
    def __init__(self, mass):
        super().__init__(mass)

    def _metric_dimensionless(self, xs):
        """
        Returns the Schwarzschild metric tensor in spherical coordinates.

        Parameters
        ----------
        xs : list
            List of coordinates [t, r, theta, phi] of spherical coordinates.
        mass : float
            Mass of the black hole in kg.
        """
        A = 1 - self.R_s / xs[1]
        B = -1 / A
        C = -xs[1] ** 2
        D = -xs[1] ** 2 * numpy.sin(xs[2]) ** 2

        metric = numpy.diag([A, B, C, D])

        return metric

    def _get_christoffel_symbols(self, xs):
        """
        Returns the Christoffel symbols of the Kerr metric.

        Parameters
        ----------
        xs : list
            List of coordinates [t, r, theta, phi] of boyer-lindquist coordinates.
        """
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
