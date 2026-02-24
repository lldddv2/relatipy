import numpy
import itertools

from .base import BaseMetric

class Schwarzschild(BaseMetric):
    def __init__(self, mass):
        super().__init__(mass, valid_coordinate="Spherical")

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
        xs = numpy.asarray(xs, dtype=float)
        
        # Single point: shape (4,)
        if xs.ndim == 1:
            A = 1 - self.R_s / xs[1]
            B = -1 / A
            C = -xs[1] ** 2
            D = -xs[1] ** 2 * numpy.sin(xs[2]) ** 2
            return numpy.diag([A, B, C, D])

        # Multiple points: shape (N, 4)
        r     = xs[:, 1]
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
        Returns the Christoffel symbols of the Schwarzschild metric.

        Parameters
        ----------
        xs : array of shape (4,) or (N, 4)
            Coordinates [x0, x1, x2, x3].
        """
        xs = numpy.asarray(xs, dtype=float)

        if xs.ndim == 2:
            N = len(xs)
            Gammas = numpy.zeros((N, 4, 4, 4))
            for i, x in enumerate(xs):
                Gammas[i] = self._get_christoffel_symbols(x)
            return Gammas

        # resto del código existente sin cambios ...
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
