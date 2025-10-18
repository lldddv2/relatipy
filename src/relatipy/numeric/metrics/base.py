import numpy
from itertools import product

from ..constants import _c

class BaseMetric:
    def __init__(self):
        pass

    def metric(self, xs):
        """
        Returns the metric tensor components g_{mu nu} as a 2D array.

        Parameters
        ----------
        xs : list
            List of coordinates [t, x, y, z].
        """
        return numpy.diag([1, -1, -1, -1])  # Minkowski metric as default
    
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

    def get_christoffel_symbols(self, xs):
        """
        Returns the Christoffel symbols of the metric.

        Parameters
        ----------
        xs : list
            List of coordinates [x0, x1, x2, x3].
        """
        raise NotImplementedError("This method should be implemented by subclasses.")