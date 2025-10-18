import numpy
from .base import CoordinateBase
from ..constants import _c


class Cartesian(CoordinateBase):
    def __init__(self, xs, vels=None, from_dxs_dt=False):
        super().__init__(
            xs, vels=vels, from_dxs_dt=from_dxs_dt, system_name="Cartesian"
        )

    def _get_vs_from_dxs_dt(self):
        x1_prime_dot = self.dxs_dt[0]
        x2_prime_dot = self.dxs_dt[1]
        x3_prime_dot = self.dxs_dt[2]

        return numpy.array([x1_prime_dot, x2_prime_dot, x3_prime_dot])

    def _get_dxs_dt_from_vs(self):
        dx1_dt = self.vs[0]
        dx2_dt = self.vs[1]
        dx3_dt = self.vs[2]
        return numpy.array([dx1_dt, dx2_dt, dx3_dt])

    def convert_to_cartesian(self):
        return self
