import numpy
from numpy import sin, cos, sqrt, arctan2
from .base import CoordinateBase


class Cylindrical(CoordinateBase):
    def __init__(self, xs, vels=None, from_dxs_dt=False):
        super().__init__(
            xs, vels=vels, from_dxs_dt=from_dxs_dt, system_name="Cylindrical"
        )

    def _get_vs_from_dxs_dt(self):
        # v_rho = rho_dot, v_phi = rho * phi_dot, v_z = z_dot
        v_rho = self.dxs_dt[0]
        v_phi = self.dxs_dt[1] * self.xs[1]
        v_z   = self.dxs_dt[2]
        return numpy.array([v_rho, v_phi, v_z])

    def _get_dxs_dt_from_vs(self):
        # drho/dt = v_rho, dphi/dt = v_phi / rho, dz/dt = v_z
        drho_dt = self.vs[0]
        dphi_dt = self.vs[1] / self.xs[1]
        dz_dt   = self.vs[2]
        return numpy.array([drho_dt, dphi_dt, dz_dt])

    @staticmethod
    def _convert_to_cartesian(xs, vs):
        xs_p = numpy.zeros_like(xs)
        vs_p = numpy.zeros_like(vs)

        rho, phi, z = xs[1], xs[2], xs[3]
        sin_phi = sin(phi)
        cos_phi = cos(phi)

        xs_p[0] = xs[0]
        xs_p[1] = rho * cos_phi
        xs_p[2] = rho * sin_phi
        xs_p[3] = z

        # vx = cos(phi)*v_rho - sin(phi)*v_phi
        vs_p[0] = cos_phi * vs[0] - sin_phi * vs[1]
        # vy = sin(phi)*v_rho + cos(phi)*v_phi
        vs_p[1] = sin_phi * vs[0] + cos_phi * vs[1]
        # vz = v_z
        vs_p[2] = vs[2]

        return xs_p, vs_p

    @staticmethod
    def _convert_from_cartesian(xs_p, vs_p):
        xs = numpy.zeros_like(xs_p)
        vs = numpy.zeros_like(vs_p)

        xs[0] = xs_p[0]
        xs[1] = sqrt(xs_p[1] ** 2 + xs_p[2] ** 2)  # rho
        xs[2] = arctan2(xs_p[2], xs_p[1])            # phi
        xs[3] = xs_p[3]                               # z

        sin_phi = sin(xs[2])
        cos_phi = cos(xs[2])

        # v_rho = cos(phi)*vx + sin(phi)*vy
        vs[0] = cos_phi * vs_p[0] + sin_phi * vs_p[1]
        # v_phi = -sin(phi)*vx + cos(phi)*vy
        vs[1] = -sin_phi * vs_p[0] + cos_phi * vs_p[1]
        # v_z = vz
        vs[2] = vs_p[2]

        coordinate = Cylindrical(xs, vels=vs, from_dxs_dt=False)
        return coordinate