import numpy
from numpy import sin, cos, sqrt, arctan2, arccos
from .base import CoordinateBase


class Spherical(CoordinateBase):
    def __init__(self, xs, vels=None, from_dxs_dt=False):
        super().__init__(
            xs, vels=vels, from_dxs_dt=from_dxs_dt, system_name="Spherical"
        )

    def _get_vs_from_dxs_dt(self):
        # [[x1_prime_dot], [x1*x2_prime_dot], [x1*x3_prime_dot*sin(x2)]]
        x1_prime_dot = self.dxs_dt[0]
        x2_prime_dot = self.dxs_dt[1] * self.xs[1]
        x3_prime_dot = self.dxs_dt[2] * self.xs[1] * sin(self.xs[2])

        return numpy.array([x1_prime_dot, x2_prime_dot, x3_prime_dot])

    def _get_dxs_dt_from_vs(self):
        # [[x1_prime_dot], [x1*x2_prime_dot], [x1*x3_prime_dot*sin(x2)]]
        dx1_dt = self.vs[0]
        dx2_dt = self.vs[1] / self.xs[1]
        dx3_dt = self.vs[2] / (self.xs[1] * sin(self.xs[2]))
        return numpy.array([dx1_dt, dx2_dt, dx3_dt])

    @staticmethod
    def _convert_to_cartesian(xs, vs):
        xs_p = numpy.zeros_like(xs)
        vs_p = numpy.zeros_like(vs)

        xs_p[0] = xs[0]
        xs_p[1] = xs[1] * sin(xs[2]) * cos(xs[3])
        xs_p[2] = xs[1] * sin(xs[2]) * sin(xs[3])
        xs_p[3] = xs[1] * cos(xs[2])

        sin_x2 = sin(xs[2])
        cos_x2 = cos(xs[2])
        sin_x3 = sin(xs[3])
        cos_x3 = cos(xs[3])

        # vs[0] = vr/c, vs[1] = r*dθ/dt/c, vs[2] = r*sin(θ)*dφ/dt/c
        # vx/c = sin(θ)cos(φ)*vr/c + cos(θ)cos(φ)*(r*dθ/dt/c) - sin(θ)sin(φ)*(r*sin(θ)*dφ/dt/c)/sin(θ)
        vs_p[0] = sin_x2 * cos_x3 * vs[0] + cos_x2 * cos_x3 * vs[1] - sin_x3 * vs[2]

        # vy/c = sin(θ)sin(φ)*vr/c + cos(θ)sin(φ)*(r*dθ/dt/c) + cos(φ)*(r*sin(θ)*dφ/dt/c)/sin(θ)
        vs_p[1] = sin_x2 * sin_x3 * vs[0] + cos_x2 * sin_x3 * vs[1] + cos_x3 * vs[2]

        # vz/c = cos(θ)*vr/c - sin(θ)*(r*dθ/dt/c)
        vs_p[2] = cos_x2 * vs[0] - sin_x2 * vs[1]

        return xs_p, vs_p

    @staticmethod
    def _convert_from_cartesian(xs_p, vs_p):
        xs = numpy.zeros_like(xs_p)
        vs = numpy.zeros_like(vs_p)

        xs[0] = xs_p[0]
        xs[1] = sqrt(xs_p[1] ** 2 + xs_p[2] ** 2 + xs_p[3] ** 2)
        xs[2] = arccos(xs_p[3] / xs[1])
        xs[3] = arctan2(xs_p[2], xs_p[1])

        sin_x2 = sin(xs[2])
        cos_x2 = cos(xs[2])
        sin_x3 = sin(xs[3])
        cos_x3 = cos(xs[3])

        # vr = sin(theta)*cos(phi)*vx + sin(theta)*sin(phi)*vy + cos(theta)*vz
        vs[0] = sin_x2 * cos_x3 * vs_p[0] + sin_x2 * sin_x3 * vs_p[1] + cos_x2 * vs_p[2]

        # vs[1] = r*dtheta/dt = cos(theta)*cos(phi)*vx + cos(theta)*sin(phi)*vy - sin(theta)*vz
        vs[1] = cos_x2 * cos_x3 * vs_p[0] + cos_x2 * sin_x3 * vs_p[1] - sin_x2 * vs_p[2]

        # vs[2] = r*sin(theta)*dphi/dt = -sin(phi)*vx + cos(phi)*vy
        vs[2] = -sin_x3 * vs_p[0] + cos_x3 * vs_p[1]

        return Spherical(xs, vels=vs, from_dxs_dt=False)
