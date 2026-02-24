import numpy
from numpy import sin, cos, sqrt, arctan2, arccos
from .base import CoordinateBase


class BoyerLindquist(CoordinateBase):
    def __init__(self, xs, vels=None, a=None, from_dxs_dt=False):
        super().__init__(
            xs, vels=vels, from_dxs_dt=from_dxs_dt, system_name="BoyerLindquist", a=a
        )
        if self.a is None:
            raise ValueError(
                "The spin parameter 'a' must be provided for Boyer-Lindquist coordinates."
            )

    def _get_dxs_dt_from_vs(self):
        # Matrix(
        # [[x1_prime_dot*sqrt(a**2*cos(x2)**2 + x1**2)/sqrt(a**2 + x1**2)],
        # [x2_prime_dot*sqrt(a**2*cos(x2)**2 + x1**2)],
        # [x3_prime_dot*sqrt(a**2 + x1**2)*sin(x2)]])
        a = self.a
        xs = self.xs
        vs = self.vs


        print(self.a)
        sqrt_cos = sqrt(a**2 * cos(xs[2]) ** 2 + xs[1] ** 2)
        sqrt_a = sqrt(a**2 + xs[1] ** 2)
        sin_x2 = sin(xs[2])

        dx1_dt = vs[0] / sqrt_cos * sqrt_a
        dx2_dt = vs[1] / sqrt_cos
        dx3_dt = vs[2] / (sqrt_a * sin_x2)

        return numpy.array([dx1_dt, dx2_dt, dx3_dt])

    def _get_vs_from_dxs_dt(self):
        # Matrix(
        # [[x1_prime_dot*sqrt(a**2*cos(x2)**2 + x1**2)/sqrt(a**2 + x1**2)],
        # [x2_prime_dot*sqrt(a**2*cos(x2)**2 + x1**2)],
        # [x3_prime_dot*sqrt(a**2 + x1**2)*sin(x2)]])
        a = self.a
        xs = self.xs
        dxs_dt = self.dxs_dt

        sqrt_cos = sqrt(a**2 * cos(xs[2]) ** 2 + xs[1] ** 2)
        sqrt_a = sqrt(a**2 + xs[1] ** 2)
        sin_x2 = sin(xs[2])

        v1 = dxs_dt[0] * sqrt_cos / sqrt_a
        v2 = dxs_dt[1] * sqrt_cos
        v3 = dxs_dt[2] * sqrt_a * sin_x2

        return numpy.array([v1, v2, v3])

    @staticmethod
    def _convert_to_cartesian(xs, vs, a):
        xs_p = numpy.zeros_like(xs)
        vs_p = numpy.zeros_like(vs)

        xs[0] = xs_p[0]

        xa = sqrt(xs[1] ** 2 + a**2)
        sin_norm = xa * sin(xs[2])
        xs_p[1] = sin_norm * cos(xs[3])
        xs_p[2] = sin_norm * sin(xs[3])
        xs_p[3] = xs[1] * cos(xs[2])

        vs_p[0] = (
            (xs[1] * vs[0] * sin(xs[2]) * cos(xs[3]) / xa)
            + (xa * cos(xs[2]) * cos(xs[3]) * vs[1])
            - (xa * sin(xs[2]) * sin(xs[3]) * vs[2])
        )
        vs_p[1] = (
            (xs[1] * vs[0] * sin(xs[2]) * sin(xs[3]) / xa)
            + (xa * cos(xs[2]) * sin(xs[3]) * vs[1])
            + (xa * sin(xs[2]) * cos(xs[3]) * vs[2])
        )
        vs_p[2] = (vs[0] * cos(xs[2])) - (xs[1] * sin(xs[2]) * vs[1])

        return xs_p, vs_p

    @staticmethod
    def _convert_from_cartesian(xs_p, vs_p, a):

        xs = numpy.zeros_like(xs_p)
        vs = numpy.zeros_like(vs_p)

        xs[0] = xs_p[0]

        w = (xs_p[1] ** 2 + xs_p[2] ** 2 + xs_p[3] ** 2) - (a**2)
        xs[1] = sqrt(0.5 * (w + sqrt((w**2) + (4 * (a**2) * (xs_p[3] ** 2)))))
        xs[2] = arccos(xs_p[3] / xs[1])
        xs[3] = arctan2(xs_p[2], xs_p[1])

        w = (xs_p[1] ** 2 + xs_p[2] ** 2 + xs_p[3] ** 2) - (a**2)
        dw_dt = 2 * (xs_p[1] * vs_p[0] + xs_p[2] * vs_p[1] + xs_p[3] * vs_p[2])

        vs[0] = (1 / (2 * xs[1])) * (
            (dw_dt / 2)
            + (
                (w * dw_dt + 4 * (a**2) * xs_p[3] * vs_p[2])
                / (2 * sqrt((w**2) + (4 * (a**2) * (xs_p[3] ** 2))))
            )
        )
        vs[1] = (-1 / sqrt(1 - (xs_p[3] / xs[1]) ** 2)) * (
            (vs_p[2] * xs[1] - vs[0] * xs_p[3]) / (xs[1] ** 2)
        )
        vs[2] = (1 / (1 + (xs_p[2] / xs_p[1]) ** 2)) * (
            (vs_p[1] * xs_p[1] - vs_p[0] * xs_p[2]) / (xs_p[1] ** 2)
        )

        coordinate = BoyerLindquist(xs, vels=vs, from_dxs_dt=False, a=a)
        return coordinate

    def _get_Q(self, metric):
        "Obtener la constante de Carter Q"
        a = self.a
        r, theta = self.xs[1], self.xs[2]

        E = self._get_E(metric)
        Lz = self._get_Lz()
        
        g = metric.metric(self.xs)  # (4, 4, N)
        
        # p_theta = g_theta_theta * theta_dot
        g_thth = g[2, 2, :]
        p_theta = g_thth * self.dxs_dt[1]
        
        return p_theta**2 + numpy.cos(theta) ** 2 * (a**2 * (1 - E**2) + Lz**2 / numpy.sin(theta) ** 2)