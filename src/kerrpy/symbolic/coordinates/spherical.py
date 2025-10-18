from .base import BaseCoordinate, xs_p
import sympy as sp

class Spherical(BaseCoordinate):
    def __init__(self):
        super().__init__()

    def cartesian_conversion_rule(self):
        x = xs_p[0]*sp.sin(xs_p[1])*sp.cos(xs_p[2])
        y = xs_p[0]*sp.sin(xs_p[1])*sp.sin(xs_p[2])
        z = xs_p[0]*sp.cos(xs_p[1])

        return sp.Matrix([x, y, z])