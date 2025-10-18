from .base import BaseCoordinate, xs_p
import sympy as sp

class BoyerLindquist(BaseCoordinate):
    def __init__(self):
        super().__init__()
        self.a = sp.symbols('a', positive=True)

    def cartesian_conversion_rule(self):
        xa = sp.sqrt(xs_p[0]**2 + self.a**2)
        sin_norm = xa * sp.sin(xs_p[1])
        x = sin_norm*sp.cos(xs_p[2])
        y = sin_norm*sp.sin(xs_p[2])
        z = xs_p[0]*sp.cos(xs_p[1])

        return sp.Matrix([x, y, z])