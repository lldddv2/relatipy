import numpy
import sympy as sp
from itertools import product

t = sp.symbols('t')
x1 = sp.Function('x^1', positive=True)(t)
x2 = sp.Function('x^2', positive=True)(t)
x3 = sp.Function('x^3', positive=True)(t)

xs_p = numpy.array([x1, x2, x3])

class BaseCoordinate:
    def __init__(self):
        pass

    def calculate_relation_vi_dxi_dt(self):
        xs = self.cartesian_conversion_rule()
        # Definimos derivadas respecto a t
        xs_prime_dot = sp.Matrix(
            [xs_p[0], xs_p[1], xs_p[2]]
        ).diff(t)

        es = [0, 0, 0]

        for i in range(3):
            es_temp = xs.diff(xs_p[i])
            es_norm = sp.sqrt((es_temp.T @ es_temp)[0,0]).simplify()
            es[i] = sp.simplify(es_temp / es_norm)

        # Velocidad cartesiana (regla de la cadena)
        v_cart = sp.zeros(1,3)
        for i, j in product(range(3), repeat=2):
                v_cart[i] += sp.diff(xs[i], xs_p[j]) * xs_prime_dot[j]

        # Proyecciones: v·e_i
        v0= sp.simplify(v_cart.dot(es[0]))
        v1= sp.simplify(v_cart.dot(es[1]))
        v2= sp.simplify(v_cart.dot(es[2]))

        vs = sp.Matrix([v0, v1, v2]).replace(sp.Abs, lambda x: x)
        vs.simplify()
        
        return vs

    def calculate_cartesian_conversion_rule_velocity(self):
        """
        Calcula la regla de conversión de velocidades generalizadas -> cartesianas.
        Retorna un vector [dx/dt, dy/dt, dz/dt] en función de (x^i, dx^i/dt).
        """
        # coordenadas cartesianas como funciones de x^i(t)
        xs = self.cartesian_conversion_rule()

        # derivadas respecto a cada coordenada generalizada
        J = sp.Matrix([[sp.diff(xs[i], xs_p[j]) for j in range(3)] for i in range(3)])

        # derivadas de x^i respecto a t (generalized velocities)
        xs_prime_dot = sp.Matrix([sp.diff(q, t) for q in xs_p])

        # velocidad cartesiana: v = J * q_dot
        v_cart = J * xs_prime_dot

        return sp.simplify(v_cart)