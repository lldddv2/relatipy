import sympy as sp
import numpy as np
import einsteinpy.symbolic as es
from itertools import product

simplify = lambda expr: expr.expand().simplify()

# coordinates
x0, x1, x2, x3 = sp.symbols("x^0 x^1 x^2 x^3")
xs = [x0, x1, x2, x3]

class Kerr:
    def __init__(self):
        self.metic_has_been_computed = False
        self.christoffel_has_been_computed = False
        self._metric = None
        self._christoffel_symbols = None

    def metric(self):
        if not self.metic_has_been_computed:
            self._metric = self._compute_metric()
            self.metic_has_been_computed = True

        return self._metric

    def christoffel_symbols(self):
        if not self.christoffel_has_been_computed:
            self._christoffel_symbols = self._compute_christoffel_symbols()
            self.christoffel_has_been_computed = True
        return self._christoffel_symbols

    def _compute_metric(self):
        # constant sp.symbols
        G, c = sp.symbols("G c")

        # parameters
        M = sp.symbols("M")
        J = sp.symbols("J")

        # Derived parameters
        R_s = 2*G*M/c**2
        a = J/(M*c)
        R_s, a = sp.symbols("R_s a")

        # Metric
        Sigma = xs[1]**2 * (1 + a**2/xs[1]**2 * sp.cos(xs[2])**2)
        Delta = xs[1]**2 * (1 - (R_s*xs[1] + a**2)/xs[1]**2)

        A = 1 - R_s*xs[1]/Sigma
        B = - Sigma/Delta
        C = - Sigma
        D = - (xs[1]**2 + a**2 + R_s*xs[1]*a**2/Sigma * sp.sin(xs[2])**2) * sp.sin(xs[2])**2
        E = 2*R_s*xs[1]*a/Sigma * sp.sin(xs[2])**2

        metric = np.diag([A, B, C, D])
        metric[0, 3] = metric[3, 0] = E/2
        metric_ = sp.Matrix(metric)

        g = es.MetricTensor(metric, xs)
        
        return g
    
    def _compute_christoffel_symbols(self):
        ch = es.ChristoffelSymbols.from_metric(self.metric())

        christoffel = np.zeros((4, 4, 4), dtype=object)

        for mu, nu, sigma in product(range(4), repeat=3):
            expr = simplify(ch.tensor()[mu, nu, sigma])#.subs(subs).subs(second_subs)
            christoffel[mu, nu, sigma] = expr

        christoffel = es.ChristoffelSymbols(christoffel, xs)

        return christoffel