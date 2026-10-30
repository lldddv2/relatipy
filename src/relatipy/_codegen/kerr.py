"""Symbolic Kerr metric and geodesic right-hand side for code generation.

Conventions follow ``native/README.md``: Boyer--Lindquist coordinates
``(t, r, theta, phi)``, signature ``(-,+,+,+)``, ``G = c = 1``, explicit mass
``M`` and spin length ``a``, ``Sigma = r^2 + a^2 cos^2 theta`` and
``Delta = r^2 - 2 M r + a^2``.
"""

from __future__ import annotations

from dataclasses import dataclass

from ._sympy import require_sympy

sp = require_sympy()

from .optimizer import contract_pythagorean, reduce_angle  # noqa: E402


@dataclass(frozen=True)
class KerrSymbols:
    """Symbols of the Kerr derivation.

    ``sigma`` and ``delta`` are opaque stand-ins for ``Sigma`` and ``Delta``
    when the derivation is ``prepared``; they appear in no expression
    otherwise.
    """

    mass: object
    spin: object
    r: object
    theta: object
    u: tuple
    sigma: object
    delta: object


def kerr_symbols() -> KerrSymbols:
    """Real symbols ``mass, spin, r, th, u0..u3, sigma, delta``."""
    mass, spin, r, theta = sp.symbols("mass spin r th", real=True)
    u = sp.symbols("u0:4", real=True)
    sigma, delta = sp.symbols("sigma delta", real=True)
    return KerrSymbols(mass, spin, r, theta, tuple(u), sigma, delta)


def sigma_expression(s: KerrSymbols):
    """``Sigma = r^2 + a^2 cos^2 theta``."""
    return s.r**2 + s.spin**2 * sp.cos(s.theta) ** 2


def delta_expression(s: KerrSymbols):
    """``Delta = r^2 - 2 M r + a^2``."""
    return s.r**2 - 2 * s.mass * s.r + s.spin**2


def metric(s: KerrSymbols, sigma=None, delta=None):
    """Covariant Kerr metric ``g_(mu nu)`` with optional opaque Sigma, Delta."""
    sigma = sigma_expression(s) if sigma is None else sigma
    delta = delta_expression(s) if delta is None else delta
    sine2 = sp.sin(s.theta) ** 2
    common = 2 * s.mass * s.r / sigma
    g = sp.zeros(4, 4)
    g[0, 0] = -(1 - common)
    g[0, 3] = g[3, 0] = -s.spin * sine2 * common
    g[1, 1] = sigma / delta
    g[2, 2] = sigma
    g[3, 3] = sine2 * (s.r**2 + s.spin**2 + s.spin**2 * sine2 * common)
    return g


def block_determinant_residual(s: KerrSymbols):
    """``g_tt g_phiphi - g_tphi^2 + Delta sin^2 theta`` simplified; zero."""
    g = metric(s)
    value = g[0, 0] * g[3, 3] - g[0, 3] ** 2 + delta_expression(s) * sp.sin(s.theta) ** 2
    return sp.simplify(sp.trigsimp(sp.cancel(value)))


def inverse_metric(g, s: KerrSymbols, delta=None):
    """Inverse of a Kerr metric from its ``t-phi`` block.

    Uses ``g_tt g_phiphi - g_tphi^2 = -Delta sin^2 theta``, checked by
    :func:`block_determinant_residual`.
    """
    delta = delta_expression(s) if delta is None else delta
    determinant = -delta * sp.sin(s.theta) ** 2
    inverse = sp.zeros(4, 4)
    inverse[0, 0] = g[3, 3] / determinant
    inverse[3, 3] = g[0, 0] / determinant
    inverse[0, 3] = inverse[3, 0] = -g[0, 3] / determinant
    inverse[1, 1] = 1 / g[1, 1]
    inverse[2, 2] = 1 / g[2, 2]
    return inverse.applyfunc(sp.cancel)


def geodesic_acceleration(s: KerrSymbols, *, prepared: bool = True) -> list:
    """``-Gamma^mu_(alpha beta) u^alpha u^beta`` of the Levi--Civita connection.

    Parameters
    ----------
    s
        Symbols from :func:`kerr_symbols`.
    prepared
        Write the result in terms of the opaque ``s.sigma`` and ``s.delta``
        (their derivatives are applied by the chain rule) instead of their
        explicit expressions.

    Returns
    -------
    list
        Four SymPy expressions, the contraction ``-g^(mu sigma)
        (d_alpha g_(sigma beta) - d_sigma g_(alpha beta) / 2) u^alpha u^beta``.
    """
    coordinates = (sp.Symbol("t"), s.r, s.theta, sp.Symbol("ph"))
    if prepared:
        sigma_f = sp.Function("Sigma")(s.r, s.theta)
        delta_f = sp.Function("Delta")(s.r)
        g = metric(s, sigma_f, delta_f)
        derivative_rules = {
            sp.Derivative(sigma_f, s.r): 2 * s.r,
            sp.Derivative(sigma_f, s.theta): -2 * s.spin**2 * sp.sin(s.theta) * sp.cos(s.theta),
            sp.Derivative(delta_f, s.r): 2 * s.r - 2 * s.mass,
        }
        opaque = {sigma_f: s.sigma, delta_f: s.delta}

        def derivative(expression, coordinate):
            value = sp.diff(expression, coordinate).subs(derivative_rules)
            return value.subs(opaque)

        inverse = inverse_metric(g.subs(opaque), s, s.delta)
    else:
        g = metric(s)

        def derivative(expression, coordinate):
            return sp.diff(expression, coordinate)

        inverse = inverse_metric(g, s)

    dg = [[[derivative(g[i, j], coordinates[k]) for k in range(4)] for j in range(4)] for i in range(4)]
    u = s.u
    acceleration = []
    for mu in range(4):
        total = sp.S.Zero
        for sigma_index in range(4):
            if inverse[mu, sigma_index] == 0:
                continue
            lower = sp.S.Zero
            for alpha in range(4):
                for beta in range(4):
                    lower += (
                        dg[sigma_index][beta][alpha]
                        - sp.Rational(1, 2) * dg[alpha][beta][sigma_index]
                    ) * u[alpha] * u[beta]
            total += inverse[mu, sigma_index] * lower
        acceleration.append(-total)
    return acceleration


def _reduce_modulo_unit_circle(polynomial, sine, cosine):
    """Remainder of a polynomial modulo ``sine^2 + cosine^2 - 1`` in ``sine``."""
    polynomial = sp.expand(polynomial)
    if not polynomial.has(sine):
        return polynomial
    divisor = sp.Poly(sine**2 + cosine**2 - 1, sine)
    return sp.expand(sp.rem(sp.Poly(polynomial, sine), divisor).as_expr())


def _prepared_factors(expr, s: KerrSymbols, sine, cosine):
    """Rewrite factors of a factored rational function without cancellation.

    ``(cos - 1)^k (cos + 1)^k`` becomes ``(-sin^2)^k``, and the factors
    ``r^2 + a^2 cos^2`` and ``r^2 - 2 M r + a^2`` become the opaque
    ``s.sigma`` and ``s.delta``.
    """
    sigma = s.r**2 + s.spin**2 * cosine**2
    delta = sp.expand(delta_expression(s))

    def rewrite(product):
        powers: dict = {}
        for factor in sp.Mul.make_args(product):
            base, exponent = factor.as_base_exp()
            powers[base] = powers.get(base, 0) + exponent
        pair = min(powers.get(cosine - 1, 0), powers.get(cosine + 1, 0))
        if pair:
            powers[cosine - 1] -= pair
            powers[cosine + 1] -= pair
            powers[sine] = powers.get(sine, 0) + 2 * pair
            powers[sp.S.NegativeOne] = powers.get(sp.S.NegativeOne, 0) + pair
        result = sp.S.One
        for base, exponent in powers.items():
            if base == sigma:
                base = s.sigma
            elif base == delta:
                base = s.delta
            result *= base**exponent
        return result

    numerator, denominator = sp.fraction(expr)
    return rewrite(numerator) / rewrite(denominator)


def simplified_geodesic_acceleration(s: KerrSymbols) -> list:
    """Levi--Civita contraction with simplified, cancellation-free coefficients.

    Starting from :func:`geodesic_acceleration` with explicit ``Sigma`` and
    ``Delta``, every coefficient of a velocity monomial ``u^alpha u^beta`` is
    reduced modulo ``sin^2 + cos^2 = 1``, cancelled and factored; then
    ``(cos - 1)(cos + 1)`` is rewritten as ``-sin^2`` and the factors
    ``Sigma`` and ``Delta`` are replaced by ``s.sigma`` and ``s.delta``, so
    the prepared values are used where they occur as factors.  Numerators
    are polynomials whose leading large-radius terms no longer cancel.

    Returns
    -------
    list
        Four expressions of ``mass, spin, r, theta, u, sigma, delta``.
    """
    sine, cosine = sp.symbols("__kerr_sin __kerr_cos")
    result = []
    for component in geodesic_acceleration(s, prepared=False):
        component = reduce_angle(component, s.theta, sine, cosine)
        numerator, denominator = sp.fraction(sp.together(component))
        numerator = _reduce_modulo_unit_circle(numerator, sine, cosine)
        denominator = _reduce_modulo_unit_circle(denominator, sine, cosine)
        polynomial = sp.Poly(numerator, *s.u)
        terms = []
        for monomial, coefficient in polynomial.terms():
            value = sp.factor(sp.cancel(coefficient.as_expr() / denominator))
            value = _prepared_factors(value, s, sine, cosine)
            value = contract_pythagorean(value, sine, cosine)
            velocity = sp.Mul(*[u**k for u, k in zip(s.u, monomial)])
            terms.append(value * velocity)
        expression = sp.Add(*terms)
        result.append(expression.xreplace({sine: sp.sin(s.theta), cosine: sp.cos(s.theta)}))
    return result
