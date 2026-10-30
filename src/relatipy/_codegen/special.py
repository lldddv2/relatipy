"""Complete elliptic integrals and run-time quadratures for generated forces.

:func:`elliptic_k` and :func:`elliptic_e` are the complete elliptic integrals
of the first and second kind in the parameter convention ``K(m)``, ``E(m)``
(DLMF 19.2.8); their derivatives are DLMF 19.4.1.  :func:`integral` is a
definite integral over one variable that :mod:`relatipy._codegen.forces`
evaluates at run time with adaptive Gauss--Kronrod quadrature; it is
differentiated under the integral sign, so a potential defined by an
integral yields an acceleration defined by integrals.

:data:`C_HELPERS` is the C source of the helpers the printed calls use.
"""

from __future__ import annotations

from ._sympy import require_sympy

sp = require_sympy()

__all__ = ["C_HELPERS", "Quad", "elliptic_e", "elliptic_k", "integral"]


class _Elliptic(sp.Function):
    """Function of ``(m, 1 - m)`` printed as a call to a C helper.

    The complement ``1 - m`` is a separate argument so that a caller can
    give it in a form free of cancellation near ``m = 1``.
    """

    nargs = 3
    _c_name = ""
    _derivative = None

    def _ccode(self, printer):
        m, complement, _ = self.args
        return f"{self._c_name}({printer._print(m)}, {printer._print(complement)})"

    def _eval_derivative(self, symbol):
        if self._derivative is None:
            raise NotImplementedError(
                "second derivatives of elliptic integrals are not supported")
        m, complement, explicit = self.args
        # dm = -d(1 - m): differentiate the form the caller wrote free of
        # cancellation when one was given.
        slope = -sp.diff(complement, symbol) if explicit else sp.diff(m, symbol)
        return self._derivative(m, complement, sp.S.Zero) * slope


class elliptic_k_prime(_Elliptic):  # noqa: N801
    """``dK/dm = (E - (1 - m) K) / (2 m (1 - m))`` (DLMF 19.4.1)."""

    _c_name = "rp_jit_ellipk_prime"


class elliptic_e_prime(_Elliptic):  # noqa: N801
    """``dE/dm = (E - K) / (2 m)`` (DLMF 19.4.1)."""

    _c_name = "rp_jit_ellipe_prime"


class _EllipticK(_Elliptic):
    _c_name = "rp_jit_ellipk"
    _derivative = elliptic_k_prime


class _EllipticE(_Elliptic):
    _c_name = "rp_jit_ellipe"
    _derivative = elliptic_e_prime


def _complement(m, complement):
    if complement is None:
        return sp.together(1 - m), sp.S.Zero
    return sp.sympify(complement), sp.S.One


def elliptic_k(m, complement=None):
    """Complete elliptic integral of the first kind ``K(m)`` (DLMF 19.2.8).

    ``complement`` is ``1 - m`` written without cancellation, for example
    ``q**2 / p**2`` for a ring; by default it is ``1 - m`` simplified.
    """
    return _EllipticK(sp.sympify(m), *_complement(m, complement))


def elliptic_e(m, complement=None):
    """Complete elliptic integral of the second kind ``E(m)`` (DLMF 19.2.8)."""
    return _EllipticE(sp.sympify(m), *_complement(m, complement))


class Quad(sp.Function):
    """``Quad(f, s, a, b, c, flags)``: integral of ``f`` over ``s`` from ``a`` to ``b``.

    ``c`` is a point where the integrand may be singular (used when bit 0 of
    ``flags`` is set): the interval is split there and the part symmetric
    about ``c`` is integrated in pairs ``f(c + t) + f(c - t)``, which gives
    the principal value of an odd ``1/(s - c)`` singularity.  Bit 1 of
    ``flags`` integrates in ``log s`` (for ``0 < a < b``).
    """

    nargs = 6

    @classmethod
    def eval(cls, integrand, variable, lower, upper, singular, flags):
        if integrand == 0:
            return sp.S.Zero
        return None

    @property
    def free_symbols(self):
        integrand, variable, *rest = self.args
        return (integrand.free_symbols - {variable}).union(*(item.free_symbols for item in rest))

    def _eval_derivative(self, symbol):
        integrand, variable, lower, upper, singular, flags = self.args
        if lower.has(symbol) or upper.has(symbol):
            raise ValueError("integration limits must not depend on the coordinates")
        return Quad(sp.diff(integrand, symbol), variable, lower, upper, singular, flags)

    def _ccode(self, printer):  # pragma: no cover - replaced before printing
        raise ValueError("integrals are evaluated by the generated quadrature")


def integral(integrand, lower, upper, *, singular=None, log=False):
    """Definite integral evaluated at run time.

    Parameters
    ----------
    integrand : callable
        Called once with the integration variable (a positive SymPy symbol);
        returns the integrand.
    lower, upper
        Limits; expressions of parameters and numbers only.
    singular : expression, optional
        Point where the integrand may be singular, for example the
        cylindrical radius of the field point for a sum of rings.
    log : bool, optional
        Integrate in ``log s``; requires ``0 < lower < upper``.

    Returns
    -------
    Quad
        Symbolic integral understood by the force code generator.
    """
    variable = sp.Dummy("s", positive=True)
    flags = (1 if singular is not None else 0) | (2 if log else 0)
    point = sp.S.Zero if singular is None else sp.sympify(singular)
    return Quad(sp.sympify(integrand(variable)), variable, sp.sympify(lower),
                sp.sympify(upper), point, sp.Integer(flags))


C_HELPERS = r"""
/* Complete elliptic integrals by the AGM (DLMF 19.8.1, 19.8.6). */
static void rp_jit_ellipke(double m, double m1, double *k, double *e)
{
    double a = 1.0, b = sqrt(m1), sum = 0.5 * m, weight = 0.5;
    for (int i = 0; i < 40 && fabs(a - b) > 1e-16 * a; ++i) {
        const double c = 0.5 * (a - b);
        const double next = 0.5 * (a + b);
        b = sqrt(a * b);
        a = next;
        weight *= 2.0;
        sum += weight * c * c;
    }
    *k = 1.5707963267948966 / a;
    *e = *k * (1.0 - sum);
}

static double rp_jit_ellipk(double m, double m1)
{
    double k, e;
    rp_jit_ellipke(m, m1, &k, &e);
    return k;
}

static double rp_jit_ellipe(double m, double m1)
{
    double k, e;
    rp_jit_ellipke(m, m1, &k, &e);
    return e;
}

/* dK/dm and dE/dm (DLMF 19.4.1); Maclaurin series (DLMF 19.5.1, 19.5.2)
 * near m = 0, where the closed forms cancel. */
static double rp_jit_ellipk_prime(double m, double m1)
{
    if (fabs(m) < 1e-4) {
        return 1.5707963267948966
            * (0.25 + m * (9.0 / 32.0 + m * (75.0 / 256.0 + m * (1225.0 / 4096.0))));
    }
    double k, e;
    rp_jit_ellipke(m, m1, &k, &e);
    return (e - m1 * k) / (2.0 * m * m1);
}

static double rp_jit_ellipe_prime(double m, double m1)
{
    if (fabs(m) < 1e-4) {
        return -1.5707963267948966
            * (0.25 + m * (3.0 / 32.0 + m * (15.0 / 256.0 + m * (175.0 / 4096.0))));
    }
    double k, e;
    rp_jit_ellipke(m, m1, &k, &e);
    return (e - k) / (2.0 * m);
}

/* Adaptive Gauss-Kronrod (7, 15) quadrature of a vector integrand. */
#define RP_JIT_QUAD_MAX 16
#define RP_JIT_QUAD_DEPTH 48
#define RP_JIT_QUAD_EPSREL 1e-10
#define RP_JIT_QUAD_BUDGET 2000

typedef void (*rp_jit_integrand)(
    double s, const double *position, const double *velocity,
    const double *parameters, double *out);

typedef struct {
    rp_jit_integrand f;
    int n;
    int log_scale;
    int paired;
    double centre;
    double centre_value;
    const double *position;
    const double *velocity;
    const double *parameters;
    double tolerance;
    double span;
    int budget;
} rp_jit_quad_context;

static const double rp_jit_xk[8] = {
    0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
    0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
    0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
    0.207784955007898467600689403773245, 0.000000000000000000000000000000000,
};
static const double rp_jit_wk[8] = {
    0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
    0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
    0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
    0.204432940075298892414161999234649, 0.209482141084727828012999174891714,
};
static const double rp_jit_wg[4] = {
    0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
    0.381830050505118944950369775488975, 0.417959183673469387755102040816327,
};

static void rp_jit_sample(const rp_jit_quad_context *q, double t, double *out, double *magnitude)
{
    double tmp[RP_JIT_QUAD_MAX];
    int n = q->n;
    for (int i = 0; i < n; ++i) {
        out[i] = 0.0;
        magnitude[i] = 0.0;
    }
    for (int side = 0; side < (q->paired ? 2 : 1); ++side) {
        /* Paired samples are placed about the exact singular point, not
         * about its logarithm, so rounding does not shift the pair. */
        double u = q->paired ? (side ? -t : t) : t;
        double s;
        if (q->paired) {
            s = q->log_scale ? q->centre_value * exp(u) : q->centre_value + u;
        } else {
            s = q->log_scale ? exp(u) : u;
        }
        q->f(s, q->position, q->velocity, q->parameters, tmp);
        for (int i = 0; i < n; ++i) {
            double value = q->log_scale ? tmp[i] * s : tmp[i];
            out[i] += value;
            magnitude[i] += fabs(value);
        }
    }
}

static void rp_jit_gk15(const rp_jit_quad_context *q, double a, double b,
                        double *kronrod, double *error)
{
    double centre = 0.5 * (a + b), half = 0.5 * (b - a);
    double gauss[RP_JIT_QUAD_MAX], f1[RP_JIT_QUAD_MAX], f2[RP_JIT_QUAD_MAX];
    double m1[RP_JIT_QUAD_MAX], m2[RP_JIT_QUAD_MAX], magnitude[RP_JIT_QUAD_MAX];
    int n = q->n;
    rp_jit_sample(q, centre, f1, m1);
    for (int i = 0; i < n; ++i) {
        kronrod[i] = rp_jit_wk[7] * f1[i];
        gauss[i] = rp_jit_wg[3] * f1[i];
        magnitude[i] = rp_jit_wk[7] * m1[i];
    }
    for (int j = 0; j < 7; ++j) {
        double dx = half * rp_jit_xk[j];
        rp_jit_sample(q, centre - dx, f1, m1);
        rp_jit_sample(q, centre + dx, f2, m2);
        for (int i = 0; i < n; ++i) {
            kronrod[i] += rp_jit_wk[j] * (f1[i] + f2[i]);
            magnitude[i] += rp_jit_wk[j] * (m1[i] + m2[i]);
            if (j % 2 == 1) {
                gauss[i] += rp_jit_wg[j / 2] * (f1[i] + f2[i]);
            }
        }
    }
    /* Error estimates below the rounding of the summed magnitudes count
     * as converged (QUADPACK's 50 eps resabs). */
    *error = 0.0;
    for (int i = 0; i < n; ++i) {
        kronrod[i] *= half;
        double e = fabs(kronrod[i] - half * gauss[i]);
        if (e > 50.0 * 2.220446049250313e-16 * half * magnitude[i] && e > *error) {
            *error = e;
        }
    }
}

static void rp_jit_adapt(rp_jit_quad_context *q, double a, double b,
                         const double *estimate, double error, int depth, double *out)
{
    double size = 0.0;
    for (int i = 0; i < q->n; ++i) {
        size = fmax(size, fabs(estimate[i]));
    }
    /* Each interval may carry tolerance / 64: the few tiny intervals at a
     * singular point, where rounding limits the attainable error, stay
     * within the budget without bisecting forever. */
    if (error <= q->tolerance * fmax((b - a) / q->span, 1.0 / 64.0) || error <= 1e-13 * size
        || depth >= RP_JIT_QUAD_DEPTH || --q->budget <= 0
        || b - a <= 1e-6 * (fabs(a) + fabs(b) + fabs(q->centre))) {
        for (int i = 0; i < q->n; ++i) {
            out[i] += estimate[i];
        }
        return;
    }
    double middle = 0.5 * (a + b);
    double left[RP_JIT_QUAD_MAX], right[RP_JIT_QUAD_MAX], left_error, right_error;
    rp_jit_gk15(q, a, middle, left, &left_error);
    rp_jit_gk15(q, middle, b, right, &right_error);
    rp_jit_adapt(q, a, middle, left, left_error, depth + 1, out);
    rp_jit_adapt(q, middle, b, right, right_error, depth + 1, out);
}

static void rp_jit_segment(rp_jit_quad_context *q, double a, double b, double *out)
{
    if (!(b > a)) {
        return;
    }
    double estimate[RP_JIT_QUAD_MAX], error, scale = 0.0;
    rp_jit_gk15(q, a, b, estimate, &error);
    for (int i = 0; i < q->n; ++i) {
        if (fabs(estimate[i]) > scale) {
            scale = fabs(estimate[i]);
        }
    }
    q->tolerance = RP_JIT_QUAD_EPSREL * (scale > 0.0 ? scale : 1.0);
    q->span = b - a;
    q->budget = RP_JIT_QUAD_BUDGET;
    rp_jit_adapt(q, a, b, estimate, error, 0, out);
}

/* Integral over [lower, upper], split at the singular point when it lies
 * inside (pairs symmetric about it); flags as in relatipy._codegen.special. */
static void rp_jit_quad(rp_jit_integrand f, int n, double lower, double upper,
                        double singular, int flags, const double *position,
                        const double *velocity, const double *parameters, double *out)
{
    rp_jit_quad_context q = {f, n, (flags & 2) != 0, 0, 0.0, singular, position, velocity,
                             parameters, 0.0, 1.0, 0};
    for (int i = 0; i < n; ++i) {
        out[i] = 0.0;
    }
    double a = lower, b = upper, c = singular;
    if (q.log_scale) {
        a = log(a);
        b = log(b);
        c = c > 0.0 ? log(c) : -INFINITY;
    }
    if (!(flags & 1) || !(a < c && c < b)) {
        rp_jit_segment(&q, a, b, out);
        return;
    }
    double h = fmin(c - a, b - c);
    rp_jit_segment(&q, a, c - h, out);
    rp_jit_segment(&q, c + h, b, out);
    q.paired = 1;
    q.centre = c;
    rp_jit_segment(&q, 0.0, h, out);
}
"""
