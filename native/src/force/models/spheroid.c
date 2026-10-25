#include "force/models/spheroid.h"

#include "force/models/quadrature.h"

#include <math.h>
#include <string.h>

#define RP_SPHEROID_PI 3.14159265358979323846

void rp_spheroid_default_params(rp_spheroid_params *params)
{
    if (params == NULL) {
        return;
    }
    params->alpha = -2.0;
    params->axis_ratio = 1.0;
    params->inner_ratio = 0.0;
    params->quadrature_points = 16U;
    params->radial_min = 1.0e-5;
    params->radial_match = 4.0;
    params->radial_far = 5.0;
    params->radial_step = 0.04;
    params->angular_step = 0.02;
    params->floor_step = 0.08;
    params->far_order = 24U;
}

/* int_rho^1 s^beta ds. */
static double power_integral(double beta, double rho)
{
    if (fabs(beta + 1.0) < 1.0e-14) {
        return -log(rho);
    }
    return (1.0 - pow(rho, beta + 1.0)) / (beta + 1.0);
}

typedef struct {
    double q;
    double q2;
    double eccentricity;       /* e = sqrt(1 - q^2) */
    double inner_ratio;
    double alpha;
    double beta;               /* alpha + 3: dM ~ s^(beta - 1) ds */
    double density;            /* rho_0 of the unit-mass spheroid */
    double zero_depth;
    rp_gauss_rule rule;
} spheroid_context;

static rp_axisymmetric_status validate_params(const rp_spheroid_params *params)
{
    if (params == NULL) {
        return RP_AXISYMMETRIC_STATUS_NULL_POINTER;
    }
    if (!isfinite(params->alpha)
        || !(params->axis_ratio > 0.0 && params->axis_ratio <= 1.0)
        || !(params->inner_ratio >= 0.0 && params->inner_ratio < 1.0)
        || (params->inner_ratio == 0.0 && !(params->alpha > -3.0))
        || params->quadrature_points < 4U
        || params->quadrature_points > RP_QUADRATURE_MAX_POINTS
        || !(params->radial_min > 0.0 && params->radial_min < 1.0)
        || !(params->radial_match > 1.0)
        || !(params->radial_far > params->radial_match)
        || !isfinite(params->radial_far)) {
        return RP_AXISYMMETRIC_STATUS_INVALID_ARGUMENT;
    }
    return RP_AXISYMMETRIC_STATUS_OK;
}

static rp_axisymmetric_status context_init(
    const rp_spheroid_params *params, spheroid_context *context
)
{
    const rp_axisymmetric_status status = validate_params(params);
    if (status != RP_AXISYMMETRIC_STATUS_OK) {
        return status;
    }
    context->q = params->axis_ratio;
    context->q2 = params->axis_ratio * params->axis_ratio;
    context->eccentricity = sqrt((1.0 - context->q) * (1.0 + context->q));
    context->inner_ratio = params->inner_ratio;
    context->alpha = params->alpha;
    context->beta = params->alpha + 3.0;
    /* Unit mass: 4 pi q rho_0 int s^(alpha + 2) ds = 1. */
    context->density = 1.0 / (4.0 * RP_SPHEROID_PI * context->q
                              * power_integral(params->alpha + 2.0, params->inner_ratio));
    /* As for the disk: the integrands are smooth in s^2 = sigma^(2/beta). */
    const double gamma = 2.0 / context->beta;
    const double target = 1.0e-13 * pow((double)params->quadrature_points, 2.0 * gamma + 2.0);
    context->zero_depth = context->inner_ratio == 0.0
            && fabs(gamma - round(gamma)) > 1.0e-12 && target < 1.0
        ? fmax(pow(target, 1.0 / context->beta), 0x1p-60) : 0.0;
    if (rp_gauss_legendre(params->quadrature_points, &context->rule) != 0) {
        return RP_AXISYMMETRIC_STATUS_INVALID_ARGUMENT;
    }
    return RP_AXISYMMETRIC_STATUS_OK;
}

/* Root u >= 0 of m^2(u) = s^2 for a shell inside the point (s <= m(0)). */
static double shell_parameter(double s2, double d2, double z2, double q2)
{
    const double b = s2 * (1.0 + q2) - d2 - z2;
    const double c = q2 * (s2 - d2) - z2;   /* <= 0 */
    const double root = sqrt(fmax(b * b - 4.0 * s2 * c, 0.0));
    return b > 0.0 ? -2.0 * c / (b + root) : (root - b) / (2.0 * s2);
}

/* I(U) = (2/e) atan(e / sqrt(q^2 + U)), continuous at e -> 0. */
static double outer_integral(double u, const spheroid_context *context)
{
    const double v = sqrt(context->q2 + u);
    const double t = context->eccentricity / v;
    const double ratio = t < 1.0e-4
        ? 1.0 - t * t / 3.0 + t * t * t * t / 5.0
        : atan(t) / t;
    return 2.0 * ratio / v;
}

static int spheroid_field(double d, double z, void *opaque, double field[3])
{
    const spheroid_context *context = opaque;
    const double q2 = context->q2;
    const double d2 = d * d;
    const double z2 = z * z;
    const double m0 = sqrt(d2 + z2 / q2);
    if (!(m0 > 0.0) || !isfinite(m0)) {
        return 1;
    }
    /* Two near-singularities in the shell label s: the kink at s = m0
     * (shells within ~(q^2 d^2 + z^2/q^2) / (2 m0) of it see u_s < q^2), and
     * the branch points s = (d +- i z) / e where the two roots u_s merge
     * (the analogue of a thin ring's R' = d +- i z). */
    const double width = fmin(m0, (q2 * d2 + z2 / q2) / (2.0 * m0));
    double breaks[RP_QUADRATURE_MAX_BREAKPOINTS];
    size_t nb = rp_quadrature_cluster_breakpoints(
        context->inner_ratio, 1.0, m0, 0.25 * width, breaks);
    if (nb < 2U) {
        return 1;
    }
    const double e = context->eccentricity;
    if (e > 0.0) {
        const double focus = d / e;
        /* Distance from the branch point to the force range s <= m0 (the
         * potential kernel is constant beyond m0). */
        const double reach = focus <= m0 ? z / e : hypot(z / e, focus - m0);
        double second[RP_QUADRATURE_MAX_BREAKPOINTS];
        const size_t n2 = reach > 0.0
            ? rp_quadrature_cluster_breakpoints(context->inner_ratio, 1.0, focus,
                                                0.25 * reach, second)
            : 0U;
        nb = rp_quadrature_merge_breakpoints(breaks, nb, second, n2);
    }
    /* Panels span a ratio of at most two in the mass coordinate. */
    nb = rp_quadrature_refine_geometric(
        breaks, nb, pow(2.0, 1.0 / fmax(1.0, context->beta)), context->zero_depth);

    const double beta = context->beta;
    const rp_gauss_rule *rule = &context->rule;
    double force_d = 0.0, force_z = 0.0, potential = 0.0;
    for (size_t k = 0; k + 1U < nb; ++k) {
        const double lo = breaks[k];
        const double hi = breaks[k + 1U];
        const int outside = hi <= m0;   /* shells inside the point: force */
        const double u_lo = lo == 0.0 ? 0.0 : rp_power_coordinate(lo, beta);
        const double u_hi = rp_power_coordinate(hi, beta);
        const double half = 0.5 * (u_hi - u_lo);
        const double middle = 0.5 * (u_hi + u_lo);
        for (size_t g = 0; g < rule->count; ++g) {
            const double s = rp_power_coordinate_inverse(middle + half * rule->node[g], beta);
            /* ds / dsigma, sigma = s^beta (or ln s). */
            const double jacobian = beta == 0.0 ? s : pow(s, 1.0 - beta) / beta;
            const double density = context->density * pow(s, context->alpha);
            const double w = rule->weight[g] * half * jacobian * density;
            const double s2 = s * s;
            const double u = outside ? shell_parameter(s2, d2, z2, q2) : 0.0;
            potential += w * s * outer_integral(u, context);
            if (outside) {
                const double one = 1.0 + u;
                const double qq = q2 + u;
                const double slope = d2 / (one * one) + z2 / (qq * qq); /* |dm^2/du| */
                const double root = sqrt(qq);
                force_d += w * 2.0 * s / (one * one * root * slope);
                force_z += w * 2.0 * s / (one * qq * root * slope);
            }
        }
    }
    const double factor = -2.0 * RP_SPHEROID_PI * context->q;
    field[0] = factor * force_d;
    field[1] = factor * force_z;
    field[2] = factor * potential;
    return 0;
}

rp_axisymmetric_status rp_spheroid_direct_field(
    const rp_spheroid_params *params, double d, double z, double field[3]
)
{
    if (params == NULL || field == NULL) {
        return RP_AXISYMMETRIC_STATUS_NULL_POINTER;
    }
    if (!(d >= 0.0) || !isfinite(d + z)) {
        return RP_AXISYMMETRIC_STATUS_INVALID_ARGUMENT;
    }
    spheroid_context context;
    const rp_axisymmetric_status status = context_init(params, &context);
    if (status != RP_AXISYMMETRIC_STATUS_OK) {
        return status;
    }
    return spheroid_field(d, fabs(z), &context, field) == 0
        ? RP_AXISYMMETRIC_STATUS_OK : RP_AXISYMMETRIC_STATUS_FIELD_FAILURE;
}

/* ------------------------------------------------------------------ */
/* Far field                                                          */
/* ------------------------------------------------------------------ */

static double legendre(size_t n, double x)
{
    double p0 = 1.0;
    double p1 = x;
    if (n == 0U) {
        return 1.0;
    }
    for (size_t k = 2; k <= n; ++k) {
        const double p2 = ((2.0 * (double)k - 1.0) * x * p1 - ((double)k - 1.0) * p0) / (double)k;
        p0 = p1;
        p1 = p2;
    }
    return p1;
}

static double binomial(size_t n, size_t k)
{
    double value = 1.0;
    for (size_t i = 1; i <= k; ++i) {
        value = value * (double)(n - k + i) / (double)i;
    }
    return value;
}

/*
 * Exterior multipoles of the unit spheroid: Phi = -sum M_k P_k(mu) / r^(k+1)
 * with M_k = int rho r^k P_k dV = 2 pi q int rho s^(k+2) ds * J_k and
 * J_k = int_-1^1 r_1^k P_k(q c / r_1) dc on the unit similar surface
 * (R, z) = (sqrt(1 - c^2), q c), r_1^2 = 1 - c^2 + q^2 c^2.  The integrand is
 * a polynomial of degree k in c, integrated exactly by Gauss-Legendre.
 */
static void set_far_coefficients(
    const rp_spheroid_params *params, rp_axisymmetric_spec *spec
)
{
    memset(spec->far_coefficients, 0, sizeof spec->far_coefficients);
    rp_gauss_rule rule;
    rp_gauss_legendre(32U, &rule);
    const double q = params->axis_ratio;
    const double normalization = power_integral(params->alpha + 2.0, params->inner_ratio);
    for (size_t kk = 0; kk <= spec->far_order / 2U; ++kk) {
        const size_t k = 2U * kk;
        double j = 0.0;
        for (size_t g = 0; g < rule.count; ++g) {
            const double c = rule.node[g];
            const double r2 = 1.0 - c * c + q * q * c * c;
            const double r1 = sqrt(r2);
            j += rule.weight[g] * pow(r2, (double)kk) * legendre(k, q * c / r1);
        }
        const double moment = 0.5 * j
            * power_integral(params->alpha + 2.0 + (double)k, params->inner_ratio) / normalization;
        /* P_k(mu) = 2^-k sum_i (-1)^i C(k, i) C(2k - 2i, k) mu^(k - 2i). */
        for (size_t i = 0; i <= kk; ++i) {
            const double coefficient = ldexp(1.0, -(int)k) * (i % 2U == 0U ? 1.0 : -1.0)
                * binomial(k, i) * binomial(2U * k - 2U * i, k);
            spec->far_coefficients[kk][kk - i] = moment * coefficient;
        }
    }
}

/* ------------------------------------------------------------------ */
/* Table and placed model                                             */
/* ------------------------------------------------------------------ */

rp_axisymmetric_status rp_spheroid_table_spec(
    const rp_spheroid_params *params, rp_axisymmetric_spec *spec
)
{
    if (params == NULL || spec == NULL) {
        return RP_AXISYMMETRIC_STATUS_NULL_POINTER;
    }
    const rp_axisymmetric_status status = validate_params(params);
    if (status != RP_AXISYMMETRIC_STATUS_OK) {
        return status;
    }
    if (params->far_order > RP_AXISYMMETRIC_MAX_FAR_ORDER || params->far_order % 2U != 0U) {
        return RP_AXISYMMETRIC_STATUS_INVALID_ARGUMENT;
    }
    memset(spec, 0, sizeof *spec);
    const double q = params->axis_ratio;
    const double inner = params->inner_ratio;
    spec->radial_min = inner > 0.0
        ? fmin(params->radial_min, 1.0e-3 * inner) : params->radial_min;
    spec->radial_match = params->radial_match;
    spec->radial_far = params->radial_far;
    spec->radial_step = params->radial_step;
    spec->angular_step = params->angular_step;
    spec->floor_step = params->floor_step;
    /* Near a surface the shells at s ~ rho see u_s ~ q^2 over
     * rho - s ~ q^2 s / 2: the radial width at an edge and, through
     * |dm^2/du| = d^2 + z^2 / q^4, the width mu ~ q^2 / 2 at the plane. */
    spec->edge_radius[spec->edge_count] = 1.0;
    spec->edge_width[spec->edge_count] = 0.5 * q * q;
    spec->edge_kink[spec->edge_count++] = 1;
    if (inner > 0.0) {
        spec->edge_radius[spec->edge_count] = inner;
        spec->edge_width[spec->edge_count] = 0.5 * q * q * inner;
        spec->edge_kink[spec->edge_count++] = 1;
    }
    spec->plane_width = 0.5 * q * q;
    spec->flattening = q;
    /* s_in = 0: self-similar rho ~ s^alpha gives A ~ r^(alpha + 3);
     * s_in > 0: empty core, A ~ r^3. */
    spec->inner_exponent = inner > 0.0 ? 3.0 : params->alpha + 3.0;
    spec->far_order = params->far_order;
    set_far_coefficients(params, spec);
    return RP_AXISYMMETRIC_STATUS_OK;
}

rp_axisymmetric_status rp_spheroid_table_build(
    const rp_spheroid_params *params, rp_axisymmetric_table **table
)
{
    if (table == NULL) {
        return RP_AXISYMMETRIC_STATUS_NULL_POINTER;
    }
    *table = NULL;
    rp_axisymmetric_spec spec;
    rp_axisymmetric_status status = rp_spheroid_table_spec(params, &spec);
    if (status != RP_AXISYMMETRIC_STATUS_OK) {
        return status;
    }
    spheroid_context context;
    status = context_init(params, &context);
    if (status != RP_AXISYMMETRIC_STATUS_OK) {
        return status;
    }
    return rp_axisymmetric_table_build(&spec, spheroid_field, &context, table);
}

void rp_spheroid_table_free(rp_axisymmetric_table *table)
{
    rp_axisymmetric_table_free(table);
}

rp_axisymmetric_status rp_spheroid_model_init(
    rp_spheroid_model *model,
    const rp_axisymmetric_table *table,
    double mass,
    double outer_radius,
    const double normal[3]
)
{
    return rp_axisymmetric_model_init(model, table, mass, outer_radius, normal);
}

int rp_spheroid_acceleration(
    const double position[3],
    const double velocity[3],
    const void *model,
    double acceleration[3]
)
{
    return rp_axisymmetric_acceleration(position, velocity, model, acceleration);
}

int rp_spheroid_potential(
    const double position[3], const void *model, double *potential
)
{
    return rp_axisymmetric_potential(position, model, potential);
}
