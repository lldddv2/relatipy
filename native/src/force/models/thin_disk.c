#include "force/models/thin_disk.h"

#include "force/models/elliptic.h"
#include "force/models/quadrature.h"

#include <math.h>
#include <string.h>

#define RP_THIN_DISK_PI 3.14159265358979323846
#define RP_THIN_DISK_SERIES_LIMIT 0.1

/*
 * Taylor coefficients of S(m) = [K - E - m E / (2 (1 - m))] / m^2
 * (mpmath, 40 digits); 17 terms reach 1e-17 at m = 0.1.
 */
static const double rp_thin_disk_s_series[] = {
    -2.9452431127404311611e-1, -3.6815538909255389513e-1,
    -4.026699568199808228e-1, -4.2280345466097986394e-1,
    -4.3601606261913548469e-1, -4.4535926396097410222e-1,
    -4.5231800246036432257e-1, -4.5770274058489246926e-1,
    -4.6199370377787583616e-1, -4.6549365607922338038e-1,
    -4.684029914297185265e-1, -4.708596504756785887e-1,
    -4.7296170248673072526e-1, -4.7478078595783353574e-1,
    -4.7637045376796021053e-1, -4.7777154333786597585e-1,
    -4.7901574006530833517e-1,
};

void rp_thin_disk_default_params(rp_thin_disk_params *params)
{
    if (params == NULL) {
        return;
    }
    params->alpha = -1.0;
    params->inner_ratio = 0.0;
    params->softening = 0.01;
    params->quadrature_points = 16U;
    params->radial_min = 1.0e-5;
    params->radial_match = 4.0;
    params->radial_far = 5.0;
    params->radial_step = 0.04;
    params->angular_step = 0.02;
    params->floor_step = 0.08;
    params->far_order = 24U;
}

int rp_thin_disk_ring_field(
    double ring_radius, double softening_length, double d, double z,
    double field[3]
)
{
    if (field == NULL || !(ring_radius >= 0.0) || !(softening_length >= 0.0)
        || !(d >= 0.0) || !isfinite(ring_radius + softening_length + d + z)) {
        return 1;
    }
    const double zeta2 = z * z + softening_length * softening_length;
    const double sum = d + ring_radius;
    const double difference = d - ring_radius;
    const double p2 = sum * sum + zeta2;
    const double q2 = difference * difference + zeta2;
    if (!(q2 > 0.0)) {
        return 1;
    }
    const double p = sqrt(p2);
    const double m = 4.0 * ring_radius * d / p2;
    const double m1 = q2 / p2;
    double k, k_minus_e;
    rp_elliptic_agm(m, m1, &k, &k_minus_e);
    const double e = k - k_minus_e;
    double s;
    if (m < RP_THIN_DISK_SERIES_LIMIT) {
        const size_t terms = sizeof rp_thin_disk_s_series / sizeof rp_thin_disk_s_series[0];
        s = 0.0;
        for (size_t i = terms; i-- > 0;) {
            s = s * m + rp_thin_disk_s_series[i];
        }
    } else {
        s = (k_minus_e / m - 0.5 * e / m1) / m;
    }
    const double inverse_pi_p = 1.0 / (RP_THIN_DISK_PI * p);
    field[0] = -inverse_pi_p
        * (16.0 * ring_radius * ring_radius * s / (p2 * p2) + 2.0 * e / q2);
    field[1] = -2.0 * e * inverse_pi_p / q2;
    field[2] = -2.0 * k * inverse_pi_p;
    return 0;
}

/* ------------------------------------------------------------------ */
/* Direct ring sum                                                    */
/* ------------------------------------------------------------------ */

typedef struct {
    double inner_ratio;
    double softening;
    double beta;               /* alpha + 2: dM ~ R^(beta - 1) dR */
    double u_inner;            /* mass coordinate at R_in */
    double inverse_mass_span;  /* 1 / (u(1) - u(R_in)) */
    /* With R_in = 0 and alpha != -1, R = u^(1/beta) is not smooth at u = 0:
     * the first panel is halved until it holds a mass fraction ~1e-13. */
    double zero_depth;
    rp_gauss_rule rule;
} disk_context;

static rp_axisymmetric_status validate_params(const rp_thin_disk_params *params)
{
    if (params == NULL) {
        return RP_AXISYMMETRIC_STATUS_NULL_POINTER;
    }
    if (!isfinite(params->alpha)
        || !(params->inner_ratio >= 0.0 && params->inner_ratio < 1.0)
        || (params->inner_ratio == 0.0 && !(params->alpha > -2.0))
        || !(params->softening > 0.0 && params->softening < 1.0)
        || params->quadrature_points < 4U
        || params->quadrature_points > RP_QUADRATURE_MAX_POINTS
        || !(params->radial_min > 0.0 && params->radial_min < 1.0)
        || !(params->radial_match > 1.0 + 4.0 * params->softening)
        || !(params->radial_far > params->radial_match)
        || !isfinite(params->radial_far)) {
        return RP_AXISYMMETRIC_STATUS_INVALID_ARGUMENT;
    }
    return RP_AXISYMMETRIC_STATUS_OK;
}

static rp_axisymmetric_status disk_context_init(
    const rp_thin_disk_params *params, disk_context *context
)
{
    const rp_axisymmetric_status status = validate_params(params);
    if (status != RP_AXISYMMETRIC_STATUS_OK) {
        return status;
    }
    context->inner_ratio = params->inner_ratio;
    context->softening = params->softening;
    context->beta = params->alpha + 2.0;
    context->u_inner = context->inner_ratio == 0.0
        ? 0.0 : rp_power_coordinate(context->inner_ratio, context->beta);
    context->inverse_mass_span =
        1.0 / (rp_power_coordinate(1.0, context->beta) - context->u_inner);
    /* The ring field is smooth in R^2 = u^gamma, gamma = 2 / beta; a
     * non-integer gamma costs Gauss-Legendre about n^-(2 gamma + 2) of the
     * first panel's mass, so halve it until that is below 1e-13. */
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

static int disk_field(double d, double z, void *opaque, double field[3])
{
    const disk_context *context = opaque;
    const double eps = context->softening;
    double width = sqrt(z * z + eps * eps * d * d);
    if (!(width > 0.0)) {
        return 1; /* the centre of the disk */
    }
    double breaks[RP_QUADRATURE_MAX_BREAKPOINTS];
    size_t nb = rp_quadrature_cluster_breakpoints(
        context->inner_ratio, 1.0, d, 0.25 * width, breaks);
    if (nb < 2U) {
        return 1;
    }
    /* Panels span a ratio of at most two in the mass coordinate. */
    nb = rp_quadrature_refine_geometric(
        breaks, nb, pow(2.0, 1.0 / fmax(1.0, context->beta)), context->zero_depth);
    double total[3] = {0.0, 0.0, 0.0};
    const rp_gauss_rule *rule = &context->rule;
    for (size_t k = 0; k + 1U < nb; ++k) {
        const double u_lo = breaks[k] == 0.0
            ? 0.0 : rp_power_coordinate(breaks[k], context->beta);
        const double u_hi = rp_power_coordinate(breaks[k + 1U], context->beta);
        const double half = 0.5 * (u_hi - u_lo);
        const double middle = 0.5 * (u_hi + u_lo);
        for (size_t q = 0; q < rule->count; ++q) {
            const double radius = rp_power_coordinate_inverse(
                middle + half * rule->node[q], context->beta);
            double ring[3];
            if (rp_thin_disk_ring_field(radius, eps * radius, d, z, ring) != 0) {
                return 1;
            }
            const double w = rule->weight[q] * half * context->inverse_mass_span;
            total[0] += w * ring[0];
            total[1] += w * ring[1];
            total[2] += w * ring[2];
        }
    }
    field[0] = total[0];
    field[1] = total[1];
    field[2] = total[2];
    return 0;
}

rp_axisymmetric_status rp_thin_disk_direct_field(
    const rp_thin_disk_params *params, double d, double z, double field[3]
)
{
    if (params == NULL || field == NULL) {
        return RP_AXISYMMETRIC_STATUS_NULL_POINTER;
    }
    if (!(d >= 0.0) || !isfinite(d + z)) {
        return RP_AXISYMMETRIC_STATUS_INVALID_ARGUMENT;
    }
    disk_context context;
    const rp_axisymmetric_status status = disk_context_init(params, &context);
    if (status != RP_AXISYMMETRIC_STATUS_OK) {
        return status;
    }
    return disk_field(d, fabs(z), &context, field) == 0
        ? RP_AXISYMMETRIC_STATUS_OK : RP_AXISYMMETRIC_STATUS_FIELD_FAILURE;
}

/* ------------------------------------------------------------------ */
/* Far field                                                          */
/* ------------------------------------------------------------------ */

/* int_rho^1 R^beta dR. */
static double power_integral(double beta, double rho)
{
    if (fabs(beta + 1.0) < 1.0e-14) {
        return -log(rho);
    }
    return (1.0 - pow(rho, beta + 1.0)) / (beta + 1.0);
}

static double pochhammer(double x, size_t n)
{
    double product = 1.0;
    for (size_t t = 0; t < n; ++t) {
        product *= x + (double)t;
    }
    return product;
}

static double factorial(size_t n)
{
    return pochhammer(1.0, n);
}

static double binomial(size_t n, size_t k)
{
    return factorial(n) / (factorial(k) * factorial(n - k));
}

/*
 * P = -r Phi = sum_k r^-k G_k(u) for the softened disk: expanding
 * (D^2 + H^2)^(-1/2) = sum_j binom(-1/2, j) H^(2j) D^-(2j+1) with
 * D^-(2j+1) = r^-(2j+1) sum_n C_n^(j+1/2)(cos gamma) (R'/r)^n and averaging
 * the Gegenbauer polynomials over the ring azimuth (cos gamma =
 * sqrt(1 - u) cos psi) gives
 * G_k(u) = <R'^k> sum_j binom(-1/2, j) eps^(2j) <C_(k-2j)^(j+1/2)>(u).
 * The j = 0 terms are the Legendre moments <R'^k> P_k(0) P_k(mu).
 */
static void set_far_coefficients(
    const rp_thin_disk_params *params, rp_axisymmetric_spec *spec
)
{
    memset(spec->far_coefficients, 0, sizeof spec->far_coefficients);
    const double rho = params->inner_ratio;
    const double eps2 = params->softening * params->softening;
    const double normalization = power_integral(params->alpha + 1.0, rho);
    for (size_t kk = 0; kk <= spec->far_order / 2U; ++kk) {
        const size_t k = 2U * kk;
        const double moment =
            power_integral(params->alpha + 1.0 + (double)k, rho) / normalization;
        double binom_half = 1.0;   /* binom(-1/2, j) */
        double eps_power = 1.0;    /* eps^(2j) */
        for (size_t j = 0; j <= kk; ++j) {
            const size_t n = k - 2U * j;
            const double lambda = (double)j + 0.5;
            const double weight = moment * binom_half * eps_power;
            for (size_t i = 0; i <= n / 2U; ++i) {
                const size_t q = n / 2U - i;
                const double term = (i % 2U == 0U ? 1.0 : -1.0)
                    * pochhammer(lambda, n - i)
                    / (factorial(i) * factorial(q) * factorial(q));
                for (size_t l = 0; l <= q; ++l) {
                    spec->far_coefficients[kk][l] += weight * term
                        * binomial(q, l) * (l % 2U == 0U ? 1.0 : -1.0);
                }
            }
            binom_half *= -((double)j + 0.5) / ((double)j + 1.0);
            eps_power *= eps2;
        }
    }
}

/* ------------------------------------------------------------------ */
/* Table and placed model                                             */
/* ------------------------------------------------------------------ */

rp_axisymmetric_status rp_thin_disk_table_spec(
    const rp_thin_disk_params *params, rp_axisymmetric_spec *spec
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
    const double eps = params->softening;
    const double inner = params->inner_ratio;
    spec->radial_min = inner > 0.0
        ? fmin(params->radial_min, 1.0e-3 * inner) : params->radial_min;
    spec->radial_match = params->radial_match;
    spec->radial_far = params->radial_far;
    spec->radial_step = params->radial_step;
    spec->angular_step = params->angular_step;
    spec->floor_step = params->floor_step;
    spec->edge_radius[spec->edge_count] = 1.0;
    spec->edge_width[spec->edge_count++] = eps;
    if (inner > 0.0) {
        spec->edge_radius[spec->edge_count] = inner;
        spec->edge_width[spec->edge_count++] = eps * inner;
    }
    spec->plane_width = eps;
    spec->flattening = 1.0;
    /* R_in = 0: self-similar Sigma ~ R^alpha gives a ~ r^alpha, A ~ r^(alpha+2).
     * R_in > 0: smooth interior, a ~ r, A ~ r^3. */
    spec->inner_exponent = inner > 0.0 ? 3.0 : params->alpha + 2.0;
    spec->far_order = params->far_order;
    set_far_coefficients(params, spec);
    return RP_AXISYMMETRIC_STATUS_OK;
}

rp_axisymmetric_status rp_thin_disk_table_build(
    const rp_thin_disk_params *params, rp_axisymmetric_table **table
)
{
    if (table == NULL) {
        return RP_AXISYMMETRIC_STATUS_NULL_POINTER;
    }
    *table = NULL;
    rp_axisymmetric_spec spec;
    rp_axisymmetric_status status = rp_thin_disk_table_spec(params, &spec);
    if (status != RP_AXISYMMETRIC_STATUS_OK) {
        return status;
    }
    disk_context context;
    status = disk_context_init(params, &context);
    if (status != RP_AXISYMMETRIC_STATUS_OK) {
        return status;
    }
    return rp_axisymmetric_table_build(&spec, disk_field, &context, table);
}

void rp_thin_disk_table_free(rp_axisymmetric_table *table)
{
    rp_axisymmetric_table_free(table);
}

rp_axisymmetric_status rp_thin_disk_rotation_from_normal(
    const double normal[3], double rotation[3][3]
)
{
    return rp_axisymmetric_rotation_from_normal(normal, rotation);
}

rp_axisymmetric_status rp_thin_disk_model_init(
    rp_thin_disk_model *model,
    const rp_axisymmetric_table *table,
    double mass,
    double outer_radius,
    const double normal[3]
)
{
    return rp_axisymmetric_model_init(model, table, mass, outer_radius, normal);
}

int rp_thin_disk_acceleration(
    const double position[3],
    const double velocity[3],
    const void *model,
    double acceleration[3]
)
{
    return rp_axisymmetric_acceleration(position, velocity, model, acceleration);
}

int rp_thin_disk_potential(
    const double position[3], const void *model, double *potential
)
{
    return rp_axisymmetric_potential(position, model, potential);
}
