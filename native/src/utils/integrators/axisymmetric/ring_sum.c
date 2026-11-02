#include "utils/integrators/axisymmetric/ring_sum.h"

#include "utils/integrators/axisymmetric/elliptic.h"

#include <math.h>
#include <string.h>

#define RP_RING_SUM_PI 3.14159265358979323846
#define RP_RING_SUM_SERIES_LIMIT 0.1
/* Points per panel of the mass and moment integrals. */
#define RP_RING_SUM_MOMENT_POINTS 32U

/*
 * Taylor coefficients of S(m) = [K - E - m E / (2 (1 - m))] / m^2
 * (mpmath, 40 digits); 17 terms reach 1e-17 at m = 0.1.
 */
static const double rp_ring_sum_s_series[] = {
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

int rp_ring_field(
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
    if (m < RP_RING_SUM_SERIES_LIMIT) {
        const size_t terms = sizeof rp_ring_sum_s_series / sizeof rp_ring_sum_s_series[0];
        s = 0.0;
        for (size_t i = terms; i-- > 0;) {
            s = s * m + rp_ring_sum_s_series[i];
        }
    } else {
        s = (k_minus_e / m - 0.5 * e / m1) / m;
    }
    const double inverse_pi_p = 1.0 / (RP_RING_SUM_PI * p);
    field[0] = -inverse_pi_p
        * (16.0 * ring_radius * ring_radius * s / (p2 * p2) + 2.0 * e / q2);
    field[1] = -2.0 * e * inverse_pi_p / q2;
    field[2] = -2.0 * k * inverse_pi_p;
    return 0;
}

/* ------------------------------------------------------------------ */
/* Mass                                                               */
/* ------------------------------------------------------------------ */

/* Mass per unit u, up to the normalization: Sigma R dR = Sigma R^(2-beta) du / beta. */
static double mass_weight(const rp_thin_mass *mass, double radius)
{
    if (mass->sigma == NULL) {
        return 1.0;
    }
    const double sigma = mass->sigma(radius, mass->context);
    if (mass->beta == 0.0) {
        return sigma * radius * radius;
    }
    return sigma * pow(radius, 2.0 - mass->beta) / mass->beta;
}

/* int_rho^1 R^beta dR. */
static double power_integral(double beta, double rho)
{
    if (fabs(beta + 1.0) < 1.0e-14) {
        return -log(rho);
    }
    return (1.0 - pow(rho, beta + 1.0)) / (beta + 1.0);
}

/* int R^power dM / inverse_mass over the whole disk, in panels of u. */
static double mass_integral(const rp_thin_ring_sum *sum, double power)
{
    const rp_thin_mass *mass = &sum->mass;
    rp_gauss_rule rule;
    if (rp_gauss_legendre(RP_RING_SUM_MOMENT_POINTS, &rule) != 0) {
        return NAN;
    }
    double breaks[RP_QUADRATURE_MAX_BREAKPOINTS];
    breaks[0] = mass->inner_ratio;
    breaks[1] = 1.0;
    const size_t nb = rp_quadrature_refine_geometric(
        breaks, 2U, pow(2.0, 1.0 / fmax(1.0, mass->beta)), sum->zero_depth);
    double total = 0.0;
    for (size_t k = 0; k + 1U < nb; ++k) {
        const double u_lo = breaks[k] == 0.0
            ? 0.0 : rp_power_coordinate(breaks[k], mass->beta);
        const double u_hi = rp_power_coordinate(breaks[k + 1U], mass->beta);
        const double half = 0.5 * (u_hi - u_lo);
        const double middle = 0.5 * (u_hi + u_lo);
        for (size_t q = 0; q < rule.count; ++q) {
            const double radius = rp_power_coordinate_inverse(
                middle + half * rule.node[q], mass->beta);
            total += rule.weight[q] * half * mass_weight(mass, radius)
                * (power == 0.0 ? 1.0 : pow(radius, power));
        }
    }
    return total;
}

rp_axisymmetric_status rp_thin_ring_sum_init(
    const rp_thin_mass *mass, rp_thin_ring_sum *sum
)
{
    if (mass == NULL || sum == NULL) {
        return RP_AXISYMMETRIC_STATUS_NULL_POINTER;
    }
    if (!isfinite(mass->beta)
        || !(mass->inner_ratio >= 0.0 && mass->inner_ratio < 1.0)
        || (mass->inner_ratio == 0.0
            && !(mass->beta > 0.0 && mass->inner_exponent > 0.0 && isfinite(mass->inner_exponent)))
        || !(mass->softening > 0.0 && mass->softening < 1.0)
        || mass->quadrature_points < 4U
        || mass->quadrature_points > RP_QUADRATURE_MAX_POINTS) {
        return RP_AXISYMMETRIC_STATUS_INVALID_ARGUMENT;
    }
    sum->mass = *mass;
    sum->u_inner = mass->inner_ratio == 0.0
        ? 0.0 : rp_power_coordinate(mass->inner_ratio, mass->beta);
    /* The ring field is smooth in R^2 = u^gamma, gamma = 2 / beta; a
     * non-integer gamma costs Gauss-Legendre about n^-(2 gamma + 2) of the
     * first panel's mass, so halve it until that is below 1e-13. */
    const double gamma = 2.0 / mass->beta;
    const double target = 1.0e-13 * pow((double)mass->quadrature_points, 2.0 * gamma + 2.0);
    sum->zero_depth = mass->inner_ratio == 0.0
            && fabs(gamma - round(gamma)) > 1.0e-12 && target < 1.0
        ? fmax(pow(target, 1.0 / mass->beta), 0x1p-60) : 0.0;
    if (rp_gauss_legendre(mass->quadrature_points, &sum->rule) != 0) {
        return RP_AXISYMMETRIC_STATUS_INVALID_ARGUMENT;
    }
    /* The power law has a uniform mass in u: exact normalization. */
    const double total = mass->sigma == NULL
        ? rp_power_coordinate(1.0, mass->beta) - sum->u_inner
        : mass_integral(sum, 0.0);
    if (!isfinite(total) || total == 0.0) {
        return RP_AXISYMMETRIC_STATUS_INVALID_ARGUMENT;
    }
    sum->inverse_mass = 1.0 / total;
    return RP_AXISYMMETRIC_STATUS_OK;
}

/* ------------------------------------------------------------------ */
/* Direct ring sum                                                    */
/* ------------------------------------------------------------------ */

int rp_thin_ring_sum_field(double d, double z, void *context, double field[3])
{
    const rp_thin_ring_sum *sum = context;
    const rp_thin_mass *mass = &sum->mass;
    const double eps = mass->softening;
    double width = sqrt(z * z + eps * eps * d * d);
    if (!(width > 0.0)) {
        return 1; /* the centre of the disk */
    }
    double breaks[RP_QUADRATURE_MAX_BREAKPOINTS];
    size_t nb = rp_quadrature_cluster_breakpoints(
        mass->inner_ratio, 1.0, d, 0.25 * width, breaks);
    if (nb < 2U) {
        return 1;
    }
    /* Panels span a ratio of at most two in the mass coordinate. */
    nb = rp_quadrature_refine_geometric(
        breaks, nb, pow(2.0, 1.0 / fmax(1.0, mass->beta)), sum->zero_depth);
    double total[3] = {0.0, 0.0, 0.0};
    const rp_gauss_rule *rule = &sum->rule;
    for (size_t k = 0; k + 1U < nb; ++k) {
        const double u_lo = breaks[k] == 0.0
            ? 0.0 : rp_power_coordinate(breaks[k], mass->beta);
        const double u_hi = rp_power_coordinate(breaks[k + 1U], mass->beta);
        const double half = 0.5 * (u_hi - u_lo);
        const double middle = 0.5 * (u_hi + u_lo);
        for (size_t q = 0; q < rule->count; ++q) {
            const double radius = rp_power_coordinate_inverse(
                middle + half * rule->node[q], mass->beta);
            double ring[3];
            if (rp_ring_field(radius, eps * radius, d, z, ring) != 0) {
                return 1;
            }
            double w = rule->weight[q] * half * sum->inverse_mass;
            if (mass->sigma != NULL) {
                w *= mass_weight(mass, radius);
                if (!isfinite(w)) {
                    return 1;
                }
            }
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

double rp_thin_ring_sum_moment(const rp_thin_ring_sum *sum, double power)
{
    if (sum->mass.sigma == NULL) {
        const double alpha = sum->mass.beta - 2.0;
        const double rho = sum->mass.inner_ratio;
        return power_integral(alpha + 1.0 + power, rho) / power_integral(alpha + 1.0, rho);
    }
    return mass_integral(sum, power) * sum->inverse_mass;
}

/* ------------------------------------------------------------------ */
/* Far field and table                                                */
/* ------------------------------------------------------------------ */

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
static int set_far_coefficients(const rp_thin_ring_sum *sum, rp_axisymmetric_spec *spec)
{
    memset(spec->far_coefficients, 0, sizeof spec->far_coefficients);
    const double eps2 = sum->mass.softening * sum->mass.softening;
    for (size_t kk = 0; kk <= spec->far_order / 2U; ++kk) {
        const size_t k = 2U * kk;
        const double moment = rp_thin_ring_sum_moment(sum, (double)k);
        if (!isfinite(moment)) {
            return 1;
        }
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
    return 0;
}

rp_axisymmetric_status rp_thin_ring_sum_table_spec(
    const rp_thin_ring_sum *sum, const rp_axisymmetric_grid *grid,
    rp_axisymmetric_spec *spec
)
{
    if (sum == NULL || grid == NULL || spec == NULL) {
        return RP_AXISYMMETRIC_STATUS_NULL_POINTER;
    }
    const double eps = sum->mass.softening;
    const double inner = sum->mass.inner_ratio;
    if (!(grid->radial_min > 0.0 && grid->radial_min < 1.0)
        || !(grid->radial_match > 1.0 + 4.0 * eps)
        || !(grid->radial_far > grid->radial_match)
        || !isfinite(grid->radial_far)
        || grid->far_order > RP_AXISYMMETRIC_MAX_FAR_ORDER
        || grid->far_order % 2U != 0U) {
        return RP_AXISYMMETRIC_STATUS_INVALID_ARGUMENT;
    }
    memset(spec, 0, sizeof *spec);
    spec->radial_min = inner > 0.0
        ? fmin(grid->radial_min, 1.0e-3 * inner) : grid->radial_min;
    spec->radial_match = grid->radial_match;
    spec->radial_far = grid->radial_far;
    spec->radial_step = grid->radial_step;
    spec->angular_step = grid->angular_step;
    spec->floor_step = grid->floor_step;
    spec->edge_radius[spec->edge_count] = 1.0;
    spec->edge_width[spec->edge_count++] = eps;
    if (inner > 0.0) {
        spec->edge_radius[spec->edge_count] = inner;
        spec->edge_width[spec->edge_count++] = eps * inner;
    }
    spec->plane_width = eps;
    /* R_in = 0: self-similar Sigma ~ R^s gives a ~ r^s, A ~ r^(s + 2).
     * R_in > 0: smooth interior, a ~ r, A ~ r^3. */
    spec->inner_exponent = inner > 0.0 ? 3.0 : sum->mass.inner_exponent;
    spec->far_order = grid->far_order;
    if (set_far_coefficients(sum, spec) != 0) {
        return RP_AXISYMMETRIC_STATUS_FIELD_FAILURE;
    }
    return RP_AXISYMMETRIC_STATUS_OK;
}

rp_axisymmetric_status rp_thin_ring_sum_table_build(
    const rp_thin_mass *mass, const rp_axisymmetric_grid *grid,
    rp_axisymmetric_table **table
)
{
    if (table == NULL) {
        return RP_AXISYMMETRIC_STATUS_NULL_POINTER;
    }
    *table = NULL;
    rp_thin_ring_sum sum;
    rp_axisymmetric_status status = rp_thin_ring_sum_init(mass, &sum);
    if (status != RP_AXISYMMETRIC_STATUS_OK) {
        return status;
    }
    rp_axisymmetric_spec spec;
    status = rp_thin_ring_sum_table_spec(&sum, grid, &spec);
    if (status != RP_AXISYMMETRIC_STATUS_OK) {
        return status;
    }
    return rp_axisymmetric_table_build(&spec, rp_thin_ring_sum_field, &sum, table);
}
