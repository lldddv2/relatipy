/**
 * @file bound.c
 * @brief Native stable bound Kerr initial-state reconstruction.
 *
 * Turning-point equations follow Schmidt (2002), Appendix B. The Mino
 * phase maps follow Fujita and Hikida (2009), equations (27) and (38).
 */
#include "bound.h"
#include "../../utils/numeric.h"

#include <float.h>
#include <math.h>
#include <stddef.h>

static void clear_state(double state[RP_INITIAL_BOUND_CANONICAL_DIM])
{
    size_t i;
    for (i = 0U; i < RP_INITIAL_BOUND_CANONICAL_DIM; ++i) {
        state[i] = 0.0;
    }
}

/* R(r) = A E^2 + B E K + C K^2 - D, L_z = x K. */
static void radial_coefficients(double r, double a, double x, double values[4])
{
    const double aa = a * a;
    const double xx = x * x;
    const double rr = r * r;
    const double delta = rr - 2.0 * r + aa;
    values[0] = (rr + aa) * (rr + aa) - delta * aa * xx;
    values[1] = -4.0 * a * x * r;
    values[2] = aa * xx - delta;
    values[3] = delta * (rr + aa * (1.0 - xx));
}

/* Exact divided differences; also equal derivatives when the radii coincide. */
static void radial_divideds(double r0, double r1, double a, double x,
                            double values[4])
{
    const double aa = a * a;
    const double xx = x * x;
    const double sum = r0 + r1;
    const double cubic = (r0 + r1) * (r0 * r0 + r1 * r1);
    values[0] = cubic + aa * (2.0 - xx) * sum
        + 2.0 * aa * xx;
    values[1] = -4.0 * a * x;
    values[2] = 2.0 - sum;
    values[3] = cubic - 2.0 * (r0 * r0 + r0 * r1 + r1 * r1)
        + aa * (2.0 - xx) * sum - 2.0 * aa * (1.0 - xx);
}

/* Invert the arithmetic-geometric mean amplitude transformation. */
static int jacobi_sn(double u, double m, double *sn, double *complete_k)
{
    double corrections[32];
    double a = 1.0;
    double b;
    double amplitude;
    size_t i;
    size_t count = 0U;

    if (!(m >= 0.0 && m < 1.0) || !isfinite(u)) {
        return 0;
    }
    b = sqrt(1.0 - m);
    for (i = 0U; i < 32U; ++i) {
        const double arithmetic = 0.5 * (a + b);
        const double geometric = sqrt(a * b);
        corrections[i] = (a - b) / (a + b);
        a = arithmetic;
        b = geometric;
        count = i + 1U;
        if (fabs(a - b) <= 8.0 * DBL_EPSILON * a) {
            break;
        }
    }
    if (count == 32U && fabs(a - b) > 8.0 * DBL_EPSILON * a) {
        return 0;
    }
    *complete_k = acos(-1.0) / (2.0 * a);
    amplitude = ldexp(a * u, (int)count);
    for (i = count; i-- > 0U;) {
        amplitude = 0.5 * (amplitude
            + asin(corrections[i] * sin(amplitude)));
    }
    *sn = sin(amplitude);
    return isfinite(*sn) && isfinite(*complete_k);
}

static int complete_elliptic_k(double m, double *value)
{
    double unused;
    return jacobi_sn(0.0, m, &unused, value);
}

static int phase_sn(double q, double m, int polar, double *sn)
{
    double k;
    double unused;
    double phase = remainder(q, 2.0 * acos(-1.0));
    if (!complete_elliptic_k(m, &k)) {
        return 0;
    }
    if (polar) {
        phase += 0.5 * acos(-1.0);
        return jacobi_sn(2.0 * k * phase / acos(-1.0), m, sn, &unused);
    }
    return jacobi_sn(k * phase / acos(-1.0), m, sn, &unused);
}

static int constants_from_turning_points(
    double a, double p, double e, double x,
    double *energy, double *angular, double *carter,
    double *r1, double *r2, double *r3, double *r4
)
{
    double inner[4];
    double outer[4];
    double qa;
    double qb;
    double qc;
    double discriminant;
    double root;
    double ratios[2];
    double chosen_e = 0.0;
    double chosen_k = 0.0;
    double chosen_q = 0.0;
    double chosen_r3 = 0.0;
    double chosen_r4 = 0.0;
    int accepted = 0;
    size_t j;

    *r2 = p / (1.0 + e);
    *r1 = p / (1.0 - e);
    radial_coefficients(*r2, a, x, inner);
    radial_divideds(*r2, *r1, a, x, outer);
    qa = inner[3] * outer[0] - outer[3] * inner[0];
    qb = inner[3] * outer[1] - outer[3] * inner[1];
    qc = inner[3] * outer[2] - outer[3] * inner[2];
    discriminant = qb * qb - 4.0 * qa * qc;
    if (!(discriminant >= 0.0) || !isfinite(discriminant) || qc == 0.0) {
        return 0;
    }
    root = -0.5 * (qb + copysign(sqrt(discriminant), qb));
    if (root == 0.0) {
        return 0;
    }
    ratios[0] = root / qc;
    ratios[1] = qa / root;
    for (j = 0U; j < 2U; ++j) {
        double ratio = ratios[j];
        double denominator;
        double esq;
        double k;
        double q;
        double sum;
        double product;
        double root_discriminant;
        double third;
        double fourth;
        if (!(ratio > 0.0) || !isfinite(ratio)) {
            continue;
        }
        denominator = inner[0] + inner[1] * ratio
            + inner[2] * ratio * ratio;
        esq = inner[3] / denominator;
        if (!(esq > 0.0 && esq < 1.0) || !isfinite(esq)) {
            continue;
        }
        k = ratio * sqrt(esq);
        q = (1.0 - x * x) * (k * k + a * a * (1.0 - esq));
        sum = 2.0 / (1.0 - esq) - *r1 - *r2;
        product = a * a * q / ((1.0 - esq) * (*r1) * (*r2));
        root_discriminant = sum * sum - 4.0 * product;
        if (!(root_discriminant >= 0.0) || !isfinite(root_discriminant)) {
            continue;
        }
        third = 0.5 * (sum + sqrt(root_discriminant));
        fourth = third == 0.0 ? 0.0 : product / third;
        if (!(*r2 - third > 256.0 * DBL_EPSILON * fmax(1.0, *r2)
              && fourth <= third && fourth >= 0.0)) {
            continue;
        }
        chosen_e = sqrt(esq);
        chosen_k = k;
        chosen_q = q;
        chosen_r3 = third;
        chosen_r4 = fourth;
        ++accepted;
    }
    if (accepted != 1 || !isfinite(chosen_q)) {
        return 0;
    }
    *energy = chosen_e;
    *angular = x * chosen_k;
    *carter = chosen_q;
    *r3 = chosen_r3;
    *r4 = chosen_r4;
    return 1;
}

rp_kerr_status rp_initial_from_bound(
    double spin,
    const double input[RP_INITIAL_BOUND_DIM],
    double canonical[RP_INITIAL_BOUND_CANONICAL_DIM]
)
{
    double p;
    double e;
    double x;
    double energy;
    double angular;
    double carter;
    double r1;
    double r2;
    double r3;
    double r4;
    double radius;
    double theta;
    double sn;
    double sn_theta;
    double m_radial;
    double m_polar;
    double sigma;
    double delta;
    double potential_p;
    double z;
    double radial_potential;
    double polar_potential;
    double sign_radial;
    double sign_polar;
    double qr;
    double qt;
    const double pi = acos(-1.0);
    int radial_turn;
    int polar_turn;
    size_t i;

    if (canonical == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    clear_state(canonical);
    if (input == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    if (!isfinite(spin)) {
        return RP_KERR_STATUS_NONFINITE_INPUT;
    }
    if (!(spin >= 0.0 && spin <= 1.0)) {
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    for (i = 0U; i < RP_INITIAL_BOUND_DIM; ++i) {
        if (!isfinite(input[i])) {
            return RP_KERR_STATUS_NONFINITE_INPUT;
        }
    }
    p = input[1];
    e = input[2];
    x = input[3];
    qr = remainder(input[4], 2.0 * pi);
    qt = remainder(input[5], 2.0 * pi);
    radial_turn = e == 0.0 || qr == 0.0 || fabs(qr) == pi;
    polar_turn = fabs(x) == 1.0 || qt == 0.0 || fabs(qt) == pi;
    if (!(p > 0.0 && e >= 0.0 && e < 1.0 && x >= -1.0 && x <= 1.0)) {
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    if (!constants_from_turning_points(
            spin, p, e, x, &energy, &angular, &carter,
            &r1, &r2, &r3, &r4)) {
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    if (!(r2 > 1.0 + sqrt((1.0 - spin) * (1.0 + spin)))) {
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    if (e == 0.0) {
        radius = p;
    } else if (qr == 0.0) {
        radius = r2;
    } else if (fabs(qr) == pi) {
        radius = r1;
    } else {
        m_radial = (r1 - r2) * (r3 - r4)
            / ((r1 - r3) * (r2 - r4));
        if (!phase_sn(input[4], m_radial, 0, &sn)) {
            return RP_KERR_STATUS_NUMERICAL_RANGE;
        }
        z = sn * sn;
        radius = (r3 * (r1 - r2) * z - r2 * (r1 - r3))
            / ((r1 - r2) * z - (r1 - r3));
    }
    if (x == 1.0 || x == -1.0) {
        theta = 0.5 * pi;
    } else if (qt == 0.0 || fabs(qt) == pi) {
        theta = atan2(fabs(x), sqrt((1.0 - x) * (1.0 + x)));
        if (fabs(qt) == pi) {
            theta = pi - theta;
        }
    } else {
        m_polar = (1.0 - x * x) * (1.0 - x * x) * spin * spin
            * (1.0 - energy * energy) / carter;
        if (!phase_sn(input[5], m_polar, 1, &sn_theta)) {
            return RP_KERR_STATUS_NUMERICAL_RANGE;
        }
        theta = acos(sqrt(1.0 - x * x) * sn_theta);
    }
    if (rp_bl_polar_axis_singular(theta)) {
        return RP_KERR_STATUS_COORDINATE_SINGULARITY;
    }
    sigma = radius * radius + spin * spin * cos(theta) * cos(theta);
    delta = radius * radius - 2.0 * radius + spin * spin;
    potential_p = energy * (radius * radius + spin * spin) - spin * angular;
    z = cos(theta) * cos(theta);
    radial_potential = (1.0 - energy * energy) * (r1 - radius)
        * (radius - r2) * (radius - r3) * (radius - r4);
    polar_potential = (1.0 - x * x - z)
        * (spin * spin * (1.0 - energy * energy)
            + (x == 0.0 ? carter - spin * spin * (1.0 - energy * energy)
                         : angular * angular / (x * x)) / (1.0 - z));
    /* Exact phase endpoints are exact turning points. Avoid amplifying
       roundoff in the potential through its square root; never snap phases
       merely close to an endpoint. */
    if (radial_turn) {
        radial_potential = 0.0;
    }
    if (polar_turn) {
        polar_potential = 0.0;
    }
    sign_radial = sin(qr);
    sign_polar = sin(qt);
    if (radial_potential < 0.0 && radial_potential > -1e-12 * radius * radius) {
        radial_potential = 0.0;
    }
    if (polar_potential < 0.0 && polar_potential > -1e-12 * (1.0 + carter)) {
        polar_potential = 0.0;
    }
    if (!(radial_potential >= 0.0 && polar_potential >= 0.0
            && delta > 0.0 && sigma > 0.0)) {
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    canonical[0] = input[0];
    canonical[1] = radius;
    canonical[2] = theta;
    canonical[3] = input[6];
    canonical[4] = ((radius * radius + spin * spin) * potential_p / delta
        + spin * (angular - spin * energy * (1.0 - z))) / sigma;
    canonical[5] = copysign(sqrt(radial_potential) / sigma, sign_radial);
    canonical[6] = copysign(sqrt(polar_potential) / sigma, sign_polar);
    canonical[7] = (spin * potential_p / delta
        + angular / (1.0 - z) - spin * energy) / sigma;
    for (i = 0U; i < RP_INITIAL_BOUND_CANONICAL_DIM; ++i) {
        if (!isfinite(canonical[i])) {
            clear_state(canonical);
            return RP_KERR_STATUS_NUMERICAL_RANGE;
        }
    }
    if (!(canonical[4] > 0.0)) {
        clear_state(canonical);
        return RP_KERR_STATUS_NON_TIMELIKE_VELOCITY;
    }
    return RP_KERR_STATUS_OK;
}
