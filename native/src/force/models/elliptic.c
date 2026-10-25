#include "force/models/elliptic.h"

#include <math.h>
#include <stddef.h>

#define RP_ELLIPTIC_PI_2 1.57079632679489661923

/* Duplication stops once 4^-n Q < |A_n|; Q uses r = 1e-16 (Carlson 1995). */
#define RP_CARLSON_RF_Q_FACTOR 386.53 /* (3 r)^(-1/6) */
#define RP_CARLSON_RD_Q_FACTOR 538.6  /* (r / 4)^(-1/6) */

static double max3(double a, double b, double c)
{
    return fmax(a, fmax(b, c));
}

double rp_carlson_rf(double x, double y, double z)
{
    if (!(x >= 0.0 && y >= 0.0 && z >= 0.0)
        || !isfinite(x + y + z)
        || (x == 0.0) + (y == 0.0) + (z == 0.0) > 2) {
        return NAN;
    }
    if ((x == 0.0) + (y == 0.0) + (z == 0.0) == 2) {
        return INFINITY;
    }
    const double x0 = x;
    const double y0 = y;
    const double a0 = (x + y + z) / 3.0;
    double a = a0;
    double scale = 1.0; /* 4^-n */
    const double q = RP_CARLSON_RF_Q_FACTOR
        * max3(fabs(a0 - x), fabs(a0 - y), fabs(a0 - z));

    while (scale * q >= fabs(a)) {
        const double sx = sqrt(x);
        const double sy = sqrt(y);
        const double sz = sqrt(z);
        const double lambda = sx * sy + sy * sz + sz * sx;
        x = 0.25 * (x + lambda);
        y = 0.25 * (y + lambda);
        z = 0.25 * (z + lambda);
        a = 0.25 * (a + lambda);
        scale *= 0.25;
    }
    /* x_n - A_n = 4^-n (x_0 - A_0) avoids a late cancellation. */
    const double dx = scale * (a0 - x0) / a;
    const double dy = scale * (a0 - y0) / a;
    const double dz = -(dx + dy);
    const double e2 = dx * dy - dz * dz;
    const double e3 = dx * dy * dz;
    return (1.0 - e2 / 10.0 + e3 / 14.0 + e2 * e2 / 24.0
            - 3.0 * e2 * e3 / 44.0) / sqrt(a);
}

double rp_carlson_rd(double x, double y, double z)
{
    if (!(x >= 0.0 && y >= 0.0 && z > 0.0)
        || !isfinite(x + y + z)
        || (x == 0.0 && y == 0.0)) {
        return NAN;
    }
    const double x0 = x;
    const double y0 = y;
    const double a0 = (x + y + 3.0 * z) / 5.0;
    double a = a0;
    double scale = 1.0; /* 4^-n */
    double sum = 0.0;
    const double q = RP_CARLSON_RD_Q_FACTOR
        * max3(fabs(a0 - x), fabs(a0 - y), fabs(a0 - z));

    while (scale * q >= fabs(a)) {
        const double sx = sqrt(x);
        const double sy = sqrt(y);
        const double sz = sqrt(z);
        const double lambda = sx * sy + sy * sz + sz * sx;
        sum += scale / (sz * (z + lambda));
        x = 0.25 * (x + lambda);
        y = 0.25 * (y + lambda);
        z = 0.25 * (z + lambda);
        a = 0.25 * (a + lambda);
        scale *= 0.25;
    }
    const double dx = scale * (a0 - x0) / a;
    const double dy = scale * (a0 - y0) / a;
    const double dz = -(dx + dy) / 3.0;
    const double xy = dx * dy;
    const double z2 = dz * dz;
    const double e2 = xy - 6.0 * z2;
    const double e3 = (3.0 * xy - 8.0 * z2) * dz;
    const double e4 = 3.0 * (xy - z2) * z2;
    const double e5 = xy * z2 * dz;
    const double series = 1.0 - 3.0 * e2 / 14.0 + e3 / 6.0
        + 9.0 * e2 * e2 / 88.0 - 3.0 * e4 / 22.0
        - 9.0 * e2 * e3 / 52.0 + 3.0 * e5 / 26.0;
    return scale * series / (a * sqrt(a)) + 3.0 * sum;
}

void rp_elliptic_agm(double m, double m1, double *k, double *k_minus_e)
{
    /* a_0 = 1, b_0 = sqrt(m1), c_0^2 = m; K = pi / (2 a_inf) and
     * K - E = K sum_n 2^(n-1) c_n^2, with c_(n+1) = c_n^2 / (4 a_(n+1)). */
    double a = 1.0;
    double b = sqrt(m1);
    double c = sqrt(m);
    double weight = 0.5;
    double sum = 0.5 * m;

    for (int iteration = 0; iteration < 40; ++iteration) {
        const double a_next = 0.5 * (a + b);
        const double c_next = 0.25 * c * c / a_next;
        b = sqrt(a * b);
        a = a_next;
        c = c_next;
        weight *= 2.0;
        const double term = weight * c * c;
        sum += term;
        if (term <= 1.0e-17 * sum && fabs(a - b) <= 1.0e-15 * a) {
            break;
        }
    }
    *k = RP_ELLIPTIC_PI_2 / a;
    *k_minus_e = *k * sum;
}

static rp_elliptic_status complete_integrals(
    double m, double m1, double *k, double *e
)
{
    if (k == NULL && e == NULL) {
        return RP_ELLIPTIC_STATUS_NULL_POINTER;
    }
    if (m1 == 0.0) {
        if (k != NULL) {
            *k = INFINITY;
        }
        if (e != NULL) {
            *e = 1.0;
        }
        return RP_ELLIPTIC_STATUS_OK;
    }
    double kv;
    double ev;
    if (m >= 0.0 && m <= 0.5) {
        double k_minus_e;
        rp_elliptic_agm(m, m1, &kv, &k_minus_e);
        ev = kv - k_minus_e;
    } else {
        /* E = 2 R_G(0, m1, 1) = m1 [R_F(0, m1, 1) + (m/3) R_D(0, 1, m1)]
         * (DLMF 19.21.10 with z = m1 between 0 and 1): for m > 0 every
         * term is positive, so m -> 1 has no cancellation. */
        kv = rp_carlson_rf(0.0, m1, 1.0);
        ev = m >= 0.0
            ? m1 * (kv + m * rp_carlson_rd(0.0, 1.0, m1) / 3.0)
            : kv - m * rp_carlson_rd(0.0, m1, 1.0) / 3.0;
    }
    if (k != NULL) {
        *k = kv;
    }
    if (e != NULL) {
        *e = ev;
    }
    return RP_ELLIPTIC_STATUS_OK;
}

rp_elliptic_status rp_elliptic_ke(double m, double *k, double *e)
{
    if (!(m <= 1.0) || !isfinite(m)) {
        return RP_ELLIPTIC_STATUS_INVALID_ARGUMENT;
    }
    return complete_integrals(m, 1.0 - m, k, e);
}

rp_elliptic_status rp_elliptic_ke_complementary(double m1, double *k, double *e)
{
    if (!(m1 >= 0.0) || !isfinite(m1)) {
        return RP_ELLIPTIC_STATUS_INVALID_ARGUMENT;
    }
    return complete_integrals(1.0 - m1, m1, k, e);
}
