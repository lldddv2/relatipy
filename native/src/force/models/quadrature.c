#include "force/models/quadrature.h"

#include <math.h>
#include <stdlib.h>

#define RP_QUADRATURE_PI 3.14159265358979323846

int rp_gauss_legendre(size_t count, rp_gauss_rule *rule)
{
    if (rule == NULL || count < 1U || count > RP_QUADRATURE_MAX_POINTS) {
        return 1;
    }
    const double n = (double)count;
    for (size_t i = 0; i < (count + 1U) / 2U; ++i) {
        double x = cos(RP_QUADRATURE_PI * ((double)i + 0.75) / (n + 0.5));
        double derivative = 1.0;
        for (int iteration = 0; iteration < 100; ++iteration) {
            double p0 = 1.0;
            double p1 = x;
            for (size_t k = 2; k <= count; ++k) {
                const double p2 = ((2.0 * (double)k - 1.0) * x * p1
                                   - ((double)k - 1.0) * p0) / (double)k;
                p0 = p1;
                p1 = p2;
            }
            derivative = n * (x * p1 - p0) / (x * x - 1.0);
            const double step = p1 / derivative;
            x -= step;
            if (fabs(step) <= 1.0e-16) {
                break;
            }
        }
        const double w = 2.0 / ((1.0 - x * x) * derivative * derivative);
        rule->node[i] = -x;
        rule->weight[i] = w;
        rule->node[count - 1U - i] = x;
        rule->weight[count - 1U - i] = w;
    }
    if (count % 2U == 1U) {
        rule->node[count / 2U] = 0.0;
    }
    rule->count = count;
    return 0;
}

static int compare_doubles(const void *a, const void *b)
{
    const double x = *(const double *)a;
    const double y = *(const double *)b;
    return (x > y) - (x < y);
}

size_t rp_quadrature_cluster_breakpoints(
    double lo, double hi, double center, double first, double *out
)
{
    if (out == NULL || !(hi > lo) || !(first > 0.0) || !isfinite(hi - lo)) {
        return 0U;
    }
    size_t count = 0;
    out[count++] = lo;
    out[count++] = hi;
    if (center > lo && center < hi) {
        out[count++] = center;
    }
    for (double offset = first; count + 2U <= RP_QUADRATURE_MAX_BREAKPOINTS; offset *= 2.0) {
        const double below = center - offset;
        const double above = center + offset;
        const int below_inside = below > lo && below < hi;
        const int above_inside = above > lo && above < hi;
        if (below_inside) {
            out[count++] = below;
        }
        if (above_inside) {
            out[count++] = above;
        }
        if (below <= lo && above >= hi) {
            break;
        }
    }
    qsort(out, count, sizeof out[0], compare_doubles);
    return count;
}

double rp_power_coordinate(double x, double beta)
{
    if (beta == 0.0) {
        return log(x);
    }
    return beta == 1.0 ? x : pow(x, beta);
}

double rp_power_coordinate_inverse(double u, double beta)
{
    if (beta == 0.0) {
        return exp(u);
    }
    return beta == 1.0 ? u : pow(u, 1.0 / beta);
}
