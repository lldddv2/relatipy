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

size_t rp_quadrature_merge_breakpoints(
    double *breaks, size_t count, const double *other, size_t other_count
)
{
    double merged[RP_QUADRATURE_MAX_BREAKPOINTS];
    size_t n = 0, i = 0, j = 0;
    while ((i < count || j < other_count) && n + 1U < RP_QUADRATURE_MAX_BREAKPOINTS) {
        const double next = j >= other_count || (i < count && breaks[i] <= other[j])
            ? breaks[i++] : other[j++];
        if (n == 0U || next > merged[n - 1U]) {
            merged[n++] = next;
        }
    }
    if (count > 0U && merged[n - 1U] < breaks[count - 1U]) {
        merged[n++] = breaks[count - 1U];
    }
    for (size_t k = 0; k < n; ++k) {
        breaks[k] = merged[k];
    }
    return n;
}

size_t rp_quadrature_refine_geometric(
    double *breaks, size_t count, double ratio, double zero_depth
)
{
    double refined[RP_QUADRATURE_MAX_BREAKPOINTS];
    size_t n = 0;
    if (count < 2U || !(ratio > 1.0)) {
        return count;
    }
    refined[n++] = breaks[0];
    for (size_t k = 0; k + 1U < count; ++k) {
        const double lo = breaks[k];
        const double hi = breaks[k + 1U];
        if (lo == 0.0 && zero_depth > 0.0 && hi > 0.0) {
            /* Halvings towards zero, written in increasing order. */
            size_t halvings = 0;
            for (double x = hi; x > zero_depth * hi && halvings < 200U; x *= 0.5) {
                ++halvings;
            }
            for (size_t h = halvings; h >= 1U; --h) {
                if (n + 2U >= RP_QUADRATURE_MAX_BREAKPOINTS) {
                    break;
                }
                refined[n++] = ldexp(hi, -(int)h);
            }
        } else if (lo > 0.0 && hi > ratio * lo) {
            const size_t pieces = (size_t)ceil(log(hi / lo) / log(ratio));
            const double factor = pow(hi / lo, 1.0 / (double)pieces);
            double x = lo;
            for (size_t p = 1; p < pieces && n + 2U < RP_QUADRATURE_MAX_BREAKPOINTS; ++p) {
                x *= factor;
                refined[n++] = x;
            }
        }
        if (n + 1U >= RP_QUADRATURE_MAX_BREAKPOINTS) {
            break;
        }
        refined[n++] = hi;
    }
    if (refined[n - 1U] != breaks[count - 1U]) {
        refined[n - 1U] = breaks[count - 1U]; /* keep the true end when full */
    }
    for (size_t k = 0; k < n; ++k) {
        breaks[k] = refined[k];
    }
    return n;
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
