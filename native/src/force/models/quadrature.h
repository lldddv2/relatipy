/**
 * @file quadrature.h
 * @brief Gauss-Legendre panels for the reference fields of the force models.
 *
 * Private to the native force models: no allocation, no retained pointer.
 */

#ifndef RELATIPY_FORCE_MODELS_QUADRATURE_H
#define RELATIPY_FORCE_MODELS_QUADRATURE_H

#include <stddef.h>

#define RP_QUADRATURE_MAX_POINTS 64U
/** Capacity that rp_quadrature_cluster_breakpoints() never exceeds. */
#define RP_QUADRATURE_MAX_BREAKPOINTS 256U

/** Nodes and weights on [-1, 1]. */
typedef struct {
    size_t count;
    double node[RP_QUADRATURE_MAX_POINTS];
    double weight[RP_QUADRATURE_MAX_POINTS];
} rp_gauss_rule;

/**
 * Fill `rule` with the `count`-point Gauss-Legendre rule (Newton on the
 * Legendre recurrence, nodes to about 1e-16).  Returns nonzero for
 * `count < 1` or `count > RP_QUADRATURE_MAX_POINTS`.
 */
int rp_gauss_legendre(size_t count, rp_gauss_rule *rule);

/**
 * Sorted breakpoints on `[lo, hi]`: both ends, `center` when inside, and
 * `center +- first 2^k` inside the interval, so that panels grow
 * geometrically away from a near-singularity of width about `4 first` at
 * `center`.  `out` holds RP_QUADRATURE_MAX_BREAKPOINTS doubles; returns
 * the count (>= 2), or 0 for invalid input.
 */
size_t rp_quadrature_cluster_breakpoints(
    double lo, double hi, double center, double first, double *out
);

/**
 * Enclosed-mass coordinate of a power law `dM ~ x^(beta - 1) dx`:
 * `u = x^beta` (or `ln x` when `beta = 0`), so that the mass is uniform
 * in `u` and the quadrature has no endpoint singularity at `x = 0`.
 */
double rp_power_coordinate(double x, double beta);

/** Inverse of rp_power_coordinate(). */
double rp_power_coordinate_inverse(double u, double beta);

#endif /* RELATIPY_FORCE_MODELS_QUADRATURE_H */
