/**
 * @file ring_sum.h
 * @brief Field of a razor-thin axisymmetric disk as a sum of softened rings.
 *
 * Private to the native core.  A model describes only its mass: a surface
 * density `Sigma(R)` on `R_in <= R <= R_out` (lengths in units of `R_out`,
 * any normalization) and the ring softening `H = eps R'`.  This helper sums
 * the closed ring field over the disk, normalizes it to unit mass, computes
 * the far multipoles and fills the generic table specification
 * (`table.h`), so the model never integrates by itself.
 *
 * Each ring of radius `R'` is softened as `|x - x'|^2 -> |x - x'|^2 + H^2`,
 * i.e. `z^2 -> z^2 + H^2` in the closed ring field; the acceleration is the
 * exact gradient of the softened potential.  With `p^2 = (d + R')^2 + zeta^2`,
 * `q^2 = (d - R')^2 + zeta^2`, `zeta^2 = z^2 + H^2` and `m = 4 R' d / p^2`,
 * a ring of mass `dM` gives
 *
 *     Phi = -2 dM K(m) / (pi p),
 *     a_z = -2 dM E(m) z / (pi p q^2),
 *     a_d = -dM d [16 R'^2 S(m) / p^4 + 2 E(m) / q^2] / (pi p),
 *
 * where `S(m) = [K - E - m E / (2 (1 - m))] / m^2` is the standard
 * `a_d = -(dM / (pi p d)) [K + E (d^2 - R'^2 - zeta^2) / q^2]` rewritten
 * without the cancellations at `d -> 0` (axis) and `R' -> 0`; `S` uses its
 * Taylor series for `m < 0.1`.
 *
 * The disk is summed with piecewise Gauss-Legendre in the coordinate
 * `u = R^beta` (`ln R` for `beta = 0`), with panels refined geometrically
 * around `R' = d`.  For `Sigma ~ R^(beta - 2)` the mass is uniform in `u`;
 * any other `Sigma` is weighted by `Sigma R^(2 - beta) / beta`, so `beta`
 * only chooses the coordinate and never changes the integral.  A density
 * that is smooth and finite at the centre wants `beta = 1`; a power-law
 * cusp `R^s` wants `beta = s + 2`.
 */

#ifndef RELATIPY_UTILS_INTEGRATORS_AXISYMMETRIC_RING_SUM_H
#define RELATIPY_UTILS_INTEGRATORS_AXISYMMETRIC_RING_SUM_H

#include "utils/integrators/axisymmetric/table.h"
#include "utils/integrators/quadrature.h"

#include <stddef.h>

/**
 * Surface density at `radius` (units of `R_out`), any normalization.  Must
 * be finite and non-negative on `[R_in, R_out]`.  `context` is the
 * `rp_thin_mass::context` pointer, passed through unchanged.
 */
typedef double (*rp_surface_density_fn)(double radius, const void *context);

/** Mass of one thin disk with `R_out = 1`. */
typedef struct {
    /**
     * Surface density, or NULL for the power law `Sigma ~ R^(beta - 2)`,
     * whose weight in `u` is then exactly constant.
     */
    rp_surface_density_fn sigma;
    const void *context;       /**< passed to `sigma`; borrowed */
    double inner_ratio;        /**< R_in / R_out in [0, 1); 0 needs beta > 0 */
    double softening;          /**< eps in (0, 1): H = eps R' */
    double beta;               /**< integration coordinate u = R^beta */
    /**
     * With `R_in = 0`: `A, C ~ r^inner_exponent` below the table, i.e. the
     * slope of `Sigma` at the centre plus two (`alpha + 2` for
     * `Sigma ~ R^alpha`).  Unused when `R_in > 0` (harmonic interior, 3).
     */
    double inner_exponent;
    size_t quadrature_points;  /**< Gauss-Legendre points per panel, 4..64 */
} rp_thin_mass;

/** Initialized ring sum; read-only after rp_thin_ring_sum_init(). */
typedef struct {
    rp_thin_mass mass;
    double u_inner;            /**< u(R_in) */
    /* With R_in = 0 and a non-integer 2 / beta, R = u^(1/beta) is not smooth
     * at u = 0: the first panel is halved until it holds a mass fraction
     * ~1e-13. */
    double zero_depth;
    double inverse_mass;       /**< 1 / int dM in the units of `sigma` */
    rp_gauss_rule rule;
} rp_thin_ring_sum;

/** Grid and asymptotics requested for a table; lengths in units of R_out. */
typedef struct {
    double radial_min;         /**< innermost node */
    double radial_match;       /**< far-field blend start, > 1 + 4 eps */
    double radial_far;         /**< far-field blend end */
    double radial_step;        /**< relative radial spacing */
    double angular_step;       /**< relative spacing in mu */
    double floor_step;         /**< spacing at edges and plane, units of eps */
    size_t far_order;          /**< even, <= RP_AXISYMMETRIC_MAX_FAR_ORDER */
} rp_axisymmetric_grid;

/**
 * Field of a unit-mass ring of radius `ring_radius >= 0` in `z = 0`,
 * softened by `softening_length >= 0`, at `d >= 0`, `z`.  Writes
 * `field = {a_d / d, a_z / z, Phi}` (finite on the axis and in the plane).
 * Returns nonzero, leaving `field` unchanged, on the unsoftened ring itself
 * or for invalid input.
 */
int rp_ring_field(
    double ring_radius, double softening_length, double d, double z,
    double field[3]
);

/**
 * Validate `mass`, build the Gauss-Legendre rule and integrate the total
 * mass.  On failure `sum` is unspecified.  `mass->context` is borrowed and
 * must outlive every use of `sum`.
 */
rp_axisymmetric_status rp_thin_ring_sum_init(
    const rp_thin_mass *mass, rp_thin_ring_sum *sum
);

/**
 * rp_axisymmetric_field_fn of the unit-mass disk (`context` is a
 * `const rp_thin_ring_sum *`): direct ring sum at `d >= 0`, `z >= 0`.
 * Returns nonzero at the centre of the disk or if `sigma` is not finite.
 */
int rp_thin_ring_sum_field(double d, double z, void *context, double field[3]);

/**
 * Moment `<R^power>` of the unit-mass disk.  Returns NaN if `sigma` is not
 * finite or the integral diverges.
 */
double rp_thin_ring_sum_moment(const rp_thin_ring_sum *sum, double power);

/**
 * Fill the table specification of the disk: grid, edges at `R_out` and
 * `R_in`, plane width, inner exponent (`inner_exponent` for `R_in = 0`, the
 * harmonic `3` otherwise) and the exact multipole series of the softened disk, from
 * the numerical moments `<R^k>`.
 */
rp_axisymmetric_status rp_thin_ring_sum_table_spec(
    const rp_thin_ring_sum *sum, const rp_axisymmetric_grid *grid,
    rp_axisymmetric_spec *spec
);

/**
 * Build the unit-disk table (`R_out = 1`, unit mass).  `*table` receives a
 * heap table owned by the caller (release with rp_axisymmetric_table_free()),
 * or NULL on failure.
 */
rp_axisymmetric_status rp_thin_ring_sum_table_build(
    const rp_thin_mass *mass, const rp_axisymmetric_grid *grid,
    rp_axisymmetric_table **table
);

#endif /* RELATIPY_UTILS_INTEGRATORS_AXISYMMETRIC_RING_SUM_H */
