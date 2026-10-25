/**
 * @file thin_disk.h
 * @brief Razor-thin power-law disk (Piran et al. 2026) as an external force.
 *
 * Surface density `Sigma(R) = Sigma_0 (R / R_out)^alpha` for
 * `R_in <= R <= R_out` (default `alpha = -1`, `R_in = 0`).  Each ring of
 * radius `R'` is softened as `|x - x'|^2 -> |x - x'|^2 + H^2` with
 * `H = eps R'` (default `eps = 0.01`), i.e. `z^2 -> z^2 + H^2` in the closed
 * ring field; the acceleration is the exact gradient of the softened
 * potential.  With `p^2 = (d + R')^2 + zeta^2`, `q^2 = (d - R')^2 + zeta^2`,
 * `zeta^2 = z^2 + H^2` and `m = 4 R' d / p^2`, a ring of mass `dM` gives
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
 * Tables are built per unit disk mass with `R_out = 1` by the generic
 * axisymmetric machinery (`axisymmetric_table.h`).  `R_out` (any length in
 * `G M / c^2`), the disk mass and the disk normal are model parameters, so
 * they can change between MCMC proposals without rebuilding the table; only
 * `alpha`, `R_in / R_out` and `eps` fix a table.
 *
 * Grid: nodes in `r / R_out` from `radial_min` to `radial_far`, relative
 * spacing `radial_step`, refined to `floor_step eps R_e` at both edges
 * `R_e`; nodes in `mu = |z| / r` refined to `floor_step eps` at the plane.
 * Below `radial_min` the field is continued as the self-similar power law
 * (`A, C ~ r^(alpha + 2)` for `R_in = 0`) or the hole's harmonic interior
 * (`~ r^3`).  The far field is the exact multipole series of the softened
 * disk to `r^-(far_order + 1)` (Legendre moments plus Gegenbauer softening
 * terms), blended in between `radial_match` and `radial_far`.
 */

#ifndef RELATIPY_FORCE_MODELS_THIN_DISK_H
#define RELATIPY_FORCE_MODELS_THIN_DISK_H

#include "force/models/axisymmetric_table.h"

#include <stddef.h>

typedef struct {
    double alpha;              /**< Sigma ~ R^alpha */
    double inner_ratio;        /**< R_in / R_out in [0, 1); 0 needs alpha > -2 */
    double softening;          /**< eps in (0, 1): H = eps R' */
    size_t quadrature_points;  /**< Gauss-Legendre points per panel, 4..64 */
    double radial_min;         /**< innermost node / R_out */
    double radial_match;       /**< far-field blend start / R_out, > 1 */
    double radial_far;         /**< far-field blend end / R_out */
    double radial_step;        /**< relative radial spacing */
    double angular_step;       /**< relative spacing in mu */
    double floor_step;         /**< spacing at edges and plane, units of eps */
    size_t far_order;          /**< even, <= RP_AXISYMMETRIC_MAX_FAR_ORDER */
} rp_thin_disk_params;

/** One placed disk; see rp_axisymmetric_model.  length_scale is R_out. */
typedef rp_axisymmetric_model rp_thin_disk_model;

/**
 * Defaults: Piran et al. fiducial disk (`alpha = -1`, `R_in = 0`,
 * `eps = 0.01`) and the grid validated by `test_thin_disk_table.c`.
 */
void rp_thin_disk_default_params(rp_thin_disk_params *params);

/**
 * Field of a unit-mass ring of radius `ring_radius >= 0` in `z = 0`,
 * softened by `softening_length >= 0`, at `d >= 0`, `z`.  Writes
 * `field = {a_d / d, a_z / z, Phi}` (finite on the axis and in the plane).
 * Returns nonzero, leaving `field` unchanged, on the unsoftened ring itself
 * or for invalid input.
 */
int rp_thin_disk_ring_field(
    double ring_radius, double softening_length, double d, double z,
    double field[3]
);

/**
 * Direct ring sum of the unit-mass disk with `R_out = 1` at `(d >= 0, z)`:
 * piecewise Gauss-Legendre in the enclosed-mass coordinate with panels
 * refined geometrically around `R' = d` down to a quarter of
 * `sqrt(z^2 + (eps d)^2)`.  Writes `field` as rp_thin_disk_ring_field().
 * This is the reference behind the table; it costs about 0.1 ms.
 */
rp_axisymmetric_status rp_thin_disk_direct_field(
    const rp_thin_disk_params *params, double d, double z, double field[3]
);

/**
 * Fill the generic table specification (grid, inner exponent, far series)
 * of a disk; exposed for tests.
 */
rp_axisymmetric_status rp_thin_disk_table_spec(
    const rp_thin_disk_params *params, rp_axisymmetric_spec *spec
);

/**
 * Build the unit-disk table.  `*table` receives a heap table owned by the
 * caller (release with rp_thin_disk_table_free()), or NULL on failure.
 */
rp_axisymmetric_status rp_thin_disk_table_build(
    const rp_thin_disk_params *params, rp_axisymmetric_table **table
);

/** Release a table from rp_thin_disk_table_build(); NULL is ignored. */
void rp_thin_disk_table_free(rp_axisymmetric_table *table);

/** Disk-frame rotation from the unit disk normal (black-hole frame). */
rp_axisymmetric_status rp_thin_disk_rotation_from_normal(
    const double normal[3], double rotation[3][3]
);

/**
 * Place a disk: borrowed table, mass `M_d` (black-hole masses, any finite
 * value), outer radius `R_out > 0` in `G M / c^2` and unit normal.
 */
rp_axisymmetric_status rp_thin_disk_model_init(
    rp_thin_disk_model *model,
    const rp_axisymmetric_table *table,
    double mass,
    double outer_radius,
    const double normal[3]
);

/** rp_external_acceleration_fn of a placed disk (model: rp_thin_disk_model). */
int rp_thin_disk_acceleration(
    const double position[3],
    const double velocity[3],
    const void *model,
    double acceleration[3]
);

/** Potential of a placed disk, consistent with its acceleration. */
int rp_thin_disk_potential(
    const double position[3], const void *model, double *potential
);

#endif /* RELATIPY_FORCE_MODELS_THIN_DISK_H */
