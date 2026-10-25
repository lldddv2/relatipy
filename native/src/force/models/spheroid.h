/**
 * @file spheroid.h
 * @brief Oblate power-law spheroid (Piran et al. 2026, Sect. 4) as a force.
 *
 * Density `rho = rho_0 s^alpha` for `s_in <= s <= s_out`, constant on the
 * similar spheroids `s^2 = R^2 + z^2 / q^2` (`q = 1` a sphere, `q -> 0` a
 * disk), at a given total mass.  The field is the homoeoid sum of Binney &
 * Tremaine (2008, Sect. 2.5).  With `m^2(u) = R^2 / (1 + u) + z^2 / (q^2 + u)`
 * and, for each shell `s`, `u_s` the root of `m^2(u_s) = s^2` (`u_s = 0` for
 * a shell enclosing the point, `m(0) < s`), the BT integrals over `u`
 * become integrals over the shells:
 *
 *     a_R / R = -2 pi q int rho(s) 2 s k_R(u_s) / |dm^2/du| ds,
 *     a_z / z = -2 pi q int rho(s) 2 s k_z(u_s) / |dm^2/du| ds,
 *     Phi     = -2 pi q int rho(s) s I(u_s) ds,
 *
 * with `k_R = (1+u)^-2 (q^2+u)^-1/2`, `k_z = (1+u)^-1 (q^2+u)^-3/2`, the force
 * integrals over the shells outside-in from the point (`s < m(0)`), and
 * `I(U) = int_U^inf du / ((1+u) sqrt(q^2+u)) = (2/e) atan(e / sqrt(q^2+U))`,
 * `e^2 = 1 - q^2`, over every shell.  The quadrature is Gauss-Legendre in
 * the enclosed-mass coordinate `s^(alpha+3)` with panels clustered at
 * `s = m(0)`; `s_in` and `s_out` are panel ends.
 *
 * Correspondence with the thin disk: as `q -> 0` at fixed mass the spheroid
 * becomes a razor-thin disk with
 * `Sigma(R) = 2 q rho_0 int_0^sqrt(s_out^2 - R^2) (R^2 + t^2)^(alpha/2) dt`,
 * which for `alpha < -1` and `s_in << R << s_out` is proportional to
 * `R^(alpha + 1)`: `alpha = -2` matches the fiducial `Sigma ~ 1/R` disk only
 * there, not at its edge (`Sigma -> 0` as `sqrt(s_out - R)`), and the
 * spheroid has no softening.  The field converges to that disk as `q^2`
 * above the plane.
 *
 * Tables use the generic machinery with flattening `q`, so the density
 * surfaces are grid lines `rho = s`, the grid is refined to the rim's
 * curvature radius `q^2` in `rho` and in `mu`, and the edges are kinks
 * (the field is only C0 across a density jump).  Unit mass with `s_out = 1`;
 * `s_out` (any length), the mass and the symmetry axis are model parameters.
 */

#ifndef RELATIPY_FORCE_MODELS_SPHEROID_H
#define RELATIPY_FORCE_MODELS_SPHEROID_H

#include "force/models/axisymmetric_table.h"

#include <stddef.h>

typedef struct {
    double alpha;              /**< rho ~ s^alpha */
    double axis_ratio;         /**< q in (0, 1] */
    double inner_ratio;        /**< s_in / s_out in [0, 1); 0 needs alpha > -3 */
    size_t quadrature_points;  /**< Gauss-Legendre points per panel, 4..64 */
    double radial_min;         /**< innermost node (in s) / s_out */
    double radial_match;       /**< far-field blend start / s_out, > 1 */
    double radial_far;         /**< far-field blend end / s_out */
    double radial_step;        /**< relative radial spacing */
    double angular_step;       /**< relative spacing in nu */
    /** Spacing at the edges (units of q^2 s_e) and the plane (units of q^2). */
    double floor_step;
    size_t far_order;          /**< even, <= RP_AXISYMMETRIC_MAX_FAR_ORDER */
} rp_spheroid_params;

/** One placed spheroid; see rp_axisymmetric_model.  length_scale is s_out. */
typedef rp_axisymmetric_model rp_spheroid_model;

/** Defaults: `alpha = -2`, `q = 1`, `s_in = 0` and the disk's grid. */
void rp_spheroid_default_params(rp_spheroid_params *params);

/**
 * Direct homoeoid sum of the unit-mass spheroid with `s_out = 1` at
 * `(d >= 0, z)`: writes `field = {a_d / d, a_z / z, Phi}`.
 */
rp_axisymmetric_status rp_spheroid_direct_field(
    const rp_spheroid_params *params, double d, double z, double field[3]
);

/** Generic table specification of a spheroid; exposed for tests. */
rp_axisymmetric_status rp_spheroid_table_spec(
    const rp_spheroid_params *params, rp_axisymmetric_spec *spec
);

/**
 * Build the unit-spheroid table.  `*table` receives a heap table owned by
 * the caller (release with rp_spheroid_table_free()), or NULL on failure.
 */
rp_axisymmetric_status rp_spheroid_table_build(
    const rp_spheroid_params *params, rp_axisymmetric_table **table
);

/** Release a table from rp_spheroid_table_build(); NULL is ignored. */
void rp_spheroid_table_free(rp_axisymmetric_table *table);

/**
 * Place a spheroid: borrowed table, mass (black-hole masses), outer label
 * `s_out > 0` in `G M / c^2` and unit symmetry axis.
 */
rp_axisymmetric_status rp_spheroid_model_init(
    rp_spheroid_model *model,
    const rp_axisymmetric_table *table,
    double mass,
    double outer_radius,
    const double normal[3]
);

/** rp_external_acceleration_fn of a placed spheroid. */
int rp_spheroid_acceleration(
    const double position[3],
    const double velocity[3],
    const void *model,
    double acceleration[3]
);

/** Potential of a placed spheroid. */
int rp_spheroid_potential(
    const double position[3], const void *model, double *potential
);

#endif /* RELATIPY_FORCE_MODELS_SPHEROID_H */
