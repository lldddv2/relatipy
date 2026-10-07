/**
 * @file kerr_null.h
 * @brief Kerr null-geodesic initial states, invariants and coordinate views.
 *
 * Internal RelatiPy interface for the scalar photon API of architecture
 * section 6.13.  It is not a stable installed ABI.
 *
 * Conventions: G = c = M = 1, dimensionless spin `a/M` in `[0, 1]`,
 * signature `(-,+,+,+)`, Boyer--Lindquist (BL) order `(t, r, theta, phi)`,
 * angles in radians.  A null state has eight components
 * `(t, r, theta, phi, k^t, k^r, k^theta, k^phi)` with
 * `k^mu = dx^mu/dlambda` for an internal affine parameter `lambda`.  States
 * built here fix the affine scale with `k^t = 1` at the initial point.
 *
 * Every function is allocation free, retains no pointer, uses no mutable
 * global state and never prints.  All buffers are caller owned; inputs and
 * outputs must not overlap.  A validated output buffer is zeroed on every
 * error path.  `rp_kerr_status` codes come from `kerr_geometry.h`.
 */

#ifndef RELATIPY_KERR_NULL_H
#define RELATIPY_KERR_NULL_H

#include <stddef.h>

#include "relatipy/kerr_geometry.h"

/** Number of components of a null state `(x^mu, k^mu)`. */
#define RP_KERR_NULL_STATE_DIM 8U

/** Number of components of a family input/view `(t, q1, q2, q3, v1, v2, v3)`. */
#define RP_KERR_NULL_FAMILY_DIM 7U

/**
 * Relative horizon margin `eps_h` of section 6.13: a null state is exterior
 * only when `r > r_+ (1 + eps_h)`.
 */
#define RP_KERR_NULL_HORIZON_MARGIN 1.0e-6

/**
 * Coordinate family of a position/velocity row `(t, q1, q2, q3, v1, v2, v3)`
 * where `v_i = d q_i / dt` (coordinate-time derivatives).
 *
 * - CARTESIAN: spin-aligned oblate Cartesian `(x, y, z)`, the same chart as
 *   `rp_initial_cartesian_to_canonical`.
 * - SPHERICAL: Euclidean spherical `(r, theta, phi)` of those Cartesian axes,
 *   as in `rp_initial_spherical_to_cartesian`.
 * - BOYER_LINDQUIST: `(R, Theta, Phi)` Boyer--Lindquist coordinates.
 */
typedef enum rp_kerr_null_family {
    RP_KERR_NULL_FAMILY_CARTESIAN = 0,
    RP_KERR_NULL_FAMILY_SPHERICAL = 1,
    RP_KERR_NULL_FAMILY_BOYER_LINDQUIST = 2
} rp_kerr_null_family;

/** Conserved and diagnostic quantities of one null state. */
typedef struct rp_kerr_null_invariants {
    /** Killing energy `E = -k_t`. */
    double energy;
    /** Axial angular momentum `L_z = k_phi`. */
    double axial_angular_momentum;
    /** Carter constant `Q = k_theta^2 + cos^2(theta) (L_z^2/sin^2(theta) - a^2 E^2)`. */
    double carter_constant;
    /** Bardeen impact parameter `b = L_z / E` (NaN when `E == 0`). */
    double impact_parameter;
    /** Bardeen `eta = Q / E^2` (NaN when `E == 0`). */
    double eta;
    /** Raw norm `g(k, k)`. */
    double norm;
    /**
     * Robust relative norm `|g(k,k)| / sum_{mu,nu} |g_{mu nu} k^mu k^nu|`.
     * The denominator never vanishes for a nonzero tangent, unlike
     * `|g_tt| (k^t)^2` on the ergosurface.
     */
    double relative_norm;
} rp_kerr_null_invariants;

/**
 * Return the exterior threshold `r_+ (1 + RP_KERR_NULL_HORIZON_MARGIN)`.
 * Returns NaN for a spin that is non-finite or outside `[0, 1]`.
 */
double rp_kerr_null_horizon_threshold(double spin);

/**
 * Convert a family row `(t, q1, q2, q3, v1, v2, v3)` to BL position and BL
 * coordinate velocity `(t, r, theta, phi, dr/dt, dtheta/dt, dphi/dt)`.
 *
 * Pure kinematics: no normalization is imposed and the velocity may be zero.
 * The BL azimuth is the principal `atan2` value for Cartesian/spherical
 * input and is copied unchanged for BL input.  Requires the exterior chart
 * `r > r_+` (not the null margin; callers apply
 * `rp_kerr_null_horizon_threshold`) and `0 < theta < pi` away from the
 * polar-axis guard of section 6.9.
 */
rp_kerr_status rp_kerr_null_family_to_bl(
    double spin,
    rp_kerr_null_family family,
    const double row[RP_KERR_NULL_FAMILY_DIM],
    double bl[RP_KERR_NULL_FAMILY_DIM]
);

/**
 * Build a null state from a coordinate direction.
 *
 * `row` is `(t, q1, q2, q3, v1, v2, v3)` in `family`; only the direction of
 * the spatial coordinate velocity is used.  After conversion to a BL
 * direction `n^i`, C solves
 * `g_tt + 2 g_ti n^i s + g_ij n^i n^j s^2 = 0` for `s > 0` and returns
 * `(t, r, theta, phi, 1, s n^r, s n^theta, s n^phi)`.
 *
 * Errors: NULL_POINTER; NONFINITE_INPUT; INVALID_PARAMETER for a spin
 * outside `[0, 1]`, `r <= rp_kerr_null_horizon_threshold(spin)`, a zero
 * direction, or two distinct positive roots (possible only where
 * `g_tt > 0`, inside the ergoregion; the root choice there is not a closed
 * decision); COORDINATE_SINGULARITY near the polar axis;
 * NO_REAL_NULL_TANGENT when no positive real root exists;
 * NUMERICAL_RANGE for non-finite intermediate values.
 */
rp_kerr_status rp_kerr_null_state_from_direction(
    double spin,
    rp_kerr_null_family family,
    const double row[RP_KERR_NULL_FAMILY_DIM],
    double state[RP_KERR_NULL_STATE_DIM]
);

/**
 * Build a null state from Bardeen constants at a BL position.
 *
 * `bl_position` is `(t, r, theta, phi)`.  Uses `rp_kerr_null_tangent` with
 * `E = 1`, `L_z = b`, `Q = eta`, then rescales to `k^t = 1` (b and eta are
 * invariant under that scale).  `radial_sign` and `polar_sign` must be -1 or
 * +1.  Other families convert their position first with
 * `rp_kerr_null_family_to_bl` (velocity columns ignored).
 *
 * Errors: NULL_POINTER; NONFINITE_INPUT; INVALID_PARAMETER for spin, signs
 * or `r <= rp_kerr_null_horizon_threshold(spin)`; COORDINATE_SINGULARITY
 * near the polar axis; NO_REAL_NULL_TANGENT when a separated potential is
 * negative at the point; NUMERICAL_RANGE when the tangent has non-positive
 * or non-finite `k^t`.
 */
rp_kerr_status rp_kerr_null_state_from_constants(
    double spin,
    const double bl_position[RP_KERR_DIM],
    double impact_parameter,
    double eta,
    int radial_sign,
    int polar_sign,
    double state[RP_KERR_NULL_STATE_DIM]
);

/**
 * Evaluate `E`, `L_z`, `Q`, `b`, `eta`, `g(k,k)` and the robust relative
 * norm of one state.  Requires the exterior chart `r > r_+` and the polar
 * guard.  On error every field of `*invariants` is zero.
 */
rp_kerr_status rp_kerr_null_invariants_evaluate(
    double spin,
    const double state[RP_KERR_NULL_STATE_DIM],
    rp_kerr_null_invariants *invariants
);

/**
 * Coordinate views of stored null states.
 *
 * Each input row is a null state; each output row is
 * `(t, q1, q2, q3, dq1/dt, dq2/dt, dq3/dt)` in `family`, with coordinate
 * velocities `k^i / k^t`.  No normalization is imposed or re-solved.  BL
 * rows copy the azimuth branch; Cartesian/spherical rows follow
 * `rp_solution_reconstruct_canonical_family_batch`.  `k^t` must be positive
 * and finite.
 *
 * Every row is cleared first; a failed row stays zero with its code in
 * `row_status`.  Returns the first failed row's code or OK; all rows are
 * processed.  `count == 0` succeeds without touching any buffer.
 */
rp_kerr_status rp_kerr_null_views_batch(
    double spin,
    rp_kerr_null_family family,
    const double *states,
    size_t count,
    double *views,
    rp_kerr_status *row_status
);

#endif /* RELATIPY_KERR_NULL_H */
