/**
 * @file reconstruct.h
 * @brief Private batch reconstruction of interpolated Kerr states.
 *
 * This interface is internal to RelatiPy. It does not define a stable C ABI.
 */

#ifndef RELATIPY_NATIVE_GEODESIC_SOLUTION_RECONSTRUCT_H
#define RELATIPY_NATIVE_GEODESIC_SOLUTION_RECONSTRUCT_H

#include <stddef.h>

#include "relatipy/kerr_geometry.h"

#define RP_SOLUTION_CARTESIAN_DIM 7U
#define RP_SOLUTION_RECONSTRUCTED_DIM 29U
#define RP_SOLUTION_ELEMENTS_DIM 6U
#define RP_SOLUTION_CONSTANTS_DIM 3U

/** Private, selected output family for canonical reconstruction. */
typedef enum rp_solution_family {
    RP_SOLUTION_FAMILY_CARTESIAN = 0,
    RP_SOLUTION_FAMILY_SPHERICAL = 1,
    RP_SOLUTION_FAMILY_ELEMENTS = 2
} rp_solution_family;

/** Column indices of one reconstructed row. */
typedef enum rp_solution_reconstructed_column {
    RP_SOL_T = 0,
    RP_SOL_BL_R = 1,
    RP_SOL_BL_THETA = 2,
    RP_SOL_BL_PHI = 3,
    RP_SOL_BL_VR = 4,
    RP_SOL_BL_VTHETA = 5,
    RP_SOL_BL_VPHI = 6,
    RP_SOL_UT = 7,
    RP_SOL_BL_UR = 8,
    RP_SOL_BL_UTHETA = 9,
    RP_SOL_BL_UPHI = 10,
    RP_SOL_SPH_R = 11,
    RP_SOL_SPH_THETA = 12,
    RP_SOL_SPH_PHI = 13,
    RP_SOL_SPH_VR = 14,
    RP_SOL_SPH_VTHETA = 15,
    RP_SOL_SPH_VPHI = 16,
    RP_SOL_UX = 17,
    RP_SOL_UY = 18,
    RP_SOL_UZ = 19,
    RP_SOL_SPH_UR = 20,
    RP_SOL_SPH_UTHETA = 21,
    RP_SOL_SPH_UPHI = 22,
    RP_SOL_SEMIMAJOR = 23,
    RP_SOL_ECCENTRICITY = 24,
    RP_SOL_INCLINATION = 25,
    RP_SOL_ASCENDING_NODE = 26,
    RP_SOL_PERIAPSIS_ARGUMENT = 27,
    RP_SOL_TRUE_ANOMALY = 28
} rp_solution_reconstructed_column;

/**
 * Reconstruct Boyer--Lindquist states from interpolated Cartesian rows.
 *
 * `spin` is the normalized Kerr parameter in [0, 1], with mass fixed to 1.
 * Each input row is `(t,x,y,z,vx,vy,vz)`, where velocities are derivatives
 * with respect to coordinate time and Cartesian axes are oblate, right-handed,
 * and aligned with the spin axis. Output columns are specified by
 * `rp_solution_reconstructed_column`: first Boyer--Lindquist coordinates,
 * coordinate velocities, and four-velocity; then Euclidean spherical
 * coordinates and coordinate velocities; Cartesian and spherical spatial
 * four-velocities; then six osculating Kepler elements `(a,e,inc,Omega,omega,f)`.
 * The elements use the instantaneous Cartesian position and coordinate-time
 * velocity with gravitational parameter 1. Exactly zero Newtonian energy
 * gives `a=+INFINITY`; any finite nonzero energy must give finite nonzero `a`.
 * Hyperbolic states have negative `a`. If the angular momentum magnitude is
 * at most `128 * DBL_EPSILON * distance * speed`, all four angles are `NAN`;
 * `a` and `e` are still computed. This includes zero coordinate velocity.
 * Otherwise the existing circular and equatorial angle conventions apply.
 * Undefined elements do not invalidate the physical state: all first 23
 * output columns and `e` remain finite. Only these declared element sentinels
 * are allowed in a successful row; other nonfinite results are range errors.
 * Angles are radians and all other values use normalized geometric units.
 * The azimuth is the principal `atan2` value; the caller may unwrap it.
 *
 * Every row is cleared before evaluation. A failed row remains zero and its
 * code is written to `row_status`. The return value is the first failed row's
 * code or `RP_KERR_STATUS_OK`. All rows are processed, even after a failure.
 * A global invalid spin clears every output row and records the same code in
 * every status. A zero count succeeds without accessing any buffer.
 *
 * All buffers are caller owned, contiguous, and mutually disjoint. This
 * function allocates no memory and retains no pointer. Cartesian input is
 * never modified. A null output or status pointer cannot be cleared.
 */
rp_kerr_status rp_solution_reconstruct_batch(
    double spin,
    const double *cartesian,
    size_t count,
    double *reconstructed,
    rp_kerr_status *row_status
);

/**
 * Reconstruct stored canonical states without renormalizing four-velocity.
 *
 * Input rows are `(t,r,theta,phi,ut,ur,utheta,uphi)`; positive finite ut and
 * the exterior BL chart are required. The supplied x and u are preserved
 * exactly in their output columns, even when integration has accumulated
 * normalization drift. Spatial four-velocities in the other charts use the
 * same supplied ut. Coordinate velocities are u^i/ut. In particular this
 * function does not solve the timelike normalization equation again.
 *
 * The outputs are the same 29-column rows as above, plus 7-column Cartesian
 * coordinate/coordinate-velocity rows. Azimuth retains the input branch.
 * Element sentinel rules, per-row errors, zero-count behavior, and ownership
 * are as above. All buffers are caller owned and mutually disjoint; no
 * allocation occurs. Both output rows are cleared on a row failure.
 */
rp_kerr_status rp_solution_reconstruct_canonical_batch(
    double spin,
    const double *canonical,
    size_t count,
    double *reconstructed,
    double *cartesian,
    rp_kerr_status *row_status
);

/**
 * Reconstruct only one requested family from canonical BL rows.
 * Cartesian and spherical rows have seven columns; element rows have six.
 * The input four-velocity is never renormalized. Failed rows are zeroed;
 * successful element rows retain the documented +INFINITY/NAN sentinels.
 * Buffers are caller-owned, contiguous, and disjoint. A zero count succeeds.
 */
rp_kerr_status rp_solution_reconstruct_canonical_family_batch(
    double spin,
    const double *canonical,
    size_t count,
    rp_solution_family family,
    double *output,
    rp_kerr_status *row_status
);

/**
 * Evaluate specific constants of motion of stored canonical BL rows.
 *
 * Input rows are `(t,r,theta,phi,ut,ur,utheta,uphi)` with mass 1 and spin in
 * [0, 1]. Each output row is `(E, Lz, Q)`: `E = -u_t`, `Lz = u_phi`, and
 * Carter `Q = u_theta^2 + cos^2(theta) (a^2 (1 - E^2) + Lz^2 / sin^2(theta))`
 * with fixed rest mass one. The four-velocity is never renormalized, so the
 * values expose integration drift. Failed rows are zeroed and report their
 * code; the return value is the first failure or `RP_KERR_STATUS_OK`. A zero
 * count succeeds. Buffers are caller owned, contiguous and disjoint; no
 * allocation occurs.
 */
rp_kerr_status rp_solution_constants_of_motion_batch(
    double spin,
    const double *canonical,
    size_t count,
    double *constants,
    rp_kerr_status *row_status
);

#endif /* RELATIPY_NATIVE_GEODESIC_SOLUTION_RECONSTRUCT_H */
