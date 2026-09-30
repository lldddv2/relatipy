#ifndef RELATIPY_NATIVE_METRIC_KERR_PROPERTIES_H
#define RELATIPY_NATIVE_METRIC_KERR_PROPERTIES_H

#include <stddef.h>
#include "relatipy/kerr_geometry.h"

/* All radii are Boyer--Lindquist radii in units GM/c^2; spin is a/M in [0, 1]. */

/*
 * Characteristic equatorial radii of a Kerr black hole.
 *
 * Writes six caller-owned doubles in this order: [0] outer horizon r_+,
 * [1] inner horizon r_-, [2] prograde ISCO, [3] retrograde ISCO,
 * [4] prograde circular photon orbit, [5] retrograde circular photon orbit.
 * Returns NULL_POINTER for radii == NULL and INVALID_PARAMETER for a
 * non-finite spin or one outside [0, 1]. The buffer is zeroed before
 * validation, so every failure with a valid buffer leaves it zero.
 * No allocation.
 */
rp_kerr_status rp_kerr_characteristic_radii(double spin, double radii[6]);

/*
 * Outer ergosurface radius r_E(theta) = 1 + sqrt(1 - spin^2 cos^2 theta).
 *
 * `theta` and `radii` are caller-owned arrays of `count` doubles and should
 * not overlap. count == 0 is valid and touches no buffer. Returns
 * NULL_POINTER when count > 0 and a buffer is NULL, and INVALID_PARAMETER
 * for a spin outside [0, 1] or a non-finite theta outside [0, pi]. Angles
 * are validated one at a time, so on INVALID_PARAMETER the entries before
 * the offending angle may already be written. No allocation.
 */
rp_kerr_status rp_kerr_ergosurface_radii(
    double spin, const double *theta, size_t count, double *radii
);

/* Reference surfaces of the Kerr exterior. */
typedef enum rp_kerr_surface {
    RP_KERR_SURFACE_OUTER_HORIZON = 0,
    RP_KERR_SURFACE_ERGOSURFACE = 1
} rp_kerr_surface;

/*
 * Cartesian meridional profile of a Kerr reference surface.
 *
 * For each polar angle theta[i] in [0, pi], evaluates the Boyer--Lindquist
 * radius r of `surface` (outer horizon r_+ or outer ergosurface r_E(theta))
 * and maps it with the core Cartesian convention:
 *     rho[i] = hypot(r, spin) * sin(theta[i]),  z[i] = r * cos(theta[i]).
 * All lengths are in GM/c^2. The caller owns `theta`, `rho` and `z`; each
 * holds `count` doubles and `rho`/`z` must not overlap `theta` or each other.
 * count == 0 is valid and touches no buffer. Returns NULL_POINTER when
 * count > 0 and a buffer is NULL, and INVALID_PARAMETER for a spin outside
 * [0, 1], an unknown surface, or a non-finite theta outside [0, pi]. All
 * angles are validated before any output is written. No allocation.
 */
rp_kerr_status rp_kerr_surface_profile(
    double spin, rp_kerr_surface surface, const double *theta, size_t count,
    double *rho, double *z
);

#endif
