#ifndef RELATIPY_GEODESIC_SOLUTION_PREVIEW_H
#define RELATIPY_GEODESIC_SOLUTION_PREVIEW_H

#include <stddef.h>
#include "relatipy/kerr_geometry.h"

/*
 * Sample the instantaneous osculating Kepler conic of one canonical state.
 *
 * `spin` is a/M in [0, 1] and `canonical` is the eight-component
 * Boyer--Lindquist state (t, r, theta, phi, u^t, u^r, u^theta, u^phi) in
 * geometric units with M = 1. The conic is computed from the oblate
 * Cartesian position and coordinate velocity of that state. `xyz` is a
 * caller-owned row-major buffer of 3 * count doubles receiving `count`
 * oblate Cartesian points in GM/c^2: one full ellipse for e < 1, or a
 * symmetric arc of the hyperbola for e > 1. `references` receives three
 * Cartesian equatorial radii hypot(r, spin) in GM/c^2 for the outer
 * horizon, prograde ISCO and retrograde ISCO, in that order.
 *
 * Requires 2 <= count <= SIZE_MAX / (3 * sizeof(double)); otherwise
 * returns INVALID_PARAMETER. Returns NULL_POINTER for a NULL buffer, the
 * status of the state conversion or reconstruction on failure,
 * INVALID_PARAMETER when the semilatus rectum is not positive and finite
 * (for example a parabolic or radial conic), and NUMERICAL_RANGE for a
 * non-finite sample. On failure the output buffers may be partially
 * written. No allocation.
 */
rp_kerr_status rp_solution_osculating_preview(
    double spin, const double canonical[8], size_t count,
    double *xyz, double references[3]
);

#endif
