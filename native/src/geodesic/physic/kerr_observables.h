/**
 * @file kerr_observables.h
 * @brief Private geometric projection of one Boyer--Lindquist stellar state.
 *
 * This is a straight-line, small-angle observation approximation. It does not
 * trace a photon, include Shapiro delay, or define the future MCMC API.
 */

#ifndef RELATIPY_NATIVE_GEODESIC_PHYSIC_KERR_OBSERVABLES_H
#define RELATIPY_NATIVE_GEODESIC_PHYSIC_KERR_OBSERVABLES_H

#include "relatipy/kerr_geometry.h"

#define RP_KERR_OBSERVABLE_STATE_DIM 8U
#define RP_KERR_OBSERVABLE_OUTPUT_DIM 4U

/**
 * Project one state into observer-frame geometric quantities.
 *
 * `state = (t,r,theta,phi,u^t,u^r,u^theta,u^phi)` uses geometrized units.
 * `spin` is the nonnegative Kerr spin length in the same length unit as `r`.
 * `rotation` maps oblate Cartesian coordinates aligned with the spin axis to
 * observer coordinates `(X,Y,Z)`, with positive Z pointing away from the
 * observer. It must be a proper orthogonal matrix (determinant +1).
 * `angle_scale` is a caller-supplied, positive angular conversion factor per
 * geometric length, such as `r_G / D_s` in radians per geometric length.
 *
 * On success, `output = (t + Z, Y * angle_scale, X * angle_scale,
 * u^t + dZ/dtau - 1)`. These are dimensionless arrival time, right-ascension
 * offset, declination offset, and redshift under the same straight-line
 * propagation approximation. The caller converts time, angle, and redshift
 * to presentation units and applies reference-frame offsets.
 *
 * All buffers are caller owned. This function allocates no memory, retains no
 * pointer, and clears a non-null output on every failure. Output must not
 * overlap either input buffer.
 */
rp_kerr_status rp_kerr_observables_project(
    const double state[RP_KERR_OBSERVABLE_STATE_DIM],
    double spin,
    const double rotation[3][3],
    double angle_scale,
    double output[RP_KERR_OBSERVABLE_OUTPUT_DIM]
);

#endif /* RELATIPY_NATIVE_GEODESIC_PHYSIC_KERR_OBSERVABLES_H */
