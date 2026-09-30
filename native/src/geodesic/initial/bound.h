/**
 * @file bound.h
 * @brief Stable bound Kerr orbital coordinates to canonical initial state.
 */
#ifndef RELATIPY_NATIVE_GEODESIC_INITIAL_BOUND_H
#define RELATIPY_NATIVE_GEODESIC_INITIAL_BOUND_H

#include "relatipy/kerr_geometry.h"

#define RP_INITIAL_BOUND_DIM 7U
#define RP_INITIAL_BOUND_CANONICAL_DIM 8U

/**
 * Convert (t,p,e,x,q_r0,q_theta0,q_phi0) to
 * (t,r,theta,phi,u^t,u^r,u^theta,u^phi), with G=c=M=1.
 * p is measured in GM/c^2; t in GM/c^3; phases are Mino angles in radians.
 * x = sign(L_z) sin(theta_min). q_r=0 is periapsis; q_theta=0 is the
 * northern polar turning point; q_phi0 is the Boyer--Lindquist phi at t.
 * The circular and equatorial limiting phases are redundant. The polar-axis
 * initial point is excluded by the Boyer--Lindquist chart policy.
 * Stability requires the pericenter and third radial root to differ by more
 * than 256*DBL_EPSILON*max(1,pericenter); unresolved near-separatrix inputs
 * are rejected. Exact radial/polar phase endpoints have zero turning speed.
 * Input and output are disjoint caller-owned buffers. No allocation occurs.
 * Output is cleared on every error when non-null.
 */
rp_kerr_status rp_initial_from_bound(
    double spin,
    const double input[RP_INITIAL_BOUND_DIM],
    double canonical[RP_INITIAL_BOUND_CANONICAL_DIM]
);

#endif /* RELATIPY_NATIVE_GEODESIC_INITIAL_BOUND_H */
