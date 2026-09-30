/**
 * @file convert.h
 * @brief Private Kerr initial-state coordinate and velocity conversion.
 */

#ifndef RELATIPY_NATIVE_GEODESIC_INITIAL_CONVERT_H
#define RELATIPY_NATIVE_GEODESIC_INITIAL_CONVERT_H

#include "relatipy/kerr_geometry.h"

#define RP_INITIAL_CARTESIAN_DIM 7U
#define RP_INITIAL_CANONICAL_DIM 8U

/**
 * Convert (t,x,y,z,vx,vy,vz) in spin-aligned oblate Cartesian coordinates
 * to (t,r,theta,phi,u^t,u^r,u^theta,u^phi) in Boyer--Lindquist coordinates.
 * Velocities are derivatives with respect to coordinate time. Spin is the
 * normalized Kerr parameter in [0,1] and mass is fixed to one.
 * phi is the principal atan2 value; callers unwrap sampled trajectories.
 * Input and output are caller-owned and disjoint. No allocation occurs.
 * A valid output is cleared on every error. The exterior domain requires
 * r > r_+ and excludes the polar chart axis.
 */
rp_kerr_status rp_initial_cartesian_to_canonical(
    double spin,
    const double cartesian[RP_INITIAL_CARTESIAN_DIM],
    double canonical[RP_INITIAL_CANONICAL_DIM]
);

/** Convert (t,R,Theta,Phi,vR,vTheta,vPhi) in Boyer--Lindquist coordinates. */
rp_kerr_status rp_initial_from_bl(
    double spin,
    const double bl[RP_INITIAL_CARTESIAN_DIM],
    double canonical[RP_INITIAL_CANONICAL_DIM]
);

/** Convert Euclidean spherical (t,r,theta,phi,vr,vtheta,vphi). */
rp_kerr_status rp_initial_from_spherical(
    double spin,
    const double spherical[RP_INITIAL_CARTESIAN_DIM],
    double canonical[RP_INITIAL_CANONICAL_DIM]
);

/** Convert Euclidean spherical state to Cartesian in the same frame. */
rp_kerr_status rp_initial_spherical_to_cartesian(
    const double spherical[RP_INITIAL_CARTESIAN_DIM],
    double cartesian[RP_INITIAL_CARTESIAN_DIM]
);

/**
 * Convert spin-aligned Kepler elements (t,a,e,inc,Omega,omega,f).
 * Mass parameter is one. Elliptic (a>0, 0<=e<1) and hyperbolic (a<0,e>1)
 * elements are supported. Parabolic e=1 lacks a finite semimajor axis.
 */
rp_kerr_status rp_initial_from_elements(
    double spin,
    const double elements[RP_INITIAL_CARTESIAN_DIM],
    double canonical[RP_INITIAL_CANONICAL_DIM]
);

/**
 * Convert Kepler elements (t,a,e,inc,Omega,omega,f) to Cartesian position
 * and coordinate-time velocity (t,x,y,z,vx,vy,vz) in the elements' frame.
 * This is the shared osculating-element convention used by Kerr.orbit and
 * the direct observable evaluator. The caller owns both disjoint buffers.
 * A non-null output is cleared on every error.
 */
rp_kerr_status rp_initial_elements_to_cartesian(
    const double elements[RP_INITIAL_CARTESIAN_DIM],
    double cartesian[RP_INITIAL_CARTESIAN_DIM]
);

/** Rotate observer-frame Cartesian position and velocity into the Kerr
 * spin-aligned frame, then form the normalized canonical state. */
rp_kerr_status rp_initial_observer_cartesian_to_canonical(
    double spin,
    const double rotation[3][3],
    const double observer[RP_INITIAL_CARTESIAN_DIM],
    double canonical[RP_INITIAL_CANONICAL_DIM]
);

/**
 * Convert one canonical Kerr state to spin-aligned oblate Cartesian position
 * and coordinate-time velocity (t,x,y,z,vx,vy,vz). u^t must be positive.
 * The exterior chart domain applies. Caller owns disjoint buffers; no
 * allocation occurs. A valid output is cleared on every error.
 */
rp_kerr_status rp_initial_canonical_to_cartesian(
    double spin,
    const double canonical[RP_INITIAL_CANONICAL_DIM],
    double cartesian[RP_INITIAL_CARTESIAN_DIM]
);

#endif /* RELATIPY_NATIVE_GEODESIC_INITIAL_CONVERT_H */
