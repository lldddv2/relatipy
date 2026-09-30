/** Private analytic Jacobian of the corrected affine Kerr RHS. */
#ifndef RELATIPY_NATIVE_GEODESIC_JACOBIAN_H
#define RELATIPY_NATIVE_GEODESIC_JACOBIAN_H

#include "relatipy/kerr_geodesic.h"

/**
 * Fill the 8-by-8 Jacobian in (x,u) order, with fixed row stride 8.
 * Inputs and output are caller owned and mutually disjoint. No allocation
 * or retained pointers. Invalid input clears the output and returns a typed
 * status; the same chart guards as the corrected RHS apply.
 */
rp_kerr_status rp_kerr_geodesic_jacobian(
    double mass,
    double spin,
    const double coordinates[RP_KERR_DIM],
    const double tangent[RP_KERR_DIM],
    double jacobian[2 * RP_KERR_DIM][2 * RP_KERR_DIM]
);

#endif
