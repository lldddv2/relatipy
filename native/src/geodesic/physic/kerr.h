/**
 * @file kerr.h
 * @brief Private boundary for the Kerr geodesic physical evaluator.
 *
 * This header is internal to the native implementation. It is not installed
 * and does not select the production canonical integration state.
 */

#ifndef RELATIPY_NATIVE_GEODESIC_PHYSIC_KERR_H
#define RELATIPY_NATIVE_GEODESIC_PHYSIC_KERR_H

#include "relatipy/kerr_geodesic.h"

/** Validated affine state and reusable Kerr point quantities. */
struct rp_kerr_geodesic_state {
    double x[RP_KERR_DIM];
    double u[RP_KERR_DIM];
    double sin_theta;
    double cos_theta;
    double sigma;
    double delta;
};

/** Complete the physical quantities of a state whose inputs were validated. */
rp_kerr_status rp_kerr_geodesic_physic_prepare_state(
    double mass,
    double spin,
    struct rp_kerr_geodesic_state *state
);

/** Evaluate the optimized affine derivatives at a prepared state. */
void rp_kerr_geodesic_physic_evaluate_rhs(
    double mass,
    double spin,
    const struct rp_kerr_geodesic_state *state,
    double coordinate_derivative[RP_KERR_DIM],
    double four_velocity_derivative[RP_KERR_DIM]
);

/** Construct a null tangent from separated Kerr constants at a prepared point. */
rp_kerr_status rp_kerr_geodesic_physic_evaluate_null_tangent(
    double spin,
    const struct rp_kerr_geodesic_state *state,
    double energy,
    double axial_angular_momentum,
    double carter_constant,
    int radial_direction,
    int polar_direction,
    double tangent[RP_KERR_DIM]
);

/** Evaluate the corrected Levi--Civita affine derivatives. */
void rp_kerr_geodesic_physic_evaluate_corrected_rhs(
    double mass,
    double spin,
    const struct rp_kerr_geodesic_state *state,
    double coordinate_derivative[RP_KERR_DIM],
    double tangent_derivative[RP_KERR_DIM]
);

/**
 * Evaluate the analytic 8-by-8 derivative of the corrected affine RHS.
 *
 * State must already be prepared and validated. The caller owns the output,
 * validates finite results, and handles domain errors. Uses geometric units
 * G = c = 1 with mass and spin length explicit; production sets mass to one.
 * No allocation, retained pointers, finite differences or mutable globals.
 */
void rp_kerr_geodesic_physic_evaluate_corrected_jacobian(
    double mass,
    double spin,
    const struct rp_kerr_geodesic_state *state,
    double jacobian[2 * RP_KERR_DIM][2 * RP_KERR_DIM]
);

#endif /* RELATIPY_NATIVE_GEODESIC_PHYSIC_KERR_H */
