/**
 * @file external.c
 * @brief Frame map, coordinate velocity and projected external acceleration.
 *
 * The prepared-state functions reuse `sin(theta)`, `cos(theta)`, `Sigma` and
 * `Delta` from `struct rp_kerr_geodesic_state`.  The inverse metric is the
 * closed-form Boyer--Lindquist expression with mass one, the same corrected
 * metric inverted in `geodesic/physic/kerr.c`; only the four components that
 * meet a covector with `F_t = 0` are formed.
 */

#include "external.h"

#include "../utils/tensor.h"

#include <math.h>
#include <stddef.h>

void rp_force_frame_prepare(
    double spin,
    const struct rp_kerr_geodesic_state *state,
    rp_force_frame *frame
)
{
    const double radius = state->x[1];
    const double sin_theta = state->sin_theta;
    const double cos_theta = state->cos_theta;
    const double sin_phi = sin(state->x[3]);
    const double cos_phi = cos(state->x[3]);
    const double rho = sqrt(radius * radius + spin * spin);
    const double radial_ratio = radius / rho;
    const double rho_sin = rho * sin_theta;
    const double rho_cos = rho * cos_theta;

    frame->position[0] = rho_sin * cos_phi;
    frame->position[1] = rho_sin * sin_phi;
    frame->position[2] = radius * cos_theta;

    frame->jacobian[0][0] = radial_ratio * sin_theta * cos_phi;
    frame->jacobian[0][1] = rho_cos * cos_phi;
    frame->jacobian[0][2] = -frame->position[1];
    frame->jacobian[1][0] = radial_ratio * sin_theta * sin_phi;
    frame->jacobian[1][1] = rho_cos * sin_phi;
    frame->jacobian[1][2] = frame->position[0];
    frame->jacobian[2][0] = cos_theta;
    frame->jacobian[2][1] = -radius * sin_theta;
    frame->jacobian[2][2] = 0.0;
}

void rp_force_frame_velocity(
    const rp_force_frame *frame,
    const double four_velocity[RP_KERR_DIM],
    double velocity[3]
)
{
    const double inverse_ut = 1.0 / four_velocity[0];
    size_t k;

    for (k = 0U; k < 3U; ++k) {
        velocity[k] = (frame->jacobian[k][0] * four_velocity[1]
            + frame->jacobian[k][1] * four_velocity[2]
            + frame->jacobian[k][2] * four_velocity[3]) * inverse_ut;
    }
}

int rp_force_evaluate_prepared(
    double spin,
    const struct rp_kerr_geodesic_state *state,
    const rp_external_force *force,
    double four_acceleration[RP_KERR_DIM]
)
{
    rp_force_frame frame;
    double velocity[3];
    double acceleration[3] = {0.0, 0.0, 0.0};
    double covector[3];
    const double *u = state->u;
    double radius;
    double spin_squared;
    double sin_squared;
    double sigma_delta;
    double u_dot_force;
    size_t i;

    rp_force_frame_prepare(spin, state, &frame);
    rp_force_frame_velocity(&frame, u, velocity);
    if (!rp_tensor_is_finite(velocity, 3U)
        || !rp_tensor_is_finite(frame.position, 3U)) {
        rp_tensor_zero(four_acceleration, RP_KERR_DIM);
        return (int)RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    if (force->evaluate(frame.position, velocity, force->model, acceleration)
        != 0) {
        rp_tensor_zero(four_acceleration, RP_KERR_DIM);
        return RP_FORCE_STATUS_CALLBACK_FAILURE;
    }
    if (!rp_tensor_is_finite(acceleration, 3U)) {
        rp_tensor_zero(four_acceleration, RP_KERR_DIM);
        return (int)RP_KERR_STATUS_NUMERICAL_RANGE;
    }

    /* F_i = A . dx/dx^i for i = r, theta, phi; F_t = 0. */
    for (i = 0U; i < 3U; ++i) {
        covector[i] = acceleration[0] * frame.jacobian[0][i]
            + acceleration[1] * frame.jacobian[1][i]
            + acceleration[2] * frame.jacobian[2][i];
    }
    u_dot_force = u[1] * covector[0] + u[2] * covector[1]
        + u[3] * covector[2];

    /* Closed-form g^(t phi), g^(rr), g^(theta theta), g^(phi phi), M = 1. */
    radius = state->x[1];
    spin_squared = spin * spin;
    sin_squared = state->sin_theta * state->sin_theta;
    sigma_delta = state->sigma * state->delta;
    four_acceleration[0] = -2.0 * radius * spin / sigma_delta * covector[2]
        + u[0] * u_dot_force;
    four_acceleration[1] = state->delta / state->sigma * covector[0]
        + u[1] * u_dot_force;
    four_acceleration[2] = covector[1] / state->sigma
        + u[2] * u_dot_force;
    four_acceleration[3] = (state->delta - spin_squared * sin_squared)
        / (sigma_delta * sin_squared) * covector[2]
        + u[3] * u_dot_force;
    if (!rp_tensor_is_finite(four_acceleration, RP_KERR_DIM)) {
        rp_tensor_zero(four_acceleration, RP_KERR_DIM);
        return (int)RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    return 0;
}
