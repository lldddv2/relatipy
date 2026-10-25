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

rp_kerr_status rp_force_bl_to_cartesian(
    double spin,
    const double coordinates[RP_KERR_DIM],
    double position[3],
    double jacobian[3][3]
)
{
    struct rp_kerr_geodesic_state state;
    rp_force_frame frame;
    size_t index;

    if (position == NULL && jacobian == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    if (coordinates != NULL) {
        for (index = 0U; index < RP_KERR_DIM; ++index) {
            state.x[index] = coordinates[index];
            state.u[index] = 0.0;
        }
    }
    if (position != NULL) {
        rp_tensor_zero(position, 3U);
    }
    if (jacobian != NULL) {
        rp_tensor_zero(&jacobian[0][0], 9U);
    }
    if (coordinates == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    if (!isfinite(spin) || !rp_tensor_is_finite(state.x, RP_KERR_DIM)) {
        return RP_KERR_STATUS_NONFINITE_INPUT;
    }
    if (fabs(spin) > 1.0) {
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    state.sin_theta = sin(state.x[2]);
    state.cos_theta = cos(state.x[2]);
    state.sigma = 0.0;
    state.delta = 0.0;
    rp_force_frame_prepare(spin, &state, &frame);
    if (!rp_tensor_is_finite(frame.position, 3U)
        || !rp_tensor_is_finite(&frame.jacobian[0][0], 9U)) {
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    for (index = 0U; index < 3U; ++index) {
        if (position != NULL) {
            position[index] = frame.position[index];
        }
        if (jacobian != NULL) {
            jacobian[index][0] = frame.jacobian[index][0];
            jacobian[index][1] = frame.jacobian[index][1];
            jacobian[index][2] = frame.jacobian[index][2];
        }
    }
    return RP_KERR_STATUS_OK;
}

int rp_force_four_acceleration(
    double spin,
    const rp_external_force *force,
    const double coordinates[RP_KERR_DIM],
    const double four_velocity[RP_KERR_DIM],
    double four_acceleration[RP_KERR_DIM]
)
{
    struct rp_kerr_geodesic_state state;
    rp_kerr_status status;
    size_t index;

    if (four_acceleration == NULL) {
        return (int)RP_KERR_STATUS_NULL_POINTER;
    }
    if (coordinates == NULL || four_velocity == NULL) {
        rp_tensor_zero(four_acceleration, RP_KERR_DIM);
        return (int)RP_KERR_STATUS_NULL_POINTER;
    }
    /* Validation mirrors prepare_geodesic_state in geodesic/kerr.c, M = 1. */
    status = RP_KERR_STATUS_OK;
    if (!isfinite(spin)) {
        status = RP_KERR_STATUS_NONFINITE_INPUT;
    }
    for (index = 0U; index < RP_KERR_DIM; ++index) {
        if (!isfinite(coordinates[index]) || !isfinite(four_velocity[index])) {
            status = RP_KERR_STATUS_NONFINITE_INPUT;
        }
        state.x[index] = coordinates[index];
        state.u[index] = four_velocity[index];
    }
    rp_tensor_zero(four_acceleration, RP_KERR_DIM);
    if (status != RP_KERR_STATUS_OK) {
        return (int)status;
    }
    if (fabs(spin) > 1.0) {
        return (int)RP_KERR_STATUS_INVALID_PARAMETER;
    }
    status = rp_kerr_geodesic_physic_prepare_state(1.0, spin, &state);
    if (status != RP_KERR_STATUS_OK) {
        return (int)status;
    }
    if (!rp_force_is_active(force)) {
        return 0;
    }
    return rp_force_evaluate_prepared(spin, &state, force, four_acceleration);
}
