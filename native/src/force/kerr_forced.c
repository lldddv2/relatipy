/**
 * @file kerr_forced.c
 * @brief Forced Kerr proper- and coordinate-time right-hand sides.
 *
 * The homogeneous term is called, never copied, through the single helper
 * `evaluate_homogeneous`.  An inactive force uses the corrected contraction
 * of `rp_kerr_null_geodesic_rhs` in `geodesic/kerr.c`, whose validation is
 * mirrored here, so it reproduces that path bit for bit.  An active force
 * uses the generated kernel of the same connection
 * (`generated/kerr_rhs.c`), which agrees with the contraction to 1e-12
 * relative and is about ten times faster.
 */

#include "kerr_forced.h"

#include "generated/kerr_rhs.h"

#include "geodesic/jacobian.h"
#include "../utils/tensor.h"

#include <math.h>
#include <stddef.h>

#define STATE_DIM (2U * RP_KERR_DIM)

/* The only calls into a homogeneous Kerr kernel (mass one). */
static void evaluate_homogeneous(
    double spin,
    const struct rp_kerr_geodesic_state *state,
    int force_active,
    double coordinate_derivative[RP_KERR_DIM],
    double four_velocity_derivative[RP_KERR_DIM]
)
{
    if (force_active) {
        rp_force_generated_kerr_rhs(
            1.0, spin, state, coordinate_derivative, four_velocity_derivative
        );
        return;
    }
    rp_kerr_geodesic_physic_evaluate_corrected_rhs(
        1.0, spin, state, coordinate_derivative, four_velocity_derivative
    );
}

/* Minimal copy of prepare_geodesic_state in geodesic/kerr.c with M = 1. */
static rp_kerr_status prepare_state(
    double spin,
    const double coordinates[RP_KERR_DIM],
    const double four_velocity[RP_KERR_DIM],
    struct rp_kerr_geodesic_state *state
)
{
    size_t index;

    if (!isfinite(spin)) {
        return RP_KERR_STATUS_NONFINITE_INPUT;
    }
    for (index = 0U; index < RP_KERR_DIM; ++index) {
        if (!isfinite(coordinates[index]) || !isfinite(four_velocity[index])) {
            return RP_KERR_STATUS_NONFINITE_INPUT;
        }
        state->x[index] = coordinates[index];
        state->u[index] = four_velocity[index];
    }
    if (fabs(spin) > 1.0) {
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    return rp_kerr_geodesic_physic_prepare_state(1.0, spin, state);
}

int rp_kerr_forced_rhs(
    double spin,
    const rp_external_force *force,
    const double coordinates[RP_KERR_DIM],
    const double four_velocity[RP_KERR_DIM],
    double coordinate_derivative[RP_KERR_DIM],
    double four_velocity_derivative[RP_KERR_DIM]
)
{
    struct rp_kerr_geodesic_state state;
    double four_acceleration[RP_KERR_DIM];
    rp_kerr_status status;
    int force_status;
    size_t index;

    if (coordinate_derivative == NULL || four_velocity_derivative == NULL) {
        if (coordinate_derivative != NULL) {
            rp_tensor_zero(coordinate_derivative, RP_KERR_DIM);
        }
        if (four_velocity_derivative != NULL) {
            rp_tensor_zero(four_velocity_derivative, RP_KERR_DIM);
        }
        return (int)RP_KERR_STATUS_NULL_POINTER;
    }
    if (coordinate_derivative == four_velocity_derivative) {
        rp_tensor_zero(coordinate_derivative, RP_KERR_DIM);
        return (int)RP_KERR_STATUS_INVALID_PARAMETER;
    }
    if (coordinates == NULL || four_velocity == NULL) {
        rp_tensor_zero(coordinate_derivative, RP_KERR_DIM);
        rp_tensor_zero(four_velocity_derivative, RP_KERR_DIM);
        return (int)RP_KERR_STATUS_NULL_POINTER;
    }

    status = prepare_state(spin, coordinates, four_velocity, &state);
    rp_tensor_zero(coordinate_derivative, RP_KERR_DIM);
    rp_tensor_zero(four_velocity_derivative, RP_KERR_DIM);
    if (status != RP_KERR_STATUS_OK) {
        return (int)status;
    }

    evaluate_homogeneous(
        spin, &state, rp_force_is_active(force),
        coordinate_derivative, four_velocity_derivative
    );
    if (rp_force_is_active(force)) {
        force_status = rp_force_evaluate_prepared(
            spin, &state, force, four_acceleration
        );
        if (force_status != 0) {
            rp_tensor_zero(coordinate_derivative, RP_KERR_DIM);
            rp_tensor_zero(four_velocity_derivative, RP_KERR_DIM);
            return force_status;
        }
        for (index = 0U; index < RP_KERR_DIM; ++index) {
            four_velocity_derivative[index] += four_acceleration[index];
        }
    }
    if (!rp_tensor_is_finite(coordinate_derivative, RP_KERR_DIM)
        || !rp_tensor_is_finite(four_velocity_derivative, RP_KERR_DIM)) {
        rp_tensor_zero(coordinate_derivative, RP_KERR_DIM);
        rp_tensor_zero(four_velocity_derivative, RP_KERR_DIM);
        return (int)RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    return 0;
}

int rp_kerr_forced_integrator_rhs(
    double proper_time,
    const double state[],
    double derivative[],
    void *opaque_context
)
{
    const rp_kerr_forced_context *context = opaque_context;

    (void)proper_time;
    if (context == NULL || state == NULL || derivative == NULL) {
        return (int)RP_KERR_STATUS_NULL_POINTER;
    }
    return rp_kerr_forced_rhs(
        context->spin, context->force,
        state, state + RP_KERR_DIM,
        derivative, derivative + RP_KERR_DIM
    );
}

int rp_kerr_forced_coordinate_time_rhs(
    double coordinate_time,
    const double state[],
    double derivative[],
    void *opaque_context
)
{
    const rp_kerr_forced_context *context = opaque_context;
    double dx[RP_KERR_DIM];
    double du[RP_KERR_DIM];
    double inverse_ut;
    int status;
    size_t index;

    (void)coordinate_time;
    if (context == NULL || state == NULL || derivative == NULL) {
        return (int)RP_KERR_STATUS_NULL_POINTER;
    }
    if (!(state[4] > 0.0) || !isfinite(state[4])) {
        return (int)RP_KERR_STATUS_NON_TIMELIKE_VELOCITY;
    }
    status = rp_kerr_forced_rhs(
        context->spin, context->force, state, state + RP_KERR_DIM, dx, du
    );
    if (status != 0) {
        return status;
    }
    /* Same operation order as coordinate_time_rhs in geodesic/kerr_mcmc.c. */
    inverse_ut = 1.0 / state[4];
    for (index = 0U; index < RP_KERR_DIM; ++index) {
        derivative[index] = dx[index] * inverse_ut;
        derivative[index + RP_KERR_DIM] = du[index] * inverse_ut;
        if (!isfinite(derivative[index])
            || !isfinite(derivative[index + RP_KERR_DIM])) {
            return (int)RP_KERR_STATUS_NUMERICAL_RANGE;
        }
    }
    return 0;
}

int rp_kerr_forced_integrator_jacobian(
    double proper_time,
    const double state[],
    double *jacobian,
    size_t dimension,
    size_t row_stride,
    void *opaque_context
)
{
    const rp_kerr_forced_context *context = opaque_context;
    double values[STATE_DIM][STATE_DIM];
    rp_kerr_status status;
    size_t row;
    size_t column;

    (void)proper_time;
    if (context == NULL || state == NULL || jacobian == NULL) {
        return (int)RP_KERR_STATUS_NULL_POINTER;
    }
    if (dimension != STATE_DIM || row_stride < dimension
        || rp_force_is_active(context->force)) {
        return (int)RP_KERR_STATUS_INVALID_PARAMETER;
    }
    /* Same body as rp_kerr_integrator_jacobian in integrators/kerr.c. */
    status = rp_kerr_geodesic_jacobian(
        1.0, context->spin, state, state + RP_KERR_DIM, values
    );
    if (status != RP_KERR_STATUS_OK) {
        return (int)status;
    }
    for (row = 0U; row < dimension; ++row) {
        for (column = 0U; column < dimension; ++column) {
            jacobian[row * row_stride + column] = values[row][column];
        }
    }
    return 0;
}

rp_integrator_status rp_kerr_forced_check_method(rp_integrator_method method)
{
    switch (method) {
    case RP_INTEGRATOR_METHOD_DOP853:
    case RP_INTEGRATOR_METHOD_DP45:
    case RP_INTEGRATOR_METHOD_RADAU:
        return RP_INTEGRATOR_STATUS_OK;
    case RP_INTEGRATOR_METHOD_PROJECTION_RADAU:
    case RP_INTEGRATOR_METHOD_COUNT:
    default:
        return RP_INTEGRATOR_STATUS_INVALID_ARGUMENT;
    }
}
