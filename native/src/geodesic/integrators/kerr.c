/**
 * @file kerr.c
 * @brief ODE callback adapters for the Kerr geodesic evaluators.
 */

#include "kerr.h"

#include "relatipy/kerr_geodesic.h"
#include "relatipy/kerr_geometry.h"
#include "geodesic/jacobian.h"

#include <float.h>
#include <math.h>
#include <stddef.h>

void rp_kerr_timelike_constants_from_momentum(
    double spin,
    double theta,
    const double momentum[4],
    double constants[3]
)
{
    const double sine = sin(theta);
    const double cosine = cos(theta);
    const double energy = -momentum[0];
    const double angular_momentum = momentum[3];

    constants[0] = energy;
    constants[1] = angular_momentum;
    constants[2] = momentum[2] * momentum[2]
        + cosine * cosine * (spin * spin * (1.0 - energy * energy)
            + angular_momentum * angular_momentum / (sine * sine));
}

int rp_kerr_timelike_constants(
    double mass,
    double spin,
    const double state[8],
    double constants[3]
)
{
    double metric[RP_KERR_DIM][RP_KERR_DIM];
    double momentum[RP_KERR_DIM] = {0.0, 0.0, 0.0, 0.0};
    rp_kerr_status status;
    size_t mu;
    size_t nu;

    if (constants == NULL) {
        return (int)RP_KERR_STATUS_NULL_POINTER;
    }
    constants[0] = constants[1] = constants[2] = 0.0;
    if (state == NULL) {
        return (int)RP_KERR_STATUS_NULL_POINTER;
    }
    for (mu = 0U; mu < 2U * RP_KERR_DIM; ++mu) {
        if (!isfinite(state[mu])) {
            return (int)RP_KERR_STATUS_NONFINITE_INPUT;
        }
    }
    status = rp_kerr_metric(mass, spin, state, metric);
    if (status != RP_KERR_STATUS_OK) {
        return (int)status;
    }
    for (mu = 0U; mu < RP_KERR_DIM; ++mu) {
        for (nu = 0U; nu < RP_KERR_DIM; ++nu) {
            momentum[mu] += metric[mu][nu] * state[RP_KERR_DIM + nu];
        }
    }
    rp_kerr_timelike_constants_from_momentum(spin, state[2], momentum, constants);
    if (!isfinite(constants[0]) || !isfinite(constants[1])
        || !isfinite(constants[2])) {
        constants[0] = constants[1] = constants[2] = 0.0;
        return (int)RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    return 0;
}

int rp_kerr_integrator_prepare_projection(
    rp_kerr_integrator_context *context,
    const double state[8]
)
{
    double metric[RP_KERR_DIM][RP_KERR_DIM];
    double momentum[RP_KERR_DIM] = {0.0, 0.0, 0.0, 0.0};
    double norm = 0.0;
    double constants[3];
    double horizon;
    rp_kerr_status status;
    size_t mu;
    size_t nu;

    if (context == NULL || state == NULL) {
        return (int)RP_KERR_STATUS_NULL_POINTER;
    }
    for (mu = 0U; mu < 2U * RP_KERR_DIM; ++mu) {
        if (!isfinite(state[mu])) {
            return (int)RP_KERR_STATUS_NONFINITE_INPUT;
        }
    }
    status = rp_kerr_metric(context->mass, context->spin, state, metric);
    if (status != RP_KERR_STATUS_OK) {
        return (int)status;
    }
    horizon = context->mass
        + sqrt((context->mass - fabs(context->spin))
            * (context->mass + fabs(context->spin)));
    if (!isfinite(horizon)) {
        return (int)RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    if (state[1] <= horizon) {
        return (int)RP_KERR_STATUS_INVALID_PARAMETER;
    }
    for (mu = 0U; mu < RP_KERR_DIM; ++mu) {
        for (nu = 0U; nu < RP_KERR_DIM; ++nu) {
            momentum[mu] += metric[mu][nu] * state[RP_KERR_DIM + nu];
        }
        norm += momentum[mu] * state[RP_KERR_DIM + mu];
    }
    if (!isfinite(norm) || fabs(norm + 1.0) > 1.0e-8) {
        return (int)RP_KERR_STATUS_NON_TIMELIKE_VELOCITY;
    }
    rp_kerr_timelike_constants_from_momentum(
        context->spin, state[2], momentum, constants
    );
    if (!isfinite(constants[0]) || !isfinite(constants[1])
        || !isfinite(constants[2])) {
        return (int)RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    context->E0 = constants[0];
    context->Lz0 = constants[1];
    context->Q0 = constants[2];
    context->projection_ready = 1;
    return 0;
}

int rp_kerr_integrator_project(double tau, double state[], void *opaque_context)
{
    const rp_kerr_integrator_context *context = opaque_context;
    double metric[RP_KERR_DIM][RP_KERR_DIM];
    double inverse[RP_KERR_DIM][RP_KERR_DIM];
    double velocity[RP_KERR_DIM];
    double momentum[RP_KERR_DIM] = {0.0, 0.0, 0.0, 0.0};
    double horizon;
    double sine;
    double cosine;
    double polar_term;
    double polar_potential;
    double radial_potential;
    double scale;
    double norm = 0.0;
    double energy;
    double angular_momentum;
    double carter_q;
    double constants[3];
    rp_kerr_status status;
    size_t mu;
    size_t nu;

    if (context == NULL || state == NULL) {
        return (int)RP_KERR_STATUS_NULL_POINTER;
    }
    if (!isfinite(tau)) {
        return (int)RP_KERR_STATUS_NONFINITE_INPUT;
    }
    if (context->projection_ready != 1 || !(context->mass > 0.0)
        || fabs(context->spin) > context->mass) {
        return (int)RP_KERR_STATUS_INVALID_PARAMETER;
    }
    if (!isfinite(context->mass) || !isfinite(context->spin)
        || !isfinite(context->E0) || !isfinite(context->Lz0)
        || !isfinite(context->Q0)) {
        return (int)RP_KERR_STATUS_NONFINITE_INPUT;
    }
    for (mu = 0U; mu < 2U * RP_KERR_DIM; ++mu) {
        if (!isfinite(state[mu])) {
            return (int)RP_KERR_STATUS_NONFINITE_INPUT;
        }
    }

    horizon = context->mass
        + sqrt((context->mass - fabs(context->spin))
            * (context->mass + fabs(context->spin)));
    if (!isfinite(horizon)) {
        return (int)RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    if (state[1] <= horizon) {
        return 0;
    }

    status = rp_kerr_metric(context->mass, context->spin, state, metric);
    if (status != RP_KERR_STATUS_OK) {
        return (int)status;
    }
    status = rp_kerr_inverse_metric(
        context->mass, context->spin, state, inverse
    );
    if (status != RP_KERR_STATUS_OK) {
        return (int)status;
    }

    /* Fix the Killing momenta, then recover polar and radial magnitudes. */
    velocity[0] = -inverse[0][0] * context->E0
        + inverse[0][3] * context->Lz0;
    velocity[3] = -inverse[3][0] * context->E0
        + inverse[3][3] * context->Lz0;
    sine = sin(state[2]);
    cosine = cos(state[2]);
    polar_term = cosine * cosine * (context->spin * context->spin
        * (1.0 - context->E0 * context->E0)
        + context->Lz0 * context->Lz0 / (sine * sine));
    polar_potential = context->Q0 - polar_term;
    scale = 1.0 + fabs(context->Q0) + fabs(polar_term);
    if (!isfinite(polar_potential) || !isfinite(scale)) {
        return (int)RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    if (polar_potential < -128.0 * DBL_EPSILON * scale) {
        return (int)RP_KERR_STATUS_NON_TIMELIKE_VELOCITY;
    }
    if (polar_potential < 0.0) {
        polar_potential = 0.0;
    }
    velocity[2] = copysign(sqrt(polar_potential) / metric[2][2], state[6]);

    radial_potential = (-1.0 + context->E0 * velocity[0]
        - context->Lz0 * velocity[3]
        - metric[2][2] * velocity[2] * velocity[2]) / metric[1][1];
    scale = (1.0 + fabs(context->E0 * velocity[0])
        + fabs(context->Lz0 * velocity[3])
        + fabs(metric[2][2] * velocity[2] * velocity[2]))
        / metric[1][1];
    if (!isfinite(radial_potential) || !isfinite(scale)) {
        return (int)RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    if (radial_potential < -128.0 * DBL_EPSILON * scale) {
        return (int)RP_KERR_STATUS_NON_TIMELIKE_VELOCITY;
    }
    if (radial_potential < 0.0) {
        radial_potential = 0.0;
    }
    velocity[1] = copysign(sqrt(radial_potential), state[5]);

    for (mu = 0U; mu < RP_KERR_DIM; ++mu) {
        if (!isfinite(velocity[mu])) {
            return (int)RP_KERR_STATUS_NUMERICAL_RANGE;
        }
        for (nu = 0U; nu < RP_KERR_DIM; ++nu) {
            momentum[mu] += metric[mu][nu] * velocity[nu];
        }
        norm += momentum[mu] * velocity[mu];
    }
    rp_kerr_timelike_constants_from_momentum(
        context->spin, state[2], momentum, constants
    );
    energy = constants[0];
    angular_momentum = constants[1];
    carter_q = constants[2];
    if (!isfinite(norm) || !isfinite(energy) || !isfinite(angular_momentum)
        || !isfinite(carter_q)) {
        return (int)RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    if (fabs(norm + 1.0) > 1.0e-10
        || fabs(energy - context->E0) > 1.0e-10 * (1.0 + fabs(context->E0))
        || fabs(angular_momentum - context->Lz0)
            > 1.0e-10 * (1.0 + fabs(context->Lz0))
        || fabs(carter_q - context->Q0)
            > 1.0e-10 * (1.0 + fabs(context->Q0))) {
        return (int)RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    for (mu = 0U; mu < RP_KERR_DIM; ++mu) {
        state[RP_KERR_DIM + mu] = velocity[mu];
    }
    return 0;
}

int rp_kerr_integrator_jacobian(
    double affine_parameter,
    const double state[],
    double *jacobian,
    size_t dimension,
    size_t row_stride,
    void *opaque_context
)
{
    const rp_kerr_integrator_context *context = opaque_context;
    double values[2 * RP_KERR_DIM][2 * RP_KERR_DIM];
    rp_kerr_status status;
    size_t row;
    size_t column;

    (void)affine_parameter;
    if (context == NULL || state == NULL || jacobian == NULL) {
        return (int)RP_KERR_STATUS_NULL_POINTER;
    }
    if (dimension != 2U * RP_KERR_DIM || row_stride < dimension) {
        return (int)RP_KERR_STATUS_INVALID_PARAMETER;
    }
    status = rp_kerr_geodesic_jacobian(
        context->mass, context->spin, state, state + RP_KERR_DIM, values
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

int rp_kerr_integrator_rhs(
    double affine_parameter,
    const double state[],
    double derivative[],
    void *opaque_context
)
{
    const rp_kerr_integrator_context *context = opaque_context;

    (void)affine_parameter;
    if (context == NULL) {
        return (int)RP_KERR_STATUS_NULL_POINTER;
    }
    /* The corrected metric contraction applies to timelike tangents too. */
    return (int)rp_kerr_null_geodesic_rhs(
        context->mass,
        context->spin,
        state,
        state + RP_KERR_DIM,
        derivative,
        derivative + RP_KERR_DIM
    );
}

int rp_kerr_null_integrator_rhs(
    double affine_parameter,
    const double state[],
    double derivative[],
    void *opaque_context
)
{
    const rp_kerr_integrator_context *context = opaque_context;

    (void)affine_parameter;
    if (context == NULL) {
        return (int)RP_KERR_STATUS_NULL_POINTER;
    }
    return (int)rp_kerr_null_geodesic_rhs(
        context->mass,
        context->spin,
        state,
        state + RP_KERR_DIM,
        derivative,
        derivative + RP_KERR_DIM
    );
}
