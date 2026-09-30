/**
 * @file kerr.c
 * @brief Orchestration layer for experimental Kerr geodesic operations.
 *
 * This file validates and copies caller inputs, enforces output and aliasing
 * contracts, propagates typed failures, and delegates legacy-compatible and
 * corrected null-geodesic physical evaluation to `physic/kerr.c`.
 *
 * The operation remains an internal local reference. It is not an integrator,
 * public C API, stable ABI, or decision about the canonical integration state.
 */

#include "relatipy/kerr_geodesic.h"
#include "jacobian.h"
#include "physic/kerr.h"
#include "../utils/tensor.h"

#include <math.h>
#include <stddef.h>

static rp_kerr_status prepare_geodesic_state(
    double mass,
    double spin,
    const double coordinates[RP_KERR_DIM],
    const double four_velocity[RP_KERR_DIM],
    struct rp_kerr_geodesic_state *state
)
{
    size_t index;

    if (coordinates == NULL || four_velocity == NULL || state == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    if (!isfinite(mass) || !isfinite(spin)) {
        return RP_KERR_STATUS_NONFINITE_INPUT;
    }
    for (index = 0; index < RP_KERR_DIM; ++index) {
        if (!isfinite(coordinates[index]) || !isfinite(four_velocity[index])) {
            return RP_KERR_STATUS_NONFINITE_INPUT;
        }
        state->x[index] = coordinates[index];
        state->u[index] = four_velocity[index];
    }
    if (!(mass > 0.0) || fabs(spin) > mass) {
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    return rp_kerr_geodesic_physic_prepare_state(mass, spin, state);
}

rp_kerr_status rp_kerr_geodesic_jacobian(
    double mass,
    double spin,
    const double coordinates[RP_KERR_DIM],
    const double tangent[RP_KERR_DIM],
    double jacobian[2 * RP_KERR_DIM][2 * RP_KERR_DIM]
)
{
    struct rp_kerr_geodesic_state state;
    rp_kerr_status status;

    if (jacobian == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    status = prepare_geodesic_state(mass, spin, coordinates, tangent, &state);
    rp_tensor_zero(&jacobian[0][0], 4U * RP_KERR_DIM * RP_KERR_DIM);
    if (status != RP_KERR_STATUS_OK) {
        return status;
    }
    rp_kerr_geodesic_physic_evaluate_corrected_jacobian(
        mass, spin, &state, jacobian
    );
    if (!rp_tensor_is_finite(&jacobian[0][0], 4U * RP_KERR_DIM * RP_KERR_DIM)) {
        rp_tensor_zero(&jacobian[0][0], 4U * RP_KERR_DIM * RP_KERR_DIM);
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    return RP_KERR_STATUS_OK;
}

rp_kerr_status rp_kerr_geodesic_rhs(
    double mass,
    double spin,
    const double coordinates[RP_KERR_DIM],
    const double four_velocity[RP_KERR_DIM],
    double coordinate_derivative[RP_KERR_DIM],
    double four_velocity_derivative[RP_KERR_DIM]
)
{
    struct rp_kerr_geodesic_state state;
    rp_kerr_status status;

    if (coordinate_derivative == NULL || four_velocity_derivative == NULL) {
        if (coordinate_derivative != NULL) {
            rp_tensor_zero(coordinate_derivative, RP_KERR_DIM);
        }
        if (four_velocity_derivative != NULL) {
            rp_tensor_zero(four_velocity_derivative, RP_KERR_DIM);
        }
        return RP_KERR_STATUS_NULL_POINTER;
    }
    if (coordinate_derivative == four_velocity_derivative) {
        rp_tensor_zero(coordinate_derivative, RP_KERR_DIM);
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    if (coordinates == NULL || four_velocity == NULL) {
        rp_tensor_zero(coordinate_derivative, RP_KERR_DIM);
        rp_tensor_zero(four_velocity_derivative, RP_KERR_DIM);
        return RP_KERR_STATUS_NULL_POINTER;
    }

    status = prepare_geodesic_state(
        mass, spin, coordinates, four_velocity, &state
    );
    rp_tensor_zero(coordinate_derivative, RP_KERR_DIM);
    rp_tensor_zero(four_velocity_derivative, RP_KERR_DIM);
    if (status != RP_KERR_STATUS_OK) {
        return status;
    }

    rp_kerr_geodesic_physic_evaluate_rhs(
        mass,
        spin,
        &state,
        coordinate_derivative,
        four_velocity_derivative
    );
    if (!rp_tensor_is_finite(coordinate_derivative, RP_KERR_DIM)
        || !rp_tensor_is_finite(four_velocity_derivative, RP_KERR_DIM)) {
        rp_tensor_zero(coordinate_derivative, RP_KERR_DIM);
        rp_tensor_zero(four_velocity_derivative, RP_KERR_DIM);
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    return RP_KERR_STATUS_OK;
}

rp_kerr_status rp_kerr_null_tangent(
    double mass,
    double spin,
    const double coordinates[RP_KERR_DIM],
    double energy,
    double axial_angular_momentum,
    double carter_constant,
    int radial_direction,
    int polar_direction,
    double tangent[RP_KERR_DIM]
)
{
    const double zero_tangent[RP_KERR_DIM] = {0.0, 0.0, 0.0, 0.0};
    struct rp_kerr_geodesic_state state;
    rp_kerr_status status;

    if (tangent == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    if (coordinates == NULL) {
        rp_tensor_zero(tangent, RP_KERR_DIM);
        return RP_KERR_STATUS_NULL_POINTER;
    }
    status = prepare_geodesic_state(
        mass, spin, coordinates, zero_tangent, &state
    );
    rp_tensor_zero(tangent, RP_KERR_DIM);
    if (status != RP_KERR_STATUS_OK) {
        return status;
    }
    if (!isfinite(energy) || !isfinite(axial_angular_momentum)
        || !isfinite(carter_constant)) {
        return RP_KERR_STATUS_NONFINITE_INPUT;
    }
    if ((radial_direction != -1 && radial_direction != 1)
        || (polar_direction != -1 && polar_direction != 1)
        || (energy == 0.0 && axial_angular_momentum == 0.0
            && carter_constant == 0.0)) {
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }

    return rp_kerr_geodesic_physic_evaluate_null_tangent(
        spin,
        &state,
        energy,
        axial_angular_momentum,
        carter_constant,
        radial_direction,
        polar_direction,
        tangent
    );
}

rp_kerr_status rp_kerr_null_geodesic_rhs(
    double mass,
    double spin,
    const double coordinates[RP_KERR_DIM],
    const double tangent[RP_KERR_DIM],
    double coordinate_derivative[RP_KERR_DIM],
    double tangent_derivative[RP_KERR_DIM]
)
{
    struct rp_kerr_geodesic_state state;
    rp_kerr_status status;

    if (coordinate_derivative == NULL || tangent_derivative == NULL) {
        if (coordinate_derivative != NULL) {
            rp_tensor_zero(coordinate_derivative, RP_KERR_DIM);
        }
        if (tangent_derivative != NULL) {
            rp_tensor_zero(tangent_derivative, RP_KERR_DIM);
        }
        return RP_KERR_STATUS_NULL_POINTER;
    }
    if (coordinate_derivative == tangent_derivative) {
        rp_tensor_zero(coordinate_derivative, RP_KERR_DIM);
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    if (coordinates == NULL || tangent == NULL) {
        rp_tensor_zero(coordinate_derivative, RP_KERR_DIM);
        rp_tensor_zero(tangent_derivative, RP_KERR_DIM);
        return RP_KERR_STATUS_NULL_POINTER;
    }

    status = prepare_geodesic_state(
        mass, spin, coordinates, tangent, &state
    );
    rp_tensor_zero(coordinate_derivative, RP_KERR_DIM);
    rp_tensor_zero(tangent_derivative, RP_KERR_DIM);
    if (status != RP_KERR_STATUS_OK) {
        return status;
    }

    rp_kerr_geodesic_physic_evaluate_corrected_rhs(
        mass,
        spin,
        &state,
        coordinate_derivative,
        tangent_derivative
    );
    if (!rp_tensor_is_finite(coordinate_derivative, RP_KERR_DIM)
        || !rp_tensor_is_finite(tangent_derivative, RP_KERR_DIM)) {
        rp_tensor_zero(coordinate_derivative, RP_KERR_DIM);
        rp_tensor_zero(tangent_derivative, RP_KERR_DIM);
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    return RP_KERR_STATUS_OK;
}
