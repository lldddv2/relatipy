/**
 * @file kerr.c
 * @brief Orchestration layer for the experimental Kerr reference.
 *
 * This file implements the operations declared in
 * `relatipy/kerr_geometry.h`.  It validates and copies caller inputs,
 * preserves the zero-on-error and input/output-aliasing contracts, and
 * delegates physical formulas to `physic/kerr.c`.
 *
 * The implementation is allocation-free and retains no caller pointer.  It
 * is not an orbit backend, integrator, Python binding, public C API, or stable
 * ABI.
 */

#include "relatipy/kerr_geometry.h"
#include "physic/kerr.h"
#include "../utils/tensor.h"

#include <math.h>
#include <stddef.h>

static rp_kerr_status prepare_point(
    double mass,
    double spin,
    const double coordinates[RP_KERR_DIM],
    struct rp_kerr_point *point
)
{
    size_t index;

    if (coordinates == NULL || point == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    if (!isfinite(mass) || !isfinite(spin)) {
        return RP_KERR_STATUS_NONFINITE_INPUT;
    }
    for (index = 0; index < RP_KERR_DIM; ++index) {
        if (!isfinite(coordinates[index])) {
            return RP_KERR_STATUS_NONFINITE_INPUT;
        }
        point->x[index] = coordinates[index];
    }
    if (!(mass > 0.0) || fabs(spin) > mass) {
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    return rp_kerr_physic_prepare_point(mass, spin, point);
}

rp_kerr_status rp_kerr_metric(
    double mass,
    double spin,
    const double coordinates[RP_KERR_DIM],
    double metric[RP_KERR_DIM][RP_KERR_DIM]
)
{
    struct rp_kerr_point point;
    rp_kerr_status status;

    if (metric == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    if (coordinates == NULL) {
        rp_tensor_zero(&metric[0][0], RP_KERR_DIM * RP_KERR_DIM);
        return RP_KERR_STATUS_NULL_POINTER;
    }
    status = prepare_point(mass, spin, coordinates, &point);
    rp_tensor_zero(&metric[0][0], RP_KERR_DIM * RP_KERR_DIM);
    if (status != RP_KERR_STATUS_OK) {
        return status;
    }

    rp_kerr_physic_evaluate_metric(mass, spin, &point, metric);
    if (!rp_tensor_is_finite(
        &metric[0][0], RP_KERR_DIM * RP_KERR_DIM
    )) {
        rp_tensor_zero(&metric[0][0], RP_KERR_DIM * RP_KERR_DIM);
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    return RP_KERR_STATUS_OK;
}

rp_kerr_status rp_kerr_inverse_metric(
    double mass,
    double spin,
    const double coordinates[RP_KERR_DIM],
    double inverse_metric[RP_KERR_DIM][RP_KERR_DIM]
)
{
    struct rp_kerr_point point;
    rp_kerr_status status;

    if (inverse_metric == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    if (coordinates == NULL) {
        rp_tensor_zero(&inverse_metric[0][0], RP_KERR_DIM * RP_KERR_DIM);
        return RP_KERR_STATUS_NULL_POINTER;
    }
    status = prepare_point(mass, spin, coordinates, &point);
    rp_tensor_zero(&inverse_metric[0][0], RP_KERR_DIM * RP_KERR_DIM);
    if (status != RP_KERR_STATUS_OK) {
        return status;
    }

    status = rp_kerr_physic_evaluate_inverse_metric(
        mass, spin, &point, inverse_metric
    );
    if (status != RP_KERR_STATUS_OK) {
        rp_tensor_zero(&inverse_metric[0][0], RP_KERR_DIM * RP_KERR_DIM);
    }
    return status;
}

rp_kerr_status rp_kerr_christoffel(
    double mass,
    double spin,
    const double coordinates[RP_KERR_DIM],
    double christoffel[RP_KERR_DIM][RP_KERR_DIM][RP_KERR_DIM]
)
{
    struct rp_kerr_point point;
    rp_kerr_status status;

    if (christoffel == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    if (coordinates == NULL) {
        rp_tensor_zero(
            &christoffel[0][0][0],
            RP_KERR_DIM * RP_KERR_DIM * RP_KERR_DIM
        );
        return RP_KERR_STATUS_NULL_POINTER;
    }
    status = prepare_point(mass, spin, coordinates, &point);
    rp_tensor_zero(
        &christoffel[0][0][0],
        RP_KERR_DIM * RP_KERR_DIM * RP_KERR_DIM
    );
    if (status != RP_KERR_STATUS_OK) {
        return status;
    }

    rp_kerr_physic_evaluate_christoffel(mass, spin, &point, christoffel);
    if (!rp_tensor_is_finite(
        &christoffel[0][0][0],
        RP_KERR_DIM * RP_KERR_DIM * RP_KERR_DIM
    )) {
        rp_tensor_zero(
            &christoffel[0][0][0],
            RP_KERR_DIM * RP_KERR_DIM * RP_KERR_DIM
        );
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    return RP_KERR_STATUS_OK;
}

rp_kerr_status rp_kerr_four_velocity(
    double mass,
    double spin,
    const double coordinates[RP_KERR_DIM],
    const double coordinate_velocity[3],
    double four_velocity[RP_KERR_DIM]
)
{
    struct rp_kerr_point point;
    double velocity[3];
    rp_kerr_status status;
    size_t index;

    if (four_velocity == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    if (coordinates == NULL || coordinate_velocity == NULL) {
        rp_tensor_zero(four_velocity, RP_KERR_DIM);
        return RP_KERR_STATUS_NULL_POINTER;
    }
    for (index = 0; index < 3; ++index) {
        velocity[index] = coordinate_velocity[index];
        if (!isfinite(velocity[index])) {
            rp_tensor_zero(four_velocity, RP_KERR_DIM);
            return RP_KERR_STATUS_NONFINITE_INPUT;
        }
    }

    status = prepare_point(mass, spin, coordinates, &point);
    rp_tensor_zero(four_velocity, RP_KERR_DIM);
    if (status != RP_KERR_STATUS_OK) {
        return status;
    }

    status = rp_kerr_physic_evaluate_four_velocity(
        mass, spin, &point, velocity, four_velocity
    );
    if (status != RP_KERR_STATUS_OK) {
        rp_tensor_zero(four_velocity, RP_KERR_DIM);
    }
    return status;
}
