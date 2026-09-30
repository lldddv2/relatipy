/**
 * @file kerr_observables.c
 * @brief Allocation-free geometric projection of a Kerr coordinate state.
 */

#include "kerr_observables.h"

#include <math.h>
#include <stddef.h>

static void clear_output(double output[RP_KERR_OBSERVABLE_OUTPUT_DIM])
{
    size_t index;

    for (index = 0U; index < RP_KERR_OBSERVABLE_OUTPUT_DIM; ++index) {
        output[index] = 0.0;
    }
}

static int proper_rotation(const double rotation[3][3])
{
    size_t row;
    size_t column;
    size_t index;
    double product;
    double determinant;

    for (row = 0U; row < 3U; ++row) {
        for (column = 0U; column < 3U; ++column) {
            if (!isfinite(rotation[row][column])) {
                return 0;
            }
        }
    }

    for (row = 0U; row < 3U; ++row) {
        for (column = 0U; column < 3U; ++column) {
            product = 0.0;
            for (index = 0U; index < 3U; ++index) {
                product += rotation[row][index] * rotation[column][index];
            }
            if (!isfinite(product)
                || fabs(product - (row == column ? 1.0 : 0.0)) > 1e-10) {
                return 0;
            }
        }
    }

    determinant = rotation[0][0] * (
        rotation[1][1] * rotation[2][2]
        - rotation[1][2] * rotation[2][1]
    ) - rotation[0][1] * (
        rotation[1][0] * rotation[2][2]
        - rotation[1][2] * rotation[2][0]
    ) + rotation[0][2] * (
        rotation[1][0] * rotation[2][1]
        - rotation[1][1] * rotation[2][0]
    );
    return isfinite(determinant) && fabs(determinant - 1.0) <= 1e-10;
}

rp_kerr_status rp_kerr_observables_project(
    const double state[RP_KERR_OBSERVABLE_STATE_DIM],
    double spin,
    const double rotation[3][3],
    double angle_scale,
    double output[RP_KERR_OBSERVABLE_OUTPUT_DIM]
)
{
    double radial_scale;
    double radial_derivative;
    double sin_theta;
    double cos_theta;
    double sin_phi;
    double cos_phi;
    double position[3];
    double velocity[3];
    double observer_position[3];
    double observer_velocity[3];
    size_t index;
    size_t component;

    if (output == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    clear_output(output);
    if (state == NULL || rotation == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    if (!isfinite(spin) || !isfinite(angle_scale)) {
        return RP_KERR_STATUS_NONFINITE_INPUT;
    }
    for (index = 0U; index < RP_KERR_OBSERVABLE_STATE_DIM; ++index) {
        if (!isfinite(state[index])) {
            return RP_KERR_STATUS_NONFINITE_INPUT;
        }
    }
    for (index = 0U; index < 3U; ++index) {
        for (component = 0U; component < 3U; ++component) {
            if (!isfinite(rotation[index][component])) {
                return RP_KERR_STATUS_NONFINITE_INPUT;
            }
        }
    }
    if (state[1] <= 0.0 || spin < 0.0 || angle_scale <= 0.0
        || !proper_rotation(rotation)) {
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }

    radial_scale = hypot(state[1], spin);
    radial_derivative = state[1] / radial_scale * state[5];
    sin_theta = sin(state[2]);
    cos_theta = cos(state[2]);
    sin_phi = sin(state[3]);
    cos_phi = cos(state[3]);

    position[0] = radial_scale * sin_theta * cos_phi;
    position[1] = radial_scale * sin_theta * sin_phi;
    position[2] = state[1] * cos_theta;

    velocity[0] = radial_derivative * sin_theta * cos_phi
        + radial_scale * cos_theta * cos_phi * state[6]
        - radial_scale * sin_theta * sin_phi * state[7];
    velocity[1] = radial_derivative * sin_theta * sin_phi
        + radial_scale * cos_theta * sin_phi * state[6]
        + radial_scale * sin_theta * cos_phi * state[7];
    velocity[2] = state[5] * cos_theta - state[1] * sin_theta * state[6];

    for (index = 0U; index < 3U; ++index) {
        observer_position[index] = 0.0;
        observer_velocity[index] = 0.0;
        for (component = 0U; component < 3U; ++component) {
            observer_position[index] += rotation[index][component]
                * position[component];
            observer_velocity[index] += rotation[index][component]
                * velocity[component];
        }
    }

    output[0] = state[0] + observer_position[2];
    output[1] = observer_position[1] * angle_scale;
    output[2] = observer_position[0] * angle_scale;
    output[3] = state[4] + observer_velocity[2] - 1.0;
    for (index = 0U; index < RP_KERR_OBSERVABLE_OUTPUT_DIM; ++index) {
        if (!isfinite(output[index])) {
            clear_output(output);
            return RP_KERR_STATUS_NUMERICAL_RANGE;
        }
    }
    return RP_KERR_STATUS_OK;
}
