#include "geodesic/physic/kerr_observables.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stddef.h>

static int close_to(double actual, double expected, double tolerance)
{
    return fabs(actual - expected) <= tolerance * (1.0 + fabs(expected));
}

static void test_schwarzschild_identity_projection(void)
{
    const double rotation[3][3] = {
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}
    };
    const double state[8] = {
        2.0, 10.0, 1.04719755119659774615, 0.52359877559829887308,
        1.1, 0.2, 0.03, 0.04
    };
    double output[4];
    const double x = state[1] * sin(state[2]) * cos(state[3]);
    const double y = state[1] * sin(state[2]) * sin(state[3]);
    const double z = state[1] * cos(state[2]);
    const double dz = state[5] * cos(state[2])
        - state[1] * sin(state[2]) * state[6];

    assert(rp_kerr_observables_project(
        state, 0.0, rotation, 0.02, output
    ) == RP_KERR_STATUS_OK);
    assert(close_to(output[0], state[0] + z, 1e-14));
    assert(close_to(output[1], 0.02 * y, 1e-14));
    assert(close_to(output[2], 0.02 * x, 1e-14));
    assert(close_to(output[3], state[4] + dz - 1.0, 1e-14));
}

static void test_rotated_velocity_matches_finite_difference(void)
{
    const double rotation[3][3] = {
        {0.0, 0.0, 1.0},
        {0.0, 1.0, 0.0},
        {-1.0, 0.0, 0.0}
    };
    const double state[8] = {
        2.0, 10.0, 1.1, 0.3, 1.15, -0.06, 0.002, 0.015
    };
    double plus_state[8];
    double minus_state[8];
    double output[4];
    double plus_output[4];
    double minus_output[4];
    const double step = 1e-5;
    size_t index;

    for (index = 0U; index < 8U; ++index) {
        plus_state[index] = state[index];
        minus_state[index] = state[index];
    }
    for (index = 0U; index < 4U; ++index) {
        plus_state[index] += step * state[index + 4U];
        minus_state[index] -= step * state[index + 4U];
    }

    assert(rp_kerr_observables_project(
        state, 0.6, rotation, 0.02, output
    ) == RP_KERR_STATUS_OK);
    assert(rp_kerr_observables_project(
        plus_state, 0.6, rotation, 0.02, plus_output
    ) == RP_KERR_STATUS_OK);
    assert(rp_kerr_observables_project(
        minus_state, 0.6, rotation, 0.02, minus_output
    ) == RP_KERR_STATUS_OK);
    assert(close_to(
        output[3] + 1.0,
        (plus_output[0] - minus_output[0]) / (2.0 * step),
        2e-10
    ));
    assert(close_to(
        output[1],
        0.02 * hypot(state[1], 0.6) * sin(state[2]) * sin(state[3]),
        1e-14
    ));
}

static void test_invalid_inputs_clear_output(void)
{
    const double identity[3][3] = {
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}
    };
    const double invalid_rotation[3][3] = {
        {2.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}
    };
    const double reflection[3][3] = {
        {-1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}
    };
    double state[8] = {0.0, 10.0, 1.0, 0.2, 1.0, 0.0, 0.0, 0.0};
    double output[4] = {9.0, 9.0, 9.0, 9.0};
    size_t index;

    assert(rp_kerr_observables_project(
        NULL, 0.0, identity, 1.0, output
    ) == RP_KERR_STATUS_NULL_POINTER);
    for (index = 0U; index < 4U; ++index) {
        assert(output[index] == 0.0);
    }
    assert(rp_kerr_observables_project(
        state, 0.0, identity, 1.0, NULL
    ) == RP_KERR_STATUS_NULL_POINTER);
    assert(rp_kerr_observables_project(
        state, -1.0, identity, 1.0, output
    ) == RP_KERR_STATUS_INVALID_PARAMETER);
    assert(rp_kerr_observables_project(
        state, 0.0, invalid_rotation, 1.0, output
    ) == RP_KERR_STATUS_INVALID_PARAMETER);
    assert(rp_kerr_observables_project(
        state, 0.0, reflection, 1.0, output
    ) == RP_KERR_STATUS_INVALID_PARAMETER);
    assert(rp_kerr_observables_project(
        state, 0.0, identity, 0.0, output
    ) == RP_KERR_STATUS_INVALID_PARAMETER);
    state[2] = NAN;
    assert(rp_kerr_observables_project(
        state, 0.0, identity, 1.0, output
    ) == RP_KERR_STATUS_NONFINITE_INPUT);
    for (index = 0U; index < 4U; ++index) {
        assert(output[index] == 0.0);
    }
    state[2] = 1.0;
    state[0] = DBL_MAX;
    state[1] = DBL_MAX;
    assert(rp_kerr_observables_project(
        state, 0.0, identity, 1.0, output
    ) == RP_KERR_STATUS_NUMERICAL_RANGE);
    for (index = 0U; index < 4U; ++index) {
        assert(output[index] == 0.0);
    }
}

int main(void)
{
    test_schwarzschild_identity_projection();
    test_rotated_velocity_matches_finite_difference();
    test_invalid_inputs_clear_output();
    return 0;
}
