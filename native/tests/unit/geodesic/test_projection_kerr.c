#include "geodesic/integrators/kerr.h"
#include "relatipy/kerr_geometry.h"

#include <assert.h>
#include <math.h>
#include <string.h>

static void make_state(
    double radius,
    double theta,
    const double coordinate_velocity[3],
    double state[8]
)
{
    state[0] = 0.0;
    state[1] = radius;
    state[2] = theta;
    state[3] = 0.0;
    assert(rp_kerr_four_velocity(
        1.0, 0.5, state, coordinate_velocity, state + 4
    ) == RP_KERR_STATUS_OK);
}

static void check_invariants(
    const rp_kerr_integrator_context *context,
    const double state[8]
)
{
    double metric[4][4];
    double momentum[4] = {0.0, 0.0, 0.0, 0.0};
    double norm = 0.0;
    double sine = sin(state[2]);
    double cosine = cos(state[2]);
    double energy;
    double angular_momentum;
    double carter_q;
    size_t mu;
    size_t nu;

    assert(rp_kerr_metric(
        context->mass, context->spin, state, metric
    ) == RP_KERR_STATUS_OK);
    for (mu = 0U; mu < 4U; ++mu) {
        for (nu = 0U; nu < 4U; ++nu) {
            momentum[mu] += metric[mu][nu] * state[4U + nu];
        }
        norm += momentum[mu] * state[4U + mu];
    }
    energy = -momentum[0];
    angular_momentum = momentum[3];
    carter_q = momentum[2] * momentum[2]
        + cosine * cosine * (context->spin * context->spin
            * (1.0 - energy * energy)
            + angular_momentum * angular_momentum / (sine * sine));
    assert(fabs(norm + 1.0) < 2.0e-11);
    assert(fabs(energy - context->E0) < 2.0e-11);
    assert(fabs(angular_momentum - context->Lz0) < 2.0e-11);
    assert(fabs(carter_q - context->Q0) < 2.0e-11);
}

static void test_general_and_signs(void)
{
    const double coordinate_velocity[3] = {-0.01, 0.002, 0.02};
    rp_kerr_integrator_context context = {.mass = 1.0, .spin = 0.5};
    double state[8];
    double position[4];

    make_state(8.0, 1.1, coordinate_velocity, state);
    assert(rp_kerr_integrator_prepare_projection(&context, state) == 0);
    memcpy(position, state, sizeof(position));
    state[4] *= 1.01;
    state[5] = fabs(state[5]) * 1.1;
    state[6] = -fabs(state[6]) * 0.9;
    state[7] *= 0.98;
    assert(rp_kerr_integrator_project(0.1, state, &context) == 0);
    assert(memcmp(position, state, sizeof(position)) == 0);
    assert(state[5] > 0.0 && state[6] < 0.0);
    check_invariants(&context, state);
}

static void test_turning_points(void)
{
    const double radial_turn[3] = {0.0, 0.002, 0.02};
    const double polar_turn[3] = {-0.01, 0.0, 0.02};
    const double circular[3] = {
        0.0, 0.0, 1.0 / (sqrt(8.0 * 8.0 * 8.0) + 0.5)
    };
    const double pi = acos(-1.0);
    const double *velocities[3] = {radial_turn, polar_turn, circular};
    const double angles[3] = {1.1, 1.1, 0.5 * pi};
    size_t index;

    for (index = 0U; index < 3U; ++index) {
        rp_kerr_integrator_context context = {.mass = 1.0, .spin = 0.5};
        double state[8];

        make_state(8.0, angles[index], velocities[index], state);
        assert(rp_kerr_integrator_prepare_projection(&context, state) == 0);
        state[4] *= 1.001;
        state[7] *= 0.999;
        assert(rp_kerr_integrator_project(0.0, state, &context) == 0);
        check_invariants(&context, state);
        if (index == 0U || index == 2U) {
            assert(fabs(state[5]) < 2.0e-7);
        }
        if (index == 1U || index == 2U) {
            assert(fabs(state[6]) < 2.0e-8);
        }
    }
}

static void test_invalid_and_horizon(void)
{
    const double coordinate_velocity[3] = {-0.01, 0.002, 0.02};
    rp_kerr_integrator_context context = {.mass = 1.0, .spin = 0.5};
    rp_kerr_integrator_context unprepared = {.mass = 1.0, .spin = 0.5};
    double state[8];
    double before[8];
    double saved_q;

    make_state(8.0, 1.1, coordinate_velocity, state);
    state[4] *= 1.2;
    assert(rp_kerr_integrator_prepare_projection(&context, state)
        == RP_KERR_STATUS_NON_TIMELIKE_VELOCITY);
    assert(context.projection_ready == 0);
    state[4] /= 1.2;
    assert(rp_kerr_integrator_prepare_projection(&context, state) == 0);
    saved_q = context.Q0;
    state[4] *= 1.2;
    assert(rp_kerr_integrator_prepare_projection(&context, state)
        == RP_KERR_STATUS_NON_TIMELIKE_VELOCITY);
    assert(context.projection_ready == 1 && context.Q0 == saved_q);
    state[4] /= 1.2;

    memcpy(before, state, sizeof(state));
    assert(rp_kerr_integrator_project(0.0, state, &unprepared)
        == RP_KERR_STATUS_INVALID_PARAMETER);
    assert(memcmp(before, state, sizeof(state)) == 0);
    saved_q = context.Q0;
    context.Q0 = -100.0;
    assert(rp_kerr_integrator_project(0.0, state, &context)
        == RP_KERR_STATUS_NON_TIMELIKE_VELOCITY);
    assert(memcmp(before, state, sizeof(state)) == 0);
    context.Q0 = saved_q;

    state[2] = 0.0;
    memcpy(before, state, sizeof(state));
    assert(rp_kerr_integrator_project(0.0, state, &context)
        == RP_KERR_STATUS_COORDINATE_SINGULARITY);
    assert(memcmp(before, state, sizeof(state)) == 0);

    state[2] = 1.1;
    state[1] = 1.5;
    memcpy(before, state, sizeof(state));
    assert(rp_kerr_integrator_project(0.0, state, &context) == 0);
    assert(memcmp(before, state, sizeof(state)) == 0);
    assert(rp_kerr_integrator_prepare_projection(&unprepared, state)
        == RP_KERR_STATUS_INVALID_PARAMETER);
    assert(unprepared.projection_ready == 0);

    state[5] = NAN;
    memcpy(before, state, sizeof(state));
    assert(rp_kerr_integrator_project(0.0, state, &context)
        == RP_KERR_STATUS_NONFINITE_INPUT);
    assert(memcmp(before, state, sizeof(state)) == 0);
}

static void test_independent_contexts(void)
{
    const double first_velocity[3] = {-0.01, 0.002, 0.02};
    const double second_velocity[3] = {0.015, -0.001, -0.012};
    rp_kerr_integrator_context first = {.mass = 1.0, .spin = 0.5};
    rp_kerr_integrator_context second = {.mass = 1.0, .spin = 0.5};
    double first_state[8];
    double second_state[8];

    make_state(8.0, 1.1, first_velocity, first_state);
    make_state(10.0, 0.9, second_velocity, second_state);
    assert(rp_kerr_integrator_prepare_projection(&first, first_state) == 0);
    assert(rp_kerr_integrator_prepare_projection(&second, second_state) == 0);
    first_state[4] *= 1.01;
    second_state[4] *= 0.99;
    assert(rp_kerr_integrator_project(0.0, second_state, &second) == 0);
    assert(rp_kerr_integrator_project(0.0, first_state, &first) == 0);
    check_invariants(&first, first_state);
    check_invariants(&second, second_state);
}

int main(void)
{
    test_general_and_signs();
    test_turning_points();
    test_invalid_and_horizon();
    test_independent_contexts();
    return 0;
}
