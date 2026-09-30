#include "geodesic/kerr_mcmc.h"
#include "geodesic/initial/convert.h"
#include "geodesic/physic/kerr_observables.h"

#include <assert.h>
#include <math.h>
#include <stddef.h>

static void check_orbit(double spin)
{
    const double identity[3][3] = {
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}
    };
    const double times[5] = {-1000.0, -100.0, 0.0, 100.0, 1000.0};
    const rp_kerr_mcmc_observation base = {
        0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    };
    const rp_kerr_mcmc_observation shifted = {
        0.01, 0.001, -0.002, 1e-6, -2e-6, 0.003, 20.0
    };
    double result[5 * 3];
    double shifted_result[5 * 3];
    rp_integrator_stats statistics;
    rp_integrator_config solver = {
        RP_INTEGRATOR_METHOD_DOP853, 8U, 1e-9, 1e-11, 0.0, 0.0, 100000U,
        NULL, NULL, NULL, NULL, NULL
    };
    size_t index;

    assert(rp_kerr_mcmc_evaluate(
        &solver, spin, identity,
        100.0, 0.2, 0.8, 0.3, 0.4,
        times, 5U, &base, result, &statistics
    ) == RP_KERR_MCMC_OK);
    assert(statistics.accepted_steps > 0U);
    assert(rp_kerr_mcmc_evaluate(
        &solver, spin, identity,
        100.0, 0.2, 0.8, 0.3, 0.4,
        times, 5U, &shifted, shifted_result, &statistics
    ) == RP_KERR_MCMC_OK);
    for (index = 0U; index < 5U; ++index) {
        assert(isfinite(result[3U * index + 0U]));
        assert(isfinite(result[3U * index + 1U]));
        assert(isfinite(result[3U * index + 2U]));
        assert(fabs(shifted_result[3U * index + 0U]
            - result[3U * index + 0U]
            - (0.001 + 1e-6 * (times[index] - 20.0))) < 1e-12);
        assert(fabs(shifted_result[3U * index + 1U]
            - result[3U * index + 1U]
            - (-0.002 - 2e-6 * (times[index] - 20.0))) < 1e-12);
        assert(fabs(shifted_result[3U * index + 2U]
            - result[3U * index + 2U] + 0.003) < 1e-12);
    }
}

static void check_osculating_initial_state(double spin)
{
    const double identity[3][3] = {
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}
    };
    const double elements[7] = {0.0, 100.0, 0.2, 0.8, 0.3, 0.4, 0.0};
    const rp_kerr_mcmc_observation observation = {
        0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    };
    rp_integrator_config solver = {
        RP_INTEGRATOR_METHOD_DOP853, 8U, 1e-9, 1e-11, 0.0, 0.0, 100000U,
        NULL, NULL, NULL, NULL, NULL
    };
    rp_integrator_stats statistics;
    double state[8];
    double projected[4];
    double output[3];

    assert(rp_initial_from_elements(spin, elements, state)
        == RP_KERR_STATUS_OK);
    assert(rp_kerr_observables_project(
        state, spin, identity, observation.angular_scale, projected
    ) == RP_KERR_STATUS_OK);
    assert(rp_kerr_mcmc_evaluate(
        &solver, spin, identity,
        elements[1], elements[2], elements[3], elements[4], elements[5],
        projected, 1U, &observation, output, &statistics
    ) == RP_KERR_MCMC_OK);
    assert(statistics.accepted_steps == 0U);
    assert(fabs(output[0] - projected[1]) < 1e-12);
    assert(fabs(output[1] - projected[2]) < 1e-12);
    assert(fabs(output[2] - projected[3]) < 1e-12);
}

static void check_invalid_times_clear_output(void)
{
    const double identity[3][3] = {
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}
    };
    const double unsorted[2] = {100.0, -100.0};
    const rp_kerr_mcmc_observation observation = {
        0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    };
    double result[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    rp_integrator_stats statistics;
    rp_integrator_config solver = {
        RP_INTEGRATOR_METHOD_DOP853, 8U, 1e-9, 1e-11, 0.0, 0.0, 100000U,
        NULL, NULL, NULL, NULL, NULL
    };
    size_t index;

    assert(rp_kerr_mcmc_evaluate(
        &solver, 0.0, identity,
        100.0, 0.2, 0.8, 0.3, 0.4,
        unsorted, 2U, &observation, result, &statistics
    ) == RP_KERR_MCMC_INVALID_INPUT);
    for (index = 0U; index < 6U; ++index) {
        assert(result[index] == 0.0);
    }
}

static void check_general_state_clocks(void)
{
    const double identity[3][3] = {
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}
    };
    const double elements[7] = {23.0, 100.0, 0.2, 0.8, 0.3, 0.4, 0.3};
    const rp_kerr_mcmc_time_kind kinds[3] = {
        RP_KERR_MCMC_TIME_COORDINATE,
        RP_KERR_MCMC_TIME_PROPER,
        RP_KERR_MCMC_TIME_ARRIVAL
    };
    rp_integrator_config solver = {
        RP_INTEGRATOR_METHOD_DOP853, 8U, 1e-9, 1e-11, 0.0, 0.0, 100000U,
        NULL, NULL, NULL, NULL, NULL
    };
    rp_integrator_stats statistics;
    double initial[8];
    double projected[4];
    double times[3];
    double output[9];
    double center;
    size_t kind;
    size_t component;

    assert(rp_initial_from_elements(0.2, elements, initial)
        == RP_KERR_STATUS_OK);
    assert(rp_kerr_observables_project(
        initial, 0.2, identity, 0.01, projected
    ) == RP_KERR_STATUS_OK);
    for (kind = 0U; kind < 3U; ++kind) {
        center = kinds[kind] == RP_KERR_MCMC_TIME_COORDINATE
            ? initial[0]
            : kinds[kind] == RP_KERR_MCMC_TIME_PROPER
                ? 7.0 : projected[0];
        times[0] = center - 100.0;
        times[1] = center;
        times[2] = center + 100.0;
        assert(rp_kerr_mcmc_evaluate_state(
            &solver, 0.2, identity, initial, 7.0, kinds[kind],
            times, 3U, 0.01, output, &statistics
        ) == RP_KERR_MCMC_OK);
        assert(statistics.accepted_steps > 0U);
        for (component = 0U; component < 3U; ++component) {
            assert(isfinite(output[component]));
            assert(isfinite(output[6U + component]));
            assert(fabs(output[3U + component]
                - projected[1U + component]) < 1e-12);
        }
        assert(fabs(output[0] - output[6]) > 1e-9);
    }
}

static void check_general_state_failures(void)
{
    const double identity[3][3] = {
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}
    };
    const double elements[7] = {0.0, 100.0, 0.2, 0.8, 0.3, 0.4, 0.0};
    const double unsorted[2] = {1.0, -1.0};
    const double centered[2] = {-1.0, 1.0};
    rp_integrator_config solver = {
        RP_INTEGRATOR_METHOD_DOP853, 8U, 1e-9, 1e-11, 0.0, 0.0, 100000U,
        NULL, NULL, NULL, NULL, NULL
    };
    rp_integrator_stats statistics;
    double initial[8];
    double output[6];
    size_t index;

    assert(rp_initial_from_elements(0.2, elements, initial)
        == RP_KERR_STATUS_OK);
    for (index = 0U; index < 6U; ++index) {
        output[index] = 1.0;
    }
    assert(rp_kerr_mcmc_evaluate_state(
        &solver, 0.2, identity, initial, 0.0,
        RP_KERR_MCMC_TIME_PROPER, unsorted, 2U, 0.01,
        output, &statistics
    ) == RP_KERR_MCMC_INVALID_INPUT);
    for (index = 0U; index < 6U; ++index) {
        assert(output[index] == 0.0);
    }
    initial[4] *= 0.9;
    for (index = 0U; index < 6U; ++index) {
        output[index] = 1.0;
    }
    assert(rp_kerr_mcmc_evaluate_state(
        &solver, 0.2, identity, initial, 0.0,
        RP_KERR_MCMC_TIME_COORDINATE, centered, 2U, 0.01,
        output, &statistics
    ) == RP_KERR_MCMC_INVALID_INPUT);
    for (index = 0U; index < 6U; ++index) {
        assert(output[index] == 0.0);
    }
}

static void check_state_arrival_matches_elements(void)
{
    const double identity[3][3] = {
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}
    };
    const double elements[7] = {0.0, 100.0, 0.2, 0.8, 0.3, 0.4, 0.0};
    const rp_kerr_mcmc_observation observation = {
        0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    };
    rp_integrator_config solver = {
        RP_INTEGRATOR_METHOD_DOP853, 8U, 1e-9, 1e-11, 0.0, 0.0, 100000U,
        NULL, NULL, NULL, NULL, NULL
    };
    rp_integrator_stats statistics;
    double initial[8];
    double projected[4];
    double times[3];
    double state_output[9];
    double element_output[9];
    size_t index;

    assert(rp_initial_from_elements(0.2, elements, initial)
        == RP_KERR_STATUS_OK);
    assert(rp_kerr_observables_project(
        initial, 0.2, identity, 0.01, projected
    ) == RP_KERR_STATUS_OK);
    times[0] = projected[0] - 100.0;
    times[1] = projected[0];
    times[2] = projected[0] + 100.0;
    assert(rp_kerr_mcmc_evaluate_state(
        &solver, 0.2, identity, initial, 0.0,
        RP_KERR_MCMC_TIME_ARRIVAL, times, 3U, 0.01,
        state_output, &statistics
    ) == RP_KERR_MCMC_OK);
    assert(rp_kerr_mcmc_evaluate(
        &solver, 0.2, identity, elements[1], elements[2],
        elements[3], elements[4], elements[5], times, 3U,
        &observation, element_output, &statistics
    ) == RP_KERR_MCMC_OK);
    for (index = 0U; index < 9U; ++index) {
        assert(fabs(state_output[index] - element_output[index]) < 1e-12);
    }
}

int main(void)
{
    check_orbit(0.0);
    check_orbit(0.2);
    check_osculating_initial_state(0.0);
    check_osculating_initial_state(0.2);
    check_invalid_times_clear_output();
    check_general_state_clocks();
    check_general_state_failures();
    check_state_arrival_matches_elements();
    return 0;
}
