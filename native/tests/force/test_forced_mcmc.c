/*
 * Forced observable model with an active force.
 *
 * Coordinate-time samples (the coordinate-time forced RHS) and proper-time
 * samples (the proper-time forced RHS) of the same forced orbit must give
 * the same observables, which validates dividing the force term by u^t.
 * The element evaluator agrees with the state evaluator at arrival times,
 * the force changes the observables, and a failing callback is reported.
 */

#include "force/kerr_forced.h"
#include "force/kerr_mcmc.h"
#include "force_test_models.h"
#include "geodesic/initial/convert.h"
#include "geodesic/physic/kerr_observables.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

#define SAMPLES 5U

static const double ELEMENTS[7] = {0.0, 80.0, 0.6, 0.9, 0.4, 2.0, 0.0};

static void check_clocks(double spin, const rp_external_force *force)
{
    const rp_integrator_config solver = {
        RP_INTEGRATOR_METHOD_DOP853, 8U, 1e-12, 1e-14, 0.0, 0.0, 1000000U,
        NULL, NULL, NULL, NULL, NULL
    };
    const double taus[SAMPLES] = {-4000.0, -900.0, 0.0, 1500.0, 6000.0};
    const double c = cos(0.5);
    const double s = sin(0.5);
    const double rotation[3][3] = {
        {1.0, 0.0, 0.0}, {0.0, c, -s}, {0.0, s, c}
    };
    double initial[8];
    double coordinate_times[SAMPLES];
    double from_proper[3 * SAMPLES];
    double from_coordinate[3 * SAMPLES];
    rp_integrator_stats stats;
    size_t k;

    assert(rp_initial_from_elements(spin, ELEMENTS, initial)
        == RP_KERR_STATUS_OK);
    /* Coordinate times of the proper-time samples along the forced orbit. */
    {
        rp_kerr_forced_context context;
        rp_integrator_config precise = solver;
        double state[8];
        size_t component;

        context.spin = spin;
        context.force = force;
        precise.relative_tolerance = 1e-13;
        for (k = 0U; k < SAMPLES; ++k) {
            for (component = 0U; component < 8U; ++component) {
                state[component] = initial[component];
            }
            if (taus[k] != 0.0) {
                assert(rp_integrator_integrate(
                    &precise, rp_kerr_forced_integrator_rhs, &context,
                    0.0, taus[k], state, NULL
                ) == RP_INTEGRATOR_STATUS_OK);
            }
            coordinate_times[k] = state[0];
        }
    }
    assert(rp_kerr_forced_mcmc_evaluate_state(
        &solver, spin, force, rotation, initial, 0.0,
        RP_KERR_MCMC_TIME_PROPER, taus, SAMPLES, 0.01, from_proper, &stats
    ) == RP_KERR_MCMC_OK);
    assert(rp_kerr_forced_mcmc_evaluate_state(
        &solver, spin, force, rotation, initial, 0.0,
        RP_KERR_MCMC_TIME_COORDINATE, coordinate_times, SAMPLES, 0.01,
        from_coordinate, &stats
    ) == RP_KERR_MCMC_OK);
    /* Integration error only; a mistransformed force term would differ
       at the force effect times (u^t - 1), far above this bound. */
    for (k = 0U; k < 3U * SAMPLES; ++k) {
        assert(fabs(from_proper[k] - from_coordinate[k]) < 1e-9);
    }
}

static void check_elements_and_effect(double spin, const rp_external_force *force)
{
    const rp_integrator_config solver = {
        RP_INTEGRATOR_METHOD_DOP853, 8U, 1e-11, 1e-13, 0.0, 0.0, 1000000U,
        NULL, NULL, NULL, NULL, NULL
    };
    const rp_kerr_mcmc_observation observation = {
        0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    };
    const double times[SAMPLES] = {-5000.0, -100.0, 300.0, 2500.0, 9000.0};
    const rp_external_force failing = {test_failing_acceleration, NULL};
    const double c = cos(1.1);
    const double s = sin(1.1);
    const double rotation[3][3] = {
        {1.0, 0.0, 0.0}, {0.0, c, -s}, {0.0, s, c}
    };
    double initial[8];
    double observer[7];
    double elements_out[3 * SAMPLES];
    double state_out[3 * SAMPLES];
    double free_out[3 * SAMPLES];
    double largest = 0.0;
    rp_integrator_stats stats;
    size_t k;

    /* The element evaluator rotates observer-frame elements into the body. */
    assert(rp_initial_elements_to_cartesian(ELEMENTS, observer)
        == RP_KERR_STATUS_OK);
    assert(rp_initial_observer_cartesian_to_canonical(
        spin, rotation, observer, initial
    ) == RP_KERR_STATUS_OK);
    assert(rp_kerr_forced_mcmc_evaluate(
        &solver, spin, force, rotation, ELEMENTS[1], ELEMENTS[2],
        ELEMENTS[3], ELEMENTS[4], ELEMENTS[5], times, SAMPLES,
        &observation, elements_out, &stats
    ) == RP_KERR_MCMC_OK);
    assert(stats.accepted_steps > 0U);
    assert(rp_kerr_forced_mcmc_evaluate_state(
        &solver, spin, force, rotation, initial, 0.0,
        RP_KERR_MCMC_TIME_ARRIVAL, times, SAMPLES, 0.01, state_out, &stats
    ) == RP_KERR_MCMC_OK);
    assert(rp_kerr_forced_mcmc_evaluate(
        &solver, spin, NULL, rotation, ELEMENTS[1], ELEMENTS[2],
        ELEMENTS[3], ELEMENTS[4], ELEMENTS[5], times, SAMPLES,
        &observation, free_out, &stats
    ) == RP_KERR_MCMC_OK);
    for (k = 0U; k < 3U * SAMPLES; ++k) {
        assert(isfinite(elements_out[k]));
        assert(fabs(elements_out[k] - state_out[k]) < 1e-10);
        largest = fmax(largest, fabs(elements_out[k] - free_out[k]));
    }
    assert(largest > 1e-3);

    elements_out[0] = 1.0;
    assert(rp_kerr_forced_mcmc_evaluate(
        &solver, spin, &failing, rotation, ELEMENTS[1], ELEMENTS[2],
        ELEMENTS[3], ELEMENTS[4], ELEMENTS[5], times, SAMPLES,
        &observation, elements_out, &stats
    ) == RP_KERR_MCMC_INTEGRATION_FAILURE);
    for (k = 0U; k < 3U * SAMPLES; ++k) {
        assert(elements_out[k] == 0.0);
    }
}

int main(void)
{
    const double axis[3] = {0.0, 0.6, 0.8};
    test_tilted tilted = {
        {0.05, 20.0, 3.0},
        {{0.0}},
        {0.02, 8.0, {12.0, -7.0, 5.0}},
        1e-6
    };
    const test_miyamoto_nagai disk = {0.05, 20.0, 3.0};
    const rp_external_force forces[2] = {
        {test_miyamoto_nagai_acceleration, &disk},
        {test_tilted_acceleration, &tilted}
    };
    const double spins[2] = {0.0, 0.9};
    size_t f;
    size_t s;

    test_axis_rotation(axis, 0.8, tilted.rotation);
    for (f = 0U; f < 2U; ++f) {
        for (s = 0U; s < 2U; ++s) {
            check_clocks(spins[s], &forces[f]);
            check_elements_and_effect(spins[s], &forces[f]);
        }
    }
    puts("test_forced_mcmc: ok");
    return 0;
}
