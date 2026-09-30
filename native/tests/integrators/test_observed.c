#include "geodesic/integrators/integrator.h"

#include <assert.h>
#include <math.h>
#include <stddef.h>

typedef struct {
    size_t count;
    double last_tau;
    double last_state[2];
    int stop_after;
    int fail_after;
} observation;

static int linear_rhs(
    double tau,
    const double state[],
    double derivative[],
    void *context
)
{
    (void)tau;
    (void)context;
    derivative[0] = state[0];
    derivative[1] = -2.0 * state[1];
    return 0;
}

static int observe(double tau, const double state[], void *context)
{
    observation *record = context;

    assert(tau > record->last_tau);
    assert(isfinite(state[0]) && isfinite(state[1]));
    ++record->count;
    record->last_tau = tau;
    record->last_state[0] = state[0];
    record->last_state[1] = state[1];
    if (record->fail_after > 0
        && record->count == (size_t)record->fail_after) {
        return -1;
    }
    if (record->stop_after > 0
        && record->count == (size_t)record->stop_after) {
        return 1;
    }
    return 0;
}

static int identity_projector(double tau, double state[], void *context)
{
    (void)tau;
    (void)state;
    (void)context;
    return 0;
}

static rp_integrator_config make_config(rp_integrator_method method)
{
    rp_integrator_config config = {0};

    config.method = method;
    config.dimension = 2U;
    config.relative_tolerance = 1.0e-8;
    config.absolute_tolerance = 1.0e-11;
    config.initial_step = 0.05;
    config.maximum_step = 0.05;
    config.maximum_steps = 100000U;
    if (method == RP_INTEGRATOR_METHOD_PROJECTION_RADAU) {
        config.projector = identity_projector;
    }
    return config;
}

static void check_observer_and_vector_tolerance(rp_integrator_method method)
{
    const double vector_atol[2] = {1.0e-11, 1.0e-7};
    rp_integrator_config config = make_config(method);
    rp_integrator_stats stats;
    observation record = {0};
    double state[2] = {1.0, 1.0};

    config.absolute_tolerances = vector_atol;
    config.step_observer = observe;
    config.step_observer_context = &record;
    assert(rp_integrator_integrate(
        &config, linear_rhs, NULL, 0.0, 1.0, state, &stats
    ) == RP_INTEGRATOR_STATUS_OK);
    assert(record.count == stats.accepted_steps);
    assert(record.count > 0U);
    assert(record.last_tau == 1.0);
    assert(stats.final_independent_variable == 1.0);
    assert(record.last_state[0] == state[0]);
    assert(record.last_state[1] == state[1]);
    assert(fabs(state[0] - exp(1.0)) < 1.0e-6);
    assert(fabs(state[1] - exp(-2.0)) < 1.0e-6);

    state[0] = 1.0;
    state[1] = 1.0;
    record = (observation){0};
    record.stop_after = 2;
    assert(rp_integrator_integrate(
        &config, linear_rhs, NULL, 0.0, 1.0, state, &stats
    ) == RP_INTEGRATOR_STATUS_OBSERVER_STOPPED);
    assert(record.count == 2U);
    assert(stats.accepted_steps == 2U);
    assert(stats.final_independent_variable == record.last_tau);
    assert(state[0] == record.last_state[0]);

    state[0] = 1.0;
    state[1] = 1.0;
    record = (observation){0};
    record.fail_after = 1;
    assert(rp_integrator_integrate(
        &config, linear_rhs, NULL, 0.0, 1.0, state, &stats
    ) == RP_INTEGRATOR_STATUS_OBSERVER_FAILURE);
    assert(record.count == 1U);
    assert(stats.accepted_steps == 1U);
    assert(stats.final_independent_variable == record.last_tau);
    assert(state[0] == record.last_state[0]);

    record = (observation){0};
    assert(rp_integrator_integrate(
        &config, linear_rhs, NULL, 0.0, 0.0, state, &stats
    ) == RP_INTEGRATOR_STATUS_OK);
    assert(record.count == 0U);
    assert(stats.accepted_steps == 0U);
}

static void check_invalid_vector_tolerance(void)
{
    double vector_atol[2] = {1.0e-10, -1.0};
    rp_integrator_config config = make_config(RP_INTEGRATOR_METHOD_DP45);
    rp_integrator_stats stats;
    double state[2] = {1.0, 1.0};

    config.absolute_tolerances = vector_atol;
    assert(rp_integrator_integrate(
        &config, linear_rhs, NULL, 0.0, 1.0, state, &stats
    ) == RP_INTEGRATOR_STATUS_INVALID_ARGUMENT);
    assert(stats.rhs_evaluations == 0U);
    vector_atol[1] = NAN;
    assert(rp_integrator_integrate(
        &config, linear_rhs, NULL, 0.0, 1.0, state, &stats
    ) == RP_INTEGRATOR_STATUS_NONFINITE_VALUE);
    vector_atol[1] = 0.0;
    assert(rp_integrator_integrate(
        &config, linear_rhs, NULL, 0.0, 0.0, state, &stats
    ) == RP_INTEGRATOR_STATUS_OK);
}

static void check_vector_changes_controller(rp_integrator_method method)
{
    const double tighter_atol[2] = {1.0e-3, 1.0e-12};
    rp_integrator_config config = make_config(method);
    rp_integrator_stats loose_stats;
    rp_integrator_stats tight_stats;
    double state[2] = {1.0, 1.0};

    config.relative_tolerance = 1.0e-12;
    config.absolute_tolerance = 1.0e-3;
    config.initial_step = 0.5;
    config.maximum_step = 0.0;
    assert(rp_integrator_integrate(
        &config, linear_rhs, NULL, 0.0, 2.0, state, &loose_stats
    ) == RP_INTEGRATOR_STATUS_OK);
    state[0] = 1.0;
    state[1] = 1.0;
    config.absolute_tolerances = tighter_atol;
    assert(rp_integrator_integrate(
        &config, linear_rhs, NULL, 0.0, 2.0, state, &tight_stats
    ) == RP_INTEGRATOR_STATUS_OK);
    assert(tight_stats.accepted_steps > loose_stats.accepted_steps);
}

int main(void)
{
    rp_integrator_method method;

    for (method = RP_INTEGRATOR_METHOD_DOP853;
         method < RP_INTEGRATOR_METHOD_COUNT;
         method = (rp_integrator_method)((int)method + 1)) {
        check_observer_and_vector_tolerance(method);
        check_vector_changes_controller(method);
    }
    check_invalid_vector_tolerance();
    return 0;
}
