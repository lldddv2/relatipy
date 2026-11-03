#include "geodesic/integrators/integrator.h"

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <string.h>

struct linear_context {
    double rate;
    double fail_after;
};

static int identity_projector(
    double independent_variable, double state[], void *context
)
{
    (void)independent_variable;
    (void)state;
    (void)context;
    return 0;
}

static int fail_projector(
    double independent_variable, double state[], void *context
)
{
    (void)state;
    (void)context;
    return independent_variable > 0.15 ? -1 : 0;
}

static int linear_rhs(
    double independent_variable,
    const double state[],
    double derivative[],
    void *opaque_context
)
{
    const struct linear_context *context = opaque_context;

    if (independent_variable > context->fail_after) {
        return 1;
    }
    derivative[0] = context->rate * state[0];
    return 0;
}

static rp_integrator_config make_config(rp_integrator_method method)
{
    const rp_integrator_config config = {
        method,
        1U,
        1.0e-8,
        1.0e-11,
        0.0,
        0.0,
        100000U,
        NULL,
        NULL,
        NULL,
        NULL,
        rp_integrator_method_projects(method) ? identity_projector : NULL
    };
    return config;
}

static double integrate_exponential(
    rp_integrator_method method,
    double relative_tolerance,
    double absolute_tolerance,
    rp_integrator_stats *stats
)
{
    struct linear_context context = {1.0, INFINITY};
    rp_integrator_config config = make_config(method);
    double state[1] = {1.0};

    config.relative_tolerance = relative_tolerance;
    config.absolute_tolerance = absolute_tolerance;
    assert(rp_integrator_integrate(
        &config, linear_rhs, &context, 0.0, 1.0, state, stats
    ) == RP_INTEGRATOR_STATUS_OK);
    return state[0];
}

static void test_registry_has_six_methods(void)
{
    assert(RP_INTEGRATOR_METHOD_COUNT == 6);
    assert(strcmp(
        rp_integrator_method_name(RP_INTEGRATOR_METHOD_DOP853), "dop853"
    ) == 0);
    assert(strcmp(
        rp_integrator_method_name(RP_INTEGRATOR_METHOD_RADAU), "radau"
    ) == 0);
    assert(strcmp(
        rp_integrator_method_name(RP_INTEGRATOR_METHOD_DP45), "dp45"
    ) == 0);
    assert(strcmp(
        rp_integrator_method_name(RP_INTEGRATOR_METHOD_PROJECTION_RADAU),
        "projection_radau"
    ) == 0);
    assert(strcmp(
        rp_integrator_method_name(RP_INTEGRATOR_METHOD_PROJECTION_DOP853),
        "projection_dop853"
    ) == 0);
    assert(strcmp(
        rp_integrator_method_name(RP_INTEGRATOR_METHOD_PROJECTION_DP45),
        "projection_dp45"
    ) == 0);
    assert(!rp_integrator_method_projects(RP_INTEGRATOR_METHOD_DOP853));
    assert(!rp_integrator_method_projects(RP_INTEGRATOR_METHOD_RADAU));
    assert(!rp_integrator_method_projects(RP_INTEGRATOR_METHOD_DP45));
    assert(rp_integrator_method_name((rp_integrator_method)99) == NULL);
}

static void test_all_methods_integrate_forward_and_backward(void)
{
    rp_integrator_method method;

    for (method = RP_INTEGRATOR_METHOD_DOP853;
         method < RP_INTEGRATOR_METHOD_COUNT;
         method = (rp_integrator_method)((int)method + 1)) {
        struct linear_context context = {1.0, INFINITY};
        rp_integrator_config config = make_config(method);
        rp_integrator_stats stats;
        double state[1] = {1.0};

        assert(rp_integrator_integrate(
            &config, linear_rhs, &context, 0.0, 1.0, state, &stats
        ) == RP_INTEGRATOR_STATUS_OK);
        assert(fabs(state[0] - exp(1.0)) < 2.0e-7);
        assert(stats.accepted_steps > 0U);
        assert(stats.rhs_evaluations > 0U);
        assert(stats.final_independent_variable == 1.0);

        assert(rp_integrator_integrate(
            &config, linear_rhs, &context, 1.0, 0.0, state, &stats
        ) == RP_INTEGRATOR_STATUS_OK);
        assert(fabs(state[0] - 1.0) < 4.0e-7);
        assert(stats.final_independent_variable == 0.0);
    }
}

static void test_tighter_tolerance_improves_each_method(void)
{
    rp_integrator_method method;

    for (method = RP_INTEGRATOR_METHOD_DOP853;
         method < RP_INTEGRATOR_METHOD_COUNT;
         method = (rp_integrator_method)((int)method + 1)) {
        const double loose = integrate_exponential(
            method, 1.0e-3, 1.0e-6, NULL
        );
        const double tight = integrate_exponential(
            method, 1.0e-10, 1.0e-13, NULL
        );
        const double loose_error = fabs(loose - exp(1.0));
        const double tight_error = fabs(tight - exp(1.0));

        assert(tight_error < loose_error || tight_error < 5.0e-13);
    }
}

static double integrate_with_fixed_cap(
    rp_integrator_method method,
    double maximum_step
)
{
    struct linear_context context = {1.0, INFINITY};
    rp_integrator_config config = make_config(method);
    double state[1] = {1.0};

    config.relative_tolerance = 1.0e6;
    config.absolute_tolerance = 1.0e6;
    config.initial_step = maximum_step;
    config.maximum_step = maximum_step;
    assert(rp_integrator_integrate(
        &config, linear_rhs, &context, 0.0, 2.0, state, NULL
    ) == RP_INTEGRATOR_STATUS_OK);
    return state[0];
}

static void test_step_refinement_reduces_global_error(void)
{
    rp_integrator_method method;

    for (method = RP_INTEGRATOR_METHOD_DOP853;
         method < RP_INTEGRATOR_METHOD_COUNT;
         method = (rp_integrator_method)((int)method + 1)) {
        const double coarse = integrate_with_fixed_cap(method, 0.5);
        const double fine = integrate_with_fixed_cap(method, 0.25);
        const double coarse_error = fabs(coarse - exp(2.0));
        const double fine_error = fabs(fine - exp(2.0));

        assert(fine_error < 0.1 * coarse_error);
    }
}

static void test_dp45_has_fifth_order_global_convergence(void)
{
    const double exact = exp(2.0);
    const double coarse_error = fabs(
        integrate_with_fixed_cap(RP_INTEGRATOR_METHOD_DP45, 0.25) - exact
    );
    const double medium_error = fabs(
        integrate_with_fixed_cap(RP_INTEGRATOR_METHOD_DP45, 0.125) - exact
    );
    const double fine_error = fabs(
        integrate_with_fixed_cap(RP_INTEGRATOR_METHOD_DP45, 0.0625) - exact
    );

    assert(coarse_error > 20.0 * medium_error);
    assert(medium_error > 20.0 * fine_error);
}

static void test_dp45_failures_and_last_accepted_state(void)
{
    struct linear_context context = {1.0, INFINITY};
    rp_integrator_config config = make_config(RP_INTEGRATOR_METHOD_DP45);
    rp_integrator_stats stats;
    double state[1] = {1.0};

    config.initial_step = 0.1;
    config.maximum_step = 0.1;
    config.maximum_steps = 1U;
    assert(rp_integrator_integrate(
        &config, linear_rhs, &context, 0.0, 1.0, state, &stats
    ) == RP_INTEGRATOR_STATUS_MAXIMUM_STEPS);
    assert(stats.accepted_steps == 1U);
    assert(stats.final_independent_variable == 0.1);
    assert(fabs(state[0] - exp(0.1)) < 1.0e-8);

    state[0] = 1.0;
    config.maximum_steps = 1000U;
    context.fail_after = 0.25;
    assert(rp_integrator_integrate(
        &config, linear_rhs, &context, 0.0, 1.0, state, &stats
    ) == RP_INTEGRATOR_STATUS_RHS_FAILURE);
    assert(stats.final_independent_variable <= 0.25);
    assert(fabs(state[0] - exp(stats.final_independent_variable)) < 1.0e-8);
}

static void test_radau_handles_a_stiff_decay(void)
{
    struct linear_context context = {-1000.0, INFINITY};
    rp_integrator_config config = make_config(RP_INTEGRATOR_METHOD_RADAU);
    rp_integrator_stats stats;
    double state[1] = {1.0};

    config.relative_tolerance = 1.0e-7;
    config.absolute_tolerance = 1.0e-10;
    assert(rp_integrator_integrate(
        &config, linear_rhs, &context, 0.0, 1.0, state, &stats
    ) == RP_INTEGRATOR_STATUS_OK);
    assert(fabs(state[0]) < 2.0e-9);
    assert(stats.jacobian_evaluations > 0U);
    assert(stats.linear_solves > 0U);
}

static void test_validation_and_rhs_failures_are_typed(void)
{
    struct linear_context context = {1.0, INFINITY};
    rp_integrator_config config = make_config(RP_INTEGRATOR_METHOD_DOP853);
    rp_integrator_stats stats;
    double state[1] = {1.0};

    config.method = (rp_integrator_method)99;
    assert(rp_integrator_integrate(
        &config, linear_rhs, &context, 0.0, 1.0, state, &stats
    ) == RP_INTEGRATOR_STATUS_INVALID_ARGUMENT);
    assert(state[0] == 1.0);

    config = make_config(RP_INTEGRATOR_METHOD_DOP853);
    context.fail_after = 0.25;
    assert(rp_integrator_integrate(
        &config, linear_rhs, &context, 0.0, 1.0, state, &stats
    ) == RP_INTEGRATOR_STATUS_RHS_FAILURE);
    assert(isfinite(state[0]));
    assert(stats.final_independent_variable <= 0.25);

    config.relative_tolerance = NAN;
    assert(rp_integrator_integrate(
        &config, linear_rhs, &context, 0.0, 1.0, state, &stats
    ) == RP_INTEGRATOR_STATUS_NONFINITE_VALUE);
}

static const rp_integrator_method PROJECTION_METHODS[3] = {
    RP_INTEGRATOR_METHOD_PROJECTION_RADAU,
    RP_INTEGRATOR_METHOD_PROJECTION_DOP853,
    RP_INTEGRATOR_METHOD_PROJECTION_DP45
};

static void test_projection_requires_callback_and_preserves_last_step(void)
{
    size_t index;

    for (index = 0U; index < 3U; ++index) {
        struct linear_context context = {1.0, INFINITY};
        rp_integrator_config config = make_config(PROJECTION_METHODS[index]);
        rp_integrator_stats stats;
        double state[1] = {1.0};

        assert(rp_integrator_method_projects(PROJECTION_METHODS[index]));
        config.projector = NULL;
        assert(rp_integrator_integrate(
            &config, linear_rhs, &context, 0.0, 1.0, state, &stats
        ) == RP_INTEGRATOR_STATUS_INVALID_ARGUMENT);
        assert(state[0] == 1.0);
        config.projector = fail_projector;
        config.relative_tolerance = 1.0e-5;
        config.absolute_tolerance = 1.0e-8;
        config.initial_step = 0.1;
        config.maximum_step = 0.1;
        assert(rp_integrator_integrate(
            &config, linear_rhs, &context, 0.0, 1.0, state, &stats
        ) == RP_INTEGRATOR_STATUS_PROJECTION_FAILURE);
        assert(stats.accepted_steps == 1U);
        assert(fabs(stats.final_independent_variable - 0.1) < 1.0e-14);
        assert(fabs(state[0] - exp(0.1)) < 1.0e-7);
    }
}

static int oscillator_rhs(
    double independent_variable,
    const double state[],
    double derivative[],
    void *context
)
{
    (void)independent_variable;
    (void)context;
    derivative[0] = state[1];
    derivative[1] = -state[0];
    return 0;
}

/* Rescale (x, v) to the unit circle, the oscillator's energy surface. */
static int circle_projector(
    double independent_variable, double state[], void *context
)
{
    const double radius = hypot(state[0], state[1]);

    (void)independent_variable;
    (void)context;
    state[0] /= radius;
    state[1] /= radius;
    return 0;
}

static void test_explicit_projection_keeps_invariant_and_refreshes_rhs(void)
{
    const rp_integrator_method plain[2] = {
        RP_INTEGRATOR_METHOD_DOP853, RP_INTEGRATOR_METHOD_DP45
    };
    const rp_integrator_method projected[2] = {
        RP_INTEGRATOR_METHOD_PROJECTION_DOP853,
        RP_INTEGRATOR_METHOD_PROJECTION_DP45
    };
    size_t index;

    for (index = 0U; index < 2U; ++index) {
        rp_integrator_config config = make_config(plain[index]);
        rp_integrator_stats plain_stats;
        rp_integrator_stats identity_stats;
        rp_integrator_stats circle_stats;
        double plain_state[2] = {1.0, 0.0};
        double identity_state[2] = {1.0, 0.0};
        double circle_state[2] = {1.0, 0.0};

        config.dimension = 2U;
        config.relative_tolerance = 1.0e-6;
        config.absolute_tolerance = 1.0e-9;
        assert(rp_integrator_integrate(
            &config, oscillator_rhs, NULL, 0.0, 20.0, plain_state,
            &plain_stats
        ) == RP_INTEGRATOR_STATUS_OK);

        /* An identity projection takes the same steps plus one RHS
         * evaluation per accepted step, at the projected state. */
        config.method = projected[index];
        config.projector = identity_projector;
        assert(rp_integrator_integrate(
            &config, oscillator_rhs, NULL, 0.0, 20.0, identity_state,
            &identity_stats
        ) == RP_INTEGRATOR_STATUS_OK);
        assert(identity_stats.accepted_steps == plain_stats.accepted_steps);
        assert(identity_stats.rejected_steps == plain_stats.rejected_steps);
        assert(identity_stats.rhs_evaluations
            == plain_stats.rhs_evaluations + plain_stats.accepted_steps);
        assert(identity_state[0] == plain_state[0]);
        assert(identity_state[1] == plain_state[1]);

        config.projector = circle_projector;
        assert(rp_integrator_integrate(
            &config, oscillator_rhs, NULL, 0.0, 20.0, circle_state,
            &circle_stats
        ) == RP_INTEGRATOR_STATUS_OK);
        assert(fabs(hypot(circle_state[0], circle_state[1]) - 1.0)
            < 4.0 * DBL_EPSILON);
        assert(fabs(circle_state[0] - cos(20.0)) < 1.0e-4);
        assert(fabs(circle_state[1] + sin(20.0)) < 1.0e-4);
    }
}

static void test_calls_keep_context_and_counters_independent(void)
{
    struct linear_context growth = {2.0, INFINITY};
    struct linear_context decay = {-2.0, INFINITY};
    rp_integrator_config config = make_config(RP_INTEGRATOR_METHOD_DOP853);
    rp_integrator_stats growth_stats;
    rp_integrator_stats decay_stats;
    double growth_state[1] = {1.0};
    double decay_state[1] = {1.0};

    assert(rp_integrator_integrate(
        &config, linear_rhs, &growth, 0.0, 0.5, growth_state, &growth_stats
    ) == RP_INTEGRATOR_STATUS_OK);
    config.method = RP_INTEGRATOR_METHOD_RADAU;
    assert(rp_integrator_integrate(
        &config, linear_rhs, &decay, 0.0, 0.5, decay_state, &decay_stats
    ) == RP_INTEGRATOR_STATUS_OK);
    assert(fabs(growth_state[0] - exp(1.0)) < 2.0e-7);
    assert(fabs(decay_state[0] - exp(-1.0)) < 2.0e-7);
    assert(growth_stats.jacobian_evaluations == 0U);
    assert(decay_stats.jacobian_evaluations > 0U);
}

static void test_dp45_calls_keep_context_and_counters_independent(void)
{
    struct linear_context growth = {2.0, INFINITY};
    struct linear_context decay = {-2.0, INFINITY};
    rp_integrator_config config = make_config(RP_INTEGRATOR_METHOD_DP45);
    rp_integrator_stats growth_stats;
    rp_integrator_stats decay_stats;
    double growth_state[1] = {1.0};
    double decay_state[1] = {1.0};

    assert(rp_integrator_integrate(
        &config, linear_rhs, &growth, 0.0, 0.5, growth_state, &growth_stats
    ) == RP_INTEGRATOR_STATUS_OK);
    assert(rp_integrator_integrate(
        &config, linear_rhs, &decay, 0.0, 0.5, decay_state, &decay_stats
    ) == RP_INTEGRATOR_STATUS_OK);
    assert(fabs(growth_state[0] - exp(1.0)) < 2.0e-7);
    assert(fabs(decay_state[0] - exp(-1.0)) < 2.0e-7);
    assert(growth_stats.accepted_steps > 0U);
    assert(decay_stats.accepted_steps > 0U);
    assert(growth_stats.jacobian_evaluations == 0U);
    assert(decay_stats.jacobian_evaluations == 0U);
    assert(growth_stats.final_independent_variable == 0.5);
    assert(decay_stats.final_independent_variable == 0.5);
}

int main(void)
{
    test_registry_has_six_methods();
    test_all_methods_integrate_forward_and_backward();
    test_tighter_tolerance_improves_each_method();
    test_step_refinement_reduces_global_error();
    test_dp45_has_fifth_order_global_convergence();
    test_dp45_failures_and_last_accepted_state();
    test_radau_handles_a_stiff_decay();
    test_validation_and_rhs_failures_are_typed();
    test_projection_requires_callback_and_preserves_last_step();
    test_explicit_projection_keeps_invariant_and_refreshes_rhs();
    test_calls_keep_context_and_counters_independent();
    test_dp45_calls_keep_context_and_counters_independent();
    return 0;
}
