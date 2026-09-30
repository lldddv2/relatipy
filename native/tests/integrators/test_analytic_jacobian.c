#include "geodesic/integrators/integrator.h"

#include <assert.h>
#include <math.h>
#include <stddef.h>

typedef enum {
    JACOBIAN_OK,
    JACOBIAN_FAILURE,
    JACOBIAN_NAN,
    JACOBIAN_INFINITY,
    JACOBIAN_INCOMPLETE
} jacobian_behavior;

typedef struct {
    size_t dimension;
    double matrix[RP_INTEGRATOR_MAX_DIMENSION][RP_INTEGRATOR_MAX_DIMENSION];
    double constant[RP_INTEGRATOR_MAX_DIMENSION];
    size_t rhs_calls;
    size_t jacobian_calls;
    size_t successful_jacobians;
    size_t observed_steps;
    size_t fail_after_accepted;
    jacobian_behavior behavior;
    double last_tau;
    double last_state[RP_INTEGRATOR_MAX_DIMENSION];
} linear_context;

static int linear_rhs(
    double tau, const double state[], double derivative[], void *opaque
)
{
    linear_context *context = opaque;
    size_t row;
    size_t column;

    assert(isfinite(tau));
    ++context->rhs_calls;
    for (row = 0U; row < context->dimension; ++row) {
        derivative[row] = context->constant[row];
        for (column = 0U; column < context->dimension; ++column) {
            derivative[row] += context->matrix[row][column] * state[column];
        }
    }
    return 0;
}

static int linear_jacobian(
    double tau, const double state[], double *jacobian, size_t dimension,
    size_t row_stride, void *opaque
)
{
    linear_context *context = opaque;
    size_t row;
    size_t column;
    const int fail_now = context->behavior != JACOBIAN_OK
        && context->observed_steps >= context->fail_after_accepted;

    assert(isfinite(tau));
    assert(dimension == context->dimension);
    assert(row_stride == RP_INTEGRATOR_MAX_DIMENSION);
    ++context->jacobian_calls;
    if (fail_now && context->behavior == JACOBIAN_FAILURE) {
        return -7;
    }
    for (row = 0U; row < dimension; ++row) {
        assert(isfinite(state[row]));
        for (column = 0U; column < dimension; ++column) {
            if (fail_now && context->behavior == JACOBIAN_INCOMPLETE
                && row == dimension - 1U && column == dimension - 1U) {
                continue;
            }
            jacobian[row * row_stride + column] = context->matrix[row][column];
        }
    }
    if (fail_now && context->behavior == JACOBIAN_NAN) {
        jacobian[(dimension - 1U) * row_stride + dimension - 1U] = NAN;
    }
    if (fail_now && context->behavior == JACOBIAN_INFINITY) {
        jacobian[(dimension - 1U) * row_stride + dimension - 1U] = INFINITY;
    }
    if (!fail_now) {
        ++context->successful_jacobians;
    }
    return 0;
}

static int observe(double tau, const double state[], void *opaque)
{
    linear_context *context = opaque;
    size_t component;

    ++context->observed_steps;
    context->last_tau = tau;
    for (component = 0U; component < context->dimension; ++component) {
        context->last_state[component] = state[component];
    }
    return 0;
}

static rp_integrator_config make_config(size_t dimension)
{
    const rp_integrator_config config = {
        .method = RP_INTEGRATOR_METHOD_RADAU,
        .dimension = dimension,
        .relative_tolerance = 1.0e-8,
        .absolute_tolerance = 1.0e-11,
        .initial_step = 0.05,
        .maximum_step = 0.05,
        .maximum_steps = 100000U,
        .jacobian = linear_jacobian
    };
    return config;
}

static void check_accounting(
    const linear_context *context, const rp_integrator_stats *stats
)
{
    assert(stats->rhs_evaluations == context->rhs_calls);
    assert(stats->jacobian_evaluations == context->successful_jacobians);
}

static void test_scalar_forward_and_backward(void)
{
    linear_context context = {.dimension = 1U};
    rp_integrator_config config = make_config(1U);
    rp_integrator_stats stats;
    double state[1] = {1.0};

    context.matrix[0][0] = 1.0;
    assert(rp_integrator_integrate(
        &config, linear_rhs, &context, 0.0, 1.0, state, &stats
    ) == RP_INTEGRATOR_STATUS_OK);
    assert(fabs(state[0] - exp(1.0)) < 2.0e-7);
    assert(stats.accepted_steps > 0U);
    assert(context.jacobian_calls == 3U * (
        stats.accepted_steps + stats.rejected_steps
    ));
    assert(stats.final_independent_variable == 1.0);
    check_accounting(&context, &stats);

    context.rhs_calls = 0U;
    context.jacobian_calls = 0U;
    context.successful_jacobians = 0U;
    assert(rp_integrator_integrate(
        &config, linear_rhs, &context, 1.0, 0.0, state, &stats
    ) == RP_INTEGRATOR_STATUS_OK);
    assert(fabs(state[0] - 1.0) < 4.0e-7);
    assert(stats.final_independent_variable == 0.0);
    check_accounting(&context, &stats);
}

static void test_coupled_linear_solution_and_row_stride(void)
{
    linear_context context = {.dimension = 2U};
    rp_integrator_config config = make_config(2U);
    rp_integrator_stats stats;
    double state[2] = {1.0, 2.0};

    context.matrix[0][0] = -2.0;
    context.matrix[0][1] = 3.0;
    context.matrix[1][1] = -5.0;
    assert(rp_integrator_integrate(
        &config, linear_rhs, &context, 0.0, 1.0, state, &stats
    ) == RP_INTEGRATOR_STATUS_OK);
    assert(fabs(state[0] - (3.0 * exp(-2.0) - 2.0 * exp(-5.0))) < 2.0e-8);
    assert(fabs(state[1] - 2.0 * exp(-5.0)) < 2.0e-8);
    assert(context.jacobian_calls > 0U);
    check_accounting(&context, &stats);
}

static void test_full_dimension_storage(void)
{
    linear_context context = {.dimension = RP_INTEGRATOR_MAX_DIMENSION};
    rp_integrator_config config = make_config(RP_INTEGRATOR_MAX_DIMENSION);
    rp_integrator_stats stats;
    double state[RP_INTEGRATOR_MAX_DIMENSION];
    size_t component;

    for (component = 0U; component < context.dimension; ++component) {
        context.matrix[component][component] = -1.0;
        state[component] = (double)(component + 1U);
    }
    assert(rp_integrator_integrate(
        &config, linear_rhs, &context, 0.0, 0.2, state, &stats
    ) == RP_INTEGRATOR_STATUS_OK);
    for (component = 0U; component < context.dimension; ++component) {
        assert(fabs(state[component] - (double)(component + 1U) * exp(-0.2))
            < 3.0e-8);
    }
    check_accounting(&context, &stats);
}

static double exponential_with_cap(double cap)
{
    linear_context context = {.dimension = 1U};
    rp_integrator_config config = make_config(1U);
    double state[1] = {1.0};

    context.matrix[0][0] = 1.0;
    config.relative_tolerance = 1.0e6;
    config.absolute_tolerance = 1.0e6;
    config.initial_step = cap;
    config.maximum_step = cap;
    assert(rp_integrator_integrate(
        &config, linear_rhs, &context, 0.0, 2.0, state, NULL
    ) == RP_INTEGRATOR_STATUS_OK);
    return state[0];
}

static void test_fifth_order_step_refinement(void)
{
    const double coarse_error = fabs(exponential_with_cap(0.5) - exp(2.0));
    const double medium_error = fabs(exponential_with_cap(0.25) - exp(2.0));
    const double fine_error = fabs(exponential_with_cap(0.125) - exp(2.0));

    assert(coarse_error > 20.0 * medium_error);
    assert(medium_error > 20.0 * fine_error);
}

static void test_stiff_decay(void)
{
    linear_context context = {.dimension = 1U};
    rp_integrator_config config = make_config(1U);
    rp_integrator_stats stats;
    double state[1] = {1.0};

    context.matrix[0][0] = -1000.0;
    config.relative_tolerance = 1.0e-7;
    config.absolute_tolerance = 1.0e-10;
    config.initial_step = 0.0;
    config.maximum_step = 0.0;
    assert(rp_integrator_integrate(
        &config, linear_rhs, &context, 0.0, 1.0, state, &stats
    ) == RP_INTEGRATOR_STATUS_OK);
    assert(fabs(state[0]) < 2.0e-9);
    assert(stats.linear_solves > 0U);
    check_accounting(&context, &stats);
}

static void test_fd_selection_and_rhs_accounting(void)
{
    linear_context analytic_context = {.dimension = 2U};
    linear_context fd_context = {.dimension = 2U};
    rp_integrator_config config = make_config(2U);
    rp_integrator_stats analytic_stats;
    rp_integrator_stats fd_stats;
    double analytic_state[2] = {1.0, 2.0};
    double fd_state[2] = {1.0, 2.0};

    /* Constant derivatives give identical stage iterations with analytic/FD J. */
    analytic_context.constant[0] = fd_context.constant[0] = 1.0;
    analytic_context.constant[1] = fd_context.constant[1] = -2.0;
    assert(rp_integrator_integrate(
        &config, linear_rhs, &analytic_context, 0.0, 0.2,
        analytic_state, &analytic_stats
    ) == RP_INTEGRATOR_STATUS_OK);
    config.jacobian = NULL;
    assert(rp_integrator_integrate(
        &config, linear_rhs, &fd_context, 0.0, 0.2, fd_state, &fd_stats
    ) == RP_INTEGRATOR_STATUS_OK);
    assert(analytic_state[0] == fd_state[0]);
    assert(analytic_state[1] == fd_state[1]);
    assert(analytic_stats.accepted_steps == fd_stats.accepted_steps);
    assert(analytic_stats.rejected_steps == fd_stats.rejected_steps);
    assert(analytic_stats.linear_solves == fd_stats.linear_solves);
    assert(analytic_stats.jacobian_evaluations == fd_stats.jacobian_evaluations);
    assert(fd_context.jacobian_calls == 0U);
    assert(fd_stats.rhs_evaluations == analytic_stats.rhs_evaluations
        + 2U * fd_stats.jacobian_evaluations);
    check_accounting(&analytic_context, &analytic_stats);
    assert(fd_context.rhs_calls == fd_stats.rhs_evaluations);
}

static void test_failure_preserves_last_accepted(
    jacobian_behavior behavior, size_t fail_after_accepted
)
{
    linear_context context = {
        .dimension = 2U,
        .fail_after_accepted = fail_after_accepted,
        .behavior = behavior
    };
    rp_integrator_config config = make_config(2U);
    rp_integrator_stats stats;
    double state[2] = {1.0, 2.0};
    const rp_integrator_status expected = behavior == JACOBIAN_FAILURE
        ? RP_INTEGRATOR_STATUS_RHS_FAILURE
        : RP_INTEGRATOR_STATUS_NONFINITE_VALUE;

    context.matrix[0][0] = 1.0;
    context.matrix[1][1] = -2.0;
    config.step_observer = observe;
    config.step_observer_context = &context;
    assert(rp_integrator_integrate(
        &config, linear_rhs, &context, 0.0, 1.0, state, &stats
    ) == expected);
    assert(stats.accepted_steps == fail_after_accepted);
    assert(stats.rejected_steps == 0U);
    assert(stats.final_independent_variable == context.last_tau);
    assert(context.jacobian_calls == context.successful_jacobians + 1U);
    if (fail_after_accepted == 0U) {
        assert(state[0] == 1.0 && state[1] == 2.0);
        assert(stats.rhs_evaluations == 2U);
        assert(stats.linear_solves == 0U);
    } else {
        assert(state[0] == context.last_state[0]);
        assert(state[1] == context.last_state[1]);
        assert(context.last_tau == 0.05);
    }
    check_accounting(&context, &stats);
}

static void test_explicit_methods_ignore_callback(void)
{
    const rp_integrator_method methods[] = {
        RP_INTEGRATOR_METHOD_DOP853, RP_INTEGRATOR_METHOD_DP45
    };
    size_t index;

    for (index = 0U; index < sizeof(methods) / sizeof(methods[0]); ++index) {
        linear_context context = {
            .dimension = 1U, .behavior = JACOBIAN_FAILURE
        };
        rp_integrator_config config = make_config(1U);
        rp_integrator_stats stats;
        double state[1] = {1.0};

        context.matrix[0][0] = 1.0;
        config.method = methods[index];
        assert(rp_integrator_integrate(
            &config, linear_rhs, &context, 0.0, 1.0, state, &stats
        ) == RP_INTEGRATOR_STATUS_OK);
        assert(fabs(state[0] - exp(1.0)) < 2.0e-7);
        assert(context.jacobian_calls == 0U);
        check_accounting(&context, &stats);
    }
}

static void test_zero_interval_does_not_call_callback(void)
{
    linear_context context = {
        .dimension = 1U, .behavior = JACOBIAN_FAILURE
    };
    rp_integrator_config config = make_config(1U);
    rp_integrator_stats stats;
    double state[1] = {1.0};

    assert(rp_integrator_integrate(
        &config, linear_rhs, &context, 0.0, 0.0, state, &stats
    ) == RP_INTEGRATOR_STATUS_OK);
    assert(state[0] == 1.0);
    assert(context.jacobian_calls == 0U);
    check_accounting(&context, &stats);
}

int main(void)
{
    jacobian_behavior behavior;

    test_scalar_forward_and_backward();
    test_coupled_linear_solution_and_row_stride();
    test_full_dimension_storage();
    test_fifth_order_step_refinement();
    test_stiff_decay();
    test_fd_selection_and_rhs_accounting();
    for (behavior = JACOBIAN_FAILURE; behavior <= JACOBIAN_INCOMPLETE;
         behavior = (jacobian_behavior)((int)behavior + 1)) {
        test_failure_preserves_last_accepted(behavior, 0U);
        test_failure_preserves_last_accepted(behavior, 1U);
    }
    test_explicit_methods_ignore_callback();
    test_zero_interval_does_not_call_callback();
    return 0;
}
