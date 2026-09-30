/* Analytic and finite-difference Radau paths use the same corrected Kerr RHS.
 * DOP853 at tighter tolerances supplies an independent integrator reference;
 * this is a Jacobian regression test, not an independent physical oracle. */
#include "geodesic/integrators/integrator.h"
#include "geodesic/integrators/kerr.h"
#include "geodesic/jacobian.h"
#include "relatipy/kerr_geometry.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static rp_integrator_config config(int analytic)
{
    const rp_integrator_config result = {
        .method = RP_INTEGRATOR_METHOD_RADAU,
        .dimension = 8U,
        .relative_tolerance = 1e-9,
        .absolute_tolerance = 1e-12,
        .initial_step = 0.01,
        .maximum_step = 0.1,
        .maximum_steps = 100000U,
        .jacobian = analytic ? rp_kerr_integrator_jacobian : NULL
    };
    return result;
}

static void invariants(const rp_kerr_integrator_context *context,
    const double state[8], double values[3])
{
    double metric[4][4];
    size_t i;
    size_t j;
    assert(rp_kerr_metric(context->mass, context->spin, state, metric)
        == RP_KERR_STATUS_OK);
    values[0] = 0.0;
    for (i = 0U; i < 4U; ++i) {
        for (j = 0U; j < 4U; ++j) {
            values[0] += metric[i][j] * state[i + 4U] * state[j + 4U];
        }
    }
    values[1] = -(metric[0][0] * state[4] + metric[0][3] * state[7]);
    values[2] = metric[3][0] * state[4] + metric[3][3] * state[7];
}

static void test_regular_endpoints_and_invariants(void)
{
    const double spins[] = {0.0, 0.5, -0.9, 1.0};
    const double initial[4] = {0.0, 8.0, 1.1, 0.3};
    const double velocity[3] = {-0.01, 0.002, 0.02};
    size_t case_index;
    for (case_index = 0U; case_index < 4U; ++case_index) {
        rp_kerr_integrator_context context = {
            .mass = 1.0, .spin = spins[case_index]
        };
        rp_integrator_config analytic_config = config(1);
        rp_integrator_config fd_config = config(0);
        rp_integrator_config reference_config = config(0);
        rp_integrator_stats analytic_stats;
        rp_integrator_stats fd_stats;
        rp_integrator_stats reference_stats;
        double analytic[8];
        double fd[8];
        double reference[8];
        double start_invariants[3];
        double end_invariants[3];
        size_t i;

        memcpy(analytic, initial, sizeof(initial));
        assert(rp_kerr_four_velocity(1.0, context.spin, initial, velocity,
            analytic + 4) == RP_KERR_STATUS_OK);
        invariants(&context, analytic, start_invariants);
        memcpy(fd, analytic, sizeof(fd));
        memcpy(reference, analytic, sizeof(reference));
        reference_config.method = RP_INTEGRATOR_METHOD_DOP853;
        reference_config.relative_tolerance = 1e-12;
        reference_config.absolute_tolerance = 1e-14;
        assert(rp_integrator_integrate(&reference_config,
            rp_kerr_integrator_rhs, &context, 0.0, 5.0,
            reference, &reference_stats) == RP_INTEGRATOR_STATUS_OK);
        assert(rp_integrator_integrate(&fd_config, rp_kerr_integrator_rhs,
            &context, 0.0, 5.0, fd, &fd_stats) == RP_INTEGRATOR_STATUS_OK);
        assert(rp_integrator_integrate(&analytic_config, rp_kerr_integrator_rhs,
            &context, 0.0, 5.0, analytic, &analytic_stats)
            == RP_INTEGRATOR_STATUS_OK);
        for (i = 0U; i < 8U; ++i) {
            assert(isfinite(analytic[i]));
            assert(fabs(analytic[i] - fd[i]) < 2e-9);
            assert(fabs(analytic[i] - reference[i]) < 2e-9);
        }
        invariants(&context, analytic, end_invariants);
        for (i = 0U; i < 3U; ++i) {
            assert(fabs(end_invariants[i] - start_invariants[i]) < 2e-9);
        }
        invariants(&context, fd, end_invariants);
        for (i = 0U; i < 3U; ++i) {
            assert(fabs(end_invariants[i] - start_invariants[i]) < 2e-9);
        }
        assert(analytic_stats.jacobian_evaluations > 0U);
        assert(analytic_stats.rhs_evaluations < fd_stats.rhs_evaluations);
    }
}

static void test_polar_valid_points_and_fd_probe_failure(void)
{
    const double pi = acos(-1.0);
    const double distances[] = {1e-4, 1e-8, 1e-10, 128.0 * DBL_EPSILON};
    rp_kerr_integrator_context context = {.mass = 1.0, .spin = 0.0};
    size_t pole;
    size_t point;
    for (pole = 0U; pole < 2U; ++pole) {
        for (point = 0U; point < 4U; ++point) {
            const double theta = pole == 0U ? distances[point]
                : pi - distances[point];
            double analytic[8] = {0.0, 8.0, theta, 0.0,
                1.0 / sqrt(0.75), 0.0, 0.0, 0.0};
            double fd[8];
            double reference[8];
            rp_integrator_config analytic_config = config(1);
            rp_integrator_config fd_config = config(0);
            rp_integrator_config reference_config = config(0);
            rp_integrator_stats stats;
            rp_integrator_status fd_status;
            size_t i;
            memcpy(fd, analytic, sizeof(fd));
            memcpy(reference, analytic, sizeof(reference));
            analytic_config.initial_step = fd_config.initial_step = 1e-4;
            analytic_config.maximum_step = fd_config.maximum_step = 1e-4;
            reference_config.method = RP_INTEGRATOR_METHOD_DOP853;
            reference_config.initial_step = reference_config.maximum_step = 1e-4;
            reference_config.relative_tolerance = 1e-12;
            reference_config.absolute_tolerance = 1e-14;
            assert(rp_integrator_integrate(&analytic_config,
                rp_kerr_integrator_rhs, &context, 0.0, 1e-3, analytic, &stats)
                == RP_INTEGRATOR_STATUS_OK);
            assert(rp_integrator_integrate(&reference_config,
                rp_kerr_integrator_rhs, &context, 0.0, 1e-3, reference, &stats)
                == RP_INTEGRATOR_STATUS_OK);
            for (i = 0U; i < 8U; ++i) {
                assert(isfinite(analytic[i]));
                assert(fabs(analytic[i] - reference[i]) < 2e-11);
            }
            assert(analytic[2] == theta);
            fd_status = rp_integrator_integrate(&fd_config,
                rp_kerr_integrator_rhs, &context, 0.0, 1e-3, fd, &stats);
            if (pole == 1U && point > 0U) {
                /* The forward theta perturbation leaves the chart even
                 * though this stationary polar coordinate never does. */
                assert(fd_status == RP_INTEGRATOR_STATUS_RHS_FAILURE);
                assert(stats.accepted_steps == 0U);
                assert(fd[2] == theta);
            } else {
                assert(fd_status == RP_INTEGRATOR_STATUS_OK);
            }
        }
    }
}

static void test_true_stage_crossing_and_nonfinite_context(void)
{
    rp_kerr_integrator_context context = {.mass = 1.0, .spin = 0.0};
    rp_integrator_config options = config(1);
    rp_integrator_stats stats;
    double state[8] = {0.0, 8.0, 0.01, 0.0,
        sqrt(65.0 / 0.75), 0.0, -1.0, 0.0};
    double original[8];
    options.initial_step = options.maximum_step = 0.02;
    memcpy(original, state, sizeof(original));
    assert(rp_integrator_integrate(&options, rp_kerr_integrator_rhs,
        &context, 0.0, 0.02, state, &stats)
        == RP_INTEGRATOR_STATUS_RHS_FAILURE);
    assert(stats.accepted_steps == 0U);
    assert(memcmp(state, original, sizeof(state)) == 0);
    context.spin = NAN;
    assert(rp_integrator_integrate(&options, rp_kerr_integrator_rhs,
        &context, 0.0, 0.02, state, &stats)
        == RP_INTEGRATOR_STATUS_RHS_FAILURE);
    assert(stats.accepted_steps == 0U);
    assert(memcmp(state, original, sizeof(state)) == 0);
}

static int record_accepted_state(double tau, const double state[], void *context)
{
    double *last = context;
    last[0] = tau;
    memcpy(last + 1U, state, 8U * sizeof(double));
    return 0;
}

static void test_failure_retains_last_accepted_state(void)
{
    rp_kerr_integrator_context context = {.mass = 1.0, .spin = 0.0};
    rp_integrator_config options = config(1);
    rp_integrator_stats stats;
    double state[8] = {0.0, 8.0, 0.03, 0.0,
        sqrt(65.0 / 0.75), 0.0, -1.0, 0.0};
    double last[9] = {0.0};
    options.initial_step = 0.005;
    options.maximum_step = 0.01;
    options.step_observer = record_accepted_state;
    options.step_observer_context = last;
    assert(rp_integrator_integrate(&options, rp_kerr_integrator_rhs,
        &context, 0.0, 0.1, state, &stats)
        == RP_INTEGRATOR_STATUS_RHS_FAILURE);
    assert(stats.accepted_steps > 0U);
    assert(last[0] == stats.final_independent_variable);
    assert(last[0] > 0.0 && last[0] < 0.1);
    assert(memcmp(state, last + 1U, sizeof(state)) == 0);
    assert(state[2] > 0.0);
}

static void test_adapter_stride_and_errors(void)
{
    rp_kerr_integrator_context context = {.mass = 1.0, .spin = 0.5};
    const double state[8] = {0.0, 8.0, 1.1, 0.3,
        1.2, -0.01, 0.002, 0.02};
    double expected[8][8];
    double padded[8][16];
    double invalid[8];
    size_t row;
    size_t column;
    assert(rp_kerr_geodesic_jacobian(1.0, 0.5, state, state + 4U, expected)
        == RP_KERR_STATUS_OK);
    for (row = 0U; row < 8U; ++row) {
        for (column = 0U; column < 16U; ++column) padded[row][column] = 12345.0;
    }
    assert(rp_kerr_integrator_jacobian(0.0, state, &padded[0][0],
        8U, 16U, &context) == 0);
    for (row = 0U; row < 8U; ++row) {
        for (column = 0U; column < 8U; ++column) {
            assert(padded[row][column] == expected[row][column]);
        }
        for (column = 8U; column < 16U; ++column) {
            assert(padded[row][column] == 12345.0);
        }
    }
    assert(rp_kerr_integrator_jacobian(0.0, state, &padded[0][0],
        7U, 16U, &context) == RP_KERR_STATUS_INVALID_PARAMETER);
    assert(rp_kerr_integrator_jacobian(0.0, state, &padded[0][0],
        8U, 7U, &context) == RP_KERR_STATUS_INVALID_PARAMETER);
    assert(rp_kerr_integrator_jacobian(0.0, state, &padded[0][0],
        8U, 16U, NULL) == RP_KERR_STATUS_NULL_POINTER);
    memcpy(invalid, state, sizeof(invalid));
    invalid[2] = 0.0;
    assert(rp_kerr_integrator_jacobian(0.0, invalid, &padded[0][0],
        8U, 16U, &context) == RP_KERR_STATUS_COORDINATE_SINGULARITY);
    invalid[2] = NAN;
    assert(rp_kerr_integrator_jacobian(0.0, invalid, &padded[0][0],
        8U, 16U, &context) == RP_KERR_STATUS_NONFINITE_INPUT);
}

int main(void)
{
    test_regular_endpoints_and_invariants();
    test_polar_valid_points_and_fd_probe_failure();
    test_true_stage_crossing_and_nonfinite_context();
    test_failure_retains_last_accepted_state();
    test_adapter_stride_and_errors();
    puts("Kerr analytic Radau integration passed");
    return 0;
}
