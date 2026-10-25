/*
 * An inactive force (NULL, or a NULL callback) reproduces the homogeneous
 * Kerr path bit for bit: the RHS and its statuses, every accepted step of
 * dop853, dp45 and radau (finite-difference Jacobian) integrations, and both
 * observable evaluators.  The forced layer rejects projection_radau.
 */

#include "force/kerr_forced.h"
#include "force/kerr_mcmc.h"
#include "force_test_models.h"
#include "geodesic/initial/convert.h"
#include "geodesic/integrators/kerr.h"
#include "geodesic/physic/kerr_observables.h"
#include "relatipy/kerr_geodesic.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static const double IDENTITY[3][3] = {
    {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}
};

static int same_bits(const double *left, const double *right, size_t count)
{
    return memcmp(left, right, count * sizeof(double)) == 0;
}

static int same_stats(const rp_integrator_stats *a, const rp_integrator_stats *b)
{
    return a->accepted_steps == b->accepted_steps
        && a->rejected_steps == b->rejected_steps
        && a->rhs_evaluations == b->rhs_evaluations
        && a->jacobian_evaluations == b->jacobian_evaluations
        && a->linear_solves == b->linear_solves
        && same_bits(&a->final_independent_variable,
            &b->final_independent_variable, 1U);
}

static void check_rhs_point(
    double spin,
    const rp_external_force *force,
    const double x[4],
    const double u[4]
)
{
    double dx_h[4] = {9.0, 9.0, 9.0, 9.0};
    double du_h[4] = {9.0, 9.0, 9.0, 9.0};
    double dx_f[4] = {7.0, 7.0, 7.0, 7.0};
    double du_f[4] = {7.0, 7.0, 7.0, 7.0};
    double state[8];
    double packed_h[8];
    double packed_f[8];
    rp_kerr_integrator_context homogeneous;
    rp_kerr_forced_context forced;
    int status_h;
    int status_f;

    status_h = (int)rp_kerr_null_geodesic_rhs(1.0, spin, x, u, dx_h, du_h);
    status_f = rp_kerr_forced_rhs(spin, force, x, u, dx_f, du_f);
    assert(status_h == status_f);
    assert(same_bits(dx_h, dx_f, 4U));
    assert(same_bits(du_h, du_f, 4U));

    if (x == NULL || u == NULL) {
        return;
    }
    memcpy(state, x, 4U * sizeof(double));
    memcpy(state + 4, u, 4U * sizeof(double));
    memset(&homogeneous, 0, sizeof(homogeneous));
    homogeneous.mass = 1.0;
    homogeneous.spin = spin;
    forced.spin = spin;
    forced.force = force;
    status_h = rp_kerr_integrator_rhs(0.0, state, packed_h, &homogeneous);
    status_f = rp_kerr_forced_integrator_rhs(0.0, state, packed_f, &forced);
    assert(status_h == status_f);
    assert(same_bits(packed_h, packed_f, 8U));
}

static void check_rhs(const rp_external_force *force)
{
    const double spins[4] = {0.0, 0.3, 0.9, 1.0};
    const double points[4][4] = {
        {0.0, 8.0, 1.1, 0.3},
        {5.0, 3.0, 0.2, -2.0},
        {0.0, 400.0, 1.5707963267948966, 1.0},
        {-2.0, 2.5, 2.9, 6.0}
    };
    const double velocities[3][4] = {
        {1.1, -0.01, 0.002, 0.02},
        {1.6, 0.3, -0.05, 0.1},
        {0.0, 0.0, 0.0, 0.0}
    };
    const double singular[4] = {0.0, 2.0, 1.1, 0.3};
    const double axis[4] = {0.0, 8.0, 0.0, 0.3};
    const double nonfinite[4] = {0.0, INFINITY, 1.1, 0.3};
    double dx[4];
    double du[4];
    double dx_h[4];
    double du_h[4];
    size_t s;
    size_t p;
    size_t v;

    for (s = 0U; s < 4U; ++s) {
        for (p = 0U; p < 4U; ++p) {
            for (v = 0U; v < 3U; ++v) {
                check_rhs_point(spins[s], force, points[p], velocities[v]);
            }
        }
        check_rhs_point(spins[s], force, axis, velocities[0]);
        check_rhs_point(spins[s], force, nonfinite, velocities[0]);
        check_rhs_point(spins[s], force, NULL, velocities[0]);
        check_rhs_point(spins[s], force, points[0], NULL);
    }
    check_rhs_point(1.0, force, singular, velocities[0]);
    check_rhs_point(1.5, force, points[0], velocities[0]);
    check_rhs_point(NAN, force, points[0], velocities[0]);
    /* Output aliasing and NULL outputs follow the homogeneous wrapper. */
    assert(rp_kerr_forced_rhs(0.5, force, points[0], velocities[0], dx, dx)
        == (int)RP_KERR_STATUS_INVALID_PARAMETER);
    assert(rp_kerr_forced_rhs(0.5, force, points[0], velocities[0], NULL, du)
        == (int)RP_KERR_STATUS_NULL_POINTER);
    /* Inputs may overlap outputs. */
    memcpy(dx, points[0], sizeof(dx));
    memcpy(du, velocities[0], sizeof(du));
    assert(rp_kerr_null_geodesic_rhs(
        1.0, 0.5, points[0], velocities[0], dx_h, du_h
    ) == RP_KERR_STATUS_OK);
    assert(rp_kerr_forced_rhs(0.5, force, dx, du, dx, du) == 0);
    assert(same_bits(dx, dx_h, 4U) && same_bits(du, du_h, 4U));
}

#define MAX_ROWS 20000U

/* Accepted steps of one integration; stops at the outer horizon. */
typedef struct {
    double horizon;
    size_t count;
    double taus[MAX_ROWS];
    double states[8U * MAX_ROWS];
} recorder;

static recorder homogeneous_steps;
static recorder forced_steps;

static int record_step(double tau, const double state[], void *opaque)
{
    recorder *steps = opaque;

    if (steps->count >= MAX_ROWS) {
        return -1;
    }
    steps->taus[steps->count] = tau;
    memcpy(steps->states + 8U * steps->count, state, 8U * sizeof(double));
    ++steps->count;
    return state[1] <= steps->horizon ? 1 : 0;
}

/* Integrate one state with the homogeneous and the forced adapters. */
static void compare_integrations(
    const rp_external_force *force,
    rp_integrator_config config,
    double spin,
    const double initial[8],
    double tau_final
)
{
    rp_kerr_integrator_context homogeneous;
    rp_kerr_forced_context forced;
    rp_integrator_stats stats_h;
    rp_integrator_stats stats_f;
    rp_integrator_status status_h;
    rp_integrator_status status_f;
    double state_h[8];
    double state_f[8];

    memset(&homogeneous, 0, sizeof(homogeneous));
    homogeneous.mass = 1.0;
    homogeneous.spin = spin;
    forced.spin = spin;
    forced.force = force;
    memcpy(state_h, initial, sizeof(state_h));
    memcpy(state_f, initial, sizeof(state_f));
    homogeneous_steps.count = forced_steps.count = 0U;
    homogeneous_steps.horizon = forced_steps.horizon
        = 1.0 + sqrt((1.0 - spin) * (1.0 + spin));

    config.step_observer = record_step;
    config.step_observer_context = &homogeneous_steps;
    status_h = rp_integrator_integrate(
        &config, rp_kerr_integrator_rhs, &homogeneous,
        0.0, tau_final, state_h, &stats_h
    );
    config.step_observer_context = &forced_steps;
    status_f = rp_integrator_integrate(
        &config, rp_kerr_forced_integrator_rhs, &forced,
        0.0, tau_final, state_f, &stats_f
    );
    assert(status_f == status_h);
    assert(homogeneous_steps.count >= 2U);
    assert(forced_steps.count == homogeneous_steps.count);
    assert(same_bits(forced_steps.taus, homogeneous_steps.taus,
        forced_steps.count));
    assert(same_bits(forced_steps.states, homogeneous_steps.states,
        8U * forced_steps.count));
    assert(same_bits(state_f, state_h, 8U));
    assert(same_stats(&stats_f, &stats_h));
}

static void check_trajectory(
    const rp_external_force *force,
    rp_integrator_method method,
    double spin,
    int vector_tolerance
)
{
    const double elements[7] = {0.0, 30.0, 0.5, 0.7, 0.2, 1.3, 0.0};
    const double vector_atol[8] = {
        1e-9, 1e-9, 1e-11, 1e-11, 1e-11, 1e-12, 1e-12, 1e-12
    };
    rp_integrator_config config = {
        RP_INTEGRATOR_METHOD_DOP853, 8U, 1e-10, 1e-12, 0.0, 0.0, 1000000U,
        NULL, NULL, NULL, NULL, NULL
    };
    double initial[8];

    config.method = method;
    config.absolute_tolerances = vector_tolerance ? vector_atol : NULL;
    assert(rp_initial_from_elements(spin, elements, initial)
        == RP_KERR_STATUS_OK);
    compare_integrations(force, config, spin, initial, 2500.0);
}

static void check_plunge_and_rejections(const rp_external_force *force)
{
    /* A plunge ends at the horizon with the homogeneous status and steps. */
    const double bl[7] = {0.0, 6.0, 1.2, 0.0, -0.2, 0.0, 0.0};
    const rp_integrator_config config = {
        RP_INTEGRATOR_METHOD_DOP853, 8U, 1e-10, 1e-12, 0.0, 0.0, 1000000U,
        NULL, NULL, NULL, NULL, NULL
    };
    double initial[8];

    assert(rp_initial_from_bl(0.6, bl, initial) == RP_KERR_STATUS_OK);
    compare_integrations(force, config, 0.6, initial, 500.0);
    assert(forced_steps.states[8U * (forced_steps.count - 1U) + 1U]
        < 1.01 * forced_steps.horizon);

    assert(rp_kerr_forced_check_method(RP_INTEGRATOR_METHOD_PROJECTION_RADAU)
        == RP_INTEGRATOR_STATUS_INVALID_ARGUMENT);
    assert(rp_kerr_forced_check_method((rp_integrator_method)17)
        == RP_INTEGRATOR_STATUS_INVALID_ARGUMENT);
    assert(rp_kerr_forced_check_method(RP_INTEGRATOR_METHOD_RADAU)
        == RP_INTEGRATOR_STATUS_OK);
}

static void check_mcmc(const rp_external_force *force, rp_integrator_method method)
{
    const double c = cos(0.4);
    const double s = sin(0.4);
    const double rotation[3][3] = {
        {c, 0.0, s}, {0.0, 1.0, 0.0}, {-s, 0.0, c}
    };
    const double times[6] = {-3000.0, -500.0, -10.0, 0.0, 800.0, 4000.0};
    const rp_kerr_mcmc_observation observation = {
        0.01, 0.001, -0.002, 1e-6, -2e-6, 0.003, 20.0
    };
    const double spins[3] = {0.0, 0.4, 0.95};
    const double elements[7] = {23.0, 60.0, 0.3, 0.8, 0.3, 0.4, 0.3};
    const rp_kerr_mcmc_time_kind kinds[3] = {
        RP_KERR_MCMC_TIME_COORDINATE,
        RP_KERR_MCMC_TIME_PROPER,
        RP_KERR_MCMC_TIME_ARRIVAL
    };
    rp_integrator_config solver = {
        RP_INTEGRATOR_METHOD_DOP853, 8U, 1e-10, 1e-12, 0.0, 0.0, 100000U,
        NULL, NULL, NULL, NULL, NULL
    };
    double out_h[18];
    double out_f[18];
    double initial[8];
    double projected[4];
    double sample[4];
    rp_integrator_stats stats_h;
    rp_integrator_stats stats_f;
    size_t k;

    solver.method = method;
    for (k = 0U; k < 3U; ++k) {
        assert(rp_kerr_mcmc_evaluate(
            &solver, spins[k], rotation, 60.0, 0.6, 0.8, 0.3, 0.4,
            times, 6U, &observation, out_h, &stats_h
        ) == RP_KERR_MCMC_OK);
        assert(rp_kerr_forced_mcmc_evaluate(
            &solver, spins[k], force, rotation, 60.0, 0.6, 0.8, 0.3, 0.4,
            times, 6U, &observation, out_f, &stats_f
        ) == RP_KERR_MCMC_OK);
        assert(stats_h.accepted_steps > 0U);
        assert(same_bits(out_h, out_f, 18U));
        assert(same_stats(&stats_h, &stats_f));
    }

    assert(rp_initial_from_elements(0.7, elements, initial)
        == RP_KERR_STATUS_OK);
    assert(rp_kerr_observables_project(
        initial, 0.7, rotation, 0.01, projected
    ) == RP_KERR_STATUS_OK);
    for (k = 0U; k < 3U; ++k) {
        const double center = kinds[k] == RP_KERR_MCMC_TIME_COORDINATE
            ? initial[0]
            : kinds[k] == RP_KERR_MCMC_TIME_PROPER ? 7.0 : projected[0];
        sample[0] = center - 900.0;
        sample[1] = center;
        sample[2] = center + 250.0;
        sample[3] = center + 2000.0;
        assert(rp_kerr_mcmc_evaluate_state(
            &solver, 0.7, rotation, initial, 7.0, kinds[k],
            sample, 4U, 0.01, out_h, &stats_h
        ) == RP_KERR_MCMC_OK);
        assert(rp_kerr_forced_mcmc_evaluate_state(
            &solver, 0.7, force, rotation, initial, 7.0, kinds[k],
            sample, 4U, 0.01, out_f, &stats_f
        ) == RP_KERR_MCMC_OK);
        assert(same_bits(out_h, out_f, 12U));
        assert(same_stats(&stats_h, &stats_f));
    }
}

static void check_mcmc_rejections(void)
{
    const double times[2] = {0.0, 100.0};
    const rp_kerr_mcmc_observation observation = {
        0.01, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    };
    const double elements[7] = {0.0, 60.0, 0.3, 0.8, 0.3, 0.4, 0.0};
    rp_integrator_config solver = {
        RP_INTEGRATOR_METHOD_PROJECTION_RADAU, 8U, 1e-10, 1e-12, 0.0, 0.0,
        100000U, NULL, NULL, NULL, NULL, NULL
    };
    double out[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double initial[8];
    rp_integrator_stats stats;
    size_t k;

    assert(rp_kerr_forced_mcmc_evaluate(
        &solver, 0.5, NULL, IDENTITY, 60.0, 0.3, 0.8, 0.3, 0.4,
        times, 2U, &observation, out, &stats
    ) == RP_KERR_MCMC_INVALID_INPUT);
    for (k = 0U; k < 6U; ++k) {
        assert(out[k] == 0.0);
    }
    assert(rp_initial_from_elements(0.5, elements, initial)
        == RP_KERR_STATUS_OK);
    assert(rp_kerr_forced_mcmc_evaluate_state(
        &solver, 0.5, NULL, IDENTITY, initial, 0.0,
        RP_KERR_MCMC_TIME_PROPER, times, 2U, 0.01, out, &stats
    ) == RP_KERR_MCMC_INVALID_INPUT);
    solver.method = RP_INTEGRATOR_METHOD_DOP853;
    assert(rp_kerr_forced_mcmc_evaluate(
        NULL, 0.5, NULL, IDENTITY, 60.0, 0.3, 0.8, 0.3, 0.4,
        times, 2U, &observation, out, &stats
    ) == RP_KERR_MCMC_INVALID_INPUT);
}

int main(void)
{
    const test_plummer model = {0.1, 5.0, {0.0, 0.0, 0.0}};
    const rp_external_force empty = {NULL, &model};
    const rp_external_force *inactive[2] = {NULL, &empty};
    const rp_integrator_method methods[3] = {
        RP_INTEGRATOR_METHOD_DOP853,
        RP_INTEGRATOR_METHOD_DP45,
        RP_INTEGRATOR_METHOD_RADAU
    };
    size_t f;
    size_t m;

    for (f = 0U; f < 2U; ++f) {
        check_rhs(inactive[f]);
        for (m = 0U; m < 3U; ++m) {
            check_trajectory(inactive[f], methods[m], 0.0, 0);
            check_trajectory(inactive[f], methods[m], 0.9, 0);
            check_trajectory(inactive[f], methods[m], 0.5, 1);
            check_mcmc(inactive[f], methods[m]);
        }
        check_plunge_and_rejections(inactive[f]);
    }
    check_mcmc_rejections();
    puts("test_null_force: ok");
    return 0;
}
