/*
 * Newtonian limit of the forced Kerr equation.
 *
 * At spin 0, large radius and weak external field, the displacement caused
 * by the force, Delta = X(forced) - X(free) at fixed coordinate time, must
 * match an independent Newtonian integration of a point mass plus the same
 * acceleration.  The residual is a 1PN correction: relative to |Delta| it
 * stays below 6 M/r (measured ~3.2 M/r) and shrinks like 1/r between radii
 * 1000 and 8000.
 */

#include "force/kerr_forced.h"
#include "force_test_models.h"
#include "geodesic/initial/convert.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static const double PI = 3.14159265358979323846;

typedef struct {
    const rp_external_force *force;
} newtonian_context;

/* Point mass M = 1 plus the external acceleration, state (x, v). */
static int newtonian_rhs(
    double time,
    const double state[],
    double derivative[],
    void *opaque_context
)
{
    const newtonian_context *context = opaque_context;
    const double r2 = state[0] * state[0] + state[1] * state[1]
        + state[2] * state[2];
    const double inverse_r3 = 1.0 / (r2 * sqrt(r2));
    double external[3] = {0.0, 0.0, 0.0};
    size_t k;

    (void)time;
    if (context->force != NULL
        && context->force->evaluate(
            state, state + 3, context->force->model, external
        ) != 0) {
        return 1;
    }
    for (k = 0U; k < 3U; ++k) {
        derivative[k] = state[3U + k];
        derivative[3U + k] = -state[k] * inverse_r3 + external[k];
    }
    return 0;
}

static void integrate_newtonian(
    const rp_external_force *force,
    const double initial[6],
    double duration,
    double final[3]
)
{
    const rp_integrator_config config = {
        RP_INTEGRATOR_METHOD_DOP853, 6U, 1e-13, 1e-13, 0.0, 0.0, 1000000U,
        NULL, NULL, NULL, NULL, NULL
    };
    newtonian_context context;
    double state[6];
    size_t k;

    context.force = force;
    for (k = 0U; k < 6U; ++k) {
        state[k] = initial[k];
    }
    assert(rp_integrator_integrate(
        &config, newtonian_rhs, &context, 0.0, duration, state, NULL
    ) == RP_INTEGRATOR_STATUS_OK);
    for (k = 0U; k < 3U; ++k) {
        final[k] = state[k];
    }
}

static void integrate_kerr(
    const rp_external_force *force,
    const double initial[6],
    double duration,
    double final[3]
)
{
    const rp_integrator_config config = {
        RP_INTEGRATOR_METHOD_DOP853, 8U, 1e-13, 1e-13, 0.0, 0.0, 1000000U,
        NULL, NULL, NULL, NULL, NULL
    };
    const double cartesian[7] = {
        0.0, initial[0], initial[1], initial[2],
        initial[3], initial[4], initial[5]
    };
    rp_kerr_forced_context context;
    double state[8];
    double out[7];
    size_t k;

    context.spin = 0.0;
    context.force = force;
    assert(rp_initial_cartesian_to_canonical(0.0, cartesian, state)
        == RP_KERR_STATUS_OK);
    assert(rp_integrator_integrate(
        &config, rp_kerr_forced_coordinate_time_rhs, &context,
        0.0, duration, state, NULL
    ) == RP_INTEGRATOR_STATUS_OK);
    state[0] = duration;
    assert(rp_initial_canonical_to_cartesian(0.0, state, out)
        == RP_KERR_STATUS_OK);
    for (k = 0U; k < 3U; ++k) {
        final[k] = out[1U + k];
    }
}

static double norm3(const double v[3])
{
    return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

/* Relative 1PN residual of the force effect over a third of an orbit. */
static double force_effect_residual(double radius)
{
    /* A distant clump pulling ~5e-3 of the black hole at every radius. */
    const test_plummer model = {
        0.04,
        0.5 * radius,
        {2.0 * radius, -1.0 * radius, 1.5 * radius}
    };
    const rp_external_force force = {test_plummer_acceleration, &model};
    const double speed = 1.1 / sqrt(radius);
    const double initial[6] = {
        radius, 0.0, 0.0,
        0.0, speed * cos(0.5), speed * sin(0.5)
    };
    const double duration = 2.0 * PI * pow(radius, 1.5) / 3.0;
    double newton_free[3];
    double newton_forced[3];
    double kerr_free[3];
    double kerr_forced[3];
    double newton_effect[3];
    double kerr_effect[3];
    double residual[3];
    double free_offset[3];
    size_t k;

    integrate_newtonian(NULL, initial, duration, newton_free);
    integrate_newtonian(&force, initial, duration, newton_forced);
    integrate_kerr(NULL, initial, duration, kerr_free);
    integrate_kerr(&force, initial, duration, kerr_forced);
    for (k = 0U; k < 3U; ++k) {
        newton_effect[k] = newton_forced[k] - newton_free[k];
        kerr_effect[k] = kerr_forced[k] - kerr_free[k];
        residual[k] = kerr_effect[k] - newton_effect[k];
        free_offset[k] = kerr_free[k] - newton_free[k];
    }
    printf("r=%7.0f |effect|/r=%.3e  free 1PN offset/r=%.3e"
        "  effect residual=%.3e (M/r=%.1e)\n",
        radius, norm3(newton_effect) / radius, norm3(free_offset) / radius,
        norm3(residual) / norm3(newton_effect), 1.0 / radius);
    fflush(stdout);
    /* The force must matter well above the integration noise. */
    assert(norm3(newton_effect) > 1e-4 * radius);
    return norm3(residual) / norm3(newton_effect);
}

int main(void)
{
    const double near = force_effect_residual(1000.0);
    const double far = force_effect_residual(8000.0);

    assert(near < 6.0 / 1000.0);
    assert(far < 6.0 / 8000.0);
    /* 1PN scaling: an eightfold radius shrinks the residual ~eightfold. */
    assert(near / far > 6.0 && near / far < 10.0);
    puts("test_newtonian_limit: ok");
    return 0;
}
