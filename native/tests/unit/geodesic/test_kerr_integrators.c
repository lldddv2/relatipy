#include "geodesic/integrators/integrator.h"
#include "geodesic/integrators/kerr.h"
#include "relatipy/kerr_geodesic.h"
#include "relatipy/kerr_geometry.h"

#include <assert.h>
#include <math.h>
#include <stddef.h>

static rp_integrator_config make_config(rp_integrator_method method)
{
    const rp_integrator_config config = {
        method,
        2U * RP_KERR_DIM,
        1.0e-9,
        1.0e-12,
        1.0e-5,
        1.0e-5,
        10000U,
        NULL,
        NULL,
        NULL,
        NULL,
        NULL
    };
    return config;
}

static double metric_norm(
    const rp_kerr_integrator_context *context,
    const double state[2U * RP_KERR_DIM]
)
{
    double metric[RP_KERR_DIM][RP_KERR_DIM];
    double norm = 0.0;
    size_t mu;
    size_t nu;

    assert(rp_kerr_metric(
        context->mass, context->spin, state, metric
    ) == RP_KERR_STATUS_OK);
    for (mu = 0U; mu < RP_KERR_DIM; ++mu) {
        for (nu = 0U; nu < RP_KERR_DIM; ++nu) {
            norm += metric[mu][nu]
                * state[RP_KERR_DIM + mu] * state[RP_KERR_DIM + nu];
        }
    }
    return norm;
}

static void timelike_constants(
    const rp_kerr_integrator_context *context,
    const double state[2U * RP_KERR_DIM],
    double *energy,
    double *angular_momentum
)
{
    double metric[RP_KERR_DIM][RP_KERR_DIM];
    const double *velocity = state + RP_KERR_DIM;

    assert(rp_kerr_metric(
        context->mass, context->spin, state, metric
    ) == RP_KERR_STATUS_OK);
    *energy = -(metric[0][0] * velocity[0]
        + metric[0][3] * velocity[3]);
    *angular_momentum = metric[3][0] * velocity[0]
        + metric[3][3] * velocity[3];
}

static void assert_timelike_constants(
    const rp_kerr_integrator_context *context,
    const double state[2U * RP_KERR_DIM],
    double expected_energy,
    double expected_angular_momentum,
    double tolerance
)
{
    double energy;
    double angular_momentum;

    timelike_constants(context, state, &energy, &angular_momentum);
    assert(fabs(metric_norm(context, state) + 1.0) < tolerance);
    assert(fabs(energy - expected_energy) < tolerance);
    assert(fabs(angular_momentum - expected_angular_momentum) < tolerance);
}

static void assert_null_constants(
    const rp_kerr_integrator_context *context,
    const double state[2U * RP_KERR_DIM],
    double expected_energy,
    double expected_angular_momentum,
    double expected_carter_constant,
    double tolerance
)
{
    const double radius = state[1];
    const double theta = state[2];
    const double sin_theta = sin(theta);
    const double cos_theta = cos(theta);
    const double sigma = radius * radius
        + context->spin * context->spin * cos_theta * cos_theta;
    const double *tangent = state + RP_KERR_DIM;
    double metric[RP_KERR_DIM][RP_KERR_DIM];
    double energy;
    double angular_momentum;
    double carter_constant;

    assert(rp_kerr_metric(
        context->mass, context->spin, state, metric
    ) == RP_KERR_STATUS_OK);
    energy = -(metric[0][0] * tangent[0] + metric[0][3] * tangent[3]);
    angular_momentum = metric[3][0] * tangent[0]
        + metric[3][3] * tangent[3];
    carter_constant = sigma * sigma * tangent[2] * tangent[2]
        + cos_theta * cos_theta
            * (angular_momentum * angular_momentum
                / (sin_theta * sin_theta)
                - context->spin * context->spin * energy * energy);

    assert(fabs(metric_norm(context, state)) < tolerance);
    assert(fabs(energy - expected_energy) < tolerance);
    assert(fabs(angular_momentum - expected_angular_momentum) < tolerance);
    assert(fabs(carter_constant - expected_carter_constant) < tolerance);
}

static void test_null_geodesic_integration(void)
{
    const double coordinates[RP_KERR_DIM] = {0.0, 8.0, 1.1, 0.3};
    const double energy = 1.0;
    const double angular_momentum = 2.0;
    const double carter_constant = 3.0;
    rp_kerr_integrator_context context = {.mass = 1.0, .spin = 0.5};
    rp_integrator_config dop853_config = make_config(
        RP_INTEGRATOR_METHOD_DOP853
    );
    rp_integrator_config radau_config = make_config(
        RP_INTEGRATOR_METHOD_RADAU
    );
    rp_integrator_stats dop853_stats;
    rp_integrator_stats radau_stats;
    double dop853_state[2U * RP_KERR_DIM];
    double radau_state[2U * RP_KERR_DIM];
    size_t component;

    for (component = 0U; component < RP_KERR_DIM; ++component) {
        dop853_state[component] = coordinates[component];
    }
    assert(rp_kerr_null_tangent(
        context.mass,
        context.spin,
        coordinates,
        energy,
        angular_momentum,
        carter_constant,
        -1,
        1,
        dop853_state + RP_KERR_DIM
    ) == RP_KERR_STATUS_OK);
    for (component = 0U; component < 2U * RP_KERR_DIM; ++component) {
        radau_state[component] = dop853_state[component];
    }

    assert(rp_integrator_integrate(
        &dop853_config,
        rp_kerr_null_integrator_rhs,
        &context,
        0.0,
        1.0e-3,
        dop853_state,
        &dop853_stats
    ) == RP_INTEGRATOR_STATUS_OK);
    assert(rp_integrator_integrate(
        &radau_config,
        rp_kerr_null_integrator_rhs,
        &context,
        0.0,
        1.0e-3,
        radau_state,
        &radau_stats
    ) == RP_INTEGRATOR_STATUS_OK);

    assert_null_constants(
        &context,
        dop853_state,
        energy,
        angular_momentum,
        carter_constant,
        2.0e-10
    );
    assert_null_constants(
        &context,
        radau_state,
        energy,
        angular_momentum,
        carter_constant,
        2.0e-10
    );
    for (component = 0U; component < 2U * RP_KERR_DIM; ++component) {
        assert(fabs(dop853_state[component] - radau_state[component])
            < 2.0e-10);
    }
    assert(dop853_stats.accepted_steps > 0U);
    assert(radau_stats.accepted_steps > 0U);
}

int main(void)
{
    const double coordinates[RP_KERR_DIM] = {0.0, 8.0, 1.1, 0.3};
    const double coordinate_velocity[3] = {-0.01, 0.002, 0.02};
    rp_kerr_integrator_context context = {.mass = 1.0, .spin = 0.5};
    rp_integrator_config dop853_config = make_config(
        RP_INTEGRATOR_METHOD_DOP853
    );
    rp_integrator_config radau_config = make_config(
        RP_INTEGRATOR_METHOD_RADAU
    );
    rp_integrator_stats dop853_stats;
    rp_integrator_stats radau_stats;
    double dop853_state[2U * RP_KERR_DIM];
    double radau_state[2U * RP_KERR_DIM];
    double initial_energy;
    double initial_angular_momentum;
    size_t component;

    for (component = 0U; component < RP_KERR_DIM; ++component) {
        dop853_state[component] = coordinates[component];
    }
    assert(rp_kerr_four_velocity(
        context.mass,
        context.spin,
        coordinates,
        coordinate_velocity,
        dop853_state + RP_KERR_DIM
    ) == RP_KERR_STATUS_OK);
    for (component = 0U; component < 2U * RP_KERR_DIM; ++component) {
        radau_state[component] = dop853_state[component];
    }
    timelike_constants(
        &context, dop853_state, &initial_energy, &initial_angular_momentum
    );

    assert(rp_integrator_integrate(
        &dop853_config,
        rp_kerr_integrator_rhs,
        &context,
        0.0,
        1.0e-2,
        dop853_state,
        &dop853_stats
    ) == RP_INTEGRATOR_STATUS_OK);
    assert(rp_integrator_integrate(
        &radau_config,
        rp_kerr_integrator_rhs,
        &context,
        0.0,
        1.0e-2,
        radau_state,
        &radau_stats
    ) == RP_INTEGRATOR_STATUS_OK);

    for (component = 0U; component < 2U * RP_KERR_DIM; ++component) {
        assert(isfinite(dop853_state[component]));
        assert(isfinite(radau_state[component]));
        assert(fabs(dop853_state[component] - radau_state[component])
            < 2.0e-10);
    }
    assert(dop853_stats.accepted_steps > 0U);
    assert(radau_stats.accepted_steps > 0U);
    assert_timelike_constants(
        &context, dop853_state, initial_energy, initial_angular_momentum,
        2.0e-10
    );
    assert_timelike_constants(
        &context, radau_state, initial_energy, initial_angular_momentum,
        2.0e-10
    );
    test_null_geodesic_integration();
    return 0;
}
