#include "geodesic/integrators/integrator.h"
#include "geodesic/integrators/kerr.h"
#include "relatipy/kerr_geodesic.h"
#include "relatipy/kerr_geometry.h"

#include <assert.h>
#include <math.h>
#include <stddef.h>

int main(void)
{
    const double coordinates[RP_KERR_DIM] = {0.0, 8.0, 1.1, 0.3};
    const double coordinate_velocity[3] = {-0.01, 0.002, 0.02};
    rp_kerr_integrator_context context = {.mass = 1.0, .spin = 0.5};
    rp_integrator_config config = {
        RP_INTEGRATOR_METHOD_DP45,
        2U * RP_KERR_DIM,
        1.0e-9,
        1.0e-12,
        1.0e-4,
        1.0e-4,
        10000U,
        NULL,
        NULL,
        NULL,
        NULL,
        NULL
    };
    rp_integrator_stats stats;
    double state[2U * RP_KERR_DIM];
    double metric[RP_KERR_DIM][RP_KERR_DIM];
    double initial_energy;
    double initial_angular_momentum;
    double energy;
    double angular_momentum;
    double norm = 0.0;
    size_t component;
    size_t mu;
    size_t nu;

    for (component = 0U; component < RP_KERR_DIM; ++component) {
        state[component] = coordinates[component];
    }
    assert(rp_kerr_four_velocity(
        context.mass, context.spin, coordinates, coordinate_velocity,
        state + RP_KERR_DIM
    ) == RP_KERR_STATUS_OK);
    assert(rp_kerr_metric(
        context.mass, context.spin, state, metric
    ) == RP_KERR_STATUS_OK);
    initial_energy = -(metric[0][0] * state[4] + metric[0][3] * state[7]);
    initial_angular_momentum = metric[3][0] * state[4]
        + metric[3][3] * state[7];

    assert(rp_integrator_integrate(
        &config, rp_kerr_integrator_rhs, &context,
        0.0, 1.0e-2, state, &stats
    ) == RP_INTEGRATOR_STATUS_OK);
    assert(stats.accepted_steps > 0U);
    assert(stats.final_independent_variable == 1.0e-2);
    assert(stats.jacobian_evaluations == 0U);
    for (component = 0U; component < 2U * RP_KERR_DIM; ++component) {
        assert(isfinite(state[component]));
    }

    assert(rp_kerr_metric(
        context.mass, context.spin, state, metric
    ) == RP_KERR_STATUS_OK);
    for (mu = 0U; mu < RP_KERR_DIM; ++mu) {
        for (nu = 0U; nu < RP_KERR_DIM; ++nu) {
            norm += metric[mu][nu] * state[RP_KERR_DIM + mu]
                * state[RP_KERR_DIM + nu];
        }
    }
    energy = -(metric[0][0] * state[4] + metric[0][3] * state[7]);
    angular_momentum = metric[3][0] * state[4]
        + metric[3][3] * state[7];
    assert(fabs(norm + 1.0) < 2.0e-10);
    assert(fabs(energy - initial_energy) < 2.0e-10);
    assert(fabs(angular_momentum - initial_angular_momentum) < 2.0e-10);
    return 0;
}
