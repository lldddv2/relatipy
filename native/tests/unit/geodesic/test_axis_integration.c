#include "geodesic/integrators/integrator.h"
#include "geodesic/integrators/kerr.h"

#include <assert.h>
#include <math.h>

int main(void)
{
    const rp_integrator_method methods[] = {
        RP_INTEGRATOR_METHOD_DOP853,
        RP_INTEGRATOR_METHOD_DP45,
        RP_INTEGRATOR_METHOD_RADAU
    };
    rp_kerr_integrator_context context = {.mass = 1.0, .spin = 0.0};
    size_t method;

    for (method = 0U; method < sizeof(methods) / sizeof(methods[0]); ++method) {
        const rp_integrator_config config = {
            methods[method], 8U, 1.0e-6, 1.0e-9,
            0.02, 0.02, 100U, NULL, NULL, NULL, NULL, NULL
        };
        double state[8] = {
            0.0, 8.0, 0.01, 0.0, sqrt(65.0 / 0.75), 0.0, -1.0, 0.0
        };
        rp_integrator_stats stats;

        assert(rp_integrator_integrate(
            &config, rp_kerr_integrator_rhs, &context,
            0.0, 0.02, state, &stats
        ) == RP_INTEGRATOR_STATUS_RHS_FAILURE);
        assert(stats.accepted_steps == 0U);
        assert(state[2] == 0.01);
    }
    return 0;
}
