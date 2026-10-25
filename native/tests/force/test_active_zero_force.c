/*
 * An active force whose acceleration is zero evaluates the same generated
 * Kerr kernel as the homogeneous path and adds a zero four-acceleration, so
 * its RHS equals the homogeneous one exactly (up to the sign of a zero),
 * and dop853 observables integrated through it equal those of the
 * homogeneous evaluator.
 */

#include "force/kerr_forced.h"
#include "force/kerr_mcmc.h"
#include "relatipy/kerr_geodesic.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static const double IDENTITY[3][3] = {
    {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}
};

static int zero_acceleration(
    const double position[3],
    const double velocity[3],
    const void *model,
    double acceleration[3]
)
{
    (void)position;
    (void)velocity;
    (void)model;
    acceleration[0] = 0.0;
    acceleration[1] = 0.0;
    acceleration[2] = 0.0;
    return 0;
}

static size_t check_rhs(const rp_external_force *force)
{
    const double spins[5] = {-0.99, -0.5, 0.0, 0.7, 0.99};
    const double radii[5] = {2.5, 7.3, 60.0, 3.0e3, 1.0e6};
    const double thetas[4] = {1e-4, 0.4, 1.2, 1.5707963267948966};
    size_t checked = 0U;
    size_t s;
    size_t i;
    size_t j;
    size_t k;

    for (s = 0U; s < 5U; ++s) {
        for (i = 0U; i < 5U; ++i) {
            for (j = 0U; j < 4U; ++j) {
                const double r = radii[i];
                const double x[4] = {0.0, r, thetas[j], 0.3};
                const double u[4] = {
                    1.2, 0.3 / sqrt(r), 0.02 / r, 0.4 / pow(r, 1.5)
                };
                double dx_h[4];
                double du_h[4];
                double dx_f[4];
                double du_f[4];

                if (rp_kerr_null_geodesic_rhs(
                        1.0, spins[s], x, u, dx_h, du_h
                    ) != RP_KERR_STATUS_OK) {
                    continue;
                }
                assert(rp_kerr_forced_rhs(
                    spins[s], force, x, u, dx_f, du_f
                ) == 0);
                for (k = 0U; k < 4U; ++k) {
                    assert(dx_h[k] == dx_f[k]);
                    assert(du_h[k] == du_f[k]);
                }
                ++checked;
            }
        }
    }
    assert(checked > 0U);
    return checked;
}

static void check_mcmc(const rp_external_force *force)
{
    const double times[6] = {-3.0e5, -4.0e4, -10.0, 0.0, 800.0, 3.0e5};
    const rp_kerr_mcmc_observation observation = {
        1e-4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    };
    const double spins[3] = {0.0, 0.5, 0.95};
    const rp_integrator_config solver = {
        RP_INTEGRATOR_METHOD_DOP853, 8U, 1e-11, 1e-13, 0.0, 0.0, 1000000U,
        NULL, NULL, NULL, NULL, NULL
    };
    double out_h[18];
    double out_f[18];
    rp_integrator_stats stats;
    size_t s;
    size_t k;

    for (s = 0U; s < 3U; ++s) {
        assert(rp_kerr_mcmc_evaluate(
            &solver, spins[s], IDENTITY, 24385.0, 0.884, 2.35, 3.98, 1.16,
            times, 6U, &observation, out_h, &stats
        ) == RP_KERR_MCMC_OK);
        assert(rp_kerr_forced_mcmc_evaluate(
            &solver, spins[s], force, IDENTITY, 24385.0, 0.884, 2.35, 3.98,
            1.16, times, 6U, &observation, out_f, &stats
        ) == RP_KERR_MCMC_OK);
        for (k = 0U; k < 18U; ++k) {
            assert(out_f[k] == out_h[k]);
        }
    }
}

int main(void)
{
    const rp_external_force zero = {zero_acceleration, NULL};
    const size_t points = check_rhs(&zero);

    check_mcmc(&zero);
    printf("%zu RHS points and 3 observable runs equal\n", points);
    printf("test_active_zero_force: ok\n");
    return 0;
}
