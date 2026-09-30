#include "geodesic/initial/bound.h"

#include <assert.h>
#include <math.h>
#include <stddef.h>

static int close_to(double a, double b, double tolerance)
{
    return fabs(a - b) <= tolerance * (1.0 + fabs(b));
}

static void check_norm(double a, const double y[8])
{
    double metric[4][4];
    double norm = 0.0;
    size_t i, j;
    assert(rp_kerr_metric(1.0, a, y, metric) == RP_KERR_STATUS_OK);
    for (i = 0U; i < 4U; ++i) {
        for (j = 0U; j < 4U; ++j) {
            norm += metric[i][j] * y[4U + i] * y[4U + j];
        }
    }
    assert(close_to(norm, -1.0, 2e-11));
    assert(y[4] > 0.0);
}

static void test_references(void)
{
    /* KerrGeoPy 0.9.3, StableOrbit.trajectory/four_velocity at lambda=0. */
    const double inputs[][7] = {
        {4.0, 9.0, 0.2, -0.6, 1.2, 0.8, 2.0},
        {4.0, 8.0, 0.4, 0.8, 2.5, 3.8, -1.0},
        {4.0, 9.0, 0.0, 0.6, 1.1, 2.2, 2.0},
        {4.0, 10.0, 0.2, 1.0, 0.0, 0.0, 0.0}
    };
    const double spins[] = {0.7, 0.5, 0.7, 0.0};
    const double expected[][8] = {
        {4.0, 8.224358688899057, 0.9794378508566266, 2.0,
         1.2711320859373498, 0.03210895662807184,
         0.03970891015040996, -0.04649950039507937},
        {4.0, 11.529727036712782, 2.0653200939088854, -1.0,
         1.150260496055867, 0.075292205515767,
         -0.01067456097095506, 0.027113742729570633},
        {4.0, 9.0, 2.0612163513526314, 2.0,
         1.2149792379518414, 0.0,
         0.03147439129690756, 0.03517934783765223},
        {4.0, 8.333333333333334, 1.5707963267948966, 0.0,
         1.2601673613391169, 0.0, 0.0, 0.054583059137681106}
    };
    double y[8];
    size_t i, j;
    for (i = 0U; i < 4U; ++i) {
        assert(rp_initial_from_bound(spins[i], inputs[i], y)
            == RP_KERR_STATUS_OK);
        for (j = 0U; j < 8U; ++j) {
            assert(close_to(y[j], expected[i][j], 4e-12));
        }
        check_norm(spins[i], y);
    }
}

static void test_phases(void)
{
    const double pi = acos(-1.0);
    double input[7] = {0.0, 9.0, 0.2, 0.6, 0.0, 0.0, 0.3};
    double y[8], repeated[8];
    assert(rp_initial_from_bound(0.7, input, y) == RP_KERR_STATUS_OK);
    assert(close_to(y[1], 7.5, 1e-13));
    assert(close_to(y[2], asin(0.6), 1e-13));
    assert(y[5] == 0.0 && y[6] == 0.0);
    input[4] = pi;
    input[5] = pi;
    assert(rp_initial_from_bound(0.7, input, y) == RP_KERR_STATUS_OK);
    assert(close_to(y[1], 11.25, 1e-13));
    assert(close_to(y[2], pi - asin(0.6), 1e-13));
    input[4] = 0.6;
    input[5] = 0.8;
    assert(rp_initial_from_bound(0.7, input, y) == RP_KERR_STATUS_OK);
    input[4] += 2.0 * pi;
    input[5] += 2.0 * pi;
    assert(rp_initial_from_bound(0.7, input, repeated)
        == RP_KERR_STATUS_OK);
    assert(close_to(y[1], repeated[1], 1e-13));
    assert(close_to(y[2], repeated[2], 1e-13));
    assert(close_to(y[5], repeated[5], 1e-13));
    assert(close_to(y[6], repeated[6], 1e-13));
    input[3] = 0.0;
    input[5] = 1.5;
    assert(rp_initial_from_bound(0.7, input, y) == RP_KERR_STATUS_OK);
    check_norm(0.7, y);
    input[3] = -1.0;
    input[1] = 15.0;
    assert(rp_initial_from_bound(1.0, input, y) == RP_KERR_STATUS_OK);
    assert(y[7] < 0.0);
    check_norm(1.0, y);
}

static void test_errors(void)
{
    double input[7] = {0.0, 9.0, 0.2, 0.6, 0.0, 0.0, 0.0};
    double y[8];
    size_t i;
    assert(rp_initial_from_bound(0.7, NULL, y)
        == RP_KERR_STATUS_NULL_POINTER);
    input[2] = 1.0;
    assert(rp_initial_from_bound(0.7, input, y)
        == RP_KERR_STATUS_INVALID_PARAMETER);
    input[2] = 0.2;
    input[3] = 1.01;
    assert(rp_initial_from_bound(0.7, input, y)
        == RP_KERR_STATUS_INVALID_PARAMETER);
    input[3] = 1.0;
    input[1] = 6.4;
    assert(rp_initial_from_bound(0.0, input, y)
        == RP_KERR_STATUS_INVALID_PARAMETER); /* p=6+2e separatrix */
    input[1] = 5.9;
    input[2] = 0.0;
    assert(rp_initial_from_bound(0.0, input, y)
        == RP_KERR_STATUS_INVALID_PARAMETER);
    input[1] = 6.0;
    assert(rp_initial_from_bound(0.0, input, y)
        == RP_KERR_STATUS_INVALID_PARAMETER); /* marginal circular orbit */
    input[1] = 9.0;
    input[2] = 0.2;
    input[3] = 0.0;
    assert(rp_initial_from_bound(0.7, input, y)
        == RP_KERR_STATUS_COORDINATE_SINGULARITY);
    input[4] = NAN;
    assert(rp_initial_from_bound(0.7, input, y)
        == RP_KERR_STATUS_NONFINITE_INPUT);
    for (i = 0U; i < 8U; ++i) {
        assert(y[i] == 0.0);
    }
}

static void test_nearly_circular_schwarzschild(void)
{
    /* Closed Schwarzschild turning-point constants avoid a second numerical
       solver losing accuracy when the two radial turning points coalesce. */
    const double eccentricities[] = {0.0, 1e-12, 1e-9, 1e-7, 0.2};
    double input[7] = {0.0, 12.0, 0.0, 0.6, 0.8, 1.2, 0.3};
    size_t i;
    for (i = 0U; i < sizeof(eccentricities) / sizeof(eccentricities[0]); ++i) {
        double y[8], metric[4][4];
        double e = eccentricities[i];
        double energy, angular;
        const double p = input[1];
        input[2] = e;
        assert(rp_initial_from_bound(0.0, input, y) == RP_KERR_STATUS_OK);
        assert(rp_kerr_metric(1.0, 0.0, y, metric) == RP_KERR_STATUS_OK);
        energy = -metric[0][0] * y[4];
        angular = metric[3][3] * y[7];
        assert(close_to(energy, sqrt(((p - 2.0) * (p - 2.0) - 4.0 * e * e)
                                    / (p * (p - 3.0 - e * e))), 2e-12));
        assert(close_to(angular, input[3] * p / sqrt(p - 3.0 - e * e), 2e-12));
        check_norm(0.0, y);
    }
}

static void test_exact_turning_velocities(void)
{
    double input[7] = {0.0, 12.0, 0.3, 0.0, 0.0, 0.0, 0.0};
    double y[8];
    const double inclinations[] = {0.5, 0.8, 0.9};
    size_t i;
    for (i = 0U; i < 3U; ++i) {
        input[3] = inclinations[i];
        input[4] = input[5] = 0.0;
        assert(rp_initial_from_bound(0.7, input, y) == RP_KERR_STATUS_OK);
        assert(y[5] == 0.0 && y[6] == 0.0);
        input[4] = input[5] = acos(-1.0);
        assert(rp_initial_from_bound(0.7, input, y) == RP_KERR_STATUS_OK);
        assert(y[5] == 0.0 && y[6] == 0.0);
        check_norm(0.7, y);
    }
}

int main(void)
{
    test_references();
    test_phases();
    test_errors();
    test_nearly_circular_schwarzschild();
    test_exact_turning_velocities();
    return 0;
}
