/*
 * Built-in Plummer model: agreement with the analytic test model, the
 * acceleration as minus the gradient of the potential (central
 * differences), a zero mass giving zero, and rejected parameters.
 */

#include "force/models/plummer.h"
#include "force_test_models.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static const double POINTS[4][3] = {
    {8.0, -3.0, 1.5}, {0.0, 0.0, 0.0}, {-120.0, 40.0, -300.0}, {0.2, 0.1, -0.3}
};
static const double VELOCITY[3] = {0.01, -0.02, 0.003};

static void check_against_test_model(void)
{
    const double model[RP_FORCE_PLUMMER_PARAMETERS] = {0.03, 4.0};
    const test_plummer reference = {0.03, 4.0, {0.0, 0.0, 0.0}};
    double value[3];
    double expected[3];
    double potential;
    size_t i;
    size_t k;

    for (i = 0U; i < 4U; ++i) {
        assert(rp_force_plummer_acceleration(POINTS[i], VELOCITY, model, value) == 0);
        test_plummer_acceleration(POINTS[i], VELOCITY, &reference, expected);
        for (k = 0U; k < 3U; ++k) {
            assert(value[k] == expected[k]);
        }
        assert(rp_force_plummer_potential(POINTS[i], model, &potential) == 0);
        assert(potential == test_plummer_potential(&reference, POINTS[i]));
    }
}

static void check_gradient(void)
{
    const double model[RP_FORCE_PLUMMER_PARAMETERS] = {0.5, 2.5};
    double value[3];
    double shifted[3];
    double plus;
    double minus;
    size_t i;
    size_t k;

    for (i = 0U; i < 4U; ++i) {
        assert(rp_force_plummer_acceleration(POINTS[i], VELOCITY, model, value) == 0);
        for (k = 0U; k < 3U; ++k) {
            const double step = 1e-5 * (1.0 + fabs(POINTS[i][k]));
            double gradient;
            shifted[0] = POINTS[i][0];
            shifted[1] = POINTS[i][1];
            shifted[2] = POINTS[i][2];
            shifted[k] = POINTS[i][k] + step;
            assert(rp_force_plummer_potential(shifted, model, &plus) == 0);
            shifted[k] = POINTS[i][k] - step;
            assert(rp_force_plummer_potential(shifted, model, &minus) == 0);
            gradient = (plus - minus) / (2.0 * step);
            assert(fabs(value[k] + gradient) <= 1e-8 * (fabs(value[k]) + 1e-6));
        }
    }
}

static void check_zero_mass_and_failures(void)
{
    const double zero[RP_FORCE_PLUMMER_PARAMETERS] = {0.0, 3.0};
    const double flat[RP_FORCE_PLUMMER_PARAMETERS] = {1.0, 0.0};
    const double nonfinite[RP_FORCE_PLUMMER_PARAMETERS] = {NAN, 3.0};
    double value[3] = {1.0, 1.0, 1.0};
    double potential = 1.0;

    assert(rp_force_plummer_acceleration(POINTS[0], VELOCITY, zero, value) == 0);
    assert(value[0] == 0.0 && value[1] == 0.0 && value[2] == 0.0);
    assert(rp_force_plummer_potential(POINTS[0], zero, &potential) == 0);
    assert(potential == 0.0);

    value[0] = 1.0;
    assert(rp_force_plummer_acceleration(POINTS[0], VELOCITY, flat, value) == 1);
    assert(value[0] == 0.0 && value[1] == 0.0 && value[2] == 0.0);
    assert(rp_force_plummer_acceleration(POINTS[0], VELOCITY, nonfinite, value) == 1);
    assert(rp_force_plummer_acceleration(POINTS[0], VELOCITY, NULL, value) == 1);
    potential = 1.0;
    assert(rp_force_plummer_potential(POINTS[0], flat, &potential) == 1);
    assert(potential == 0.0);
}

int main(void)
{
    check_against_test_model();
    check_gradient();
    check_zero_mass_and_failures();
    puts("test_plummer_model: ok");
    return 0;
}
