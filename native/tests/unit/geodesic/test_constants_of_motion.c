#include "geodesic/solution/reconstruct.h"

#include <assert.h>
#include <math.h>

/* Bardeen, Press & Teukolsky (1972), Eqs. (2.12)-(2.13), prograde, M = 1. */
static void circular_equatorial(
    double spin, double radius, double polar_velocity, double state[8]
)
{
    const double root = sqrt(radius);
    const double denominator = pow(radius, 0.75)
        * sqrt(radius * root - 3.0 * root + 2.0 * spin);
    const double ut = (radius * root + spin) / denominator;

    state[0] = 0.0;
    state[1] = radius;
    state[2] = acos(0.0);
    state[3] = 0.4;
    state[4] = ut;
    state[5] = 0.0;
    state[6] = polar_velocity;
    state[7] = ut / (radius * root + spin);
}

static int close(double actual, double expected)
{
    return fabs(actual - expected) <= 1.0e-13 * (1.0 + fabs(expected));
}

int main(void)
{
    const double spin = 0.7;
    const double radius = 9.0;
    const double root = sqrt(radius);
    const double denominator = pow(radius, 0.75)
        * sqrt(radius * root - 3.0 * root + 2.0 * spin);
    const double energy = (radius * root - 2.0 * root + spin) / denominator;
    const double angular_momentum
        = (radius * radius - 2.0 * spin * root + spin * spin) / denominator;
    double states[4][8];
    double constants[4][3];
    rp_kerr_status statuses[4];
    size_t column;

    circular_equatorial(spin, radius, 0.0, states[0]);
    /* At theta = pi/2, Q = u_theta^2 = (r^2 u^theta)^2 and E, Lz are unchanged. */
    circular_equatorial(spin, radius, 0.01, states[1]);
    circular_equatorial(spin, radius, 0.0, states[2]);
    states[2][2] = 0.0;
    circular_equatorial(spin, radius, 0.0, states[3]);
    states[3][5] = NAN;

    assert(rp_solution_constants_of_motion_batch(spin, &states[0][0], 4U,
        &constants[0][0], statuses) != RP_KERR_STATUS_OK);
    assert(statuses[0] == RP_KERR_STATUS_OK);
    assert(statuses[1] == RP_KERR_STATUS_OK);
    assert(statuses[2] != RP_KERR_STATUS_OK);
    assert(statuses[3] == RP_KERR_STATUS_NONFINITE_INPUT);
    assert(close(constants[0][0], energy));
    assert(close(constants[0][1], angular_momentum));
    assert(fabs(constants[0][2]) <= 1.0e-24);
    assert(close(constants[1][0], energy));
    assert(close(constants[1][1], angular_momentum));
    assert(close(constants[1][2], pow(radius * radius * 0.01, 2.0)));
    for (column = 0U; column < 3U; ++column) {
        assert(constants[2][column] == 0.0);
        assert(constants[3][column] == 0.0);
    }

    assert(rp_solution_constants_of_motion_batch(1.5, &states[0][0], 1U,
        &constants[0][0], statuses) == RP_KERR_STATUS_INVALID_PARAMETER);
    assert(constants[0][0] == 0.0);
    assert(rp_solution_constants_of_motion_batch(spin, NULL, 0U, NULL, NULL)
        == RP_KERR_STATUS_OK);
    return 0;
}
