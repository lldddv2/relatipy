#include "geodesic/solution/reconstruct.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stddef.h>

static int close_to(double actual, double expected, double tolerance)
{
    return fabs(actual - expected) <= tolerance * (1.0 + fabs(expected));
}

static void assert_physical_columns_finite(const double *reconstructed)
{
    size_t column;

    for (column = 0U; column < RP_SOL_SEMIMAJOR; ++column) {
        assert(isfinite(reconstructed[column]));
    }
}

static void assert_angles_undefined(const double *reconstructed)
{
    size_t column;

    for (column = RP_SOL_INCLINATION;
         column < RP_SOLUTION_RECONSTRUCTED_DIM; ++column) {
        assert(isnan(reconstructed[column]));
    }
}

static void forward_cartesian(
    const double state[7],
    double spin,
    double cartesian[RP_SOLUTION_CARTESIAN_DIM]
)
{
    const double r = state[1];
    const double theta = state[2];
    const double phi = state[3];
    const double vr = state[4];
    const double vtheta = state[5];
    const double vphi = state[6];
    const double radial_scale = hypot(r, spin);
    const double radial_scale_dot = r / radial_scale * vr;
    const double st = sin(theta);
    const double ct = cos(theta);
    const double sp = sin(phi);
    const double cp = cos(phi);

    cartesian[0] = state[0];
    cartesian[1] = radial_scale * st * cp;
    cartesian[2] = radial_scale * st * sp;
    cartesian[3] = r * ct;
    cartesian[4] = radial_scale_dot * st * cp
        + radial_scale * ct * cp * vtheta
        - radial_scale * st * sp * vphi;
    cartesian[5] = radial_scale_dot * st * sp
        + radial_scale * ct * sp * vtheta
        + radial_scale * st * cp * vphi;
    cartesian[6] = vr * ct - r * st * vtheta;
}

static void test_schwarzschild_circular_state(void)
{
    const double cartesian[7] = {
        3.0, 10.0, 0.0, 0.0, 0.0, 0.31622776601683794, 0.0
    };
    double reconstructed[RP_SOLUTION_RECONSTRUCTED_DIM];
    double expected_four_velocity[RP_KERR_DIM];
    const double coordinates[RP_KERR_DIM] = {
        3.0, 10.0, 1.57079632679489661923, 0.0
    };
    const double velocity[3] = {0.0, 0.0, 0.031622776601683794};
    rp_kerr_status row_status;

    assert(rp_solution_reconstruct_batch(
        0.0, cartesian, 1U, reconstructed, &row_status
    ) == RP_KERR_STATUS_OK);
    assert(row_status == RP_KERR_STATUS_OK);
    assert(close_to(reconstructed[RP_SOL_BL_R], 10.0, 1e-14));
    assert(close_to(reconstructed[RP_SOL_BL_THETA], coordinates[2], 1e-14));
    assert(close_to(reconstructed[RP_SOL_BL_VPHI], velocity[2], 1e-14));
    assert(close_to(reconstructed[RP_SOL_SPH_R], 10.0, 1e-14));
    assert(close_to(reconstructed[RP_SOL_SPH_VPHI], velocity[2], 1e-14));
    assert(close_to(reconstructed[RP_SOL_SEMIMAJOR], 10.0, 1e-13));
    assert(reconstructed[RP_SOL_ECCENTRICITY] == 0.0);
    assert(reconstructed[RP_SOL_INCLINATION] == 0.0);
    assert(reconstructed[RP_SOL_ASCENDING_NODE] == 0.0);
    assert(reconstructed[RP_SOL_PERIAPSIS_ARGUMENT] == 0.0);
    assert(reconstructed[RP_SOL_TRUE_ANOMALY] == 0.0);

    assert(rp_kerr_four_velocity(
        1.0, 0.0, coordinates, velocity, expected_four_velocity
    ) == RP_KERR_STATUS_OK);
    assert(close_to(reconstructed[RP_SOL_UT], expected_four_velocity[0], 1e-14));
    assert(close_to(reconstructed[RP_SOL_BL_UPHI], expected_four_velocity[3], 1e-14));
    assert(close_to(reconstructed[RP_SOL_UY],
        expected_four_velocity[0] * cartesian[5], 1e-14));
}

static void test_spinning_inverse_and_four_velocity(void)
{
    const double state[7] = {2.0, 8.0, 1.2, 0.7, 0.012, -0.001, 0.008};
    const double velocity[3] = {0.012, -0.001, 0.008};
    double cartesian[RP_SOLUTION_CARTESIAN_DIM];
    double original[RP_SOLUTION_CARTESIAN_DIM];
    double reconstructed[RP_SOLUTION_RECONSTRUCTED_DIM];
    double four_velocity[RP_KERR_DIM];
    rp_kerr_status row_status;
    size_t column;

    forward_cartesian(state, 1.0, cartesian);
    for (column = 0U; column < RP_SOLUTION_CARTESIAN_DIM; ++column) {
        original[column] = cartesian[column];
    }
    assert(rp_solution_reconstruct_batch(
        1.0, cartesian, 1U, reconstructed, &row_status
    ) == RP_KERR_STATUS_OK);
    assert(row_status == RP_KERR_STATUS_OK);
    for (column = 0U; column < RP_SOLUTION_CARTESIAN_DIM; ++column) {
        assert(cartesian[column] == original[column]);
    }
    for (column = 0U; column < 4U; ++column) {
        assert(close_to(reconstructed[column], state[column], 2e-14));
    }
    for (column = 0U; column < 3U; ++column) {
        assert(close_to(reconstructed[RP_SOL_BL_VR + column],
            velocity[column], 2e-14));
    }
    assert(rp_kerr_four_velocity(
        1.0, 1.0, state, velocity, four_velocity
    ) == RP_KERR_STATUS_OK);
    for (column = 0U; column < 4U; ++column) {
        assert(close_to(reconstructed[RP_SOL_UT + column],
            four_velocity[column], 2e-14));
    }
    for (column = 0U; column < 3U; ++column) {
        assert(close_to(reconstructed[RP_SOL_UX + column],
            four_velocity[0] * cartesian[4U + column], 2e-14));
        assert(close_to(reconstructed[RP_SOL_SPH_UR + column],
            four_velocity[0] * reconstructed[RP_SOL_SPH_VR + column], 2e-14));
    }
}

static void test_batch_row_failures(void)
{
    const double cartesian[4][7] = {
        {0.0, 10.0, 0.0, 0.0, 0.0, 0.1, 0.0},
        {0.0, 0.0, 0.0, 10.0, 0.0, 0.1, 0.0},
        {0.0, 10.0, 0.0, 0.0, 0.0, 2.0, 0.0},
        {0.0, 8.0, 0.0, 0.0, 0.0, 0.1, 0.0}
    };
    double reconstructed[4][RP_SOLUTION_RECONSTRUCTED_DIM];
    rp_kerr_status row_status[4];
    size_t column;

    assert(rp_solution_reconstruct_batch(
        0.0, &cartesian[0][0], 4U, &reconstructed[0][0], row_status
    ) == RP_KERR_STATUS_COORDINATE_SINGULARITY);
    assert(row_status[0] == RP_KERR_STATUS_OK);
    assert(row_status[1] == RP_KERR_STATUS_COORDINATE_SINGULARITY);
    assert(row_status[2] == RP_KERR_STATUS_NON_TIMELIKE_VELOCITY);
    assert(row_status[3] == RP_KERR_STATUS_OK);
    for (column = 0U; column < RP_SOLUTION_RECONSTRUCTED_DIM; ++column) {
        assert(reconstructed[1][column] == 0.0);
        assert(reconstructed[2][column] == 0.0);
    }
}

static void test_osculating_element_conventions(void)
{
    const double semimajor = 12.0;
    const double eccentricity = 0.25;
    const double inclination = 0.6;
    const double ascending_node = 0.3;
    const double periapsis_argument = 0.5;
    const double radius = semimajor * (1.0 - eccentricity);
    const double speed = sqrt((1.0 + eccentricity) / radius);
    const double cp = cos(ascending_node);
    const double sp = sin(ascending_node);
    const double ci = cos(inclination);
    const double si = sin(inclination);
    const double cw = cos(periapsis_argument);
    const double sw = sin(periapsis_argument);
    const double peri_p[3] = {
        cp * cw - sp * sw * ci,
        sp * cw + cp * sw * ci,
        sw * si
    };
    const double peri_q[3] = {
        -cp * sw - sp * cw * ci,
        -sp * sw + cp * cw * ci,
        cw * si
    };
    double cartesian[7] = {0.0};
    double reconstructed[RP_SOLUTION_RECONSTRUCTED_DIM];
    rp_kerr_status row_status;
    size_t column;

    for (column = 0U; column < 3U; ++column) {
        cartesian[1U + column] = radius * peri_p[column];
        cartesian[4U + column] = speed * peri_q[column];
    }
    assert(rp_solution_reconstruct_batch(
        0.0, cartesian, 1U, reconstructed, &row_status
    ) == RP_KERR_STATUS_OK);
    assert(close_to(reconstructed[RP_SOL_SEMIMAJOR], semimajor, 3e-14));
    assert(close_to(reconstructed[RP_SOL_ECCENTRICITY], eccentricity, 3e-14));
    assert(close_to(reconstructed[RP_SOL_INCLINATION], inclination, 3e-14));
    assert(close_to(reconstructed[RP_SOL_ASCENDING_NODE], ascending_node, 3e-14));
    assert(close_to(reconstructed[RP_SOL_PERIAPSIS_ARGUMENT],
        periapsis_argument, 3e-14));
    assert(close_to(reconstructed[RP_SOL_TRUE_ANOMALY], 0.0, 3e-14));

    /* Equatorial eccentric: Omega=0 and omega is longitude of periapsis. */
    cartesian[1] = radius * cw;
    cartesian[2] = radius * sw;
    cartesian[3] = 0.0;
    cartesian[4] = -speed * sw;
    cartesian[5] = speed * cw;
    cartesian[6] = 0.0;
    assert(rp_solution_reconstruct_batch(
        0.0, cartesian, 1U, reconstructed, &row_status
    ) == RP_KERR_STATUS_OK);
    assert(reconstructed[RP_SOL_ASCENDING_NODE] == 0.0);
    assert(close_to(reconstructed[RP_SOL_PERIAPSIS_ARGUMENT],
        periapsis_argument, 3e-14));

    /* Circular inclined: omega=0 and f is argument of latitude. */
    {
        const double phase = 0.4;
        const double circular_radius = 12.0;
        const double circular_speed = sqrt(1.0 / circular_radius);
        const double node_p[3] = {cp, sp, 0.0};
        const double node_q[3] = {-sp * ci, cp * ci, si};
        for (column = 0U; column < 3U; ++column) {
            cartesian[1U + column] = circular_radius
                * (cos(phase) * node_p[column] + sin(phase) * node_q[column]);
            cartesian[4U + column] = circular_speed
                * (-sin(phase) * node_p[column] + cos(phase) * node_q[column]);
        }
        assert(rp_solution_reconstruct_batch(
            0.0, cartesian, 1U, reconstructed, &row_status
        ) == RP_KERR_STATUS_OK);
        assert(reconstructed[RP_SOL_ECCENTRICITY] == 0.0);
        assert(reconstructed[RP_SOL_PERIAPSIS_ARGUMENT] == 0.0);
        assert(close_to(reconstructed[RP_SOL_TRUE_ANOMALY], phase, 3e-14));
    }
}

static void test_degenerate_elements_preserve_physical_states(void)
{
    const double cartesian[5][RP_SOLUTION_CARTESIAN_DIM] = {
        {2.0, 10.0, 0.0, 0.0, -0.01, 0.0, 0.0},
        {3.0, 10.0, 0.0, 0.0, 0.01, 0.0, 0.0},
        {4.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {5.0, 32.0, 0.0, 0.0, 0.0, 0.25, 0.0},
        {6.0, 32.0, 0.0, 0.0, 0.25, 0.0, 0.0}
    };
    double reconstructed[5][RP_SOLUTION_RECONSTRUCTED_DIM];
    rp_kerr_status row_status[5];
    size_t row;

    assert(rp_solution_reconstruct_batch(
        0.0, &cartesian[0][0], 5U, &reconstructed[0][0], row_status
    ) == RP_KERR_STATUS_OK);
    for (row = 0U; row < 5U; ++row) {
        assert(row_status[row] == RP_KERR_STATUS_OK);
        assert_physical_columns_finite(reconstructed[row]);
        assert(reconstructed[row][RP_SOL_T] == cartesian[row][0]);
        assert(close_to(reconstructed[row][RP_SOL_BL_R],
            cartesian[row][1], 1e-14));
        assert(close_to(reconstructed[row][RP_SOL_BL_VR],
            cartesian[row][4], 1e-14));
        assert(isfinite(reconstructed[row][RP_SOL_ECCENTRICITY]));
        assert(close_to(reconstructed[row][RP_SOL_ECCENTRICITY], 1.0, 1e-14));
    }
    assert_angles_undefined(reconstructed[0]);
    assert_angles_undefined(reconstructed[1]);
    assert_angles_undefined(reconstructed[2]);
    assert_angles_undefined(reconstructed[4]);
    assert(close_to(reconstructed[0][RP_SOL_SEMIMAJOR],
        -0.5 / (0.5 * 0.01 * 0.01 - 0.1), 1e-14));
    assert(reconstructed[0][RP_SOL_SEMIMAJOR]
        == reconstructed[1][RP_SOL_SEMIMAJOR]);
    assert(reconstructed[2][RP_SOL_SEMIMAJOR] == 5.0);
    assert(reconstructed[2][RP_SOL_BL_VR] == 0.0);
    assert(reconstructed[2][RP_SOL_UX] == 0.0);
    assert(reconstructed[3][RP_SOL_SEMIMAJOR] == INFINITY);
    assert(reconstructed[4][RP_SOL_SEMIMAJOR] == INFINITY);
    assert(reconstructed[3][RP_SOL_INCLINATION] == 0.0);
    assert(reconstructed[3][RP_SOL_ASCENDING_NODE] == 0.0);
    assert(reconstructed[3][RP_SOL_PERIAPSIS_ARGUMENT] == 0.0);
    assert(reconstructed[3][RP_SOL_TRUE_ANOMALY] == 0.0);
}

static void test_unresolved_angular_momentum_threshold(void)
{
    double cartesian[7] = {0.0, 10.0, 0.0, 0.0, 0.1, 0.0, 0.0};
    const double threshold_factors[3] = {64.0, 128.0, 256.0};
    double reconstructed[RP_SOLUTION_RECONSTRUCTED_DIM];
    rp_kerr_status row_status;
    size_t row;
    size_t column;

    for (row = 0U; row < 3U; ++row) {
        cartesian[5] = threshold_factors[row] * DBL_EPSILON * cartesian[4];
        assert(rp_solution_reconstruct_batch(
            0.0, cartesian, 1U, reconstructed, &row_status
        ) == RP_KERR_STATUS_OK);
        assert(row_status == RP_KERR_STATUS_OK);
        assert_physical_columns_finite(reconstructed);
        assert(isfinite(reconstructed[RP_SOL_SEMIMAJOR]));
        assert(isfinite(reconstructed[RP_SOL_ECCENTRICITY]));
        if (row < 2U) {
            assert_angles_undefined(reconstructed);
        } else {
            for (column = RP_SOL_INCLINATION;
                 column < RP_SOLUTION_RECONSTRUCTED_DIM; ++column) {
                assert(isfinite(reconstructed[column]));
            }
        }
    }
}

static void test_near_parabolic_energy_remains_finite(void)
{
    double cartesian[7] = {0.0, 32.0, 0.0, 0.0, 0.0, 0.25, 0.0};
    double reconstructed[RP_SOLUTION_RECONSTRUCTED_DIM];
    rp_kerr_status row_status;
    const double directions[2] = {0.0, INFINITY};
    size_t row;

    for (row = 0U; row < 2U; ++row) {
        double energy;
        cartesian[5] = nextafter(0.25, directions[row]);
        energy = 0.5 * cartesian[5] * cartesian[5] - 1.0 / cartesian[1];
        assert(energy != 0.0);
        assert(rp_solution_reconstruct_batch(
            0.0, cartesian, 1U, reconstructed, &row_status
        ) == RP_KERR_STATUS_OK);
        assert(row_status == RP_KERR_STATUS_OK);
        assert(isfinite(reconstructed[RP_SOL_SEMIMAJOR]));
        assert(reconstructed[RP_SOL_SEMIMAJOR] == -0.5 / energy);
        assert_physical_columns_finite(reconstructed);
    }
}

static void test_mixed_degenerate_rows_and_failures(void)
{
    const double cartesian[8][RP_SOLUTION_CARTESIAN_DIM] = {
        {0.0, 10.0, 0.0, 0.0, 0.0, 0.1, 0.0},
        {1.0, 10.0, 0.0, 0.0, -0.01, 0.0, 0.0},
        {2.0, 0.0, 0.0, 10.0, 0.0, 0.1, 0.0},
        {3.0, 32.0, 0.0, 0.0, 0.0, 0.25, 0.0},
        {4.0, 10.0, 0.0, 0.0, 0.0, 2.0, 0.0},
        {5.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {6.0, DBL_MAX, 0.0, 0.0, 0.0, 0.1, 0.0},
        {7.0, 10.0, 0.0, 0.0, NAN, 0.1, 0.0}
    };
    const rp_kerr_status expected[8] = {
        RP_KERR_STATUS_OK,
        RP_KERR_STATUS_OK,
        RP_KERR_STATUS_COORDINATE_SINGULARITY,
        RP_KERR_STATUS_OK,
        RP_KERR_STATUS_NON_TIMELIKE_VELOCITY,
        RP_KERR_STATUS_OK,
        RP_KERR_STATUS_NUMERICAL_RANGE,
        RP_KERR_STATUS_NONFINITE_INPUT
    };
    double reconstructed[8][RP_SOLUTION_RECONSTRUCTED_DIM];
    rp_kerr_status row_status[8];
    size_t row;
    size_t column;

    for (row = 0U; row < 8U; ++row) {
        for (column = 0U; column < RP_SOLUTION_RECONSTRUCTED_DIM; ++column) {
            reconstructed[row][column] = 9.0;
        }
    }
    assert(rp_solution_reconstruct_batch(
        0.0, &cartesian[0][0], 8U, &reconstructed[0][0], row_status
    ) == RP_KERR_STATUS_COORDINATE_SINGULARITY);
    for (row = 0U; row < 8U; ++row) {
        assert(row_status[row] == expected[row]);
        if (expected[row] == RP_KERR_STATUS_OK) {
            assert_physical_columns_finite(reconstructed[row]);
            assert(reconstructed[row][RP_SOL_T] == cartesian[row][0]);
        } else {
            for (column = 0U; column < RP_SOLUTION_RECONSTRUCTED_DIM; ++column) {
                assert(reconstructed[row][column] == 0.0);
            }
        }
    }
    assert_angles_undefined(reconstructed[1]);
    assert_angles_undefined(reconstructed[5]);
    assert(reconstructed[3][RP_SOL_SEMIMAJOR] == INFINITY);
}

static void test_global_errors(void)
{
    const double cartesian[7] = {0.0, 10.0, 0.0, 0.0, 0.0, 0.1, 0.0};
    double reconstructed[RP_SOLUTION_RECONSTRUCTED_DIM];
    rp_kerr_status row_status;
    size_t column;

    for (column = 0U; column < RP_SOLUTION_RECONSTRUCTED_DIM; ++column) {
        reconstructed[column] = 9.0;
    }
    assert(rp_solution_reconstruct_batch(
        -0.1, cartesian, 1U, reconstructed, &row_status
    ) == RP_KERR_STATUS_INVALID_PARAMETER);
    assert(row_status == RP_KERR_STATUS_INVALID_PARAMETER);
    for (column = 0U; column < RP_SOLUTION_RECONSTRUCTED_DIM; ++column) {
        assert(reconstructed[column] == 0.0);
    }
    assert(rp_solution_reconstruct_batch(
        0.0, NULL, 1U, reconstructed, &row_status
    ) == RP_KERR_STATUS_NULL_POINTER);
    assert(row_status == RP_KERR_STATUS_NULL_POINTER);
    assert(rp_solution_reconstruct_batch(
        0.0, cartesian, 1U, NULL, &row_status
    ) == RP_KERR_STATUS_NULL_POINTER);
}

static void test_domain_errors(void)
{
    double cartesian[7] = {0.0, 10.0, 0.0, 0.0, 0.0, 0.1, 0.0};
    double reconstructed[RP_SOLUTION_RECONSTRUCTED_DIM];
    rp_kerr_status row_status;
    size_t column;

    cartesian[4] = NAN;
    assert(rp_solution_reconstruct_batch(
        0.0, cartesian, 1U, reconstructed, &row_status
    ) == RP_KERR_STATUS_NONFINITE_INPUT);
    cartesian[4] = 0.0;
    cartesian[1] = 1.5;
    assert(rp_solution_reconstruct_batch(
        0.0, cartesian, 1U, reconstructed, &row_status
    ) == RP_KERR_STATUS_INVALID_PARAMETER);
    cartesian[1] = 10.0;
    cartesian[4] = 0.0;
    cartesian[5] = 0.5;
    assert(rp_solution_reconstruct_batch(
        0.0, cartesian, 1U, reconstructed, &row_status
    ) == RP_KERR_STATUS_OK);
    assert(close_to(reconstructed[RP_SOL_SEMIMAJOR], -20.0, 2e-14));
    assert(close_to(reconstructed[RP_SOL_ECCENTRICITY], 1.5, 2e-14));
    cartesian[4] = NAN;
    assert(rp_solution_reconstruct_batch(
        0.0, cartesian, 1U, reconstructed, &row_status
    ) == RP_KERR_STATUS_NONFINITE_INPUT);
    for (column = 0U; column < RP_SOLUTION_RECONSTRUCTED_DIM; ++column) {
        assert(reconstructed[column] == 0.0);
    }
    cartesian[4] = 0.0;
    cartesian[1] = DBL_MAX;
    assert(rp_solution_reconstruct_batch(
        0.0, cartesian, 1U, reconstructed, &row_status
    ) == RP_KERR_STATUS_NUMERICAL_RANGE);
    assert(row_status == RP_KERR_STATUS_NUMERICAL_RANGE);
    for (column = 0U; column < RP_SOLUTION_RECONSTRUCTED_DIM; ++column) {
        assert(reconstructed[column] == 0.0);
    }
}

int main(void)
{
    test_schwarzschild_circular_state();
    test_spinning_inverse_and_four_velocity();
    test_batch_row_failures();
    test_osculating_element_conventions();
    test_degenerate_elements_preserve_physical_states();
    test_unresolved_angular_momentum_threshold();
    test_near_parabolic_energy_remains_finite();
    test_mixed_degenerate_rows_and_failures();
    test_global_errors();
    test_domain_errors();
    return 0;
}
