#include "relatipy/kerr_geodesic.h"

#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>

static int failures = 0;

#define CHECK(condition)                                                        \
    do {                                                                        \
        if (!(condition)) {                                                     \
            (void)fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__,     \
                #condition);                                                    \
            ++failures;                                                         \
        }                                                                       \
    } while (0)

static int close_enough(double actual, double expected, double relative, double absolute)
{
    const double scale = fmax(fabs(actual), fabs(expected));
    return fabs(actual - expected) <= absolute + relative * scale;
}

#define CHECK_CLOSE(actual, expected, relative, absolute)                       \
    do {                                                                        \
        const double actual_value = (actual);                                   \
        const double expected_value = (expected);                               \
        if (!close_enough(actual_value, expected_value, relative, absolute)) {  \
            (void)fprintf(stderr,                                               \
                "FAIL %s:%d: %.17g != %.17g (%s)\n",                          \
                __FILE__, __LINE__, actual_value, expected_value, #actual);     \
            ++failures;                                                         \
        }                                                                       \
    } while (0)

static double tangent_norm(
    double mass,
    double spin,
    const double coordinates[RP_KERR_DIM],
    const double tangent[RP_KERR_DIM]
)
{
    double metric[RP_KERR_DIM][RP_KERR_DIM];
    double norm = 0.0;
    size_t mu;
    size_t nu;

    CHECK(rp_kerr_metric(mass, spin, coordinates, metric)
        == RP_KERR_STATUS_OK);
    for (mu = 0U; mu < RP_KERR_DIM; ++mu) {
        for (nu = 0U; nu < RP_KERR_DIM; ++nu) {
            norm += metric[mu][nu] * tangent[mu] * tangent[nu];
        }
    }
    return norm;
}

static void test_null_tangent_from_constants(void)
{
    const double mass = 1.0;
    const double spin = 0.5;
    const double coordinates[RP_KERR_DIM] = {0.0, 8.0, 1.1, 0.3};
    const double energy = 1.0;
    const double angular_momentum = 2.0;
    const double carter_constant = 3.0;
    double metric[RP_KERR_DIM][RP_KERR_DIM];
    double tangent[RP_KERR_DIM];
    double recovered_energy;
    double recovered_angular_momentum;
    double recovered_carter_constant;
    double sigma;
    double sin_theta;
    double cos_theta;

    CHECK(rp_kerr_null_tangent(
        mass,
        spin,
        coordinates,
        energy,
        angular_momentum,
        carter_constant,
        -1,
        1,
        tangent
    ) == RP_KERR_STATUS_OK);
    CHECK(rp_kerr_metric(mass, spin, coordinates, metric)
        == RP_KERR_STATUS_OK);
    CHECK_CLOSE(tangent_norm(mass, spin, coordinates, tangent), 0.0,
        0.0, 2e-15);
    CHECK(tangent[0] > 0.0);
    CHECK(tangent[1] < 0.0);
    CHECK(tangent[2] > 0.0);

    recovered_energy = -(metric[0][0] * tangent[0]
        + metric[0][3] * tangent[3]);
    recovered_angular_momentum = metric[3][0] * tangent[0]
        + metric[3][3] * tangent[3];
    sin_theta = sin(coordinates[2]);
    cos_theta = cos(coordinates[2]);
    sigma = coordinates[1] * coordinates[1]
        + spin * spin * cos_theta * cos_theta;
    recovered_carter_constant = sigma * sigma * tangent[2] * tangent[2]
        + cos_theta * cos_theta
            * (recovered_angular_momentum * recovered_angular_momentum
                / (sin_theta * sin_theta)
                - spin * spin * recovered_energy * recovered_energy);

    CHECK_CLOSE(recovered_energy, energy, 3e-14, 3e-14);
    CHECK_CLOSE(recovered_angular_momentum, angular_momentum, 3e-14, 3e-14);
    CHECK_CLOSE(recovered_carter_constant, carter_constant, 3e-14, 3e-14);
}

static void test_schwarzschild_radial_null_tangent(void)
{
    const double radius = 10.0;
    const double coordinates[RP_KERR_DIM] = {
        0.0, radius, 1.5707963267948966, 0.0
    };
    double tangent[RP_KERR_DIM];

    CHECK(rp_kerr_null_tangent(
        1.0, 0.0, coordinates, 1.0, 0.0, 0.0, 1, -1, tangent
    ) == RP_KERR_STATUS_OK);
    CHECK_CLOSE(tangent[0], 1.0 / (1.0 - 2.0 / radius), 3e-14, 3e-14);
    CHECK_CLOSE(tangent[1], 1.0, 3e-14, 3e-14);
    CHECK_CLOSE(tangent[2], 0.0, 0.0, 0.0);
    CHECK_CLOSE(tangent[3], 0.0, 0.0, 0.0);
    CHECK_CLOSE(tangent_norm(1.0, 0.0, coordinates, tangent), 0.0,
        0.0, 2e-15);
}

static void test_equatorial_kerr_null_tangent(void)
{
    const double coordinates[RP_KERR_DIM] = {
        0.0, 20.0, 1.5707963267948966, 0.0
    };
    double tangent[RP_KERR_DIM];

    CHECK(rp_kerr_null_tangent(
        1.0, 0.5, coordinates, 1.0, 6.0, 0.0, -1, 1, tangent
    ) == RP_KERR_STATUS_OK);
    CHECK_CLOSE(tangent[2], 0.0, 0.0, 0.0);
    CHECK_CLOSE(tangent_norm(1.0, 0.5, coordinates, tangent), 0.0,
        0.0, 3e-15);
}

static void test_null_rhs_matches_numerical_metric_connection(void)
{
    const double mass = 1.0;
    const double spin = 0.5;
    const double coordinates[RP_KERR_DIM] = {0.0, 8.0, 1.1, 0.3};
    const double steps[RP_KERR_DIM] = {0.0, 1e-5, 1e-6, 0.0};
    double tangent[RP_KERR_DIM];
    double dx[RP_KERR_DIM];
    double du[RP_KERR_DIM];
    double inverse[RP_KERR_DIM][RP_KERR_DIM];
    double metric_plus[RP_KERR_DIM][RP_KERR_DIM];
    double metric_minus[RP_KERR_DIM][RP_KERR_DIM];
    double derivative[RP_KERR_DIM][RP_KERR_DIM][RP_KERR_DIM] = {{{0.0}}};
    double point_plus[RP_KERR_DIM];
    double point_minus[RP_KERR_DIM];
    size_t alpha;
    size_t beta;
    size_t lambda;
    size_t mu;
    size_t sigma;

    CHECK(rp_kerr_null_tangent(
        mass, spin, coordinates, 1.0, 2.0, 3.0, -1, 1, tangent
    ) == RP_KERR_STATUS_OK);
    CHECK(rp_kerr_null_geodesic_rhs(
        mass, spin, coordinates, tangent, dx, du
    ) == RP_KERR_STATUS_OK);
    CHECK(rp_kerr_inverse_metric(mass, spin, coordinates, inverse)
        == RP_KERR_STATUS_OK);

    for (alpha = 1U; alpha <= 2U; ++alpha) {
        for (mu = 0U; mu < RP_KERR_DIM; ++mu) {
            point_plus[mu] = coordinates[mu];
            point_minus[mu] = coordinates[mu];
        }
        point_plus[alpha] += steps[alpha];
        point_minus[alpha] -= steps[alpha];
        CHECK(rp_kerr_metric(mass, spin, point_plus, metric_plus)
            == RP_KERR_STATUS_OK);
        CHECK(rp_kerr_metric(mass, spin, point_minus, metric_minus)
            == RP_KERR_STATUS_OK);
        for (mu = 0U; mu < RP_KERR_DIM; ++mu) {
            for (beta = 0U; beta < RP_KERR_DIM; ++beta) {
                derivative[alpha][mu][beta] =
                    (metric_plus[mu][beta] - metric_minus[mu][beta])
                    / (2.0 * steps[alpha]);
            }
        }
    }

    for (lambda = 0U; lambda < RP_KERR_DIM; ++lambda) {
        double expected = 0.0;

        CHECK_CLOSE(dx[lambda], tangent[lambda], 0.0, 0.0);
        for (sigma = 0U; sigma < RP_KERR_DIM; ++sigma) {
            for (alpha = 0U; alpha < RP_KERR_DIM; ++alpha) {
                for (beta = 0U; beta < RP_KERR_DIM; ++beta) {
                    expected -= inverse[lambda][sigma]
                        * (derivative[alpha][sigma][beta]
                            - 0.5 * derivative[sigma][alpha][beta])
                        * tangent[alpha] * tangent[beta];
                }
            }
        }
        CHECK_CLOSE(du[lambda], expected, 2e-9, 2e-11);
    }
}

static void test_optimized_rhs_matches_christoffel_contraction(void)
{
    const double masses[3] = {1.0, 2.0, 0.75};
    const double spins[3] = {0.5, -1.4, 0.0};
    const double coordinates[3][RP_KERR_DIM] = {
        {0.0, 8.0, 1.1, 0.3},
        {0.2, 15.0, 0.8, -0.4},
        {-1.0, 6.0, 1.35, 2.0}
    };
    const double coordinate_velocities[3][3] = {
        {-0.01, 0.002, 0.02},
        {0.015, -0.001, -0.015},
        {-0.02, 0.003, 0.01}
    };
    double gamma[RP_KERR_DIM][RP_KERR_DIM][RP_KERR_DIM];
    double four_velocity[RP_KERR_DIM];
    double dx[RP_KERR_DIM];
    double du[RP_KERR_DIM];
    double expected_du;
    size_t point;
    size_t lambda;
    size_t mu;
    size_t nu;

    for (point = 0; point < 3; ++point) {
        CHECK(rp_kerr_four_velocity(
            masses[point],
            spins[point],
            coordinates[point],
            coordinate_velocities[point],
            four_velocity
        ) == RP_KERR_STATUS_OK);
        CHECK(rp_kerr_christoffel(
            masses[point], spins[point], coordinates[point], gamma)
            == RP_KERR_STATUS_OK);
        CHECK(rp_kerr_geodesic_rhs(
            masses[point],
            spins[point],
            coordinates[point],
            four_velocity,
            dx,
            du
        ) == RP_KERR_STATUS_OK);

        for (lambda = 0; lambda < RP_KERR_DIM; ++lambda) {
            expected_du = 0.0;
            for (mu = 0; mu < RP_KERR_DIM; ++mu) {
                for (nu = 0; nu < RP_KERR_DIM; ++nu) {
                    expected_du -= gamma[lambda][mu][nu]
                        * four_velocity[mu] * four_velocity[nu];
                }
            }
            CHECK_CLOSE(dx[lambda], four_velocity[lambda], 0.0, 0.0);
            CHECK_CLOSE(du[lambda], expected_du, 5e-13, 5e-14);
        }
    }
}

static void check_zero_vector(const double vector[RP_KERR_DIM])
{
    size_t index;

    for (index = 0; index < RP_KERR_DIM; ++index) {
        CHECK_CLOSE(vector[index], 0.0, 0.0, 0.0);
    }
}

static void fill_outputs(
    double coordinate_derivative[RP_KERR_DIM],
    double four_velocity_derivative[RP_KERR_DIM]
)
{
    size_t index;

    for (index = 0; index < RP_KERR_DIM; ++index) {
        coordinate_derivative[index] = 7.0;
        four_velocity_derivative[index] = 7.0;
    }
}

static void test_invalid_inputs_and_zeroed_outputs(void)
{
    const double regular[RP_KERR_DIM] = {0.0, 8.0, 1.1, 0.3};
    const double axis[RP_KERR_DIM] = {0.0, 8.0, 0.0, 0.3};
    const double horizon[RP_KERR_DIM] = {
        0.0, 1.0 + sqrt(0.75), 1.1, 0.3
    };
    const double ring[RP_KERR_DIM] = {0.0, 0.0, 1.5707963267948966, 0.3};
    const double overflow[RP_KERR_DIM] = {0.0, DBL_MAX, 1.1, 0.3};
    const double regular_u[RP_KERR_DIM] = {1.1, -0.01, 0.002, 0.02};
    const double nonfinite_u[RP_KERR_DIM] = {1.1, -0.01, NAN, 0.02};
    const double overflow_u[RP_KERR_DIM] = {DBL_MAX, DBL_MAX, 0.0, 0.0};
    double dx[RP_KERR_DIM];
    double du[RP_KERR_DIM];

    fill_outputs(dx, du);
    CHECK(rp_kerr_geodesic_rhs(1.0, 0.5, regular, regular_u, NULL, du)
        == RP_KERR_STATUS_NULL_POINTER);
    check_zero_vector(du);

    fill_outputs(dx, du);
    CHECK(rp_kerr_geodesic_rhs(1.0, 0.5, regular, regular_u, dx, NULL)
        == RP_KERR_STATUS_NULL_POINTER);
    check_zero_vector(dx);

#define CHECK_RHS_FAILURE(expected_status, mass, spin, point, velocity)         \
    do {                                                                        \
        fill_outputs(dx, du);                                                   \
        CHECK(rp_kerr_geodesic_rhs(                                             \
            (mass), (spin), (point), (velocity), dx, du)                       \
            == (expected_status));                                              \
        check_zero_vector(dx);                                                  \
        check_zero_vector(du);                                                  \
    } while (0)

    CHECK_RHS_FAILURE(
        RP_KERR_STATUS_NULL_POINTER, 1.0, 0.5, NULL, regular_u);
    CHECK_RHS_FAILURE(
        RP_KERR_STATUS_NULL_POINTER, 1.0, 0.5, regular, NULL);
    CHECK_RHS_FAILURE(
        RP_KERR_STATUS_INVALID_PARAMETER, 0.0, 0.0, regular, regular_u);
    CHECK_RHS_FAILURE(
        RP_KERR_STATUS_INVALID_PARAMETER, 1.0, 1.1, regular, regular_u);
    CHECK_RHS_FAILURE(
        RP_KERR_STATUS_NONFINITE_INPUT, 1.0, 0.5, regular, nonfinite_u);
    CHECK_RHS_FAILURE(
        RP_KERR_STATUS_COORDINATE_SINGULARITY, 1.0, 0.5, axis, regular_u);
    CHECK_RHS_FAILURE(
        RP_KERR_STATUS_COORDINATE_SINGULARITY, 1.0, 0.5, horizon, regular_u);
    CHECK_RHS_FAILURE(
        RP_KERR_STATUS_PHYSICAL_SINGULARITY, 1.0, 0.5, ring, regular_u);
    CHECK_RHS_FAILURE(
        RP_KERR_STATUS_NUMERICAL_RANGE, 1.0, 0.5, overflow, regular_u);
    CHECK_RHS_FAILURE(
        RP_KERR_STATUS_NUMERICAL_RANGE, 1.0, 0.5, regular, overflow_u);

#undef CHECK_RHS_FAILURE

    fill_outputs(dx, du);
    CHECK(rp_kerr_geodesic_rhs(1.0, 0.5, regular, regular_u, dx, dx)
        == RP_KERR_STATUS_INVALID_PARAMETER);
    check_zero_vector(dx);
}

static void test_invalid_null_inputs_and_zeroed_outputs(void)
{
    const double regular[RP_KERR_DIM] = {0.0, 8.0, 1.1, 0.3};
    const double regular_tangent[RP_KERR_DIM] = {1.0, -0.5, 0.01, 0.02};
    double tangent[RP_KERR_DIM];
    double dx[RP_KERR_DIM];
    double dk[RP_KERR_DIM];

    fill_outputs(tangent, dx);
    CHECK(rp_kerr_null_tangent(
        1.0, 0.5, regular, 1.0, 20.0, 0.0, 1, 1, tangent
    ) == RP_KERR_STATUS_NO_REAL_NULL_TANGENT);
    check_zero_vector(tangent);

    fill_outputs(tangent, dx);
    CHECK(rp_kerr_null_tangent(
        1.0, 0.5, regular, 1.0, 2.0, 3.0, 0, 1, tangent
    ) == RP_KERR_STATUS_INVALID_PARAMETER);
    check_zero_vector(tangent);

    fill_outputs(tangent, dx);
    CHECK(rp_kerr_null_tangent(
        1.0, 0.5, regular, NAN, 2.0, 3.0, 1, 1, tangent
    ) == RP_KERR_STATUS_NONFINITE_INPUT);
    check_zero_vector(tangent);

    fill_outputs(dx, dk);
    CHECK(rp_kerr_null_geodesic_rhs(
        1.0, 0.5, regular, regular_tangent, dx, dx
    ) == RP_KERR_STATUS_INVALID_PARAMETER);
    check_zero_vector(dx);

    fill_outputs(dx, dk);
    CHECK(rp_kerr_null_geodesic_rhs(
        1.0, 0.5, regular, regular_tangent, dx, NULL
    ) == RP_KERR_STATUS_NULL_POINTER);
    check_zero_vector(dx);
}

static void test_input_output_overlap(void)
{
    const double coordinates[RP_KERR_DIM] = {0.0, 8.0, 1.1, 0.3};
    const double coordinate_velocity[3] = {-0.01, 0.002, 0.02};
    double four_velocity[RP_KERR_DIM];
    double expected_dx[RP_KERR_DIM];
    double expected_du[RP_KERR_DIM];
    double expected_null_tangent[RP_KERR_DIM];
    double expected_null_dx[RP_KERR_DIM];
    double expected_null_dk[RP_KERR_DIM];
    double coordinate_dx[RP_KERR_DIM];
    double velocity_du[RP_KERR_DIM];
    double null_coordinate_dx[RP_KERR_DIM];
    double null_tangent_dk[RP_KERR_DIM];
    size_t index;

    CHECK(rp_kerr_four_velocity(
        1.0, 0.5, coordinates, coordinate_velocity, four_velocity)
        == RP_KERR_STATUS_OK);
    CHECK(rp_kerr_geodesic_rhs(
        1.0,
        0.5,
        coordinates,
        four_velocity,
        expected_dx,
        expected_du
    ) == RP_KERR_STATUS_OK);

    for (index = 0; index < RP_KERR_DIM; ++index) {
        coordinate_dx[index] = coordinates[index];
        velocity_du[index] = four_velocity[index];
    }
    CHECK(rp_kerr_geodesic_rhs(
        1.0,
        0.5,
        coordinate_dx,
        velocity_du,
        coordinate_dx,
        velocity_du
    ) == RP_KERR_STATUS_OK);

    for (index = 0; index < RP_KERR_DIM; ++index) {
        CHECK_CLOSE(coordinate_dx[index], expected_dx[index], 0.0, 0.0);
        CHECK_CLOSE(velocity_du[index], expected_du[index], 2e-14, 2e-14);
    }

    for (index = 0; index < RP_KERR_DIM; ++index) {
        null_coordinate_dx[index] = coordinates[index];
    }
    CHECK(rp_kerr_null_tangent(
        1.0, 0.5, coordinates, 1.0, 2.0, 3.0, -1, 1,
        expected_null_tangent
    ) == RP_KERR_STATUS_OK);
    CHECK(rp_kerr_null_tangent(
        1.0, 0.5, null_coordinate_dx, 1.0, 2.0, 3.0, -1, 1,
        null_coordinate_dx
    ) == RP_KERR_STATUS_OK);
    for (index = 0; index < RP_KERR_DIM; ++index) {
        CHECK_CLOSE(
            null_coordinate_dx[index], expected_null_tangent[index], 0.0, 0.0
        );
        null_coordinate_dx[index] = coordinates[index];
        null_tangent_dk[index] = expected_null_tangent[index];
    }
    CHECK(rp_kerr_null_geodesic_rhs(
        1.0,
        0.5,
        coordinates,
        expected_null_tangent,
        expected_null_dx,
        expected_null_dk
    ) == RP_KERR_STATUS_OK);
    CHECK(rp_kerr_null_geodesic_rhs(
        1.0,
        0.5,
        null_coordinate_dx,
        null_tangent_dk,
        null_coordinate_dx,
        null_tangent_dk
    ) == RP_KERR_STATUS_OK);
    for (index = 0; index < RP_KERR_DIM; ++index) {
        CHECK_CLOSE(
            null_coordinate_dx[index], expected_null_dx[index], 0.0, 0.0
        );
        CHECK_CLOSE(
            null_tangent_dk[index], expected_null_dk[index], 2e-14, 2e-14
        );
    }
}

static void test_polar_axis_rhs_domain(void)
{
    const double pi = acos(-1.0);
    const double epsilon = 64.0 * DBL_EPSILON;
    const double invalid[] = {
        -0.1, 0.0, epsilon, pi - epsilon, pi, pi + 0.1
    };
    const double tangent[RP_KERR_DIM] = {1.0, 0.0, 0.0, 0.0};
    double coordinates[RP_KERR_DIM] = {0.0, 8.0, 1.0, 0.0};
    double dx[RP_KERR_DIM];
    double du[RP_KERR_DIM];
    size_t index;

    for (index = 0U; index < sizeof(invalid) / sizeof(invalid[0]); ++index) {
        coordinates[2] = invalid[index];
        fill_outputs(dx, du);
        CHECK(rp_kerr_null_geodesic_rhs(
            1.0, 0.5, coordinates, tangent, dx, du
        ) == RP_KERR_STATUS_COORDINATE_SINGULARITY);
        check_zero_vector(dx);
        check_zero_vector(du);
    }
    coordinates[2] = 2.0 * epsilon;
    CHECK(rp_kerr_null_geodesic_rhs(
        1.0, 0.5, coordinates, tangent, dx, du
    ) == RP_KERR_STATUS_OK);
    coordinates[2] = pi - 2.0 * epsilon;
    CHECK(rp_kerr_null_geodesic_rhs(
        1.0, 0.5, coordinates, tangent, dx, du
    ) == RP_KERR_STATUS_OK);
}

int main(void)
{
    test_null_tangent_from_constants();
    test_schwarzschild_radial_null_tangent();
    test_equatorial_kerr_null_tangent();
    test_null_rhs_matches_numerical_metric_connection();
    test_optimized_rhs_matches_christoffel_contraction();
    test_invalid_inputs_and_zeroed_outputs();
    test_invalid_null_inputs_and_zeroed_outputs();
    test_polar_axis_rhs_domain();
    test_input_output_overlap();

    if (failures != 0) {
        (void)fprintf(stderr, "%d Kerr geodesic assertion(s) failed\n", failures);
        return 1;
    }
    (void)puts("Kerr geodesic tests passed");
    return 0;
}
