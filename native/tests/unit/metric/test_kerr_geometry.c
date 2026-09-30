#include "relatipy/kerr_geometry.h"

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

static void test_schwarzschild_metric_and_inverse(void)
{
    const double coordinates[RP_KERR_DIM] = {0.0, 10.0, 1.1, 0.2};
    double metric[RP_KERR_DIM][RP_KERR_DIM];
    double inverse[RP_KERR_DIM][RP_KERR_DIM];
    size_t mu;
    size_t nu;
    size_t sigma;
    double product;

    CHECK(rp_kerr_metric(1.0, 0.0, coordinates, metric) == RP_KERR_STATUS_OK);
    CHECK_CLOSE(metric[0][0], -0.8, 1e-14, 1e-14);
    CHECK_CLOSE(metric[1][1], 1.25, 1e-14, 1e-14);
    CHECK_CLOSE(metric[2][2], 100.0, 1e-14, 1e-14);
    CHECK_CLOSE(metric[3][3], 100.0 * sin(1.1) * sin(1.1), 1e-14, 1e-14);

    CHECK(rp_kerr_inverse_metric(1.0, 0.0, coordinates, inverse)
        == RP_KERR_STATUS_OK);
    for (mu = 0; mu < RP_KERR_DIM; ++mu) {
        for (nu = 0; nu < RP_KERR_DIM; ++nu) {
            product = 0.0;
            for (sigma = 0; sigma < RP_KERR_DIM; ++sigma) {
                product += metric[mu][sigma] * inverse[sigma][nu];
            }
            CHECK_CLOSE(product, mu == nu ? 1.0 : 0.0, 2e-14, 2e-14);
        }
    }
}

static void test_schwarzschild_christoffel(void)
{
    const double radius = 10.0;
    const double theta = 1.1;
    const double coordinates[RP_KERR_DIM] = {0.0, radius, theta, 0.2};
    double gamma[RP_KERR_DIM][RP_KERR_DIM][RP_KERR_DIM];
    size_t lambda;
    size_t mu;
    size_t nu;

    CHECK(rp_kerr_christoffel(1.0, 0.0, coordinates, gamma)
        == RP_KERR_STATUS_OK);
    CHECK_CLOSE(gamma[0][0][1], 1.0 / (radius * (radius - 2.0)), 2e-14, 2e-14);
    CHECK_CLOSE(gamma[0][1][0], gamma[0][0][1], 0.0, 0.0);
    CHECK_CLOSE(gamma[1][0][0], (radius - 2.0) / (radius * radius * radius),
        2e-14, 2e-14);
    CHECK_CLOSE(gamma[1][1][1], -1.0 / (radius * (radius - 2.0)),
        2e-14, 2e-14);
    CHECK_CLOSE(gamma[1][2][2], -(radius - 2.0), 2e-14, 2e-14);
    CHECK_CLOSE(gamma[1][3][3], -(radius - 2.0) * sin(theta) * sin(theta),
        2e-14, 2e-14);
    CHECK_CLOSE(gamma[2][1][2], 1.0 / radius, 2e-14, 2e-14);
    CHECK_CLOSE(gamma[2][3][3], -sin(theta) * cos(theta), 2e-14, 2e-14);
    CHECK_CLOSE(gamma[3][1][3], 1.0 / radius, 2e-14, 2e-14);
    CHECK_CLOSE(gamma[3][2][3], cos(theta) / sin(theta), 2e-14, 2e-14);

    for (lambda = 0; lambda < RP_KERR_DIM; ++lambda) {
        for (mu = 0; mu < RP_KERR_DIM; ++mu) {
            for (nu = 0; nu < RP_KERR_DIM; ++nu) {
                CHECK_CLOSE(gamma[lambda][mu][nu], gamma[lambda][nu][mu],
                    0.0, 0.0);
            }
        }
    }
}

static void test_kerr_symmetries_and_inverse(void)
{
    const double masses[3] = {1.0, 2.0, 0.75};
    const double spins[3] = {0.5, -1.7, 0.0};
    const double points[3][RP_KERR_DIM] = {
        {0.0, 8.0, 1.1, 0.3},
        {4.0, 12.0, 0.7, -2.0},
        {-1.0, 5.0, 2.2, 1.5}
    };
    double metric[RP_KERR_DIM][RP_KERR_DIM];
    double inverse[RP_KERR_DIM][RP_KERR_DIM];
    double gamma[RP_KERR_DIM][RP_KERR_DIM][RP_KERR_DIM];
    size_t point;
    size_t lambda;
    size_t mu;
    size_t nu;
    size_t sigma;
    double product;

    for (point = 0; point < 3; ++point) {
        CHECK(rp_kerr_metric(masses[point], spins[point], points[point], metric)
            == RP_KERR_STATUS_OK);
        CHECK(rp_kerr_inverse_metric(
            masses[point], spins[point], points[point], inverse)
            == RP_KERR_STATUS_OK);
        CHECK(rp_kerr_christoffel(masses[point], spins[point], points[point], gamma)
            == RP_KERR_STATUS_OK);
        for (mu = 0; mu < RP_KERR_DIM; ++mu) {
            for (nu = 0; nu < RP_KERR_DIM; ++nu) {
                CHECK_CLOSE(metric[mu][nu], metric[nu][mu], 0.0, 0.0);
                CHECK_CLOSE(inverse[mu][nu], inverse[nu][mu], 0.0, 0.0);
                product = 0.0;
                for (sigma = 0; sigma < RP_KERR_DIM; ++sigma) {
                    product += metric[mu][sigma] * inverse[sigma][nu];
                }
                CHECK_CLOSE(product, mu == nu ? 1.0 : 0.0, 3e-14, 3e-14);
            }
        }
        for (lambda = 0; lambda < RP_KERR_DIM; ++lambda) {
            for (mu = 0; mu < RP_KERR_DIM; ++mu) {
                for (nu = 0; nu < RP_KERR_DIM; ++nu) {
                    CHECK_CLOSE(gamma[lambda][mu][nu], gamma[lambda][nu][mu],
                        0.0, 0.0);
                }
            }
        }
    }
}

static void test_literal_legacy_christoffel_regression(void)
{
    /*
     * Frozen from Kerr._get_christoffel_symbols in the 2026-09-01 legacy
     * snapshot.  The C API receives the legacy internal spin length, so the
     * Python constructor used dimensionless spin a/M in each case.  These are
     * literal compatibility values: this test intentionally does not require
     * them to be the Levi-Civita connection of the separately tested metric.
     */
    const double masses[3] = {1.0, 2.0, 0.75};
    const double spins[3] = {0.5, -1.4, 0.675};
    const double coordinates[3][RP_KERR_DIM] = {
        {0.0, 8.0, 1.1, 0.3},
        {0.2, 15.0, 0.8, -0.4},
        {-1.0, 6.0, 1.35, 2.0}
    };
    const double expected[3][RP_KERR_DIM][RP_KERR_DIM][RP_KERR_DIM] = {
        {
            {{0.0, 0.020756247543677078, -0.00039413983281485126, 0.0},
             {0.020756247543677078, 0.0, 0.0, -0.02469076273038043},
             {-0.00039413983281485126, 0.0, 0.0, 0.00015652289119530666},
             {0.0, -0.02469076273038043, 0.00015652289119530666, 0.0}},
            {{0.011620304889560888, 0.0, 0.0, -0.004614716824978808},
             {0.0, -0.02169724155453523, -0.001577826425797284, 0.0},
             {0.0, -0.001577826425797284, -5.9639567157709426, 0.0},
             {-0.004614716824978808, 0.0, 0.0, -4.735043332424259}},
            {{-6.153489274526096e-06, 0.0, 0.0, 0.0007907233717766035},
             {0.0, 3.304348535701119e-05, 0.12489961708420821, 0.0},
             {0.0, 0.12489961708420821, -0.001577826425797284, 0.0},
             {0.0007907233717766035, 0.0, 0.0, -0.40612845345216936}},
            {{0.0, 0.0001615272182387321, -0.000992482355935277, 0.0},
             {0.0001615272182387321, 0.0, 0.0, 0.12432147266357679},
             {-0.000992482355935277, 0.0, 0.0, 0.5093622450718792},
             {0.0, 0.12432147266357679, 0.5093622450718792, 0.0}}
        },
        {
            {{0.0, 0.0119310628525823, -0.0011512299811903667, 0.0},
             {0.0119310628525823, 0.0, 0.0, 0.025783109802546084},
             {-0.0011512299811903667, 0.0, 0.0, -0.000829391742690033},
             {0.0, 0.025783109802546084, -0.000829391742690033, 0.0}},
            {{0.00633317380482608, 0.0, 0.0, 0.0045626696182046265},
             {0.0, -0.013349072449853517, -0.004335366801519993, 0.0},
             {0.0, -0.004335366801519993, -10.823567227775882, 0.0},
             {0.0045626696182046265, 0.0, 0.0, -5.5665179820373805}},
            {{-5.09503397777441e-06, 0.0, 0.0, -0.0008259777939969145},
             {0.0, 2.6590816986751676e-05, 0.06638596189754589, 0.0},
             {0.0, 0.06638596189754589, -0.004335366801519993, 0.0},
             {-0.0008259777939969145, 0.0, 0.0, -0.5032052700803743}},
            {{0.0, -7.359661611568215e-05, 0.001597954743669163, 0.0},
             {-7.359661611568215e-05, 0.0, 0.0, 0.06593189833572628},
             {0.001597954743669163, 0.0, 0.0, 0.9723658306316648},
             {0.0, 0.06593189833572628, 0.9723658306316648, 0.0}}
        },
        {
            {{0.0, 0.027612208840615205, -0.0006753081944310627, 0.0},
             {0.027612208840615205, 0.0, 0.0, -0.052831868149926374},
             {-0.0006753081944310627, 0.0, 0.0, 0.0004339694880985594},
             {0.0, -0.052831868149926374, 0.0004339694880985594, 0.0}},
            {{0.01532407904006055, 0.0, 0.0, -0.009847626300758312},
             {0.0, -0.031216461830180325, -0.0027028725434599314, 0.0},
             {0.0, -0.0027028725434599314, -4.421378531006508, 0.0},
             {-0.009847626300758312, 0.0, 0.0, -4.202983520670714}},
            {{-1.874718060273453e-05, 0.0, 0.0, 0.0010125039790526878},
             {0.0, 0.00010182468200739069, 0.16656555413365384, 0.0},
             {0.0, 0.16656555413365384, -0.0027028725434599314, 0.0},
             {0.0010125039790526878, 0.0, 0.0, -0.21755674974308947}},
            {{0.0, 0.0005112583028658887, -0.0010508599566847194, 0.0},
             {0.0005112583028658887, 0.0, 0.0, 0.16360543781649062},
             {-0.0010508599566847194, 0.0, 0.0, 0.22513102644315708},
             {0.0, 0.16360543781649062, 0.22513102644315708, 0.0}}
        }
    };
    double gamma[RP_KERR_DIM][RP_KERR_DIM][RP_KERR_DIM];
    size_t point;
    size_t lambda;
    size_t mu;
    size_t nu;

    for (point = 0; point < 3; ++point) {
        CHECK(rp_kerr_christoffel(
            masses[point], spins[point], coordinates[point], gamma)
            == RP_KERR_STATUS_OK);
        for (lambda = 0; lambda < RP_KERR_DIM; ++lambda) {
            for (mu = 0; mu < RP_KERR_DIM; ++mu) {
                for (nu = 0; nu < RP_KERR_DIM; ++nu) {
                    CHECK_CLOSE(gamma[lambda][mu][nu],
                        expected[point][lambda][mu][nu], 5e-13, 5e-14);
                }
            }
        }
    }
}

static void test_four_velocity_regression_and_normalization(void)
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
    const double expected[3][RP_KERR_DIM] = {
        {1.167733652142157, -0.01167733652142157,
            0.0023354673042843143, 0.02335467304284314},
        {1.1839872523784072, 0.017759808785676106,
            -0.0011839872523784073, -0.017759808785676106},
        {1.158013073513648, -0.02316026147027296,
            0.003474039220540944, 0.01158013073513648}
    };
    double metric[RP_KERR_DIM][RP_KERR_DIM];
    double four_velocity[RP_KERR_DIM];
    double norm;
    size_t point;
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
        CHECK(rp_kerr_metric(
            masses[point], spins[point], coordinates[point], metric)
            == RP_KERR_STATUS_OK);

        norm = 0.0;
        for (mu = 0; mu < RP_KERR_DIM; ++mu) {
            CHECK_CLOSE(four_velocity[mu], expected[point][mu], 3e-14, 3e-14);
            for (nu = 0; nu < RP_KERR_DIM; ++nu) {
                norm += metric[mu][nu]
                    * four_velocity[mu] * four_velocity[nu];
            }
        }
        CHECK(four_velocity[0] > 0.0);
        CHECK_CLOSE(norm, -1.0, 5e-14, 5e-14);
        CHECK_CLOSE(four_velocity[1] / four_velocity[0],
            coordinate_velocities[point][0], 2e-14, 2e-14);
        CHECK_CLOSE(four_velocity[2] / four_velocity[0],
            coordinate_velocities[point][1], 2e-14, 2e-14);
        CHECK_CLOSE(four_velocity[3] / four_velocity[0],
            coordinate_velocities[point][2], 2e-14, 2e-14);
    }
}

static void test_invalid_inputs_and_zeroed_outputs(void)
{
    const double regular[RP_KERR_DIM] = {0.0, 8.0, 1.1, 0.0};
    const double axis[RP_KERR_DIM] = {0.0, 8.0, 0.0, 0.0};
    const double horizon[RP_KERR_DIM] = {
        0.0, 1.0 + sqrt(0.75), 1.1, 0.0
    };
    const double nonfinite[RP_KERR_DIM] = {0.0, NAN, 1.1, 0.0};
    const double ring[RP_KERR_DIM] = {0.0, 0.0, 1.5707963267948966, 0.0};
    const double overflow[RP_KERR_DIM] = {0.0, DBL_MAX, 1.1, 0.0};
    const double regular_velocity[3] = {-0.01, 0.002, 0.02};
    const double nonfinite_velocity[3] = {-0.01, NAN, 0.02};
    const double nontimelike_velocity[3] = {0.0, 0.0, 1.0};
    const double overflow_velocity[3] = {DBL_MAX, 0.0, 0.0};
    double metric[RP_KERR_DIM][RP_KERR_DIM];
    double inverse[RP_KERR_DIM][RP_KERR_DIM];
    double gamma[RP_KERR_DIM][RP_KERR_DIM][RP_KERR_DIM];
    double four_velocity[RP_KERR_DIM];
    size_t lambda;
    size_t mu;
    size_t nu;

    CHECK(rp_kerr_metric(1.0, 0.5, NULL, metric)
        == RP_KERR_STATUS_NULL_POINTER);
    for (mu = 0; mu < RP_KERR_DIM; ++mu) {
        for (nu = 0; nu < RP_KERR_DIM; ++nu) {
            inverse[mu][nu] = 7.0;
        }
    }
    CHECK(rp_kerr_inverse_metric(1.0, 0.5, NULL, inverse)
        == RP_KERR_STATUS_NULL_POINTER);
    for (mu = 0; mu < RP_KERR_DIM; ++mu) {
        for (nu = 0; nu < RP_KERR_DIM; ++nu) {
            CHECK_CLOSE(metric[mu][nu], 0.0, 0.0, 0.0);
            CHECK_CLOSE(inverse[mu][nu], 0.0, 0.0, 0.0);
        }
    }
    CHECK(rp_kerr_metric(1.0, 0.5, regular, NULL)
        == RP_KERR_STATUS_NULL_POINTER);
    CHECK(rp_kerr_metric(0.0, 0.0, regular, metric)
        == RP_KERR_STATUS_INVALID_PARAMETER);
    CHECK(rp_kerr_metric(1.0, 1.1, regular, metric)
        == RP_KERR_STATUS_INVALID_PARAMETER);
    CHECK(rp_kerr_metric(1.0, 0.5, nonfinite, metric)
        == RP_KERR_STATUS_NONFINITE_INPUT);
    CHECK(rp_kerr_metric(1.0, 0.5, axis, metric)
        == RP_KERR_STATUS_COORDINATE_SINGULARITY);
    CHECK(rp_kerr_metric(1.0, 0.5, horizon, metric)
        == RP_KERR_STATUS_COORDINATE_SINGULARITY);
    CHECK(rp_kerr_metric(1.0, 0.5, ring, metric)
        == RP_KERR_STATUS_PHYSICAL_SINGULARITY);
    CHECK(rp_kerr_metric(1.0, 0.5, overflow, metric)
        == RP_KERR_STATUS_NUMERICAL_RANGE);
    for (mu = 0; mu < RP_KERR_DIM; ++mu) {
        for (nu = 0; nu < RP_KERR_DIM; ++nu) {
            CHECK_CLOSE(metric[mu][nu], 0.0, 0.0, 0.0);
        }
    }

    for (lambda = 0; lambda < RP_KERR_DIM; ++lambda) {
        for (mu = 0; mu < RP_KERR_DIM; ++mu) {
            for (nu = 0; nu < RP_KERR_DIM; ++nu) {
                gamma[lambda][mu][nu] = 7.0;
            }
        }
    }
    CHECK(rp_kerr_christoffel(1.0, 0.5, NULL, gamma)
        == RP_KERR_STATUS_NULL_POINTER);
    for (lambda = 0; lambda < RP_KERR_DIM; ++lambda) {
        for (mu = 0; mu < RP_KERR_DIM; ++mu) {
            for (nu = 0; nu < RP_KERR_DIM; ++nu) {
                CHECK_CLOSE(gamma[lambda][mu][nu], 0.0, 0.0, 0.0);
            }
        }
    }
    CHECK(rp_kerr_christoffel(1.0, 0.5, regular, NULL)
        == RP_KERR_STATUS_NULL_POINTER);

#define CHECK_GAMMA_FAILURE(expected_status, mass, spin, point)                 \
    do {                                                                        \
        for (lambda = 0; lambda < RP_KERR_DIM; ++lambda) {                     \
            for (mu = 0; mu < RP_KERR_DIM; ++mu) {                             \
                for (nu = 0; nu < RP_KERR_DIM; ++nu) {                         \
                    gamma[lambda][mu][nu] = 7.0;                               \
                }                                                               \
            }                                                                   \
        }                                                                       \
        CHECK(rp_kerr_christoffel((mass), (spin), (point), gamma)              \
            == (expected_status));                                              \
        for (lambda = 0; lambda < RP_KERR_DIM; ++lambda) {                     \
            for (mu = 0; mu < RP_KERR_DIM; ++mu) {                             \
                for (nu = 0; nu < RP_KERR_DIM; ++nu) {                         \
                    CHECK_CLOSE(gamma[lambda][mu][nu], 0.0, 0.0, 0.0);         \
                }                                                               \
            }                                                                   \
        }                                                                       \
    } while (0)

    CHECK_GAMMA_FAILURE(RP_KERR_STATUS_INVALID_PARAMETER, 0.0, 0.0, regular);
    CHECK_GAMMA_FAILURE(RP_KERR_STATUS_INVALID_PARAMETER, 1.0, 1.1, regular);
    CHECK_GAMMA_FAILURE(RP_KERR_STATUS_NONFINITE_INPUT, 1.0, 0.5, nonfinite);
    CHECK_GAMMA_FAILURE(
        RP_KERR_STATUS_COORDINATE_SINGULARITY, 1.0, 0.5, axis);
    CHECK_GAMMA_FAILURE(
        RP_KERR_STATUS_COORDINATE_SINGULARITY, 1.0, 0.5, horizon);
    CHECK_GAMMA_FAILURE(RP_KERR_STATUS_PHYSICAL_SINGULARITY, 1.0, 0.5, ring);
    CHECK_GAMMA_FAILURE(RP_KERR_STATUS_NUMERICAL_RANGE, 1.0, 0.5, overflow);

#undef CHECK_GAMMA_FAILURE

    CHECK(rp_kerr_four_velocity(
        1.0, 0.5, regular, regular_velocity, NULL)
        == RP_KERR_STATUS_NULL_POINTER);

#define CHECK_FOUR_VELOCITY_FAILURE(                                          \
    expected_status, mass, spin, point, velocity)                             \
    do {                                                                       \
        for (mu = 0; mu < RP_KERR_DIM; ++mu) {                                \
            four_velocity[mu] = 7.0;                                           \
        }                                                                      \
        CHECK(rp_kerr_four_velocity(                                           \
            (mass), (spin), (point), (velocity), four_velocity)                \
            == (expected_status));                                             \
        for (mu = 0; mu < RP_KERR_DIM; ++mu) {                                \
            CHECK_CLOSE(four_velocity[mu], 0.0, 0.0, 0.0);                    \
        }                                                                      \
    } while (0)

    CHECK_FOUR_VELOCITY_FAILURE(
        RP_KERR_STATUS_NULL_POINTER, 1.0, 0.5, NULL, regular_velocity);
    CHECK_FOUR_VELOCITY_FAILURE(
        RP_KERR_STATUS_NULL_POINTER, 1.0, 0.5, regular, NULL);
    CHECK_FOUR_VELOCITY_FAILURE(
        RP_KERR_STATUS_INVALID_PARAMETER, 0.0, 0.0, regular, regular_velocity);
    CHECK_FOUR_VELOCITY_FAILURE(
        RP_KERR_STATUS_NONFINITE_INPUT, 1.0, 0.5, regular, nonfinite_velocity);
    CHECK_FOUR_VELOCITY_FAILURE(
        RP_KERR_STATUS_COORDINATE_SINGULARITY,
        1.0,
        0.5,
        horizon,
        regular_velocity
    );
    CHECK_FOUR_VELOCITY_FAILURE(
        RP_KERR_STATUS_NON_TIMELIKE_VELOCITY,
        1.0,
        0.5,
        regular,
        nontimelike_velocity
    );
    CHECK_FOUR_VELOCITY_FAILURE(
        RP_KERR_STATUS_NUMERICAL_RANGE,
        1.0,
        0.5,
        regular,
        overflow_velocity
    );

#undef CHECK_FOUR_VELOCITY_FAILURE
}

static void test_polar_axis_domain(void)
{
    const double pi = acos(-1.0);
    const double epsilon = 64.0 * DBL_EPSILON;
    const double angles[] = {
        -0.1, 0.0, epsilon, pi - epsilon, pi, pi + 0.1
    };
    double coordinates[RP_KERR_DIM] = {0.0, 8.0, 1.0, 0.0};
    double metric[RP_KERR_DIM][RP_KERR_DIM];
    size_t index;

    for (index = 0U; index < sizeof(angles) / sizeof(angles[0]); ++index) {
        coordinates[2] = angles[index];
        CHECK(rp_kerr_metric(1.0, 0.5, coordinates, metric)
            == RP_KERR_STATUS_COORDINATE_SINGULARITY);
    }
    coordinates[2] = 2.0 * epsilon;
    CHECK(rp_kerr_metric(1.0, 0.5, coordinates, metric)
        == RP_KERR_STATUS_OK);
    coordinates[2] = pi - 2.0 * epsilon;
    CHECK(rp_kerr_metric(1.0, 0.5, coordinates, metric)
        == RP_KERR_STATUS_OK);
}

static void test_christoffel_overlapping_buffers(void)
{
    const double regular[RP_KERR_DIM] = {0.0, 8.0, 1.1, 0.3};
    double expected[RP_KERR_DIM][RP_KERR_DIM][RP_KERR_DIM];
    double aliased[RP_KERR_DIM * RP_KERR_DIM * RP_KERR_DIM] = {0.0};
    size_t lambda;
    size_t mu;
    size_t nu;

    CHECK(rp_kerr_christoffel(1.0, 0.5, regular, expected)
        == RP_KERR_STATUS_OK);
    for (mu = 0; mu < RP_KERR_DIM; ++mu) {
        aliased[mu] = regular[mu];
    }
    CHECK(rp_kerr_christoffel(
        1.0,
        0.5,
        aliased,
        (double (*)[RP_KERR_DIM][RP_KERR_DIM])aliased
    ) == RP_KERR_STATUS_OK);
    for (lambda = 0; lambda < RP_KERR_DIM; ++lambda) {
        for (mu = 0; mu < RP_KERR_DIM; ++mu) {
            for (nu = 0; nu < RP_KERR_DIM; ++nu) {
                const size_t flat = (lambda * RP_KERR_DIM + mu)
                    * RP_KERR_DIM + nu;
                CHECK_CLOSE(
                    aliased[flat], expected[lambda][mu][nu], 2e-14, 2e-14);
            }
        }
    }
}

static void test_four_velocity_overlapping_buffers(void)
{
    const double coordinates[RP_KERR_DIM] = {0.0, 8.0, 1.1, 0.3};
    const double coordinate_velocity[3] = {-0.01, 0.002, 0.02};
    double expected[RP_KERR_DIM];
    double coordinate_output[RP_KERR_DIM];
    double velocity_output[RP_KERR_DIM];
    size_t index;

    CHECK(rp_kerr_four_velocity(
        1.0, 0.5, coordinates, coordinate_velocity, expected)
        == RP_KERR_STATUS_OK);

    for (index = 0; index < RP_KERR_DIM; ++index) {
        coordinate_output[index] = coordinates[index];
    }
    CHECK(rp_kerr_four_velocity(
        1.0,
        0.5,
        coordinate_output,
        coordinate_velocity,
        coordinate_output
    ) == RP_KERR_STATUS_OK);

    for (index = 0; index < 3; ++index) {
        velocity_output[index] = coordinate_velocity[index];
    }
    velocity_output[3] = 0.0;
    CHECK(rp_kerr_four_velocity(
        1.0,
        0.5,
        coordinates,
        velocity_output,
        velocity_output
    ) == RP_KERR_STATUS_OK);

    for (index = 0; index < RP_KERR_DIM; ++index) {
        CHECK_CLOSE(coordinate_output[index], expected[index], 2e-14, 2e-14);
        CHECK_CLOSE(velocity_output[index], expected[index], 2e-14, 2e-14);
    }
}

int main(void)
{
    test_schwarzschild_metric_and_inverse();
    test_schwarzschild_christoffel();
    test_kerr_symmetries_and_inverse();
    test_literal_legacy_christoffel_regression();
    test_four_velocity_regression_and_normalization();
    test_invalid_inputs_and_zeroed_outputs();
    test_polar_axis_domain();
    test_christoffel_overlapping_buffers();
    test_four_velocity_overlapping_buffers();

    if (failures != 0) {
        (void)fprintf(stderr, "%d Kerr geometry assertion(s) failed\n", failures);
        return 1;
    }
    (void)puts("Kerr geometry tests passed");
    return 0;
}
