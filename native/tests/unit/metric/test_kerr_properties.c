#include "metric/kerr_properties.h"

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

#define CHECK_CLOSE(actual, expected, tolerance)                                \
    do {                                                                        \
        const double actual_value = (actual);                                   \
        const double expected_value = (expected);                               \
        if (!(fabs(actual_value - expected_value) <= (tolerance))) {            \
            (void)fprintf(stderr, "FAIL %s:%d: %.17g != %.17g\n", __FILE__,  \
                __LINE__, actual_value, expected_value);                        \
            ++failures;                                                         \
        }                                                                       \
    } while (0)

#define PI 3.14159265358979323846
#define COUNT 5U

/* Reference: r_+ = 1 + sqrt(1 - a^2), r_E(theta) = 1 + sqrt(1 - a^2 cos^2),
 * rho = hypot(r, a) sin(theta), z = r cos(theta); units GM/c^2. */

static void test_schwarzschild_horizon_is_sphere(void)
{
    const double theta[COUNT] = {0.0, PI / 6.0, PI / 2.0, 2.0, PI};
    double rho[COUNT];
    double z[COUNT];
    size_t i;

    CHECK(rp_kerr_surface_profile(
        0.0, RP_KERR_SURFACE_OUTER_HORIZON, theta, COUNT, rho, z
    ) == RP_KERR_STATUS_OK);
    for (i = 0U; i < COUNT; ++i) {
        CHECK_CLOSE(hypot(rho[i], z[i]), 2.0, 1e-14);
        CHECK_CLOSE(rho[i], 2.0 * sin(theta[i]), 1e-14);
        CHECK_CLOSE(z[i], 2.0 * cos(theta[i]), 1e-14);
    }
}

static void test_poles_and_equator(void)
{
    const double spin = 0.9;
    const double horizon = 1.0 + sqrt(1.0 - spin * spin);
    const double theta[3] = {0.0, PI / 2.0, PI};
    double hr[3];
    double hz[3];
    double er[3];
    double ez[3];

    CHECK(rp_kerr_surface_profile(
        spin, RP_KERR_SURFACE_OUTER_HORIZON, theta, 3U, hr, hz
    ) == RP_KERR_STATUS_OK);
    CHECK(rp_kerr_surface_profile(
        spin, RP_KERR_SURFACE_ERGOSURFACE, theta, 3U, er, ez
    ) == RP_KERR_STATUS_OK);

    /* Poles: both surfaces meet on the axis at z = +/- r_+. */
    CHECK_CLOSE(hr[0], 0.0, 1e-15);
    CHECK_CLOSE(hz[0], horizon, 1e-14);
    CHECK_CLOSE(hz[2], -horizon, 1e-14);
    CHECK_CLOSE(er[0], hr[0], 1e-15);
    CHECK_CLOSE(ez[0], hz[0], 1e-14);
    CHECK_CLOSE(er[2], hr[2], 1e-15);
    CHECK_CLOSE(ez[2], hz[2], 1e-14);

    /* Equator: horizon at hypot(r_+, a); ergosurface at hypot(2, a). */
    CHECK_CLOSE(hr[1], hypot(horizon, spin), 1e-14);
    CHECK_CLOSE(hz[1], 0.0, 1e-15);
    CHECK_CLOSE(er[1], hypot(2.0, spin), 1e-14);
    CHECK_CLOSE(ez[1], 0.0, 1e-15);
    CHECK(er[1] > hr[1]);
}

static void test_ergosurface_matches_radii(void)
{
    const double spin = 0.6;
    const double theta[COUNT] = {0.1, 0.7, 1.2, 2.2, 3.0};
    double radii[COUNT];
    double rho[COUNT];
    double z[COUNT];
    size_t i;

    CHECK(rp_kerr_ergosurface_radii(spin, theta, COUNT, radii)
        == RP_KERR_STATUS_OK);
    CHECK(rp_kerr_surface_profile(
        spin, RP_KERR_SURFACE_ERGOSURFACE, theta, COUNT, rho, z
    ) == RP_KERR_STATUS_OK);
    for (i = 0U; i < COUNT; ++i) {
        CHECK_CLOSE(rho[i], hypot(radii[i], spin) * sin(theta[i]), 1e-14);
        CHECK_CLOSE(z[i], radii[i] * cos(theta[i]), 1e-14);
    }
}

static void test_invalid_inputs(void)
{
    const double good[2] = {0.0, 1.0};
    const double negative[2] = {0.5, -1e-3};
    const double beyond[2] = {0.5, PI + 1e-9};
    double not_finite[2];
    double rho[2] = {-7.0, -7.0};
    double z[2] = {-7.0, -7.0};

    not_finite[0] = 0.5;
    not_finite[1] = NAN;

    CHECK(rp_kerr_surface_profile(
        0.5, RP_KERR_SURFACE_ERGOSURFACE, NULL, 0U, NULL, NULL
    ) == RP_KERR_STATUS_OK);
    CHECK(rp_kerr_surface_profile(
        0.5, RP_KERR_SURFACE_ERGOSURFACE, NULL, 2U, rho, z
    ) == RP_KERR_STATUS_NULL_POINTER);
    CHECK(rp_kerr_surface_profile(
        0.5, RP_KERR_SURFACE_ERGOSURFACE, good, 2U, NULL, z
    ) == RP_KERR_STATUS_NULL_POINTER);
    CHECK(rp_kerr_surface_profile(
        0.5, RP_KERR_SURFACE_ERGOSURFACE, good, 2U, rho, NULL
    ) == RP_KERR_STATUS_NULL_POINTER);
    CHECK(rp_kerr_surface_profile(
        1.5, RP_KERR_SURFACE_OUTER_HORIZON, good, 2U, rho, z
    ) == RP_KERR_STATUS_INVALID_PARAMETER);
    CHECK(rp_kerr_surface_profile(
        -0.1, RP_KERR_SURFACE_OUTER_HORIZON, good, 2U, rho, z
    ) == RP_KERR_STATUS_INVALID_PARAMETER);
    CHECK(rp_kerr_surface_profile(
        NAN, RP_KERR_SURFACE_OUTER_HORIZON, good, 2U, rho, z
    ) == RP_KERR_STATUS_INVALID_PARAMETER);
    CHECK(rp_kerr_surface_profile(
        0.5, (rp_kerr_surface)7, good, 2U, rho, z
    ) == RP_KERR_STATUS_INVALID_PARAMETER);
    CHECK(rp_kerr_surface_profile(
        0.5, RP_KERR_SURFACE_ERGOSURFACE, negative, 2U, rho, z
    ) == RP_KERR_STATUS_INVALID_PARAMETER);
    CHECK(rp_kerr_surface_profile(
        0.5, RP_KERR_SURFACE_OUTER_HORIZON, beyond, 2U, rho, z
    ) == RP_KERR_STATUS_INVALID_PARAMETER);
    CHECK(rp_kerr_surface_profile(
        0.5, RP_KERR_SURFACE_ERGOSURFACE, not_finite, 2U, rho, z
    ) == RP_KERR_STATUS_INVALID_PARAMETER);
    /* Angles are validated before any output is written. */
    CHECK(rho[0] == -7.0 && z[0] == -7.0);
}

static void test_extremal_spin(void)
{
    const double theta[1] = {PI / 2.0};
    double rho[1];
    double z[1];

    CHECK(rp_kerr_surface_profile(
        1.0, RP_KERR_SURFACE_OUTER_HORIZON, theta, 1U, rho, z
    ) == RP_KERR_STATUS_OK);
    CHECK_CLOSE(rho[0], sqrt(2.0), 1e-14);
    CHECK_CLOSE(z[0], 0.0, 1e-15);
}

int main(void)
{
    test_schwarzschild_horizon_is_sphere();
    test_poles_and_equator();
    test_ergosurface_matches_radii();
    test_invalid_inputs();
    test_extremal_spin();

    if (failures != 0) {
        (void)fprintf(stderr, "%d Kerr properties assertion(s) failed\n", failures);
        return 1;
    }
    (void)puts("Kerr properties tests passed");
    return 0;
}
