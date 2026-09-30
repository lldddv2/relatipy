#include "geodesic/initial/convert.h"
#include "geodesic/solution/reconstruct.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stddef.h>

static int close_to(double actual, double expected, double tolerance)
{
    return fabs(actual - expected) <= tolerance * (1.0 + fabs(expected));
}

static void assert_zero(const double *values, size_t count)
{
    size_t index;
    for (index = 0U; index < count; ++index) {
        assert(values[index] == 0.0);
    }
}

static void assert_normalized(double spin, const double canonical[8])
{
    double metric[RP_KERR_DIM][RP_KERR_DIM];
    double norm = 0.0;
    size_t mu;
    size_t nu;

    assert(rp_kerr_metric(1.0, spin, canonical, metric)
        == RP_KERR_STATUS_OK);
    for (mu = 0U; mu < RP_KERR_DIM; ++mu) {
        for (nu = 0U; nu < RP_KERR_DIM; ++nu) {
            norm += metric[mu][nu]
                * canonical[4U + mu] * canonical[4U + nu];
        }
    }
    assert(close_to(norm, -1.0, 5e-13));
}

static void test_schwarzschild_cartesian(void)
{
    const double input[7] = {
        3.0, 10.0, 0.0, 0.0, 0.0, 0.31622776601683794, 0.0
    };
    double canonical[8];
    double back[7];
    size_t index;

    assert(rp_initial_cartesian_to_canonical(
        0.0, input, canonical
    ) == RP_KERR_STATUS_OK);
    assert(close_to(canonical[0], 3.0, 1e-14));
    assert(close_to(canonical[1], 10.0, 1e-14));
    assert(close_to(canonical[2], acos(-1.0) / 2.0, 1e-14));
    assert(close_to(canonical[3], 0.0, 1e-14));
    assert(canonical[4] > 1.0);
    assert_normalized(0.0, canonical);
    assert(rp_initial_canonical_to_cartesian(
        0.0, canonical, back
    ) == RP_KERR_STATUS_OK);
    for (index = 0U; index < 7U; ++index) {
        assert(close_to(back[index], input[index], 5e-14));
    }
}

static void test_spinning_family_agreement(void)
{
    const double bl[7] = {2.0, 8.0, 1.2, 0.7, 0.012, -0.001, 0.008};
    double canonical[8];
    double from_cart[8];
    double from_sph[8];
    double cart[7];
    double spherical[7];
    double radial;
    double transverse;
    double transverse_velocity;
    size_t index;

    assert(rp_initial_from_bl(
        1.0, bl, canonical
    ) == RP_KERR_STATUS_OK);
    assert_normalized(1.0, canonical);
    assert(rp_initial_canonical_to_cartesian(
        1.0, canonical, cart
    ) == RP_KERR_STATUS_OK);
    assert(rp_initial_cartesian_to_canonical(
        1.0, cart, from_cart
    ) == RP_KERR_STATUS_OK);
    for (index = 0U; index < 8U; ++index) {
        assert(close_to(from_cart[index], canonical[index], 3e-13));
    }

    transverse = hypot(cart[1], cart[2]);
    radial = hypot(transverse, cart[3]);
    transverse_velocity = (cart[1] * cart[4] + cart[2] * cart[5])
        / transverse;
    spherical[0] = cart[0];
    spherical[1] = radial;
    spherical[2] = atan2(transverse, cart[3]);
    spherical[3] = atan2(cart[2], cart[1]);
    spherical[4] = (transverse * transverse_velocity
        + cart[3] * cart[6]) / radial;
    spherical[5] = (cart[3] * transverse_velocity
        - transverse * cart[6]) / (radial * radial);
    spherical[6] = (cart[1] * cart[5] - cart[2] * cart[4])
        / (transverse * transverse);
    assert(rp_initial_from_spherical(
        1.0, spherical, from_sph
    ) == RP_KERR_STATUS_OK);
    for (index = 0U; index < 8U; ++index) {
        assert(close_to(from_sph[index], canonical[index], 3e-13));
    }
}

static void test_elements_and_retrograde(void)
{
    const double elliptic[7] = {0.0, 12.0, 0.25, 0.6, 0.3, 0.5, 0.4};
    const double retrograde[7] = {
        0.0, 12.0, 0.0, 3.14159265358979323846, 0.0, 0.0, 0.0
    };
    const double hyperbolic[7] = {0.0, -10.0, 1.5, 0.0, 0.0, 0.0, 0.0};
    double canonical[8];
    double cart[7];

    assert(rp_initial_from_elements(
        0.5, elliptic, canonical
    ) == RP_KERR_STATUS_OK);
    assert_normalized(0.5, canonical);
    assert(rp_initial_from_elements(
        0.0, retrograde, canonical
    ) == RP_KERR_STATUS_OK);
    assert_normalized(0.0, canonical);
    assert(rp_initial_canonical_to_cartesian(
        0.0, canonical, cart
    ) == RP_KERR_STATUS_OK);
    assert(cart[5] < 0.0);
    assert(rp_initial_from_elements(
        0.0, hyperbolic, canonical
    ) == RP_KERR_STATUS_OK);
    assert_normalized(0.0, canonical);
}

static void test_errors(void)
{
    double cart[7] = {0.0, 10.0, 0.0, 0.0, 0.0, 0.1, 0.0};
    double bl[7] = {0.0, 8.0, 1.2, 0.0, 0.0, 0.0, 0.0};
    double elements[7] = {0.0, 10.0, 1.0, 0.0, 0.0, 0.0, 0.0};
    double canonical[8];
    size_t index;

    for (index = 0U; index < 8U; ++index) {
        canonical[index] = 9.0;
    }
    assert(rp_initial_cartesian_to_canonical(
        0.0, NULL, canonical
    ) == RP_KERR_STATUS_NULL_POINTER);
    assert_zero(canonical, 8U);
    assert(rp_initial_cartesian_to_canonical(
        NAN, cart, canonical
    ) == RP_KERR_STATUS_NONFINITE_INPUT);
    assert_zero(canonical, 8U);
    cart[1] = 2.0;
    assert(rp_initial_cartesian_to_canonical(
        0.0, cart, canonical
    ) == RP_KERR_STATUS_INVALID_PARAMETER);
    assert_zero(canonical, 8U);
    cart[1] = 10.0;
    cart[2] = 0.0;
    cart[3] = 10.0;
    cart[1] = 0.0;
    assert(rp_initial_cartesian_to_canonical(
        0.0, cart, canonical
    ) == RP_KERR_STATUS_COORDINATE_SINGULARITY);
    assert_zero(canonical, 8U);
    cart[1] = 10.0;
    cart[3] = 0.0;
    cart[5] = 2.0;
    assert(rp_initial_cartesian_to_canonical(
        0.0, cart, canonical
    ) == RP_KERR_STATUS_NON_TIMELIKE_VELOCITY);
    assert_zero(canonical, 8U);
    bl[1] = 2.0;
    assert(rp_initial_from_bl(
        0.0, bl, canonical
    ) == RP_KERR_STATUS_INVALID_PARAMETER);
    assert_zero(canonical, 8U);
    assert(rp_initial_from_elements(
        0.0, elements, canonical
    ) == RP_KERR_STATUS_INVALID_PARAMETER);
    assert_zero(canonical, 8U);
    elements[2] = 0.0;
    elements[1] = -10.0;
    assert(rp_initial_from_elements(
        0.0, elements, canonical
    ) == RP_KERR_STATUS_INVALID_PARAMETER);
    assert_zero(canonical, 8U);
    bl[1] = 8.0;
    bl[2] = 0.0;
    assert(rp_initial_from_bl(
        0.0, bl, canonical
    ) == RP_KERR_STATUS_COORDINATE_SINGULARITY);
    assert_zero(canonical, 8U);
    canonical[1] = 8.0;
    canonical[2] = 0.0;
    canonical[4] = 1.0;
    assert(rp_initial_canonical_to_cartesian(
        0.0, canonical, cart
    ) == RP_KERR_STATUS_COORDINATE_SINGULARITY);
    assert_zero(cart, 7U);
    bl[1] = 1.0;
    bl[2] = 1.2;
    assert(rp_initial_from_bl(
        1.0, bl, canonical
    ) == RP_KERR_STATUS_INVALID_PARAMETER);
    assert_zero(canonical, 8U);
}


static void test_elements_reconstruct_roundtrip(void)
{
    const double elements[7] = {0.0, 12.0, 0.25, 0.6, 0.3, 0.5, 0.4};
    const double hyperbolic[7] = {0.0, -10.0, 1.5, 0.2, 0.3, 0.4, 0.0};
    double canonical[8];
    double cartesian[7];
    double row[RP_SOLUTION_RECONSTRUCTED_DIM];
    rp_kerr_status row_status;

    assert(rp_initial_from_elements(
        0.5, elements, canonical
    ) == RP_KERR_STATUS_OK);
    assert(rp_initial_canonical_to_cartesian(
        0.5, canonical, cartesian
    ) == RP_KERR_STATUS_OK);
    assert(rp_solution_reconstruct_batch(
        0.5, cartesian, 1U, row, &row_status
    ) == RP_KERR_STATUS_OK);
    assert(row_status == RP_KERR_STATUS_OK);
    assert(close_to(row[RP_SOL_SEMIMAJOR], elements[1], 4e-13));
    assert(close_to(row[RP_SOL_ECCENTRICITY], elements[2], 4e-13));
    assert(close_to(row[RP_SOL_INCLINATION], elements[3], 4e-13));
    assert(close_to(row[RP_SOL_ASCENDING_NODE], elements[4], 4e-13));
    assert(close_to(row[RP_SOL_PERIAPSIS_ARGUMENT], elements[5], 4e-13));
    assert(close_to(row[RP_SOL_TRUE_ANOMALY], elements[6], 4e-13));
    assert(rp_initial_from_elements(
        0.0, hyperbolic, canonical
    ) == RP_KERR_STATUS_OK);
    assert(rp_initial_canonical_to_cartesian(
        0.0, canonical, cartesian
    ) == RP_KERR_STATUS_OK);
    assert(rp_solution_reconstruct_batch(
        0.0, cartesian, 1U, row, &row_status
    ) == RP_KERR_STATUS_OK);
    assert(close_to(row[RP_SOL_SEMIMAJOR], hyperbolic[1], 4e-13));
    assert(close_to(row[RP_SOL_ECCENTRICITY], hyperbolic[2], 4e-13));
}

static void test_polar_axis_conversions(void)
{
    const double pi = acos(-1.0);
    const double epsilon = 64.0 * DBL_EPSILON;
    const double rejected[] = {
        -0.1, 0.0, epsilon, pi - epsilon, pi, pi + 0.1
    };
    double bl[7] = {0.0, 8.0, 1.0, 0.0, 0.0, 0.0, 0.0};
    double spherical[7] = {0.0, 8.0, 1.0, 0.0, 0.0, 0.0, 0.0};
    double canonical[8] = {0.0, 8.0, 1.0, 0.0, 1.2, 0.0, 0.0, 0.0};
    double cartesian[7] = {0.0, 8.0, 0.0, 8.0, 0.0, 0.0, 0.0};
    double output[8];
    double back[7];
    size_t index;

    for (index = 0U; index < sizeof(rejected) / sizeof(rejected[0]); ++index) {
        bl[2] = rejected[index];
        spherical[2] = rejected[index];
        canonical[2] = rejected[index];
        assert(rp_initial_from_bl(0.0, bl, output)
            == RP_KERR_STATUS_COORDINATE_SINGULARITY);
        assert_zero(output, 8U);
        assert(rp_initial_from_spherical(0.0, spherical, output)
            == RP_KERR_STATUS_COORDINATE_SINGULARITY);
        assert_zero(output, 8U);
        assert(rp_initial_canonical_to_cartesian(0.0, canonical, back)
            == RP_KERR_STATUS_COORDINATE_SINGULARITY);
        assert_zero(back, 7U);
    }

    bl[2] = 2.0 * epsilon;
    spherical[2] = bl[2];
    canonical[2] = bl[2];
    assert(rp_initial_from_bl(0.0, bl, output) == RP_KERR_STATUS_OK);
    assert(rp_initial_from_spherical(0.0, spherical, output)
        == RP_KERR_STATUS_OK);
    assert(rp_initial_canonical_to_cartesian(0.0, canonical, back)
        == RP_KERR_STATUS_OK);
    bl[2] = pi - 2.0 * epsilon;
    assert(rp_initial_from_bl(0.0, bl, output) == RP_KERR_STATUS_OK);

    cartesian[1] = 8.0 * (epsilon / 2.0);
    assert(rp_initial_cartesian_to_canonical(0.0, cartesian, output)
        == RP_KERR_STATUS_COORDINATE_SINGULARITY);
    assert_zero(output, 8U);
    cartesian[1] = 8.0 * (2.0 * epsilon);
    assert(rp_initial_cartesian_to_canonical(0.0, cartesian, output)
        == RP_KERR_STATUS_OK);
}

int main(void)
{
    test_schwarzschild_cartesian();
    test_spinning_family_agreement();
    test_elements_and_retrograde();
    test_elements_reconstruct_roundtrip();
    test_errors();
    test_polar_axis_conversions();
    return 0;
}
