/** Contract tests for null initial states; G = c = M = 1, angles in radians. */
#include "relatipy/kerr_null.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

static int failures;

#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "%s:%d: FAIL: %s\n", __FILE__, __LINE__, #condition); \
        ++failures; \
    } \
} while (0)
#define REQUIRE(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "%s:%d: FAIL: %s\n", __FILE__, __LINE__, #condition); \
        ++failures; \
        return; \
    } \
} while (0)

/* Regular points: 128 eps allows accumulated trig/Jacobian roundoff only. */
static int close_roundoff(double actual, double expected)
{
    return isfinite(actual) && isfinite(expected)
        && fabs(actual - expected) <= 128.0 * DBL_EPSILON
            * fmax(1.0, fabs(expected));
}

static void check_zero(const double *values, size_t count)
{
    size_t i;
    for (i = 0U; i < count; ++i) {
        CHECK(values[i] == 0.0);
    }
}

/* Independent chart reference: x=hypot(r,a) sin(theta) cos(phi), etc.
 * Differentiate explicitly with respect to coordinate time. */
static void bl_to_cartesian(double spin, const double bl[7], double cart[7])
{
    const double r = bl[1];
    const double rho = hypot(r, spin);
    const double st = sin(bl[2]);
    const double ct = cos(bl[2]);
    const double sp = sin(bl[3]);
    const double cp = cos(bl[3]);
    const double rho_dot = r / rho * bl[4];
    cart[0] = bl[0];
    cart[1] = rho * st * cp;
    cart[2] = rho * st * sp;
    cart[3] = r * ct;
    cart[4] = rho_dot * st * cp + rho * ct * cp * bl[5]
        - rho * st * sp * bl[6];
    cart[5] = rho_dot * st * sp + rho * ct * sp * bl[5]
        + rho * st * cp * bl[6];
    cart[6] = bl[4] * ct - r * st * bl[5];
}

static void cartesian_to_spherical(const double cart[7], double spherical[7])
{
    const double rho = hypot(cart[1], cart[2]);
    const double r = hypot(rho, cart[3]);
    const double rho_dot = (cart[1] * cart[4] + cart[2] * cart[5]) / rho;
    spherical[0] = cart[0];
    spherical[1] = r;
    spherical[2] = atan2(rho, cart[3]);
    spherical[3] = atan2(cart[2], cart[1]);
    spherical[4] = (rho * rho_dot + cart[3] * cart[6]) / r;
    spherical[5] = (cart[3] * rho_dot - rho * cart[6]) / (r * r);
    spherical[6] = (cart[1] * cart[5] - cart[2] * cart[4]) / (rho * rho);
}

static void test_threshold(void)
{
    CHECK(rp_kerr_null_horizon_threshold(0.0) == 2.0 * (1.0 + 1.0e-6));
    CHECK(close_roundoff(rp_kerr_null_horizon_threshold(0.9),
        (1.0 + sqrt(1.0 - 0.9 * 0.9)) * (1.0 + 1.0e-6)));
    CHECK(isnan(rp_kerr_null_horizon_threshold(-0.1)));
    CHECK(isnan(rp_kerr_null_horizon_threshold(1.1)));
    CHECK(isnan(rp_kerr_null_horizon_threshold(NAN)));
}

static void test_family_conversion(void)
{
    const double bl[7] = {2.5, 10.0, 1.1, 0.4, -0.3, 0.02, 0.04};
    double cart[7];
    double converted[7];
    size_t i;
    REQUIRE(rp_kerr_null_family_to_bl(0.9,
        RP_KERR_NULL_FAMILY_BOYER_LINDQUIST, bl, converted) == RP_KERR_STATUS_OK);
    CHECK(memcmp(bl, converted, sizeof(bl)) == 0);
    REQUIRE(rp_kerr_null_family_to_bl(0.0,
        RP_KERR_NULL_FAMILY_SPHERICAL, bl, converted) == RP_KERR_STATUS_OK);
    for (i = 0U; i < 7U; ++i) {
        CHECK(close_roundoff(converted[i], bl[i]));
    }
    bl_to_cartesian(0.0, bl, cart);
    REQUIRE(rp_kerr_null_family_to_bl(0.0,
        RP_KERR_NULL_FAMILY_CARTESIAN, cart, converted) == RP_KERR_STATUS_OK);
    for (i = 0U; i < 7U; ++i) {
        CHECK(close_roundoff(converted[i], bl[i]));
    }
    /* Conversion is kinematic: a zero velocity remains admissible. */
    memcpy(cart, bl, sizeof(cart));
    cart[4] = cart[5] = cart[6] = 0.0;
    REQUIRE(rp_kerr_null_family_to_bl(0.9,
        RP_KERR_NULL_FAMILY_BOYER_LINDQUIST, cart, converted) == RP_KERR_STATUS_OK);
    CHECK(memcmp(cart, converted, sizeof(cart)) == 0);
}

static void test_direction_families(void)
{
    const double spins[2] = {0.0, 0.9};
    const double bl[7] = {2.5, 10.0, 1.1, 0.4, -0.3, 0.02, 0.04};
    const double scales[2] = {2.0, 1.0e-3};
    size_t a;
    int f;
    for (a = 0U; a < 2U; ++a) {
        double rows[3][7];
        double reference[8];
        bl_to_cartesian(spins[a], bl, rows[0]);
        cartesian_to_spherical(rows[0], rows[1]);
        memcpy(rows[2], bl, sizeof(bl));
        REQUIRE(rp_kerr_null_state_from_direction(spins[a],
            RP_KERR_NULL_FAMILY_BOYER_LINDQUIST, bl, reference) == RP_KERR_STATUS_OK);
        for (f = 0; f < 3; ++f) {
            const rp_kerr_null_family family = (rp_kerr_null_family)f;
            rp_kerr_null_invariants invariants;
            rp_kerr_status row_status;
            double converted[7];
            double state[8];
            double view[7];
            double factor;
            size_t i;
            size_t s;
            REQUIRE(rp_kerr_null_family_to_bl(spins[a], family,
                rows[f], converted) == RP_KERR_STATUS_OK);
            REQUIRE(rp_kerr_null_state_from_direction(spins[a], family,
                rows[f], state) == RP_KERR_STATUS_OK);
            CHECK(state[4] == 1.0);
            REQUIRE(rp_kerr_null_invariants_evaluate(spins[a], state,
                &invariants) == RP_KERR_STATUS_OK);
            CHECK(invariants.relative_norm < 1.0e-13);
            factor = state[5] / converted[4];
            CHECK(factor > 0.0);
            for (i = 0U; i < 3U; ++i) {
                CHECK(close_roundoff(state[5U + i], factor * converted[4U + i]));
            }
            for (i = 0U; i < 8U; ++i) {
                CHECK(close_roundoff(state[i], reference[i]));
            }
            REQUIRE(rp_kerr_null_views_batch(spins[a], family, state, 1U,
                view, &row_status) == RP_KERR_STATUS_OK);
            CHECK(row_status == RP_KERR_STATUS_OK);
            for (i = 0U; i < 4U; ++i) {
                CHECK(close_roundoff(view[i], rows[f][i]));
            }
            for (i = 4U; i < 7U; ++i) {
                CHECK(close_roundoff(view[i], factor * rows[f][i]));
            }
            for (s = 0U; s < 2U; ++s) {
                double scaled[7];
                double scaled_state[8];
                memcpy(scaled, rows[f], sizeof(scaled));
                for (i = 4U; i < 7U; ++i) {
                    scaled[i] *= scales[s];
                }
                REQUIRE(rp_kerr_null_state_from_direction(spins[a], family,
                    scaled, scaled_state) == RP_KERR_STATUS_OK);
                for (i = 0U; i < 8U; ++i) {
                    CHECK(close_roundoff(scaled_state[i], state[i]));
                }
            }
        }
    }
}

static void test_radial_direction(void)
{
    const double row[7] = {0.0, 10.0, 1.57079632679489661923, 0.0, -1.0, 0.0, 0.0};
    double state[8];
    REQUIRE(rp_kerr_null_state_from_direction(0.0,
        RP_KERR_NULL_FAMILY_BOYER_LINDQUIST, row, state) == RP_KERR_STATUS_OK);
    CHECK(state[4] == 1.0);
    CHECK(close_roundoff(state[5], -(1.0 - 2.0 / 10.0)));
    CHECK(state[6] == 0.0 && state[7] == 0.0);
}

static void expect_direction_error(double spin, const double row[7],
    rp_kerr_status expected)
{
    double state[8] = {9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0};
    const rp_kerr_status actual = rp_kerr_null_state_from_direction(spin,
        RP_KERR_NULL_FAMILY_BOYER_LINDQUIST, row, state);
    if (actual != expected) {
        fprintf(stderr, "direction status: actual=%d expected=%d\n", (int)actual, (int)expected);
    }
    CHECK(actual == expected);
    check_zero(state, 8U);
}

static void test_direction_errors(void)
{
    const double valid[7] = {0.0, 10.0, 1.1, 0.4, -0.3, 0.02, 0.04};
    double row[7];
    double state[8];
    size_t i;
    memcpy(row, valid, sizeof(row));
    row[4] = row[5] = row[6] = 0.0;
    expect_direction_error(0.0, row, RP_KERR_STATUS_INVALID_PARAMETER);
    for (i = 0U; i < 2U; ++i) {
        const double spin = i == 0U ? 0.0 : 0.9;
        memcpy(row, valid, sizeof(row));
        row[1] = rp_kerr_null_horizon_threshold(spin);
        expect_direction_error(spin, row, RP_KERR_STATUS_INVALID_PARAMETER);
        row[1] -= 1.0e-7;
        expect_direction_error(spin, row, RP_KERR_STATUS_INVALID_PARAMETER);
    }
    for (i = 0U; i < 7U; ++i) {
        memcpy(row, valid, sizeof(row));
        row[i] = NAN;
        expect_direction_error(0.9, row, RP_KERR_STATUS_NONFINITE_INPUT);
    }
    expect_direction_error(0.9, NULL, RP_KERR_STATUS_NULL_POINTER);
    CHECK(rp_kerr_null_state_from_direction(0.9,
        RP_KERR_NULL_FAMILY_BOYER_LINDQUIST, valid, NULL) == RP_KERR_STATUS_NULL_POINTER);
    memcpy(row, valid, sizeof(row));
    row[2] = 0.0;
    expect_direction_error(0.9, row, RP_KERR_STATUS_COORDINATE_SINGULARITY);
    row[2] = 32.0 * DBL_EPSILON;
    expect_direction_error(0.9, row, RP_KERR_STATUS_COORDINATE_SINGULARITY);
    expect_direction_error(-0.1, valid, RP_KERR_STATUS_INVALID_PARAMETER);
    expect_direction_error(1.1, valid, RP_KERR_STATUS_INVALID_PARAMETER);
    expect_direction_error(NAN, valid, RP_KERR_STATUS_NONFINITE_INPUT);
    CHECK(rp_kerr_null_state_from_direction(0.9, (rp_kerr_null_family)99,
        valid, state) == RP_KERR_STATUS_INVALID_PARAMETER);
    check_zero(state, 8U);
}

static void test_ergoregion_roots(void)
{
    double row[7] = {0.0, 1.8, 1.57079632679489661923, 0.0, 0.0, 0.0, 1.0};
    double metric[4][4];
    double discriminant;
    double root_minus;
    double root_plus;
    REQUIRE(rp_kerr_metric(1.0, 0.9, row, metric) == RP_KERR_STATUS_OK);
    /* A=g_phiphi > 0, B=2*g_tphi < 0, C=g_tt > 0.
     * At r=1.8, a=0.9: both corotating angular roots are positive;
     * reversing n^phi makes both scale roots negative. No root is chosen. */
    discriminant = metric[0][3] * metric[0][3] - metric[3][3] * metric[0][0];
    REQUIRE(discriminant > 0.0);
    root_minus = (-metric[0][3] - sqrt(discriminant)) / metric[3][3];
    root_plus = (-metric[0][3] + sqrt(discriminant)) / metric[3][3];
    CHECK(root_minus > 0.0 && root_plus > root_minus);
    printf("ergoregion angular roots: %.17g %.17g\n", root_minus, root_plus);
    expect_direction_error(0.9, row, RP_KERR_STATUS_INVALID_PARAMETER);
    row[6] = -1.0;
    expect_direction_error(0.9, row, RP_KERR_STATUS_NO_REAL_NULL_TANGENT);
}

static void test_constants(void)
{
    const double position[4] = {2.5, 8.0, 1.0, 0.4};
    const double spins[2] = {0.0, 0.9};
    size_t a;
    int radial;
    int polar;
    for (a = 0U; a < 2U; ++a) {
        for (radial = -1; radial <= 1; radial += 2) {
            for (polar = -1; polar <= 1; polar += 2) {
                double state[8];
                rp_kerr_null_invariants inv;
                REQUIRE(rp_kerr_null_state_from_constants(spins[a], position,
                    2.0, 20.0, radial, polar, state) == RP_KERR_STATUS_OK);
                CHECK(state[4] == 1.0);
                CHECK(state[5] * radial > 0.0);
                CHECK(state[6] * polar > 0.0);
                REQUIRE(rp_kerr_null_invariants_evaluate(spins[a], state,
                    &inv) == RP_KERR_STATUS_OK);
                /* Constants reconstructed after affine scaling: 1e-12 relative
                 * is the contract accuracy, well above regular-point roundoff. */
                CHECK(fabs(inv.impact_parameter - 2.0) / 2.0 < 1.0e-12);
                CHECK(fabs(inv.eta - 20.0) / 20.0 < 1.0e-12);
                CHECK(close_roundoff(inv.axial_angular_momentum, 2.0 * inv.energy));
                CHECK(close_roundoff(inv.carter_constant, 20.0 * inv.energy * inv.energy));
                CHECK(inv.relative_norm < 1.0e-13);
            }
        }
    }
}

static void test_constants_errors(void)
{
    const double position[4] = {0.0, 8.0, 1.0, 0.4};
    double state[8];
    const int invalid[2] = {0, 2};
    size_t i;
    CHECK(rp_kerr_null_state_from_constants(0.9, position,
        2.0, -100.0, -1, 1, state) == RP_KERR_STATUS_NO_REAL_NULL_TANGENT);
    check_zero(state, 8U);
    for (i = 0U; i < 2U; ++i) {
        CHECK(rp_kerr_null_state_from_constants(0.9, position,
            2.0, 20.0, invalid[i], 1, state) == RP_KERR_STATUS_INVALID_PARAMETER);
        check_zero(state, 8U);
        CHECK(rp_kerr_null_state_from_constants(0.9, position,
            2.0, 20.0, -1, invalid[i], state) == RP_KERR_STATUS_INVALID_PARAMETER);
        check_zero(state, 8U);
    }
    CHECK(rp_kerr_null_state_from_constants(0.9, NULL,
        2.0, 20.0, -1, 1, state) == RP_KERR_STATUS_NULL_POINTER);
    check_zero(state, 8U);
    CHECK(rp_kerr_null_state_from_constants(0.9, position,
        2.0, 20.0, -1, 1, NULL) == RP_KERR_STATUS_NULL_POINTER);
    CHECK(rp_kerr_null_state_from_constants(0.9, position,
        NAN, 20.0, -1, 1, state) == RP_KERR_STATUS_NONFINITE_INPUT);
    check_zero(state, 8U);
}

static void test_invariants_definition(void)
{
    const double position[4] = {0.0, 8.0, 1.0, 0.4};
    double state[8];
    double metric[4][4];
    double norm = 0.0;
    double denominator = 0.0;
    rp_kerr_null_invariants inv;
    size_t mu;
    size_t nu;
    REQUIRE(rp_kerr_null_state_from_constants(0.9, position,
        2.0, 20.0, -1, 1, state) == RP_KERR_STATUS_OK);
    state[5] *= 1.01; /* Deliberately off-null: make denominator observable. */
    REQUIRE(rp_kerr_metric(1.0, 0.9, state, metric) == RP_KERR_STATUS_OK);
    for (mu = 0U; mu < 4U; ++mu) {
        for (nu = 0U; nu < 4U; ++nu) {
            const double term = metric[mu][nu] * state[4U + mu] * state[4U + nu];
            norm += term;
            denominator += fabs(term);
        }
    }
    REQUIRE(rp_kerr_null_invariants_evaluate(0.9, state, &inv) == RP_KERR_STATUS_OK);
    CHECK(denominator > 0.0);
    CHECK(close_roundoff(inv.norm, norm));
    CHECK(close_roundoff(inv.relative_norm, fabs(norm) / denominator));
    CHECK(close_roundoff(inv.energy,
        -(metric[0][0] * state[4] + metric[0][3] * state[7])));
    CHECK(close_roundoff(inv.axial_angular_momentum,
        metric[3][0] * state[4] + metric[3][3] * state[7]));
    {
        const double ptheta = metric[2][2] * state[6];
        const double q = ptheta * ptheta + cos(state[2]) * cos(state[2])
            * (inv.axial_angular_momentum * inv.axial_angular_momentum
                / (sin(state[2]) * sin(state[2])) - 0.9 * 0.9 * inv.energy * inv.energy);
        CHECK(close_roundoff(inv.carter_constant, q));
    }
    CHECK(rp_kerr_null_invariants_evaluate(0.9, NULL, &inv) == RP_KERR_STATUS_NULL_POINTER);
    CHECK(inv.energy == 0.0 && inv.axial_angular_momentum == 0.0);
    CHECK(inv.carter_constant == 0.0 && inv.impact_parameter == 0.0 && inv.eta == 0.0);
    CHECK(inv.norm == 0.0 && inv.relative_norm == 0.0);
    CHECK(rp_kerr_null_invariants_evaluate(0.9, state, NULL) == RP_KERR_STATUS_NULL_POINTER);
}

static void test_views_batch(void)
{
    const double position[4] = {2.5, 8.0, 1.0, 0.4};
    double states[4][8];
    double views[4][7];
    rp_kerr_status statuses[4];
    rp_kerr_status result;
    size_t i;
    int f;
    REQUIRE(rp_kerr_null_state_from_constants(0.9, position,
        2.0, 20.0, -1, 1, states[0]) == RP_KERR_STATUS_OK);
    memcpy(states[1], states[0], sizeof(states[0]));
    memcpy(states[2], states[0], sizeof(states[0]));
    memcpy(states[3], states[0], sizeof(states[0]));
    states[1][4] = 0.0;
    states[2][4] = -1.0;
    for (i = 4U; i < 8U; ++i) {
        states[3][i] *= 2.0; /* Views preserve scale ratios without normalization. */
    }
    CHECK(rp_kerr_null_views_batch(NAN, (rp_kerr_null_family)99,
        NULL, 0U, NULL, NULL) == RP_KERR_STATUS_OK);
    for (f = 0; f < 3; ++f) {
        size_t row;
        for (row = 0U; row < 4U; ++row) {
            for (i = 0U; i < 7U; ++i) {
                views[row][i] = 9.0;
            }
        }
        result = rp_kerr_null_views_batch(0.9, (rp_kerr_null_family)f,
            &states[0][0], 4U, &views[0][0], statuses);
        CHECK(result != RP_KERR_STATUS_OK);
        CHECK(result == statuses[1]);
        CHECK(statuses[0] == RP_KERR_STATUS_OK && statuses[3] == RP_KERR_STATUS_OK);
        CHECK(statuses[1] != RP_KERR_STATUS_OK && statuses[2] != RP_KERR_STATUS_OK);
        check_zero(views[1], 7U);
        check_zero(views[2], 7U);
        for (i = 0U; i < 7U; ++i) {
            CHECK(close_roundoff(views[0][i], views[3][i]));
        }
        if (f == RP_KERR_NULL_FAMILY_BOYER_LINDQUIST) {
            for (i = 0U; i < 4U; ++i) {
                CHECK(views[3][i] == states[3][i]);
            }
            for (i = 4U; i < 7U; ++i) {
                CHECK(views[3][i] == states[3][i + 1U] / states[3][4]);
            }
        }
    }
}

static void run_case(const char *name, void (*test)(void))
{
    const int before = failures;
    test();
    printf("%s: %s\n", name, failures == before ? "PASS" : "FAIL");
}

int main(void)
{
    run_case("threshold", test_threshold);
    run_case("family conversion", test_family_conversion);
    run_case("direction families / scale / Cartesian round trip", test_direction_families);
    run_case("Schwarzschild radial direction", test_radial_direction);
    run_case("direction errors / cleared output", test_direction_errors);
    run_case("ergoregion roots", test_ergoregion_roots);
    run_case("constants / signs / invariants", test_constants);
    run_case("constants errors", test_constants_errors);
    run_case("relative norm definition / invariants errors", test_invariants_definition);
    run_case("views batch / failed rows / affine scale", test_views_batch);
    printf("null initial: %d failure(s)\n", failures);
    return failures == 0 ? 0 : 1;
}
