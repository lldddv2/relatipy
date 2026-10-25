/*
 * Frame map, coordinate velocity and projected four-acceleration.
 *
 * The map and velocity are checked against the independent conversion in
 * geodesic/initial/convert.c, the Jacobian against central differences of
 * the map, and f^mu against (g^(mu nu) + u^mu u^nu) F_nu with F_nu taken as
 * central differences of a Plummer potential composed with the map.
 */

#include "force/external.h"
#include "force_test_models.h"
#include "geodesic/initial/convert.h"
#include "relatipy/kerr_geometry.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

#define POINT_COUNT 4

static const double SPINS[3] = {0.0, 0.5, 0.9};
static const double POINTS[POINT_COUNT][7] = {
    /* t, r, theta, phi, dr/dt, dtheta/dt, dphi/dt */
    {0.0, 8.0, 1.1, 0.3, -0.05, 0.004, 0.03},
    {3.0, 25.0, 0.4, 2.7, 0.02, -0.002, 0.006},
    {-1.0, 120.0, 2.6, -1.9, 0.01, 0.0005, -0.0007},
    {0.0, 6.0, 1.5707963267948966, 4.0, 0.0, 0.0, 0.06}
};

static double relative_error(double value, double reference, double scale)
{
    return fabs(value - reference) / (fabs(reference) + scale);
}

static void state_at(double spin, size_t point, double state[8])
{
    double norm = 0.0;
    double metric[RP_KERR_DIM][RP_KERR_DIM];
    size_t mu;
    size_t nu;

    assert(rp_initial_from_bl(spin, POINTS[point], state)
        == RP_KERR_STATUS_OK);
    assert(rp_kerr_metric(1.0, spin, state, metric) == RP_KERR_STATUS_OK);
    for (mu = 0U; mu < 4U; ++mu) {
        for (nu = 0U; nu < 4U; ++nu) {
            norm += metric[mu][nu] * state[4U + mu] * state[4U + nu];
        }
    }
    assert(fabs(norm + 1.0) < 1e-12);
}

static void check_map_and_velocity(double spin, size_t point)
{
    struct rp_kerr_geodesic_state prepared;
    rp_force_frame frame;
    double state[8];
    double cartesian[7];
    double position[3];
    double jacobian[3][3];
    double velocity[3];
    size_t k;
    size_t j;

    state_at(spin, point, state);
    assert(rp_initial_canonical_to_cartesian(spin, state, cartesian)
        == RP_KERR_STATUS_OK);
    assert(rp_force_bl_to_cartesian(spin, state, position, jacobian)
        == RP_KERR_STATUS_OK);
    for (k = 0U; k < RP_KERR_DIM; ++k) {
        prepared.x[k] = state[k];
        prepared.u[k] = state[4U + k];
    }
    assert(rp_kerr_geodesic_physic_prepare_state(1.0, spin, &prepared)
        == RP_KERR_STATUS_OK);
    rp_force_frame_prepare(spin, &prepared, &frame);
    rp_force_frame_velocity(&frame, prepared.u, velocity);
    for (k = 0U; k < 3U; ++k) {
        assert(relative_error(position[k], cartesian[1U + k], 1e-12) < 1e-14);
        assert(frame.position[k] == position[k]);
        assert(relative_error(velocity[k], cartesian[4U + k], 1e-12) < 1e-13);
    }

    /* Central differences of the map in r, theta and phi. */
    for (j = 0U; j < 3U; ++j) {
        double plus[4];
        double minus[4];
        double position_plus[3];
        double position_minus[3];
        const double step = 1e-5 * (j == 0U ? state[1] : 1.0);

        for (k = 0U; k < 4U; ++k) {
            plus[k] = minus[k] = state[k];
        }
        plus[1U + j] += step;
        minus[1U + j] -= step;
        assert(rp_force_bl_to_cartesian(spin, plus, position_plus, NULL)
            == RP_KERR_STATUS_OK);
        assert(rp_force_bl_to_cartesian(spin, minus, position_minus, NULL)
            == RP_KERR_STATUS_OK);
        for (k = 0U; k < 3U; ++k) {
            const double numeric =
                (position_plus[k] - position_minus[k]) / (2.0 * step);
            assert(frame.jacobian[k][j] == jacobian[k][j]);
            assert(fabs(jacobian[k][j] - numeric)
                < 1e-8 * (1.0 + fabs(numeric)) * (j == 0U ? 1.0 : state[1]));
        }
    }
}

static double composed_potential(
    const test_plummer *model,
    double spin,
    const double coordinates[4]
)
{
    double position[3];

    assert(rp_force_bl_to_cartesian(spin, coordinates, position, NULL)
        == RP_KERR_STATUS_OK);
    return test_plummer_potential(model, position);
}

static void check_projection(double spin, size_t point)
{
    const test_plummer model = {0.3, 4.0, {1.5, -2.0, 0.7}};
    const rp_external_force force = {test_plummer_acceleration, &model};
    double state[8];
    double inverse[RP_KERR_DIM][RP_KERR_DIM];
    double metric[RP_KERR_DIM][RP_KERR_DIM];
    double covector[RP_KERR_DIM] = {0.0, 0.0, 0.0, 0.0};
    double expected[RP_KERR_DIM];
    double acceleration[RP_KERR_DIM];
    double u_dot_force = 0.0;
    double orthogonality = 0.0;
    double scale = 0.0;
    size_t mu;
    size_t nu;

    state_at(spin, point, state);
    assert(rp_force_four_acceleration(
        spin, &force, state, state + 4, acceleration
    ) == 0);
    assert(rp_kerr_inverse_metric(1.0, spin, state, inverse)
        == RP_KERR_STATUS_OK);
    assert(rp_kerr_metric(1.0, spin, state, metric) == RP_KERR_STATUS_OK);
    /* F_nu = -d_nu Phi; F_t = 0 because the potential is static. */
    for (mu = 1U; mu < 4U; ++mu) {
        double plus[4];
        double minus[4];
        const double step = 1e-5 * (mu == 1U ? state[1] : 1.0);

        for (nu = 0U; nu < 4U; ++nu) {
            plus[nu] = minus[nu] = state[nu];
        }
        plus[mu] += step;
        minus[mu] -= step;
        covector[mu] = -(composed_potential(&model, spin, plus)
            - composed_potential(&model, spin, minus)) / (2.0 * step);
        u_dot_force += state[4U + mu] * covector[mu];
    }
    for (mu = 0U; mu < 4U; ++mu) {
        expected[mu] = state[4U + mu] * u_dot_force;
        for (nu = 0U; nu < 4U; ++nu) {
            expected[mu] += inverse[mu][nu] * covector[nu];
        }
        scale += fabs(expected[mu]);
    }
    for (mu = 0U; mu < 4U; ++mu) {
        assert(fabs(acceleration[mu] - expected[mu]) < 1e-7 * scale);
        for (nu = 0U; nu < 4U; ++nu) {
            orthogonality += metric[mu][nu] * state[4U + mu]
                * acceleration[nu];
        }
    }
    assert(fabs(orthogonality) < 1e-14 * scale * state[4]);
}

static void check_statuses(void)
{
    const test_plummer model = {0.3, 4.0, {0.0, 0.0, 0.0}};
    const rp_external_force failing = {test_failing_acceleration, NULL};
    const rp_external_force empty = {NULL, &model};
    const double x[4] = {0.0, 8.0, 1.1, 0.3};
    const double u[4] = {1.2, 0.01, 0.002, 0.03};
    const double bad[4] = {0.0, NAN, 1.1, 0.3};
    const double axis[4] = {0.0, 8.0, 0.0, 0.3};
    double f[4] = {1.0, 1.0, 1.0, 1.0};
    double position[3] = {1.0, 1.0, 1.0};
    size_t i;

    assert(rp_force_four_acceleration(0.5, NULL, x, u, f) == 0);
    for (i = 0U; i < 4U; ++i) {
        assert(f[i] == 0.0);
        f[i] = 1.0;
    }
    assert(rp_force_four_acceleration(0.5, &empty, x, u, f) == 0);
    assert(f[0] == 0.0 && f[3] == 0.0);
    f[0] = 1.0;
    assert(rp_force_four_acceleration(0.5, &failing, x, u, f)
        == RP_FORCE_STATUS_CALLBACK_FAILURE);
    assert(f[0] == 0.0 && f[1] == 0.0 && f[2] == 0.0 && f[3] == 0.0);
    assert(rp_force_four_acceleration(0.5, &failing, bad, u, f)
        == (int)RP_KERR_STATUS_NONFINITE_INPUT);
    assert(rp_force_four_acceleration(1.5, &failing, x, u, f)
        == (int)RP_KERR_STATUS_INVALID_PARAMETER);
    assert(rp_force_four_acceleration(0.5, &failing, axis, u, f)
        == (int)RP_KERR_STATUS_COORDINATE_SINGULARITY);
    assert(rp_force_four_acceleration(0.5, &failing, NULL, u, f)
        == (int)RP_KERR_STATUS_NULL_POINTER);
    assert(rp_force_four_acceleration(0.5, &failing, x, u, NULL)
        == (int)RP_KERR_STATUS_NULL_POINTER);
    assert(rp_force_bl_to_cartesian(0.5, x, NULL, NULL)
        == RP_KERR_STATUS_NULL_POINTER);
    assert(rp_force_bl_to_cartesian(0.5, bad, position, NULL)
        == RP_KERR_STATUS_NONFINITE_INPUT);
    assert(position[0] == 0.0 && position[1] == 0.0 && position[2] == 0.0);
    assert(rp_force_bl_to_cartesian(-1.5, x, position, NULL)
        == RP_KERR_STATUS_INVALID_PARAMETER);
    /* The map itself is regular on the axis. */
    assert(rp_force_bl_to_cartesian(0.5, axis, position, NULL)
        == RP_KERR_STATUS_OK);
    assert(fabs(position[2] - 8.0) < 1e-15);
}

int main(void)
{
    size_t spin;
    size_t point;

    for (spin = 0U; spin < 3U; ++spin) {
        for (point = 0U; point < POINT_COUNT; ++point) {
            check_map_and_velocity(SPINS[spin], point);
            check_projection(SPINS[spin], point);
        }
    }
    check_statuses();
    puts("test_external: ok");
    return 0;
}
