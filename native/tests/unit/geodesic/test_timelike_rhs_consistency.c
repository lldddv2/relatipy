#include "relatipy/kerr_geodesic.h"
#include "geodesic/integrators/kerr.h"

#include <math.h>
#include <stddef.h>
#include <stdio.h>

/* Independent Boyer--Lindquist metric and finite-difference connection.
 * The production metric, inverse, derivative and RHS evaluators are not used
 * by this oracle. All coordinates are normalized with G = c = M = 1. */
static void reference_metric(double spin, const double x[4], double g[4][4])
{
    const double r = x[1];
    const double sine = sin(x[2]);
    const double sigma = r * r + spin * spin * cos(x[2]) * cos(x[2]);
    const double delta = r * r - 2.0 * r + spin * spin;
    const double s2 = sine * sine;
    const double common = 2.0 * r / sigma;
    size_t i;
    size_t j;

    for (i = 0; i < 4; ++i) {
        for (j = 0; j < 4; ++j) {
            g[i][j] = 0.0;
        }
    }
    g[0][0] = -(1.0 - common);
    g[0][3] = -spin * s2 * common;
    g[3][0] = g[0][3];
    g[1][1] = sigma / delta;
    g[2][2] = sigma;
    g[3][3] = s2 * (r * r + spin * spin + spin * spin * s2 * common);
}

static void reference_inverse(double g[4][4], double inverse[4][4])
{
    const double determinant = g[0][0] * g[3][3] - g[0][3] * g[0][3];
    size_t i;
    size_t j;

    for (i = 0; i < 4; ++i) {
        for (j = 0; j < 4; ++j) {
            inverse[i][j] = 0.0;
        }
    }
    inverse[0][0] = g[3][3] / determinant;
    inverse[0][3] = -g[0][3] / determinant;
    inverse[3][0] = inverse[0][3];
    inverse[3][3] = g[0][0] / determinant;
    inverse[1][1] = 1.0 / g[1][1];
    inverse[2][2] = 1.0 / g[2][2];
}

static void reference_rhs(double spin, const double y[8], double derivative[8])
{
    double g[4][4];
    double inverse[4][4];
    double dg[4][4][4] = {{{0.0}}};
    double plus[4][4];
    double minus[4][4];
    double point[4];
    size_t i;
    size_t j;
    size_t k;
    size_t alpha;
    size_t beta;
    size_t sigma;

    reference_metric(spin, y, g);
    reference_inverse(g, inverse);
    for (k = 1; k <= 2; ++k) {
        const double h = 1e-5 * fmax(1.0, fabs(y[k]));
        for (i = 0; i < 4; ++i) {
            point[i] = y[i];
        }
        point[k] += h;
        reference_metric(spin, point, plus);
        point[k] -= 2.0 * h;
        reference_metric(spin, point, minus);
        for (i = 0; i < 4; ++i) {
            for (j = 0; j < 4; ++j) {
                dg[k][i][j] = (plus[i][j] - minus[i][j]) / (2.0 * h);
            }
        }
    }
    for (i = 0; i < 4; ++i) {
        derivative[i] = y[i + 4];
        derivative[i + 4] = 0.0;
        for (sigma = 0; sigma < 4; ++sigma) {
            for (alpha = 0; alpha < 4; ++alpha) {
                for (beta = 0; beta < 4; ++beta) {
                    derivative[i + 4] -= 0.5 * inverse[i][sigma]
                        * (dg[alpha][sigma][beta]
                            + dg[beta][sigma][alpha]
                            - dg[sigma][alpha][beta])
                        * y[alpha + 4] * y[beta + 4];
                }
            }
        }
    }
}

static int corrected_rhs(double spin, const double y[8], double derivative[8])
{
    return rp_kerr_null_geodesic_rhs(
        1.0, spin, y, y + 4, derivative, derivative + 4
    ) == RP_KERR_STATUS_OK;
}

static int adapter_rhs(double spin, const double y[8], double derivative[8])
{
    rp_kerr_integrator_context context = {.mass = 1.0, .spin = spin};
    return rp_kerr_integrator_rhs(0.0, y, derivative, &context)
        == RP_KERR_STATUS_OK;
}

static double norm(double spin, const double y[8])
{
    double g[4][4];
    double result = 0.0;
    size_t i;
    size_t j;

    reference_metric(spin, y, g);
    for (i = 0; i < 4; ++i) {
        for (j = 0; j < 4; ++j) {
            result += g[i][j] * y[i + 4] * y[j + 4];
        }
    }
    return result;
}

static void invariants(double spin, const double y[8], double result[3])
{
    double g[4][4];
    reference_metric(spin, y, g);
    result[0] = norm(spin, y);
    result[1] = -(g[0][0] * y[4] + g[0][3] * y[7]);
    result[2] = g[3][0] * y[4] + g[3][3] * y[7];
}

static void initial_state(double spin, const double x[4],
    const double velocity[3], double y[8])
{
    double g[4][4];
    double w[4] = {1.0, velocity[0], velocity[1], velocity[2]};
    double squared = 0.0;
    size_t i;
    size_t j;
    reference_metric(spin, x, g);
    for (i = 0; i < 4; ++i) {
        for (j = 0; j < 4; ++j) {
            squared += g[i][j] * w[i] * w[j];
        }
    }
    for (i = 0; i < 4; ++i) {
        y[i] = x[i];
        y[i + 4] = w[i] / sqrt(-squared);
    }
}

static int close_value(double actual, double expected, double tolerance,
    const char *label)
{
    const double difference = fabs(actual - expected);
    if (!isfinite(actual) || difference > tolerance) {
        (void)fprintf(stderr, "%s: actual %.17g expected %.17g difference %.3g tolerance %.3g\n",
            label, actual, expected, difference, tolerance);
        return 0;
    }
    return 1;
}

static int rk4_step(double spin, double y[8], double h,
    int (*rhs)(double, const double[8], double[8]))
{
    double k1[8];
    double k2[8];
    double k3[8];
    double k4[8];
    double stage[8];
    size_t i;
    if (!rhs(spin, y, k1)) return 0;
    for (i = 0; i < 8; ++i) stage[i] = y[i] + 0.5 * h * k1[i];
    if (!rhs(spin, stage, k2)) return 0;
    for (i = 0; i < 8; ++i) stage[i] = y[i] + 0.5 * h * k2[i];
    if (!rhs(spin, stage, k3)) return 0;
    for (i = 0; i < 8; ++i) stage[i] = y[i] + h * k3[i];
    if (!rhs(spin, stage, k4)) return 0;
    for (i = 0; i < 8; ++i) {
        y[i] += h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
    }
    return 1;
}

static int reference_wrapper(double spin, const double y[8], double derivative[8])
{
    reference_rhs(spin, y, derivative);
    return 1;
}

int main(void)
{
    const double spins[3] = {0.0, 0.5, -0.9};
    const double positions[3][4] = {
        {0.0, 8.0, 1.1, 0.3},
        {0.2, 12.0, 0.9, -0.4},
        {-1.0, 6.0, 1.35, 2.0}
    };
    const double velocities[3][3] = {
        {-0.01, 0.002, 0.02},
        {0.015, -0.001, -0.015},
        {-0.02, 0.003, 0.01}
    };
    size_t point;
    int failures = 0;
    double largest_rhs_difference = 0.0;
    double largest_trajectory_difference = 0.0;
    double largest_norm_drift = 0.0;
    double largest_energy_drift = 0.0;
    double largest_angular_momentum_drift = 0.0;
    double largest_legacy_difference = 0.0;
    {
        const double radii[2] = {8.0, 15.0};
        size_t case_index;
        for (case_index = 0; case_index < 2; ++case_index) {
            const double r = radii[case_index];
            const double y[8] = {
                0.0, r, 1.5707963267948966, 0.0,
                1.0 / sqrt(1.0 - 2.0 / r), 0.0, 0.0, 0.0
            };
            double derivative[8];
            if (!corrected_rhs(0.0, y, derivative)) return 2;
            failures += !close_value(derivative[5], -1.0 / (r * r),
                2e-14, "Schwarzschild stationary radial acceleration");
            failures += !close_value(derivative[4], 0.0, 2e-14,
                "Schwarzschild stationary time acceleration");
            failures += !close_value(derivative[6], 0.0, 2e-14,
                "Schwarzschild stationary polar acceleration");
            failures += !close_value(derivative[7], 0.0, 2e-14,
                "Schwarzschild stationary azimuth acceleration");
        }
    }
    for (point = 0; point < 3; ++point) {
        double y[8];
        double numerical[8];
        double adapted[8];
        double production[8];
        double independent[8];
        double constructed[4];
        double legacy_dx[4];
        double legacy_du[4];
        double initial[3];
        double final[3];
        size_t i;
        size_t step;
        const double spin = spins[point];
        initial_state(spin, positions[point], velocities[point], y);
        if (rp_kerr_four_velocity(1.0, spin, y, velocities[point], constructed)
            != RP_KERR_STATUS_OK) return 2;
        for (i = 0; i < 4; ++i) {
            failures += !close_value(constructed[i], y[i + 4], 2e-14,
                "initial four-velocity order");
        }
        failures += !close_value(norm(spin, y), -1.0, 2e-14,
            "initial timelike norm");
        if (!corrected_rhs(spin, y, numerical)) return 2;
        if (!adapter_rhs(spin, y, adapted)) return 2;
        reference_rhs(spin, y, independent);
        for (i = 0; i < 8; ++i) {
            largest_rhs_difference = fmax(largest_rhs_difference,
                fabs(numerical[i] - independent[i]));
            failures += !close_value(numerical[i], independent[i], 4e-9,
                "corrected timelike RHS versus independent connection");
            failures += !close_value(adapted[i], independent[i], 4e-9,
                "timelike adapter versus independent connection");
        }
        if (rp_kerr_geodesic_rhs(1.0, spin, y, y + 4,
            legacy_dx, legacy_du) != RP_KERR_STATUS_OK) return 2;
        for (i = 0; i < 4; ++i) {
            largest_legacy_difference = fmax(largest_legacy_difference,
                fabs(legacy_du[i] - numerical[i + 4]));
        }
        invariants(spin, y, initial);
        for (i = 0; i < 8; ++i) production[i] = independent[i] = y[i];
        /* Both trajectories advance in s = tau/T0; their RHS implementations
           are independent, while the fixed RK4 tableau is shared. */
        for (step = 0; step < 200; ++step) {
            if (!rk4_step(spin, production, 0.005, adapter_rhs)
                || !rk4_step(spin, independent, 0.005, reference_wrapper)) {
                return 2;
            }
        }
        for (i = 0; i < 8; ++i) {
            largest_trajectory_difference = fmax(largest_trajectory_difference,
                fabs(production[i] - independent[i]));
            failures += !close_value(production[i], independent[i], 2e-8,
                "short proper-time trajectory");
        }
        invariants(spin, production, final);
        largest_norm_drift = fmax(largest_norm_drift, fabs(final[0] + 1.0));
        largest_energy_drift = fmax(largest_energy_drift,
            fabs(final[1] - initial[1]));
        largest_angular_momentum_drift = fmax(
            largest_angular_momentum_drift, fabs(final[2] - initial[2]));
        failures += !close_value(final[0], -1.0, 2e-10,
            "trajectory timelike norm");
        failures += !close_value(final[1], initial[1], 2e-10,
            "trajectory energy");
        failures += !close_value(final[2], initial[2], 2e-10,
            "trajectory axial angular momentum");
    }
    (void)printf("max |RHS correction - oracle| %.3g, max |state correction - oracle| %.3g\n",
        largest_rhs_difference, largest_trajectory_difference);
    (void)printf("max drifts: norm %.3g, E %.3g, Lz %.3g; max |legacy - correction| %.3g\n",
        largest_norm_drift, largest_energy_drift,
        largest_angular_momentum_drift, largest_legacy_difference);
    if (failures != 0) {
        (void)fprintf(stderr, "%d timelike consistency checks failed\n", failures);
        return 1;
    }
    (void)puts("Timelike corrected RHS consistency passed");
    return 0;
}
