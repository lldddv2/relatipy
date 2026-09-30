/* Analytic corrected Kerr Jacobian, in geometric units G = c = 1.
 * Central differences below are test oracles only. The production evaluator
 * contains no perturbations. These checks do not prove global chart accuracy. */
#include "geodesic/jacobian.h"
#include "geodesic/physic/kerr.h"

#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

static int close_value(double actual, double expected, double tolerance,
    const char *label, size_t row, size_t column)
{
    if (!isfinite(actual) || !isfinite(expected)
        || fabs(actual - expected) > tolerance) {
        (void)fprintf(stderr, "%s[%zu,%zu]: %.17g expected %.17g tolerance %.3g\n",
            label, row, column, actual, expected, tolerance);
        return 0;
    }
    return 1;
}

static int rhs(double mass, double spin, const double y[8], double f[8])
{
    return rp_kerr_null_geodesic_rhs(mass, spin, y, y + 4, f, f + 4)
        == RP_KERR_STATUS_OK;
}

static int check_regular(double mass, double spin, const double y[8])
{
    const double relative_steps[3] = {1e-4, 3e-5, 1e-5};
    double jacobian[8][8];
    double f[8];
    size_t row;
    size_t column;
    int failures = 0;

    if (rp_kerr_geodesic_jacobian(mass, spin, y, y + 4, jacobian)
        != RP_KERR_STATUS_OK || !rhs(mass, spin, y, f)) return 1;
    for (column = 0; column < 8; ++column) {
        double estimates[3][8];
        size_t probe;
        for (probe = 0; probe < 3; ++probe) {
            const double h = relative_steps[probe] * fmax(1.0, fabs(y[column]));
            double plus[8];
            double minus[8];
            double f_plus[8];
            double f_minus[8];
            (void)memcpy(plus, y, sizeof(plus));
            (void)memcpy(minus, y, sizeof(minus));
            plus[column] += h;
            minus[column] -= h;
            if (!rhs(mass, spin, plus, f_plus)
                || !rhs(mass, spin, minus, f_minus)) return failures + 1;
            for (row = 0; row < 8; ++row) {
                estimates[probe][row] = (f_plus[row] - f_minus[row])
                    / (plus[column] - minus[column]);
            }
        }
        for (row = 0; row < 8; ++row) {
            double best = estimates[0][row];
            size_t probe;
            for (probe = 1; probe < 3; ++probe) {
                if (fabs(estimates[probe][row] - jacobian[row][column])
                    < fabs(best - jacobian[row][column])) {
                    best = estimates[probe][row];
                }
            }
            failures += !close_value(jacobian[row][column], best,
                3e-10 + 3e-7 * fabs(jacobian[row][column]),
                "central differences, three scales", row, column);
            if (row < 4) {
                failures += !close_value(jacobian[row][column],
                    column == row + 4 ? 1.0 : 0.0, 0.0,
                    "exact upper blocks", row, column);
            }
            if (column == 0 || column == 3) {
                failures += !close_value(jacobian[row][column], 0.0, 0.0,
                    "stationarity and axial symmetry", row, column);
            }
        }
    }

    /* Polarization is an exact reference for the quadratic velocity RHS:
     * (a(u + e_j) - a(u - e_j))/2 = d a / d u_j, for unit e_j.
     * There is no step-size/truncation approximation in this identity. */
    for (column = 0; column < 4; ++column) {
        double plus[8];
        double minus[8];
        double f_plus[8];
        double f_minus[8];
        (void)memcpy(plus, y, sizeof(plus));
        (void)memcpy(minus, y, sizeof(minus));
        plus[column + 4] += 1.0;
        minus[column + 4] -= 1.0;
        if (!rhs(mass, spin, plus, f_plus)
            || !rhs(mass, spin, minus, f_minus)) return failures + 1;
        for (row = 4; row < 8; ++row) {
            const double expected = 0.5 * (f_plus[row] - f_minus[row]);
            failures += !close_value(jacobian[row][column + 4], expected,
                5e-13 * (1.0 + fabs(expected)),
                "quadratic polarization", row, column + 4);
        }
    }
    for (row = 4; row < 8; ++row) {
        double contraction = 0.0;
        for (column = 0; column < 4; ++column) {
            contraction += jacobian[row][column + 4] * y[column + 4];
        }
        failures += !close_value(contraction, 2.0 * f[row],
            5e-13 * (1.0 + fabs(f[row])), "quadratic homogeneity", row, 0);
    }
    return failures;
}

static int check_reflection(void)
{
    const double signs[8] = {1.0, 1.0, -1.0, 1.0, 1.0, 1.0, -1.0, 1.0};
    const double y[8] = {0.0, 7.0, 0.7, 0.2, 1.2, -0.03, 0.009, 0.02};
    double reflected[8];
    double first[8][8];
    double second[8][8];
    size_t row;
    size_t column;
    int failures = 0;
    (void)memcpy(reflected, y, sizeof(reflected));
    reflected[2] = acos(-1.0) - y[2];
    reflected[6] = -y[6];
    if (rp_kerr_geodesic_jacobian(1.0, 0.5, y, y + 4, first)
            != RP_KERR_STATUS_OK
        || rp_kerr_geodesic_jacobian(1.0, 0.5, reflected, reflected + 4, second)
            != RP_KERR_STATUS_OK) return 1;
    for (row = 0; row < 8; ++row) {
        for (column = 0; column < 8; ++column) {
            const double expected = signs[row] * signs[column] * first[row][column];
            failures += !close_value(second[row][column], expected,
                5e-13 * (1.0 + fabs(expected)), "equatorial reflection", row, column);
        }
    }
    return failures;
}

static int check_schwarzschild(void)
{
    const double mass = 2.0;
    const double radius = 14.0;
    const double y[8] = {0.0, 14.0, 1.1, 0.0, 1.3, 0.0, 0.0, 0.0};
    const double delta = radius * radius - 2.0 * mass * radius;
    double jacobian[8][8];
    int failures = 0;
    if (rp_kerr_geodesic_jacobian(mass, 0.0, y, y + 4, jacobian)
        != RP_KERR_STATUS_OK) return 1;
    /* d_r [-M (1 - 2M/r) u_t^2/r^2], with u fixed. */
    failures += !close_value(jacobian[5][1],
        2.0 * mass * y[4] * y[4] * (radius - 3.0 * mass)
            / (radius * radius * radius * radius),
        1e-15, "Schwarzschild radial coordinate derivative", 5, 1);
    failures += !close_value(jacobian[5][4],
        -2.0 * mass * delta * y[4]
            / (radius * radius * radius * radius),
        1e-15, "Schwarzschild time tangent derivative", 5, 4);
    failures += !close_value(jacobian[4][5],
        -2.0 * mass * y[4] / delta,
        1e-15, "Schwarzschild radial tangent derivative", 4, 5);
    return failures;
}

static int check_near_pole(double theta)
{
    const double pi = acos(-1.0);
    const double guard = 64.0 * DBL_EPSILON;
    const double relative_steps[3] = {1e-3, 3e-4, 1e-4};
    const double distance = fmin(theta, pi - theta);
    const double y[8] = {0.0, 8.0, theta, 0.0, 1.2, -0.01, 0.002, 0.03};
    double jacobian[8][8];
    double best[4] = {INFINITY, INFINITY, INFINITY, INFINITY};
    double best_roundoff[4] = {0.0};
    size_t probe;
    size_t row;
    int failures = 0;
    if (rp_kerr_geodesic_jacobian(1.0, 0.5, y, y + 4, jacobian)
        != RP_KERR_STATUS_OK) return 1;
    for (probe = 0; probe < 3; ++probe) {
        const double h = relative_steps[probe] * distance;
        double plus[8];
        double minus[8];
        double f_plus[8];
        double f_minus[8];
        double separation;
        (void)memcpy(plus, y, sizeof(plus));
        (void)memcpy(minus, y, sizeof(minus));
        plus[2] += h;
        minus[2] -= h;
        /* Both probes stay inside the original chart. Use represented probe
         * separation, essential near pi where angular ULP is much larger. */
        if (!(minus[2] > guard && pi - plus[2] > guard)
            || plus[2] == minus[2]
            || !rhs(1.0, 0.5, plus, f_plus)
            || !rhs(1.0, 0.5, minus, f_minus)) return failures + 1;
        separation = plus[2] - minus[2];
        for (row = 4; row < 8; ++row) {
            const double estimate = (f_plus[row] - f_minus[row]) / separation;
            const double error = fabs(estimate - jacobian[row][2]);
            const double roundoff = 128.0 * DBL_EPSILON
                * (fabs(f_plus[row]) + fabs(f_minus[row])) / separation;
            if (error < best[row - 4]) {
                best[row - 4] = error;
                best_roundoff[row - 4] = roundoff;
            }
        }
    }
    for (row = 4; row < 8; ++row) {
        failures += !close_value(best[row - 4], 0.0,
            1e-9 + 3e-6 * fabs(jacobian[row][2]) + best_roundoff[row - 4],
            "local polar coordinate probes", row, 2);
    }
    /* Guard-adjacent north angles permit relative probes in double. These
     * tests demonstrate finite local agreement, not uniform precision at an
     * axis. Cancellation and angular representation still bound accuracy. */
    return failures;
}

static int check_domain(void)
{
    const double pi = acos(-1.0);
    const double guard = 64.0 * DBL_EPSILON;
    const double angles[6] = {0.0, guard, -0.1, pi, pi - guard, pi + 0.1};
    const double u[4] = {1.0, -0.01, 0.002, 0.01};
    size_t index;
    int failures = 0;
    for (index = 0; index < 6; ++index) {
        const double x[4] = {0.0, 8.0, angles[index], 0.0};
        double jacobian[8][8];
        double dx[4];
        double du[4];
        rp_kerr_status status;
        rp_kerr_status rhs_status;
        size_t row;
        size_t column;
        for (row = 0; row < 8; ++row) {
            for (column = 0; column < 8; ++column) jacobian[row][column] = 17.0;
        }
        status = rp_kerr_geodesic_jacobian(1.0, 0.5, x, u, jacobian);
        rhs_status = rp_kerr_null_geodesic_rhs(1.0, 0.5, x, u, dx, du);
        failures += status != RP_KERR_STATUS_COORDINATE_SINGULARITY;
        failures += status != rhs_status;
        for (row = 0; row < 8; ++row) {
            for (column = 0; column < 8; ++column) {
                failures += !close_value(jacobian[row][column], 0.0, 0.0,
                    "domain failure clears output", row, column);
            }
        }
    }
    return failures;
}

int main(void)
{
    const double spins[3] = {0.0, 0.5, 1.0};
    const double states[4][8] = {
        {0.0, 8.0, 1.1, 0.2, 1.3, -0.02, 0.003, 0.02},
        {-0.3, 3.0, 0.6, 4.0, 1.8, 0.06, -0.02, 0.12},
        {2.0, 23.0, 2.1, -0.6, 1.05, 0.015, 0.004, -0.006},
        {1.0, 6.0, 1.5707963267948966, 0.0, 1.2, -0.03, 0.01, 0.02}
    };
    const double explicit_mass_state[8] = {
        0.0, 16.0, 0.8, 0.0, 1.4, 0.03, -0.007, 0.01
    };
    const double zero_tangent_state[8] = {0.0, 8.0, 0.8, 0.0, 0.0, 0.0, 0.0, 0.0};
    size_t spin;
    size_t state;
    int failures = 0;
    for (spin = 0; spin < 3; ++spin) {
        for (state = 0; state < 4; ++state) {
            failures += check_regular(1.0, spins[spin], states[state]);
        }
    }
    failures += check_regular(2.3, 0.75, explicit_mass_state);
    failures += check_regular(1.0, 0.5, zero_tangent_state);
    failures += check_reflection();
    failures += check_schwarzschild();
    failures += check_near_pole(1e-7);
    failures += check_near_pole(128.0 * DBL_EPSILON);
    failures += check_near_pole(acos(-1.0) - 1e-7);
    failures += check_domain();
    if (failures != 0) {
        (void)fprintf(stderr, "Kerr Jacobian failures: %d\n", failures);
        return 1;
    }
    (void)puts("Kerr analytic Jacobian: PASS");
    return 0;
}
