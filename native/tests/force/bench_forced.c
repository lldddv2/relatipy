/*
 * Benchmark of the forced Kerr layer (not a test; prints timings).
 *
 * 1. Nanoseconds per call of the homogeneous corrected RHS adapter
 *    (rp_kerr_integrator_rhs), the legacy optimized RHS for reference, and
 *    the forced adapter with no force, a Plummer sphere and a
 *    Miyamoto--Nagai disk.
 * 2. Time per observable evaluation of an S2-like orbit (a = 0.125 arcsec,
 *    e = 0.884, 145 arrival epochs over 30 yr, Sgr A* at 4.3e6 Msun and
 *    8.28 kpc) with rp_kerr_mcmc_evaluate and rp_kerr_forced_mcmc_evaluate.
 *
 * Build like the force tests, with -std=c11 (timespec_get) and -O3, and run
 * without arguments.  RHS variants are interleaved within each repetition
 * and every figure is the best repetition, to reduce noise on a shared
 * machine.
 */

#include "force/kerr_forced.h"
#include "force/kerr_mcmc.h"
#include "force_test_models.h"
#include "geodesic/integrators/kerr.h"
#include "geodesic/kerr_mcmc.h"
#include "relatipy/kerr_geodesic.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#define POINTS 512U
#define RHS_SWEEPS 200U
#define REPEATS 41U
#define EPOCHS 145U

static const double PI = 3.14159265358979323846;

static double now_seconds(void)
{
    struct timespec spec;
    timespec_get(&spec, TIME_UTC);
    return (double)spec.tv_sec + 1e-9 * (double)spec.tv_nsec;
}

static double STATES[POINTS][8];
static volatile double SINK;

static void fill_states(void)
{
    size_t i;

    for (i = 0U; i < POINTS; ++i) {
        const double phase = (double)i / (double)POINTS;
        STATES[i][0] = 0.0;
        STATES[i][1] = 10.0 + 4000.0 * phase;
        STATES[i][2] = 0.3 + 2.5 * phase;
        STATES[i][3] = 6.0 * phase;
        STATES[i][4] = 1.05;
        STATES[i][5] = 0.01 * (phase - 0.5);
        STATES[i][6] = 1e-4 * phase;
        STATES[i][7] = 1e-3 / STATES[i][1];
    }
}

typedef int (*packed_rhs)(double, const double[], double[], void *);

/* One sweep set of a packed adapter; returns seconds. */
static double sweep_packed(packed_rhs rhs, void *context)
{
    double derivative[8];
    double sum = 0.0;
    double start;
    size_t sweep;
    size_t i;

    start = now_seconds();
    for (sweep = 0U; sweep < RHS_SWEEPS; ++sweep) {
        for (i = 0U; i < POINTS; ++i) {
            if (rhs(0.0, STATES[i], derivative, context) != 0) {
                return NAN;
            }
            sum += derivative[5];
        }
    }
    SINK = sum;
    return now_seconds() - start;
}

static double sweep_legacy(double spin)
{
    double dx[4];
    double du[4];
    double sum = 0.0;
    double start;
    size_t sweep;
    size_t i;

    start = now_seconds();
    for (sweep = 0U; sweep < RHS_SWEEPS; ++sweep) {
        for (i = 0U; i < POINTS; ++i) {
            if (rp_kerr_geodesic_rhs(
                    1.0, spin, STATES[i], STATES[i] + 4, dx, du
                ) != RP_KERR_STATUS_OK) {
                return NAN;
            }
            sum += du[1];
        }
    }
    SINK = sum;
    return now_seconds() - start;
}

/* The force term alone at prepared states: frame, callback, projection. */
static double sweep_force_term(double spin, const rp_external_force *force)
{
    static struct rp_kerr_geodesic_state prepared[POINTS];
    double acceleration[4];
    double sum = 0.0;
    double start;
    size_t sweep;
    size_t i;
    size_t k;

    for (i = 0U; i < POINTS; ++i) {
        for (k = 0U; k < 4U; ++k) {
            prepared[i].x[k] = STATES[i][k];
            prepared[i].u[k] = STATES[i][4U + k];
        }
        assert(rp_kerr_geodesic_physic_prepare_state(1.0, spin, &prepared[i])
            == RP_KERR_STATUS_OK);
    }
    start = now_seconds();
    for (sweep = 0U; sweep < RHS_SWEEPS; ++sweep) {
        for (i = 0U; i < POINTS; ++i) {
            if (rp_force_evaluate_prepared(
                    spin, &prepared[i], force, acceleration
                ) != 0) {
                return NAN;
            }
            sum += acceleration[1];
        }
    }
    SINK = sum;
    return now_seconds() - start;
}

#define RHS_VARIANTS 8U

static void bench_rhs(void)
{
    const double spin = 0.7;
    const test_plummer plummer = {1e-3, 12000.0, {0.0, 0.0, 0.0}};
    const test_miyamoto_nagai disk = {1e-3, 8000.0, 80.0};
    const rp_external_force plummer_force = {test_plummer_acceleration, &plummer};
    const rp_external_force disk_force = {test_miyamoto_nagai_acceleration, &disk};
    const char *names[RHS_VARIANTS] = {
        "homogeneous corrected (rp_kerr_integrator_rhs)",
        "legacy optimized (rp_kerr_geodesic_rhs)",
        "forced, NULL force",
        "forced, Plummer sphere",
        "forced, Miyamoto-Nagai disk",
        "forced coordinate time, Miyamoto-Nagai disk",
        "force term only, Plummer sphere",
        "force term only, Miyamoto-Nagai disk"
    };
    rp_kerr_integrator_context homogeneous;
    rp_kerr_forced_context none;
    rp_kerr_forced_context with_plummer;
    rp_kerr_forced_context with_disk;
    double best[RHS_VARIANTS];
    double scale;
    size_t repeat;
    size_t v;

    fill_states();
    memset(&homogeneous, 0, sizeof(homogeneous));
    homogeneous.mass = 1.0;
    homogeneous.spin = spin;
    none.spin = with_plummer.spin = with_disk.spin = spin;
    none.force = NULL;
    with_plummer.force = &plummer_force;
    with_disk.force = &disk_force;
    for (v = 0U; v < RHS_VARIANTS; ++v) {
        best[v] = INFINITY;
    }
    /* Interleave variants so that machine load affects all of them alike. */
    for (repeat = 0U; repeat < REPEATS; ++repeat) {
        best[0] = fmin(best[0], sweep_packed(rp_kerr_integrator_rhs, &homogeneous));
        best[1] = fmin(best[1], sweep_legacy(spin));
        best[2] = fmin(best[2], sweep_packed(rp_kerr_forced_integrator_rhs, &none));
        best[3] = fmin(best[3],
            sweep_packed(rp_kerr_forced_integrator_rhs, &with_plummer));
        best[4] = fmin(best[4],
            sweep_packed(rp_kerr_forced_integrator_rhs, &with_disk));
        best[5] = fmin(best[5],
            sweep_packed(rp_kerr_forced_coordinate_time_rhs, &with_disk));
        best[6] = fmin(best[6], sweep_force_term(spin, &plummer_force));
        best[7] = fmin(best[7], sweep_force_term(spin, &disk_force));
    }
    scale = 1e9 / (double)(RHS_SWEEPS * POINTS);
    printf("RHS ns/call (spin %.1f, %u points x %u sweeps, best of %u)\n",
        spin, POINTS, RHS_SWEEPS, REPEATS);
    for (v = 0U; v < RHS_VARIANTS; ++v) {
        printf("  %-47s %7.1f  x%.2f\n", names[v], best[v] * scale,
            best[v] / best[0]);
    }
}

static void bench_s2(void)
{
    /* Sgr A*: GM/c^2 = 6.3495e9 m, D = 8.28 kpc, GM/c^3 = 21.18 s. */
    const double gravitational_radius = 1476.625 * 4.3e6;
    const double distance = 8.28e3 * 3.0856775814913673e16;
    const double angular_scale = gravitational_radius / distance;
    const double arcsec = PI / (180.0 * 3600.0);
    const double semi_major_axis = 0.125 * arcsec / angular_scale;
    const double year = 3.15576e7 / (gravitational_radius / 299792458.0);
    const double degree = PI / 180.0;
    const double identity[3][3] = {
        {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}
    };
    const rp_kerr_mcmc_observation observation = {
        angular_scale, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    };
    const rp_integrator_config solver = {
        RP_INTEGRATOR_METHOD_DOP853, 8U, 1e-9, 1e-12, 0.0, 0.0, 100000U,
        NULL, NULL, NULL, NULL, NULL
    };
    /* An extended mass of 1e-3 M inside a few S2 semi-major axes. */
    const test_plummer plummer = {1e-3, 2.0 * semi_major_axis, {0.0, 0.0, 0.0}};
    const test_miyamoto_nagai disk = {1e-3, semi_major_axis, 0.01 * semi_major_axis};
    const rp_external_force plummer_force = {test_plummer_acceleration, &plummer};
    const rp_external_force disk_force = {test_miyamoto_nagai_acceleration, &disk};
    const rp_external_force *forces[3] = {NULL, &plummer_force, &disk_force};
    const char *names[3] = {"NULL force", "Plummer", "Miyamoto-Nagai"};
    double times[EPOCHS];
    double output_h[3 * EPOCHS];
    double output_f[3 * EPOCHS];
    rp_integrator_stats stats;
    double best;
    double base;
    size_t repeat;
    size_t k;
    size_t f;

    for (k = 0U; k < EPOCHS; ++k) {
        times[k] = (-26.0 + 30.0 * (double)k / (double)(EPOCHS - 1U)) * year;
    }
    printf("S2-like evaluation: a = %.0f M, e = 0.884, %u epochs over 30 yr,"
        " dop853 rtol 1e-9 atol 1e-12, spin 0.5\n",
        semi_major_axis, EPOCHS);
    best = INFINITY;
    for (repeat = 0U; repeat < REPEATS; ++repeat) {
        const double start = now_seconds();
        assert(rp_kerr_mcmc_evaluate(
            &solver, 0.5, identity, semi_major_axis, 0.884,
            134.0 * degree, 228.0 * degree, 66.0 * degree,
            times, EPOCHS, &observation, output_h, &stats
        ) == RP_KERR_MCMC_OK);
        best = fmin(best, now_seconds() - start);
    }
    base = best;
    printf("  rp_kerr_mcmc_evaluate          %8.3f ms  x1.00  (%zu steps, %zu RHS)\n",
        1e3 * best, stats.accepted_steps, stats.rhs_evaluations);
    for (f = 0U; f < 3U; ++f) {
        double largest = 0.0;
        best = INFINITY;
        for (repeat = 0U; repeat < REPEATS; ++repeat) {
            const double start = now_seconds();
            assert(rp_kerr_forced_mcmc_evaluate(
                &solver, 0.5, forces[f], identity, semi_major_axis, 0.884,
                134.0 * degree, 228.0 * degree, 66.0 * degree,
                times, EPOCHS, &observation, output_f, &stats
            ) == RP_KERR_MCMC_OK);
            best = fmin(best, now_seconds() - start);
        }
        for (k = 0U; k < 3U * EPOCHS; ++k) {
            largest = fmax(largest, fabs(output_f[k] - output_h[k])
                / (k % 3U == 2U ? 1.0 : arcsec));
        }
        printf("  forced, %-15s        %8.3f ms  x%.2f  (%zu steps, %zu RHS;"
            " max |diff| %.2e arcsec or v/c)\n",
            names[f], 1e3 * best, best / base, stats.accepted_steps,
            stats.rhs_evaluations, largest);
    }
}

int main(void)
{
    bench_rhs();
    bench_s2();
    return 0;
}
