/*
 * Timing of the thin-disk callback and table build (compile with -O3).
 * The homogeneous Kerr RHS is timed on the same path as a reference.
 * Reports the minimum over repetitions; fails only on a gross regression.
 */
#define _POSIX_C_SOURCE 199309L
#include "force/external_force.h"
#include "force/models/thin_disk.h"
#include "geodesic/physic/kerr.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define PI 3.14159265358979323846
#define POINTS 100000
#define REPEATS 15

static double now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + 1.0e-9 * (double)t.tv_nsec;
}

static double time_force(const rp_external_force *force, const double *x, double *sink)
{
    const double v[3] = {0.0, 0.0, 0.0};
    double best = INFINITY;
    for (int rep = 0; rep < REPEATS; ++rep) {
        double a[3];
        const double t0 = now();
        for (int i = 0; i < POINTS; ++i) {
            force->evaluate(x + 3 * i, v, force->model, a);
            *sink += a[0];
        }
        best = fmin(best, (now() - t0) / POINTS * 1.0e9);
    }
    return best;
}

int main(void)
{
    rp_thin_disk_params params;
    rp_thin_disk_default_params(&params);
    rp_axisymmetric_table *table = NULL;
    const double t0 = now();
    if (rp_thin_disk_table_build(&params, &table) != RP_AXISYMMETRIC_STATUS_OK) {
        printf("FAIL table build\n");
        return EXIT_FAILURE;
    }
    const double build = now() - t0;
    printf("table build: %.2f s, %zu x %zu nodes, %.1f MB\n", build,
           table->radius.count, table->cosine.count, (double)table->bytes / 1.0e6);

    /* Piran et al. fiducial placement: R_out = 31740, M_d = 1e3 / 4.3e6. */
    const double normal[3] = {0.3 / 1.0, -0.5, 0.8};
    double unit[3];
    const double n = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
    for (int k = 0; k < 3; ++k) {
        unit[k] = normal[k] / n;
    }
    rp_thin_disk_model model;
    rp_thin_disk_model_init(&model, table, 1.0e3 / 4.3e6, 31740.0, unit);
    const rp_external_force force = {rp_thin_disk_acceleration, &model};

    double *x = malloc(3 * POINTS * sizeof *x);
    struct rp_kerr_geodesic_state *states = malloc(POINTS * sizeof *states);
    if (x == NULL || states == NULL) {
        printf("FAIL allocation\n");
        return EXIT_FAILURE;
    }
    /* An S301-like ellipse (a = 16010, e = 0.98), sampled in eccentric
     * anomaly as an integrator would visit it. */
    for (int i = 0; i < POINTS; ++i) {
        const double ecc = 2.0 * PI * i / POINTS;
        const double a = 16010.0, e = 0.98;
        const double p = a * (cos(ecc) - e);
        const double q = a * sqrt(1.0 - e * e) * sin(ecc);
        x[3 * i] = 0.8 * p;
        x[3 * i + 1] = q;
        x[3 * i + 2] = 0.6 * p;
    }
    double sink = 0.0;
    const double orbit = time_force(&force, x, &sink);

    /* Kerr reference on the same path: prepare + corrected RHS. */
    for (int i = 0; i < POINTS; ++i) {
        const double *p = x + 3 * i;
        const double r = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
        states[i].x[0] = 0.0;
        states[i].x[1] = r;
        states[i].x[2] = acos(p[2] / r);
        states[i].x[3] = atan2(p[1], p[0]);
        states[i].u[0] = 1.0;
        states[i].u[1] = 1.0e-3;
        states[i].u[2] = 1.0e-6;
        states[i].u[3] = 1.0e-5;
    }
    double kerr = INFINITY;
    for (int rep = 0; rep < REPEATS; ++rep) {
        double dx[4], du[4];
        const double t1 = now();
        for (int i = 0; i < POINTS; ++i) {
            rp_kerr_geodesic_physic_prepare_state(1.0, 0.7, &states[i]);
            rp_kerr_geodesic_physic_evaluate_corrected_rhs(1.0, 0.7, &states[i], dx, du);
            sink += du[1];
        }
        kerr = fmin(kerr, (now() - t1) / POINTS * 1.0e9);
    }
    double legacy = INFINITY;
    for (int rep = 0; rep < REPEATS; ++rep) {
        double dx[4], du[4];
        const double t1 = now();
        for (int i = 0; i < POINTS; ++i) {
            rp_kerr_geodesic_physic_prepare_state(1.0, 0.7, &states[i]);
            rp_kerr_geodesic_physic_evaluate_rhs(1.0, 0.7, &states[i], dx, du);
            sink += du[1];
        }
        legacy = fmin(legacy, (now() - t1) / POINTS * 1.0e9);
    }

    /* Random points from 100 r_g to 3 R_out: cache-cold worst case. */
    unsigned long long s = 88172645463325252ULL;
    for (int i = 0; i < POINTS; ++i) {
        double u[3];
        for (int k = 0; k < 3; ++k) {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            u[k] = (double)(s >> 11) * 0x1.0p-53;
        }
        const double r = 100.0 * exp(u[0] * log(3.0 * 31740.0 / 100.0));
        const double c = 2.0 * u[1] - 1.0, st = sqrt(1.0 - c * c), ph = 2.0 * PI * u[2];
        x[3 * i] = r * st * cos(ph);
        x[3 * i + 1] = r * st * sin(ph);
        x[3 * i + 2] = r * c;
    }
    const double random = time_force(&force, x, &sink);

    printf("thin disk callback: %.1f ns/call along the orbit, %.1f ns/call at random points\n",
           orbit, random);
    printf("Kerr prepare + RHS on the same orbit: corrected %.1f ns/call, legacy optimized %.1f ns/call\n",
           kerr, legacy);
    printf("(checksum %.3e)\n", sink);
    free(x);
    free(states);
    rp_thin_disk_table_free(table);
    /* Generous bounds: a shared machine, but a 4x regression is a bug. */
    if (orbit > 4.0 * 40.0 || build > 30.0) {
        printf("FAIL timing regression\n");
        return EXIT_FAILURE;
    }
    printf("test_thin_disk_bench: ok\n");
    return EXIT_SUCCESS;
}
