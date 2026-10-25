/* Tabulated thin disk against the direct ring sum, and the placed model. */
#define _POSIX_C_SOURCE 199309L
#include "force/external_force.h"
#include "force/models/thin_disk.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int failures = 0;
static unsigned long long rng_state = 88172645463325252ULL;

static double uniform(void)
{
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return (double)(rng_state >> 11) * 0x1.0p-53;
}

static double now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + 1.0e-9 * (double)t.tv_nsec;
}

static void check(const char *label, double error, double tolerance)
{
    if (!(error <= tolerance)) {
        printf("FAIL %s: %.3e > %.1e\n", label, error, tolerance);
        ++failures;
    }
}

static int compare_doubles(const void *a, const void *b)
{
    const double x = *(const double *)a;
    const double y = *(const double *)b;
    return (x > y) - (x < y);
}

static double norm3(const double v[3])
{
    return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

static double distance3(const double a[3], const double b[3])
{
    const double d[3] = {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
    return norm3(d);
}

/* ------------------------------------------------------------------ */

enum { SET_GLOBAL, SET_PLANE, SET_EDGE, SET_INNER_EDGE, SET_CORE, SET_FAR, SET_COUNT };
static const char *set_names[SET_COUNT] = {
    "global 1e-4..6", "near plane mu<0.05", "outer edge", "inner edge",
    "inside r_min", "blend and far 3.5..9",
};

/* A random disk-frame point of one test region (lengths in R_out). */
static void sample(int set, double inner, double r_min, double *r, double *mu)
{
    switch (set) {
    case SET_GLOBAL:
        *r = exp(log(1.0e-4) + uniform() * log(6.0e4));
        *mu = uniform();
        break;
    case SET_PLANE:
        *r = exp(log(1.0e-4) + uniform() * log(1.0e4));
        *mu = 0.05 * uniform();
        break;
    case SET_EDGE:
        *r = 0.9 + 0.2 * uniform();
        *mu = 0.05 * uniform();
        break;
    case SET_INNER_EDGE:
        *r = inner * (0.9 + 0.2 * uniform());
        *mu = 0.05 * uniform();
        break;
    case SET_CORE:
        *r = r_min * exp(-uniform() * log(100.0));
        *mu = uniform();
        break;
    default:
        *r = 3.5 + 5.5 * uniform();
        *mu = uniform();
        break;
    }
}

/*
 * Table versus direct sum through the placed model (unit mass and radius,
 * normal e_z, random azimuth).  Reports median / p99 / max relative errors
 * of |a| and Phi per region and returns the worst acceleration error.
 */
static double accuracy(
    const char *label, const rp_thin_disk_params *params,
    const rp_axisymmetric_table *table, int points, const double tolerance[SET_COUNT]
)
{
    rp_thin_disk_model model;
    const double normal[3] = {0.0, 0.0, 1.0};
    rp_thin_disk_model_init(&model, table, 1.0, 1.0, normal);
    double *ea = malloc((size_t)points * sizeof *ea);
    double *ep = malloc((size_t)points * sizeof *ep);
    double worst = 0.0;
    printf("%s (%zu x %zu nodes):\n", label, table->radius.count, table->cosine.count);
    for (int set = 0; set < SET_COUNT; ++set) {
        if (set == SET_INNER_EDGE && params->inner_ratio == 0.0) {
            continue;
        }
        for (int k = 0; k < points; ++k) {
            double r, mu;
            sample(set, params->inner_ratio, table->radial_min, &r, &mu);
            const double phi = 6.283185307179586 * uniform();
            const double d = r * sqrt((1.0 - mu) * (1.0 + mu));
            const double sign = uniform() < 0.5 ? -1.0 : 1.0;
            const double x[3] = {d * cos(phi), d * sin(phi), sign * r * mu};
            double f[3];
            rp_thin_disk_direct_field(params, d, x[2], f);
            const double want[3] = {f[0] * x[0], f[0] * x[1], f[1] * x[2]};
            double got[3], potential;
            const double v[3] = {0.0, 0.0, 0.0};
            rp_thin_disk_acceleration(x, v, &model, got);
            rp_thin_disk_potential(x, &model, &potential);
            ea[k] = distance3(got, want) / norm3(want);
            ep[k] = fabs(potential - f[2]) / fabs(f[2]);
        }
        qsort(ea, (size_t)points, sizeof *ea, compare_doubles);
        qsort(ep, (size_t)points, sizeof *ep, compare_doubles);
        const int p99 = (int)(0.99 * points);
        printf("  %-22s |a|: median %.1e p99 %.1e max %.1e   Phi: median %.1e max %.1e\n",
               set_names[set], ea[points / 2], ea[p99], ea[points - 1],
               ep[points / 2], ep[points - 1]);
        char name[96];
        snprintf(name, sizeof name, "%s, %s", label, set_names[set]);
        check(name, ea[points - 1], tolerance[set]);
        check(name, ep[points - 1], tolerance[set]);
        worst = fmax(worst, ea[points - 1]);
    }
    free(ea);
    free(ep);
    return worst;
}

/* ------------------------------------------------------------------ */

static void random_unit(double n[3])
{
    const double c = 2.0 * uniform() - 1.0;
    const double s = sqrt(1.0 - c * c);
    const double phi = 6.283185307179586 * uniform();
    n[0] = s * cos(phi);
    n[1] = s * sin(phi);
    n[2] = c;
}

static void random_rotation(double q[3][3])
{
    double n[3];
    random_unit(n);
    rp_thin_disk_rotation_from_normal(n, q);
    /* Add a random spin about the normal so that q is generic. */
    const double t = 6.283185307179586 * uniform();
    double rows[2][3];
    for (int k = 0; k < 3; ++k) {
        rows[0][k] = cos(t) * q[0][k] + sin(t) * q[1][k];
        rows[1][k] = -sin(t) * q[0][k] + cos(t) * q[1][k];
    }
    memcpy(q[0], rows[0], sizeof rows[0]);
    memcpy(q[1], rows[1], sizeof rows[1]);
}

static void apply(double q[3][3], const double x[3], double out[3])
{
    for (int i = 0; i < 3; ++i) {
        out[i] = q[i][0] * x[0] + q[i][1] * x[1] + q[i][2] * x[2];
    }
}

/* A random black-hole-frame point for a disk of outer radius L. */
static void random_point(double length, double x[3])
{
    double n[3];
    random_unit(n);
    const double r = length * exp(log(1.0e-4) + uniform() * log(1.0e5));
    for (int k = 0; k < 3; ++k) {
        x[k] = r * n[k];
    }
}

static void placed_model_properties(const rp_axisymmetric_table *table)
{
    const double v[3] = {0.0, 0.0, 0.0};
    const double length = 31740.0;
    const double mass = 2.3e-4;
    double normal[3];
    random_unit(normal);
    rp_thin_disk_model model, unit, heavy;
    rp_thin_disk_model_init(&model, table, mass, length, normal);
    rp_thin_disk_model_init(&unit, table, 1.0, 1.0, normal);
    rp_thin_disk_model_init(&heavy, table, 3.7 * mass, length, normal);
    const rp_external_force force = {rp_thin_disk_acceleration, &model};

    double symmetry = 0.0, scaling = 0.0, length_scaling = 0.0, rotation = 0.0;
    for (int k = 0; k < 2000; ++k) {
        double x[3], a[3], b[3];
        random_point(length, x);
        force.evaluate(x, v, force.model, a);
        /* z -> -z: reflect through the disk plane. */
        const double h = x[0] * normal[0] + x[1] * normal[1] + x[2] * normal[2];
        double xr[3], ar[3];
        for (int i = 0; i < 3; ++i) {
            xr[i] = x[i] - 2.0 * h * normal[i];
        }
        rp_thin_disk_acceleration(xr, v, &model, b);
        const double hb = b[0] * normal[0] + b[1] * normal[1] + b[2] * normal[2];
        for (int i = 0; i < 3; ++i) {
            ar[i] = b[i] - 2.0 * hb * normal[i];
        }
        symmetry = fmax(symmetry, distance3(a, ar) / norm3(a));
        /* Linear mass scaling. */
        rp_thin_disk_acceleration(x, v, &heavy, b);
        for (int i = 0; i < 3; ++i) {
            b[i] /= 3.7;
        }
        scaling = fmax(scaling, distance3(a, b) / norm3(a));
        /* a(x; M, L) = M / L^2 a_unit(x / L). */
        const double xu[3] = {x[0] / length, x[1] / length, x[2] / length};
        rp_thin_disk_acceleration(xu, v, &unit, b);
        for (int i = 0; i < 3; ++i) {
            b[i] *= mass / (length * length);
        }
        length_scaling = fmax(length_scaling, distance3(a, b) / norm3(a));
        /* Rotation equivariance: a(Qx; Qn) = Q a(x; n). */
        double q[3][3], qx[3], qn[3], qa[3];
        random_rotation(q);
        apply(q, x, qx);
        apply(q, normal, qn);
        rp_thin_disk_model rotated;
        rp_thin_disk_model_init(&rotated, table, mass, length, qn);
        rp_thin_disk_acceleration(qx, v, &rotated, b);
        apply(q, a, qa);
        rotation = fmax(rotation, distance3(b, qa) / norm3(qa));
    }
    printf("z -> -z %.1e, mass scaling %.1e, length scaling %.1e, rotation %.1e\n",
           symmetry, scaling, length_scaling, rotation);
    check("z -> -z symmetry", symmetry, 1.0e-13);
    check("mass scaling", scaling, 1.0e-15);
    check("length scaling", length_scaling, 1.0e-13);
    check("rotation equivariance", rotation, 1.0e-12);

    /* Monopole far field: |a - a_mono| / |a| ~ <R^2> / r^2 = 1/3 (L / r)^2. */
    double worst_mono = 0.0;
    for (int k = 0; k < 200; ++k) {
        double n[3], x[3], a[3];
        random_unit(n);
        const double r = 1.0e3 * length;
        for (int i = 0; i < 3; ++i) {
            x[i] = r * n[i];
        }
        rp_thin_disk_acceleration(x, v, &model, a);
        const double mono[3] = {-mass * n[0] / (r * r), -mass * n[1] / (r * r), -mass * n[2] / (r * r)};
        worst_mono = fmax(worst_mono, distance3(a, mono) / norm3(mono));
    }
    printf("monopole at 1e3 R_out: max rel %.2e (quadrupole bound 1e-6)\n", worst_mono);
    check("monopole far field", worst_mono, 1.0e-6);

    /* -grad Phi against the interpolated acceleration. */
    double *errors = malloc(2000 * sizeof *errors);
    for (int k = 0; k < 2000; ++k) {
        double x[3], a[3], g[3];
        random_point(length, x);
        rp_thin_disk_acceleration(x, v, &model, a);
        const double h = 1.0e-6 * norm3(x);
        for (int i = 0; i < 3; ++i) {
            double xp[3] = {x[0], x[1], x[2]}, xm[3] = {x[0], x[1], x[2]};
            double pp, pm;
            xp[i] += h;
            xm[i] -= h;
            rp_thin_disk_potential(xp, &model, &pp);
            rp_thin_disk_potential(xm, &model, &pm);
            g[i] = -(pp - pm) / (2.0 * h);
        }
        errors[k] = distance3(g, a) / norm3(a);
    }
    qsort(errors, 2000, sizeof *errors, compare_doubles);
    /* Phi and a are separate C2 interpolants: their mismatch is the O(h^3)
     * derivative error of the Phi spline, largest inside the softened
     * plane (mu < eps). */
    printf("-grad Phi (central differences) vs a: median %.1e max %.1e\n",
           errors[1000], errors[1999]);
    check("potential gradient consistency", errors[1999], 2.0e-5);
    check("potential gradient consistency (median)", errors[1000], 2.0e-6);
    free(errors);

    /* Invalid input. */
    double a[3], p;
    const double origin[3] = {0.0, 0.0, 0.0};
    const double bad[3] = {NAN, 1.0, 0.0};
    const double tilted[3] = {0.0, 0.0, 1.1};
    if (rp_thin_disk_acceleration(origin, v, &model, a) == 0
        || rp_thin_disk_acceleration(bad, v, &model, a) == 0
        || rp_thin_disk_acceleration(bad, v, NULL, a) == 0
        || rp_thin_disk_potential(origin, &model, &p) == 0
        || rp_thin_disk_model_init(&unit, table, 1.0, 1.0, tilted) != RP_AXISYMMETRIC_STATUS_INVALID_ARGUMENT
        || rp_thin_disk_model_init(&unit, table, 1.0, -1.0, normal) != RP_AXISYMMETRIC_STATUS_INVALID_ARGUMENT
        || rp_thin_disk_model_init(&unit, NULL, 1.0, 1.0, normal) != RP_AXISYMMETRIC_STATUS_NULL_POINTER) {
        printf("FAIL invalid input accepted\n");
        ++failures;
    }
}

/* Continuity across the region joins and C1/C2 across grid nodes. */
static void smoothness(const rp_axisymmetric_table *table)
{
    double jump = 0.0;
    const double joins[3] = {table->radial_min, table->radial_match, table->radial_far};
    for (int j = 0; j < 3; ++j) {
        for (int k = 0; k <= 20; ++k) {
            const double mu = k / 20.0;
            double lo[3], hi[3];
            rp_axisymmetric_table_fields(table, joins[j] * (1.0 - 1.0e-13), mu, lo);
            rp_axisymmetric_table_fields(table, joins[j] * (1.0 + 1.0e-13), mu, hi);
            for (int q = 0; q < 3; ++q) {
                jump = fmax(jump, fabs(hi[q] - lo[q]) / fmax(fabs(hi[q]), 1.0e-300));
            }
        }
    }
    printf("jumps at r_min, r_match, r_far: max rel %.1e\n", jump);
    check("continuity at region joins", jump, 1.0e-11);

    /* One-sided first and second differences on both sides of interior
     * nodes, in r and in mu, relative to the local scale of A and C. */
    double c1 = 0.0, c2 = 0.0;
    const rp_axisymmetric_axis *ar = &table->radius;
    const rp_axisymmetric_axis *am = &table->cosine;
    for (int k = 0; k < 400; ++k) {
        const int along_r = k % 2 == 0;
        const rp_axisymmetric_axis *axis = along_r ? ar : am;
        const size_t i = 2U + (size_t)(uniform() * (double)(axis->count - 4U));
        const double h = fmin(axis->node[i + 1U] - axis->node[i], axis->node[i] - axis->node[i - 1U]);
        const double delta = 1.0e-3 * h;
        double r = along_r ? ar->node[i] : exp(log(1.0e-3) + uniform() * log(2.0e3));
        double mu = along_r ? uniform() : am->node[i];
        if (r >= table->radial_match) {
            r = 0.5 * table->radial_match;
        }
        double f[5][3];
        for (int s = -2; s <= 2; ++s) {
            rp_axisymmetric_table_fields(table, along_r ? r + s * delta : r,
                                         along_r ? mu : mu + s * delta, f[s + 2]);
        }
        for (int q = 0; q < 2; ++q) {
            const double scale = fabs(f[2][q]);
            const double left1 = (f[2][q] - f[1][q]) / delta;
            const double right1 = (f[3][q] - f[2][q]) / delta;
            const double left2 = (f[2][q] - 2.0 * f[1][q] + f[0][q]) / (delta * delta);
            const double right2 = (f[4][q] - 2.0 * f[3][q] + f[2][q]) / (delta * delta);
            /* A C1 jump would be O(scale / h); a C2 jump O(scale / h^2). */
            c1 = fmax(c1, fabs(right1 - left1) * h / scale);
            c2 = fmax(c2, fabs(right2 - left2) * h * h / scale);
        }
    }
    printf("across nodes: first-derivative jump %.1e, second-derivative jump %.1e (scaled)\n", c1, c2);
    check("C1 across nodes", c1, 1.0e-3);
    check("C2 across nodes", c2, 1.0e-2);
}

int main(void)
{
    rp_thin_disk_params params;
    rp_thin_disk_default_params(&params);
    rp_axisymmetric_table *table = NULL;
    const double t0 = now();
    if (rp_thin_disk_table_build(&params, &table) != RP_AXISYMMETRIC_STATUS_OK) {
        printf("FAIL default table build\n");
        return EXIT_FAILURE;
    }
    printf("default table: build %.2f s, %.1f MB\n", now() - t0, (double)table->bytes / 1.0e6);
    const double tolerance[SET_COUNT] = {1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-7};
    accuracy("alpha=-1, R_in=0, eps=0.01", &params, table, 3000, tolerance);
    placed_model_properties(table);
    smoothness(table);
    rp_thin_disk_table_free(table);

    /* Other profiles of Piran et al.: steeper, flat with a hole. */
    static const double configs[][3] = {{-1.5, 0.0, 0.01}, {0.0, 0.1, 0.01}, {-1.0, 0.0, 0.03}};
    for (size_t c = 0; c < sizeof configs / sizeof configs[0]; ++c) {
        rp_thin_disk_default_params(&params);
        params.alpha = configs[c][0];
        params.inner_ratio = configs[c][1];
        params.softening = configs[c][2];
        const double t1 = now();
        if (rp_thin_disk_table_build(&params, &table) != RP_AXISYMMETRIC_STATUS_OK) {
            printf("FAIL table build for alpha=%g\n", params.alpha);
            ++failures;
            continue;
        }
        char label[64];
        snprintf(label, sizeof label, "alpha=%g, R_in=%g, eps=%g", params.alpha,
                 params.inner_ratio, params.softening);
        printf("build %.2f s\n", now() - t1);
        accuracy(label, &params, table, 500, tolerance);
        rp_thin_disk_table_free(table);
    }

    /* Invalid parameters. */
    rp_thin_disk_default_params(&params);
    params.softening = 0.0;
    if (rp_thin_disk_table_build(&params, &table) != RP_AXISYMMETRIC_STATUS_INVALID_ARGUMENT
        || table != NULL) {
        printf("FAIL invalid parameters accepted\n");
        ++failures;
    }

    if (failures != 0) {
        printf("test_thin_disk_table: %d failure(s)\n", failures);
        return EXIT_FAILURE;
    }
    printf("test_thin_disk_table: ok\n");
    return EXIT_SUCCESS;
}
