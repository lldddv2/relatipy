/* Oblate power-law spheroid: direct homoeoid sum and its table. */
#define _POSIX_C_SOURCE 199309L
#include "force/models/spheroid.h"
#include "force/models/quadrature.h"
#include "force/models/thin_disk.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define PI 3.14159265358979323846

static int failures = 0;
static unsigned long long rng_state = 2463534242ULL;

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

/* Relative error of the field vectors (a_d, a_z) and of Phi. */
static double field_error(const double f[3], const double g[3], double d, double z)
{
    const double norm = hypot(g[0] * d, g[1] * z);
    const double ea = norm > 0.0 ? hypot((f[0] - g[0]) * d, (f[1] - g[1]) * z) / norm : 0.0;
    return fmax(ea, fabs(f[2] - g[2]) / fabs(g[2]));
}

/* q = 1: the spherical power law, G M(<r) / r^2 and its potential. */
static void sphere(void)
{
    static const double configs[][2] = {{-2.0, 0.0}, {-1.5, 0.2}, {0.0, 0.0}, {-2.9, 0.0}};
    double worst = 0.0;
    for (size_t c = 0; c < sizeof configs / sizeof configs[0]; ++c) {
        rp_spheroid_params params;
        rp_spheroid_default_params(&params);
        params.alpha = configs[c][0];
        params.inner_ratio = configs[c][1];
        const double b = params.alpha + 3.0;
        const double s_in = params.inner_ratio;
        const double total = 1.0 - pow(s_in, b);
        for (int k = 0; k < 200; ++k) {
            const double r = exp(log(1.0e-4) + uniform() * log(3.0e4));
            const double mu = uniform();
            const double d = r * sqrt(1.0 - mu * mu), z = r * mu;
            double f[3];
            rp_spheroid_direct_field(&params, d, z, f);
            const double inside = fmin(fmax(r, s_in), 1.0);
            const double enclosed = (pow(inside, b) - pow(s_in, b)) / total;
            /* Phi = -M(<r)/r - 4 pi rho_0 int_max(r, s_in)^1 s^(alpha+1) ds,
             * with 4 pi rho_0 = b / total. */
            const double lower = fmax(r, s_in);
            const double c2 = params.alpha + 2.0;
            const double shell = lower < 1.0
                ? (b / total) * (fabs(c2) < 1e-14 ? -log(lower) : -expm1(c2 * log(lower)) / c2)
                : 0.0;
            const double want[3] = {-enclosed / (r * r * r), -enclosed / (r * r * r),
                                    -enclosed / r - shell};
            if (enclosed > 0.0) {
                worst = fmax(worst, field_error(f, want, d, z));
            } else {
                worst = fmax(worst, fabs(f[2] / want[2] - 1.0) + fabs(f[0]) + fabs(f[1]));
            }
        }
    }
    printf("q = 1 vs spherical power law G M(<r)/r^2 and Phi: max rel %.1e\n", worst);
    check("sphere: direct sum vs analytic", worst, 2.0e-14);
}

static void quadrature_and_mass(void)
{
    static const double qs[] = {1.0, 0.5, 0.1, 0.01};
    static const double points[][2] = {
        {0.3, 0.0}, {0.3, 0.001}, {0.3, 0.05}, {1.0, 0.0}, {1.001, 0.0},
        {0.5, 0.003}, {2.0, 0.1}, {0.0, 0.05}, {0.9999, 5e-4}, {1e-4, 1e-6},
    };
    double worst = 0.0, mass = 0.0;
    for (size_t iq = 0; iq < sizeof qs / sizeof qs[0]; ++iq) {
        for (int profile = 0; profile < 2; ++profile) {
            rp_spheroid_params base, fine;
            rp_spheroid_default_params(&base);
            base.axis_ratio = qs[iq];
            base.alpha = profile == 0 ? -2.0 : -1.3;
            base.inner_ratio = profile == 0 ? 0.0 : 0.05;
            fine = base;
            fine.quadrature_points = 48U;
            for (size_t k = 0; k < sizeof points / sizeof points[0]; ++k) {
                double f[3], g[3];
                rp_spheroid_direct_field(&base, points[k][0], points[k][1], f);
                rp_spheroid_direct_field(&fine, points[k][0], points[k][1], g);
                worst = fmax(worst, field_error(f, g, points[k][0], points[k][1]));
            }
            /* Unit mass: monopole at 1e5 s_out, quadrupole ~ 1e-10. */
            double f[3];
            const double far = 1.0e5;
            rp_spheroid_direct_field(&base, 0.6 * far, 0.8 * far, f);
            mass = fmax(mass, fabs(-f[2] * far - 1.0));
        }
    }
    printf("direct sum, 16 vs 48 Gauss points: max rel %.1e; unit mass to %.1e\n", worst, mass);
    check("spheroid quadrature convergence", worst, 1.0e-13);
    check("spheroid unit mass", mass, 1.0e-9);
}

/*
 * q -> 0 at fixed mass: a razor-thin disk with the projected density
 * Sigma(R) = 2 q rho_0 int_0^sqrt(1 - R^2) (R^2 + t^2)^(alpha/2) dt, which for
 * alpha = -2 is arccos(R) / (2 pi R).  Its field from unsoftened rings,
 * R = cos(phi): int_0^(pi/2) phi sin(phi) ring(cos phi) dphi.
 */
static void projected_disk(double d, double z, double out[3])
{
    rp_gauss_rule rule;
    rp_gauss_legendre(16U, &rule);
    const int panels = 64;
    out[0] = out[1] = out[2] = 0.0;
    for (int p = 0; p < panels; ++p) {
        const double lo = 0.5 * PI * p / panels, hi = 0.5 * PI * (p + 1) / panels;
        for (size_t g = 0; g < rule.count; ++g) {
            const double phi = 0.5 * (lo + hi) + 0.5 * (hi - lo) * rule.node[g];
            const double w = 0.5 * (hi - lo) * rule.weight[g] * phi * sin(phi);
            double ring[3];
            rp_thin_disk_ring_field(cos(phi), 0.0, d, z, ring);
            for (int k = 0; k < 3; ++k) {
                out[k] += w * ring[k];
            }
        }
    }
}

static void thin_limit(void)
{
    static const double points[][2] = {{0.3, 0.2}, {0.7, 0.1}, {1.2, 0.3}, {0.0, 0.5}, {0.5, 0.5}};
    static const double qs[] = {0.03, 0.01, 0.003};
    double errors[3] = {0.0, 0.0, 0.0};
    double power_law = 0.0;
    rp_thin_disk_params disk;
    rp_thin_disk_default_params(&disk);
    for (size_t k = 0; k < sizeof points / sizeof points[0]; ++k) {
        const double d = points[k][0], z = points[k][1];
        double sigma_disk[3], fiducial[3];
        projected_disk(d, z, sigma_disk);
        rp_thin_disk_direct_field(&disk, d, z, fiducial);
        power_law = fmax(power_law, field_error(fiducial, sigma_disk, d, z));
        for (size_t iq = 0; iq < 3; ++iq) {
            rp_spheroid_params params;
            rp_spheroid_default_params(&params);
            params.axis_ratio = qs[iq];
            double f[3];
            rp_spheroid_direct_field(&params, d, z, f);
            errors[iq] = fmax(errors[iq], field_error(f, sigma_disk, d, z));
        }
    }
    printf("q -> 0 vs the projected thin disk Sigma = acos(R)/(2 pi R): q=0.03 %.1e, "
           "0.01 %.1e, 0.003 %.1e (rate ~q^2)\n", errors[0], errors[1], errors[2]);
    printf("  (projected disk vs the fiducial Sigma ~ 1/R disk at the same mass: %.1e)\n",
           power_law);
    check("thin limit at q = 0.003", errors[2], 5.0e-5);
    const double rate = log(errors[1] / errors[2]) / log(0.01 / 0.003);
    printf("  convergence order %.2f\n", rate);
    check("thin limit is second order", fabs(rate - 2.0), 0.3);
}

static double table_accuracy(const rp_spheroid_params *params, int points)
{
    rp_axisymmetric_table *table = NULL;
    const double t0 = now();
    if (rp_spheroid_table_build(params, &table) != RP_AXISYMMETRIC_STATUS_OK) {
        printf("FAIL spheroid table build\n");
        ++failures;
        return INFINITY;
    }
    printf("q=%g alpha=%g s_in=%g: %zu x %zu nodes, %.1f MB, build %.2f s\n",
           params->axis_ratio, params->alpha, params->inner_ratio, table->radius.count,
           table->cosine.count, (double)table->bytes / 1.0e6, now() - t0);
    rp_spheroid_model model;
    const double normal[3] = {0.0, 0.0, 1.0};
    rp_spheroid_model_init(&model, table, 1.0, 1.0, normal);
    const double q = params->axis_ratio;
    static const char *names[] = {"global", "near plane", "near surface", "far"};
    double *errors = malloc((size_t)points * sizeof *errors);
    double worst = 0.0;
    for (int set = 0; set < 4; ++set) {
        for (int k = 0; k < points; ++k) {
            /* Sample in the grid coordinates rho, nu. */
            double rho, nu;
            if (set == 0) {
                rho = exp(log(1.0e-4) + uniform() * log(4.0e4));
                nu = uniform();
            } else if (set == 1) {
                rho = exp(log(1.0e-3) + uniform() * log(2.0e3));
                nu = 0.05 * uniform();
            } else if (set == 2) {
                rho = 0.95 + 0.1 * uniform();
                nu = uniform();
            } else {
                rho = (3.5 + 5.0 * uniform()) / q;
                nu = uniform();
            }
            /* nu is the grid's mu = |z| / r; r = rho / sqrt(1 + (1/q^2 - 1) mu^2). */
            const double r = rho / sqrt(1.0 + (1.0 / (q * q) - 1.0) * nu * nu);
            const double d = r * sqrt(1.0 - nu * nu);
            const double z = r * nu;
            const double x[3] = {d, 0.0, uniform() < 0.5 ? z : -z};
            const double v[3] = {0.0, 0.0, 0.0};
            double g[3], a[3], potential;
            rp_spheroid_direct_field(params, d, z, g);
            rp_spheroid_acceleration(x, v, &model, a);
            rp_spheroid_potential(x, &model, &potential);
            const double want[3] = {g[0] * d, 0.0, g[1] * x[2]};
            /* Inside an empty core the force vanishes: use |Phi| / r as scale. */
            const double norm = fmax(sqrt(want[0] * want[0] + want[2] * want[2]),
                                     1.0e-6 * fabs(g[2]) / r);
            errors[k] = fmax(sqrt((a[0] - want[0]) * (a[0] - want[0]) + a[1] * a[1]
                                  + (a[2] - want[2]) * (a[2] - want[2])) / norm,
                             fabs(potential / g[2] - 1.0));
        }
        qsort(errors, (size_t)points, sizeof *errors, compare_doubles);
        printf("  %-13s |a|, Phi: median %.1e p99 %.1e max %.1e\n", names[set],
               errors[points / 2], errors[(int)(0.99 * points)], errors[points - 1]);
        worst = fmax(worst, errors[points - 1]);
    }
    free(errors);

    /* z -> -z and the C0 kink: no jump across the surface. */
    double symmetry = 0.0, jump = 0.0;
    for (int k = 0; k < 200; ++k) {
        const double nu = uniform();
        const double v[3] = {0.0, 0.0, 0.0};
        double a[3], b[3];
        const double rho = exp(log(1.0e-3) + uniform() * log(3.0e3));
        const double rr = rho / sqrt(1.0 + (1.0 / (q * q) - 1.0) * nu * nu);
        const double x[3] = {rr * sqrt(1.0 - nu * nu), 0.0, rr * nu};
        const double y[3] = {x[0], 0.0, -x[2]};
        rp_spheroid_acceleration(x, v, &model, a);
        rp_spheroid_acceleration(y, v, &model, b);
        symmetry = fmax(symmetry, hypot(a[0] - b[0], a[2] + b[2]) / hypot(a[0], a[2]));
        double lo[3], hi[3];
        /* The surface s = 1 along the direction mu = nu. */
        const double rs = 1.0 / sqrt(1.0 + (1.0 / (q * q) - 1.0) * nu * nu);
        const double s = sqrt(1.0 - nu * nu);
        const double xi[3] = {(1.0 - 1e-13) * rs * s, 0.0, (1.0 - 1e-13) * rs * nu};
        const double xo[3] = {(1.0 + 1e-13) * rs * s, 0.0, (1.0 + 1e-13) * rs * nu};
        rp_spheroid_acceleration(xi, v, &model, lo);
        rp_spheroid_acceleration(xo, v, &model, hi);
        jump = fmax(jump, hypot(lo[0] - hi[0], lo[2] - hi[2]) / hypot(lo[0], lo[2]));
    }
    printf("  z -> -z %.1e, jump across the surface %.1e\n", symmetry, jump);
    check("spheroid z -> -z symmetry", symmetry, 1.0e-15);
    check("spheroid continuity across the surface", jump, 1.0e-9);
    rp_spheroid_table_free(table);
    return worst;
}

int main(void)
{
    sphere();
    quadrature_and_mass();
    thin_limit();

    rp_spheroid_params params;
    rp_spheroid_default_params(&params);
    check("table q=1 alpha=-2", table_accuracy(&params, 1000), 3.0e-7);
    params.axis_ratio = 0.3;
    params.alpha = -1.5;
    params.inner_ratio = 0.05;
    /* The worst points lie within ~1e-3 s_in outside the empty core, where
     * the force vanishes linearly; elsewhere this table is below 1e-6. */
    check("table q=0.3 alpha=-1.5 s_in=0.05", table_accuracy(&params, 1000), 2.0e-4);
    rp_spheroid_default_params(&params);
    params.axis_ratio = 0.05;
    check("table q=0.05 alpha=-2", table_accuracy(&params, 1000), 5.0e-6);

    rp_spheroid_default_params(&params);
    params.alpha = -3.5;
    double f[3];
    if (rp_spheroid_direct_field(&params, 0.5, 0.0, f) != RP_AXISYMMETRIC_STATUS_INVALID_ARGUMENT) {
        printf("FAIL alpha <= -3 with s_in = 0 accepted\n");
        ++failures;
    }
    if (failures != 0) {
        printf("test_thin_disk_spheroid: %d failure(s)\n", failures);
        return EXIT_FAILURE;
    }
    printf("test_thin_disk_spheroid: ok\n");
    return EXIT_SUCCESS;
}
