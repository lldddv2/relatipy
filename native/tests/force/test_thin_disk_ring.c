/* Ring field and direct disk sum of the thin-disk model. */
#include "force/models/thin_disk.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159265358979323846

static int failures = 0;

static void check(const char *label, double error, double tolerance)
{
    if (!(error <= tolerance)) {
        printf("FAIL %s: %.3e > %.1e\n", label, error, tolerance);
        ++failures;
    }
}

/* Softened point masses on the ring: a periodic trapezoidal sum, which
 * converges exponentially at the rate H / R used here. */
static void brute_ring(double radius, double h, double d, double z, double out[3])
{
    const int n = 40000;
    double ax = 0.0, az = 0.0, phi = 0.0;
    for (int i = 0; i < n; ++i) {
        const double t = 2.0 * PI * (i + 0.5) / n;
        const double dx = d - radius * cos(t);
        const double dy = -radius * sin(t);
        const double s2 = dx * dx + dy * dy + z * z + h * h;
        const double s = sqrt(s2);
        ax -= dx / (s2 * s);
        az -= z / (s2 * s);
        phi -= 1.0 / s;
    }
    out[0] = ax / n;
    out[1] = az / n;
    out[2] = phi / n;
}

static void ring_versus_brute_force(void)
{
    /* (R, H, d, z) */
    static const double cases[][4] = {
        {0.7, 0.007, 0.3, 0.05}, {0.7, 0.007, 0.7, 0.0}, {0.7, 0.007, 0.7, 0.003},
        {0.7, 0.007, 0.69, 0.001}, {0.7, 0.007, 2.0, 1.0}, {1.0, 0.01, 0.05, 0.0},
        {1e-3, 1e-5, 0.4, 0.2}, {0.2, 0.0, 0.5, 0.1}, {0.5, 0.005, 0.1, 3.0},
    };
    double worst_a = 0.0, worst_phi = 0.0;
    for (size_t k = 0; k < sizeof cases / sizeof cases[0]; ++k) {
        const double *c = cases[k];
        double f[3], b[3];
        if (rp_thin_disk_ring_field(c[0], c[1], c[2], c[3], f) != 0) {
            printf("FAIL ring status\n");
            ++failures;
            continue;
        }
        brute_ring(c[0], c[1], c[2], c[3], b);
        const double ad = f[0] * c[2];
        const double az = f[1] * c[3];
        const double err_a = hypot(ad - b[0], az - b[1]) / hypot(b[0], b[1]);
        const double err_phi = fabs(f[2] - b[2]) / fabs(b[2]);
        worst_a = fmax(worst_a, err_a);
        worst_phi = fmax(worst_phi, err_phi);
    }
    printf("ring vs brute-force point masses: max rel |a| %.1e, Phi %.1e\n",
           worst_a, worst_phi);
    check("ring acceleration vs brute force", worst_a, 1.0e-12);
    check("ring potential vs brute force", worst_phi, 1.0e-13);
}

static void ring_limits(void)
{
    double f[3], g[3];
    /* d -> 0: a_d/d = (R^2 - 2 zeta^2) / (2 s^5), a_z/z = -1/s^3, Phi = -1/s. */
    const double radius = 0.6, h = 0.006, z = 0.2;
    const double zeta2 = z * z + h * h;
    const double s = sqrt(radius * radius + zeta2);
    rp_thin_disk_ring_field(radius, h, 0.0, z, f);
    check("axis a_d/d", fabs(f[0] - (radius * radius - 2.0 * zeta2) / (2.0 * pow(s, 5))) * pow(s, 3), 1.0e-15);
    check("axis a_z/z", fabs(f[1] + 1.0 / (s * s * s)) * pow(s, 3), 1.0e-15);
    check("axis Phi", fabs(f[2] + 1.0 / s) * s, 1.0e-15);
    /* Smooth approach to the axis (no 0/0 at small d). */
    rp_thin_disk_ring_field(radius, h, 1.0e-9, z, g);
    check("axis continuity", fabs(g[0] - f[0]) / fabs(f[0]), 1.0e-14);
    /* R -> 0: a softened point mass at the origin. */
    const double d = 0.3, zz = 0.4, hh = 1.0e-3;
    const double s2 = d * d + zz * zz + hh * hh;
    rp_thin_disk_ring_field(0.0, hh, d, zz, f);
    check("point a_d/d", fabs(f[0] * pow(s2, 1.5) + 1.0), 1.0e-15);
    check("point a_z/z", fabs(f[1] * pow(s2, 1.5) + 1.0), 1.0e-15);
    rp_thin_disk_ring_field(1.0e-9, hh, d, zz, g);
    check("small-ring continuity", fabs(g[0] - f[0]) / fabs(f[0]), 1.0e-14);
    /* Far field: monopole with an O(R^2 / r^2) quadrupole. */
    const double far = 1.0e4;
    rp_thin_disk_ring_field(1.0, 0.01, far * 0.6, far * 0.8, f);
    const double ad = f[0] * far * 0.6, az = f[1] * far * 0.8;
    const double mono = 1.0 / (far * far);
    check("ring monopole", hypot(ad + 0.6 * mono, az + 0.8 * mono) / mono, 2.0e-8);
    /* Invalid input. */
    if (rp_thin_disk_ring_field(0.5, 0.0, 0.5, 0.0, f) == 0
        || rp_thin_disk_ring_field(-1.0, 0.0, 0.5, 0.0, f) == 0
        || rp_thin_disk_ring_field(0.5, 0.0, NAN, 0.0, f) == 0) {
        printf("FAIL ring invalid input accepted\n");
        ++failures;
    }
}

/* Central differences of Phi against the closed-form field. */
static void ring_gradient(void)
{
    static const double points[][2] = {{0.3, 0.05}, {0.71, 0.004}, {1.5, 0.3}, {0.05, 0.6}};
    const double radius = 0.7, h = 0.007;
    double worst = 0.0;
    for (size_t k = 0; k < sizeof points / sizeof points[0]; ++k) {
        const double d = points[k][0], z = points[k][1];
        const double step = 1.0e-5 * fmax(fabs(d - radius), 0.01);
        double f[3], p[3], m[3];
        rp_thin_disk_ring_field(radius, h, d, z, f);
        rp_thin_disk_ring_field(radius, h, d + step, z, p);
        rp_thin_disk_ring_field(radius, h, d - step, z, m);
        const double ad = -(p[2] - m[2]) / (2.0 * step);
        rp_thin_disk_ring_field(radius, h, d, z + step, p);
        rp_thin_disk_ring_field(radius, h, d, z - step, m);
        const double az = -(p[2] - m[2]) / (2.0 * step);
        const double norm = hypot(f[0] * d, f[1] * z);
        worst = fmax(worst, hypot(ad - f[0] * d, az - f[1] * z) / norm);
    }
    printf("ring -grad Phi vs field: max rel %.1e\n", worst);
    check("ring potential gradient", worst, 1.0e-8);
}

/* The direct disk sum: quadrature convergence and the total mass. */
static void direct_sum(void)
{
    rp_thin_disk_params base, fine;
    rp_thin_disk_default_params(&base);
    static const double points[][2] = {
        {0.3, 0.0}, {0.3, 0.003}, {1.0, 0.0}, {1.0, 0.01}, {1e-4, 0.0},
        {2.0, 0.01}, {0.0, 0.5}, {0.999, 2e-4}, {1e-6, 1e-8}, {0.5, 0.5},
    };
    static const double configs[][2] = {{-1.0, 0.0}, {-1.7, 0.0}, {0.5, 0.0}, {-1.0, 0.2}, {-3.0, 0.05}};
    double worst = 0.0;
    for (size_t c = 0; c < sizeof configs / sizeof configs[0]; ++c) {
        base.alpha = configs[c][0];
        base.inner_ratio = configs[c][1];
        fine = base;
        fine.quadrature_points = 40U;
        for (size_t k = 0; k < sizeof points / sizeof points[0]; ++k) {
            double a[3], b[3];
            if (rp_thin_disk_direct_field(&base, points[k][0], points[k][1], a) != RP_AXISYMMETRIC_STATUS_OK
                || rp_thin_disk_direct_field(&fine, points[k][0], points[k][1], b) != RP_AXISYMMETRIC_STATUS_OK) {
                printf("FAIL direct status\n");
                ++failures;
                continue;
            }
            const double d = points[k][0], z = points[k][1];
            const double norm = hypot(b[0] * d, b[1] * z);
            if (norm > 0.0) {
                worst = fmax(worst, hypot((a[0] - b[0]) * d, (a[1] - b[1]) * z) / norm);
            }
            worst = fmax(worst, fabs(a[2] - b[2]) / fabs(b[2]));
        }
        /* Total mass: far away the field is the unit monopole. */
        double f[3];
        const double far = 1.0e5;
        rp_thin_disk_direct_field(&base, 0.0, far, f);
        check("direct unit mass (axis)", fabs(-f[1] * far * far * far - 1.0), 1.0e-9);
        rp_thin_disk_direct_field(&base, far, 0.0, f);
        check("direct unit mass (plane)", fabs(-f[0] * far * far * far - 1.0), 1.0e-9);
        check("direct unit potential", fabs(-f[2] * far - 1.0), 1.0e-9);
    }
    printf("direct sum, 16 vs 40 Gauss points: max rel %.1e\n", worst);
    check("direct quadrature convergence", worst, 1.0e-13);

    rp_thin_disk_default_params(&base);
    base.alpha = -2.5; /* diverging mass without an inner edge */
    double f[3];
    if (rp_thin_disk_direct_field(&base, 0.5, 0.0, f) != RP_AXISYMMETRIC_STATUS_INVALID_ARGUMENT) {
        printf("FAIL alpha <= -2 with R_in = 0 accepted\n");
        ++failures;
    }
}

int main(void)
{
    ring_versus_brute_force();
    ring_limits();
    ring_gradient();
    direct_sum();
    if (failures != 0) {
        printf("test_thin_disk_ring: %d failure(s)\n", failures);
        return EXIT_FAILURE;
    }
    printf("test_thin_disk_ring: ok\n");
    return EXIT_SUCCESS;
}
