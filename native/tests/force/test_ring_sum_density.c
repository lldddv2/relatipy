/* Ring sum of a thin disk with a surface density given as a callback. */
#include "utils/integrators/axisymmetric/ring_sum.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures = 0;

static void check(const char *label, double error, double tolerance)
{
    if (!(error <= tolerance)) {
        printf("FAIL %s: %.3e > %.1e\n", label, error, tolerance);
        ++failures;
    }
}

static double power_law(double radius, const void *context)
{
    return pow(radius, *(const double *)context);
}

/* Truncated exponential disk; context = scale length / R_out. */
static double exponential(double radius, const void *context)
{
    return exp(-radius / *(const double *)context);
}

static double field_error(const double a[3], const double b[3], double d, double z)
{
    const double norm = hypot(b[0] * d, b[1] * z);
    double error = fabs(a[2] - b[2]) / fabs(b[2]);
    if (norm > 0.0) {
        error = fmax(error, hypot((a[0] - b[0]) * d, (a[1] - b[1]) * z) / norm);
    }
    return error;
}

static const double points[][2] = {
    {0.3, 0.0}, {0.3, 0.003}, {1.0, 0.0}, {1.0, 0.01}, {1e-4, 0.0},
    {2.0, 0.01}, {0.0, 0.5}, {0.999, 2e-4}, {0.5, 0.5},
};

/* A power-law callback reproduces the exact power-law path. */
static void callback_matches_power_law(void)
{
    static const double configs[][2] = {{-1.0, 0.0}, {-1.7, 0.0}, {0.5, 0.0}, {-1.0, 0.2}, {-3.0, 0.05}};
    double worst_field = 0.0, worst_moment = 0.0;
    for (size_t c = 0; c < sizeof configs / sizeof configs[0]; ++c) {
        const double alpha = configs[c][0];
        const rp_thin_mass exact = {
            .sigma = NULL, .context = NULL, .inner_ratio = configs[c][1],
            .softening = 0.01, .beta = alpha + 2.0, .inner_exponent = alpha + 2.0,
            .quadrature_points = 16U,
        };
        rp_thin_mass callback = exact;
        callback.sigma = power_law;
        callback.context = &alpha;
        rp_thin_ring_sum a, b;
        if (rp_thin_ring_sum_init(&exact, &a) != RP_AXISYMMETRIC_STATUS_OK
            || rp_thin_ring_sum_init(&callback, &b) != RP_AXISYMMETRIC_STATUS_OK) {
            printf("FAIL init status (alpha %.1f)\n", alpha);
            ++failures;
            continue;
        }
        for (size_t k = 0; k < sizeof points / sizeof points[0]; ++k) {
            double fa[3], fb[3];
            if (rp_thin_ring_sum_field(points[k][0], points[k][1], &a, fa) != 0
                || rp_thin_ring_sum_field(points[k][0], points[k][1], &b, fb) != 0) {
                printf("FAIL field status\n");
                ++failures;
                continue;
            }
            worst_field = fmax(worst_field, field_error(fb, fa, points[k][0], points[k][1]));
        }
        for (int k = 0; k <= 24; k += 2) {
            const double ma = rp_thin_ring_sum_moment(&a, k);
            const double mb = rp_thin_ring_sum_moment(&b, k);
            worst_moment = fmax(worst_moment, fabs(mb - ma) / fabs(ma));
        }
    }
    printf("power-law callback vs exact path: field %.1e, moments %.1e\n", worst_field, worst_moment);
    check("callback field", worst_field, 1.0e-12);
    check("callback moments", worst_moment, 1.0e-12);
}

/* Composite Simpson reference of int_0^1 exp(-R/h) R^(k+1) dR. */
static double simpson_moment(double h, int k)
{
    const int n = 200000;
    double total = 0.0;
    for (int i = 0; i <= n; ++i) {
        const double r = (double)i / n;
        const double w = (i == 0 || i == n) ? 1.0 : (i % 2 ? 4.0 : 2.0);
        total += w * exp(-r / h) * pow(r, k + 1);
    }
    return total / (3.0 * n);
}

/* A non-power-law density: unit mass, exact moments and convergence. */
static void exponential_disk(void)
{
    const double h = 0.25;
    const rp_thin_mass mass = {
        .sigma = exponential, .context = &h, .inner_ratio = 0.0,
        .softening = 0.01, .beta = 1.0, .inner_exponent = 2.0, .quadrature_points = 16U,
    };
    rp_thin_mass fine_mass = mass;
    fine_mass.quadrature_points = 40U;
    rp_thin_ring_sum sum, fine;
    if (rp_thin_ring_sum_init(&mass, &sum) != RP_AXISYMMETRIC_STATUS_OK
        || rp_thin_ring_sum_init(&fine_mass, &fine) != RP_AXISYMMETRIC_STATUS_OK) {
        printf("FAIL exponential init\n");
        ++failures;
        return;
    }
    double worst_moment = 0.0;
    const double norm = simpson_moment(h, 0);
    for (int k = 2; k <= 24; k += 2) {
        const double reference = simpson_moment(h, k) / norm;
        worst_moment = fmax(worst_moment,
                            fabs(rp_thin_ring_sum_moment(&sum, k) - reference) / reference);
    }
    check("exponential moments vs Simpson", worst_moment, 1.0e-10);

    double worst = 0.0;
    for (size_t k = 0; k < sizeof points / sizeof points[0]; ++k) {
        double a[3], b[3];
        if (rp_thin_ring_sum_field(points[k][0], points[k][1], &sum, a) != 0
            || rp_thin_ring_sum_field(points[k][0], points[k][1], &fine, b) != 0) {
            printf("FAIL exponential field status\n");
            ++failures;
            continue;
        }
        worst = fmax(worst, field_error(a, b, points[k][0], points[k][1]));
    }
    check("exponential quadrature convergence", worst, 1.0e-12);

    /* Far away the field is the unit monopole. */
    const double far = 1.0e5;
    double f[3];
    rp_thin_ring_sum_field(0.0, far, &sum, f);
    check("exponential unit mass (axis)", fabs(-f[1] * far * far * far - 1.0), 1.0e-9);
    rp_thin_ring_sum_field(far, 0.0, &sum, f);
    check("exponential unit mass (plane)", fabs(-f[0] * far * far * far - 1.0), 1.0e-9);
    printf("exponential disk: moments %.1e, 16 vs 40 points %.1e\n", worst_moment, worst);

    /* The table built from the callback interpolates the direct sum. */
    const rp_axisymmetric_grid grid = {
        .radial_min = 1.0e-5, .radial_match = 4.0, .radial_far = 5.0,
        .radial_step = 0.04, .angular_step = 0.02, .floor_step = 0.08, .far_order = 24U,
    };
    rp_axisymmetric_table *table = NULL;
    if (rp_thin_ring_sum_table_build(&mass, &grid, &table) != RP_AXISYMMETRIC_STATUS_OK) {
        printf("FAIL exponential table build\n");
        ++failures;
        return;
    }
    double worst_table = 0.0;
    for (size_t k = 0; k < sizeof points / sizeof points[0]; ++k) {
        const double d = points[k][0], z = points[k][1];
        const double r = hypot(d, z);
        if (r < 1.0e-3) {
            continue;
        }
        double exact[3], interpolated[3];
        rp_thin_ring_sum_field(d, z, &sum, exact);
        rp_axisymmetric_table_fields(table, r, z / r, interpolated);
        /* table stores A = -r^3 a_d / d, C = -r^3 a_z / z, P = -r Phi */
        const double r3 = r * r * r;
        const double a_d = -interpolated[0] / r3, a_z = -interpolated[1] / r3;
        const double reference = hypot(exact[0] * d, exact[1] * z);
        worst_table = fmax(worst_table, hypot((a_d - exact[0]) * d, (a_z - exact[1]) * z) / reference);
        worst_table = fmax(worst_table, fabs(-interpolated[2] / r - exact[2]) / fabs(exact[2]));
    }
    rp_axisymmetric_table_free(table);
    printf("exponential table vs direct sum: %.1e\n", worst_table);
    check("exponential table", worst_table, 1.0e-6);
}

static void invalid_input(void)
{
    const double alpha = -2.5;
    rp_thin_mass mass = {
        .sigma = power_law, .context = &alpha, .inner_ratio = 0.0,
        .softening = 0.01, .beta = -0.5, .inner_exponent = -0.5, .quadrature_points = 16U,
    };
    rp_thin_ring_sum sum;
    if (rp_thin_ring_sum_init(&mass, &sum) != RP_AXISYMMETRIC_STATUS_INVALID_ARGUMENT) {
        printf("FAIL beta <= 0 with R_in = 0 accepted\n");
        ++failures;
    }
    mass.beta = 1.0;
    mass.softening = 0.0;
    if (rp_thin_ring_sum_init(&mass, &sum) != RP_AXISYMMETRIC_STATUS_INVALID_ARGUMENT) {
        printf("FAIL zero softening accepted\n");
        ++failures;
    }
}

int main(void)
{
    callback_matches_power_law();
    exponential_disk();
    invalid_input();
    if (failures != 0) {
        printf("test_ring_sum_density: %d failure(s)\n", failures);
        return EXIT_FAILURE;
    }
    printf("test_ring_sum_density: ok\n");
    return EXIT_SUCCESS;
}
