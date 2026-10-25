#include "force/models/axisymmetric_table.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Guide keys keep the exponent and 7 mantissa bits: 128 bins per octave. */
#define RP_AXISYMMETRIC_GUIDE_SHIFT 45U
/* Sub-steps per local spacing when integrating the node density. */
#define RP_AXISYMMETRIC_PLACEMENT_SUBSTEPS 16.0
#define RP_AXISYMMETRIC_MAX_NODES 200000U
#define RP_AXISYMMETRIC_ALIGNMENT 64U

/* ------------------------------------------------------------------ */
/* Axis construction                                                  */
/* ------------------------------------------------------------------ */

static double radial_spacing(const rp_axisymmetric_spec *spec, double r)
{
    double h = spec->radial_step * r;
    for (size_t e = 0; e < spec->edge_count; ++e) {
        const double near_edge = fmax(
            spec->floor_step * spec->edge_width[e],
            spec->radial_step * fabs(r - spec->edge_radius[e])
        );
        h = fmin(h, near_edge);
    }
    return h;
}

static double angular_spacing(const rp_axisymmetric_spec *spec, double mu)
{
    return fmax(spec->floor_step * spec->plane_width, spec->angular_step * mu);
}

typedef double (*spacing_fn)(const rp_axisymmetric_spec *spec, double x);

/*
 * Nodes on [a, b) with local spacing close to spacing(x): the node density
 * 1 / spacing is integrated by the midpoint rule and the nodes are put at
 * equal steps of that integral, so the spacing varies smoothly and the
 * segment ends exactly at b.  With out == NULL the nodes are only counted.
 * Returns the number of nodes written (b itself is not written).
 */
static size_t place_segment(
    const rp_axisymmetric_spec *spec, spacing_fn spacing,
    double a, double b, double *out
)
{
    double x = a;
    double s = 0.0;
    while (x < b) {
        const double sub = fmin(
            spacing(spec, x) / RP_AXISYMMETRIC_PLACEMENT_SUBSTEPS, b - x);
        s += sub / spacing(spec, x + 0.5 * sub);
        x = sub == b - x ? b : x + sub;
    }
    size_t count = (size_t)ceil(s - 1.0e-9);
    if (count < 1U) {
        count = 1U;
    }
    if (out == NULL) {
        return count;
    }
    const double target_step = s / (double)count;
    size_t written = 0;
    out[written++] = a;
    x = a;
    s = 0.0;
    while (x < b && written < count) {
        const double sub = fmin(
            spacing(spec, x) / RP_AXISYMMETRIC_PLACEMENT_SUBSTEPS, b - x);
        const double ds = sub / spacing(spec, x + 0.5 * sub);
        while (written < count && s + ds >= target_step * (double)written) {
            const double fraction = (target_step * (double)written - s) / ds;
            out[written++] = x + fraction * sub;
        }
        s += ds;
        x = sub == b - x ? b : x + sub;
    }
    while (written < count) { /* rounding at the very end */
        out[written] = 0.5 * (out[written - 1] + b);
        ++written;
    }
    return count;
}

static uint64_t guide_key(double x)
{
    uint64_t bits;
    memcpy(&bits, &x, sizeof bits);
    return bits >> RP_AXISYMMETRIC_GUIDE_SHIFT;
}

static double guide_lower_bound(uint64_t key)
{
    const uint64_t bits = key << RP_AXISYMMETRIC_GUIDE_SHIFT;
    double x;
    memcpy(&x, &bits, sizeof x);
    return x;
}

static void free_axis(rp_axisymmetric_axis *axis)
{
    free(axis->node);
    free(axis->inverse_step);
    free(axis->guide);
    memset(axis, 0, sizeof *axis);
}

/* Fill inverse steps and the guide of an axis whose nodes are set. */
static rp_axisymmetric_status finish_axis(rp_axisymmetric_axis *axis, size_t *bytes)
{
    const size_t n = axis->count;
    axis->inverse_step = malloc((n - 1U) * sizeof *axis->inverse_step);
    if (axis->inverse_step == NULL) {
        return RP_AXISYMMETRIC_STATUS_ALLOCATION_FAILURE;
    }
    for (size_t i = 0; i + 1U < n; ++i) {
        const double h = axis->node[i + 1U] - axis->node[i];
        if (!(h > 0.0)) {
            return RP_AXISYMMETRIC_STATUS_INVALID_ARGUMENT;
        }
        axis->inverse_step[i] = 1.0 / h;
    }
    /* Keys below node[1] all belong to cell 0. */
    axis->guide_base = guide_key(axis->node[1]);
    axis->guide_count = (size_t)(guide_key(axis->node[n - 1U]) - axis->guide_base) + 2U;
    axis->guide = malloc(axis->guide_count * sizeof *axis->guide);
    if (axis->guide == NULL) {
        return RP_AXISYMMETRIC_STATUS_ALLOCATION_FAILURE;
    }
    size_t cell = 0;
    for (size_t g = 0; g < axis->guide_count; ++g) {
        const double lower = guide_lower_bound(axis->guide_base + g);
        while (cell + 2U < n && axis->node[cell + 1U] <= lower) {
            ++cell;
        }
        axis->guide[g] = (uint32_t)cell;
    }
    *bytes += n * sizeof *axis->node + (n - 1U) * sizeof *axis->inverse_step
        + axis->guide_count * sizeof *axis->guide;
    return RP_AXISYMMETRIC_STATUS_OK;
}

static int compare_doubles(const void *a, const void *b)
{
    const double x = *(const double *)a;
    const double y = *(const double *)b;
    return (x > y) - (x < y);
}

static rp_axisymmetric_status build_radial_axis(
    const rp_axisymmetric_spec *spec, rp_axisymmetric_axis *axis, size_t *bytes
)
{
    double breaks[RP_AXISYMMETRIC_MAX_EDGES + 2U];
    size_t nb = 0;
    breaks[nb++] = spec->radial_min;
    for (size_t e = 0; e < spec->edge_count; ++e) {
        const double r = spec->edge_radius[e];
        if (r > spec->radial_min * (1.0 + 1.0e-9) && r < spec->radial_far * (1.0 - 1.0e-9)) {
            breaks[nb++] = r;
        }
    }
    breaks[nb++] = spec->radial_far;
    qsort(breaks, nb, sizeof breaks[0], compare_doubles);

    size_t count = 1U;
    for (size_t k = 0; k + 1U < nb; ++k) {
        if (breaks[k + 1U] > breaks[k]) {
            count += place_segment(spec, radial_spacing, breaks[k], breaks[k + 1U], NULL);
        }
    }
    if (count > RP_AXISYMMETRIC_MAX_NODES) {
        return RP_AXISYMMETRIC_STATUS_INVALID_ARGUMENT;
    }
    axis->node = malloc(count * sizeof *axis->node);
    if (axis->node == NULL) {
        return RP_AXISYMMETRIC_STATUS_ALLOCATION_FAILURE;
    }
    size_t at = 0;
    for (size_t k = 0; k + 1U < nb; ++k) {
        if (breaks[k + 1U] > breaks[k]) {
            at += place_segment(spec, radial_spacing, breaks[k], breaks[k + 1U],
                                axis->node + at);
        }
    }
    axis->node[at++] = spec->radial_far;
    axis->count = at;
    if (axis->count < 4U) {
        return RP_AXISYMMETRIC_STATUS_INVALID_ARGUMENT;
    }
    return finish_axis(axis, bytes);
}

static rp_axisymmetric_status build_angular_axis(
    const rp_axisymmetric_spec *spec, rp_axisymmetric_axis *axis, size_t *bytes
)
{
    const size_t count = place_segment(spec, angular_spacing, 0.0, 1.0, NULL) + 1U;
    if (count > RP_AXISYMMETRIC_MAX_NODES) {
        return RP_AXISYMMETRIC_STATUS_INVALID_ARGUMENT;
    }
    axis->node = malloc(count * sizeof *axis->node);
    if (axis->node == NULL) {
        return RP_AXISYMMETRIC_STATUS_ALLOCATION_FAILURE;
    }
    const size_t at = place_segment(spec, angular_spacing, 0.0, 1.0, axis->node);
    axis->node[at] = 1.0;
    axis->count = at + 1U;
    if (axis->count < 4U) {
        return RP_AXISYMMETRIC_STATUS_INVALID_ARGUMENT;
    }
    return finish_axis(axis, bytes);
}

/* ------------------------------------------------------------------ */
/* Cubic spline slopes                                                */
/* ------------------------------------------------------------------ */

/*
 * Slopes of the C2 cubic spline through (x[i], y[i * y_stride]), written to
 * slope[i * s_stride].  The right end is not-a-knot; the left end is
 * not-a-knot, or a zero slope when zero_left_slope (plane symmetry).
 * work holds 2 n doubles.  Requires n >= 4.
 */
static void spline_slopes(
    size_t n, const double *x, const double *y, size_t y_stride,
    double *slope, size_t s_stride, int zero_left_slope, double *work
)
{
    double *cp = work;      /* modified super-diagonal */
    double *dp = work + n;  /* modified right-hand side */
#define H(i) (x[(i) + 1] - x[(i)])
#define DELTA(i) ((y[((i) + 1) * y_stride] - y[(i) * y_stride]) / H(i))
    double b, c, d;
    if (zero_left_slope) {
        b = 1.0;
        c = 0.0;
        d = 0.0;
    } else {
        const double h0 = H(0);
        const double h1 = H(1);
        b = h1;
        c = h0 + h1;
        d = ((h0 + 2.0 * (h0 + h1)) * h1 * DELTA(0) + h0 * h0 * DELTA(1)) / (h0 + h1);
    }
    cp[0] = c / b;
    dp[0] = d / b;
    for (size_t i = 1; i + 1U < n; ++i) {
        const double a_i = H(i);
        const double b_i = 2.0 * (H(i - 1) + H(i));
        const double c_i = H(i - 1);
        const double d_i = 3.0 * (H(i) * DELTA(i - 1) + H(i - 1) * DELTA(i));
        const double m = b_i - a_i * cp[i - 1];
        cp[i] = c_i / m;
        dp[i] = (d_i - a_i * dp[i - 1]) / m;
    }
    {
        const size_t i = n - 1U;
        const double hl = H(n - 2U);  /* last interval */
        const double hp = H(n - 3U);  /* the one before */
        const double a_i = hl + hp;
        const double b_i = hp;
        const double d_i = (hl * hl * DELTA(n - 3U)
                            + (2.0 * (hp + hl) + hl) * hp * DELTA(n - 2U)) / (hp + hl);
        const double m = b_i - a_i * cp[i - 1];
        dp[i] = (d_i - a_i * dp[i - 1]) / m;
        cp[i] = 0.0;
    }
#undef H
#undef DELTA
    slope[(n - 1U) * s_stride] = dp[n - 1U];
    for (size_t i = n - 1U; i-- > 0;) {
        slope[i * s_stride] = dp[i] - cp[i] * slope[(i + 1U) * s_stride];
    }
}

/*
 * Tensor-spline Hermite data of one quantity on the grid: value, d/dr,
 * d/dmu and d2/drdmu, stored at out[stride * node + offset + {0, 2, 4, 6}
 * / pitch], see pack_quantity().
 */
static void spline_derivatives(
    size_t nr, size_t nm, const double *r, const double *mu,
    const double *value, double *d_r, double *d_mu, double *d_rmu, double *work
)
{
    for (size_t j = 0; j < nm; ++j) {
        spline_slopes(nr, r, value + j * nr, 1U, d_r + j * nr, 1U, 0, work);
    }
    for (size_t i = 0; i < nr; ++i) {
        spline_slopes(nm, mu, value + i, nr, d_mu + i, nr, 1, work);
        spline_slopes(nm, mu, d_r + i, nr, d_rmu + i, nr, 1, work);
    }
}

/* ------------------------------------------------------------------ */
/* Far field                                                          */
/* ------------------------------------------------------------------ */

static double far_series(
    const double coefficient[][RP_AXISYMMETRIC_FAR_TERMS], size_t terms,
    double rho2, double u
)
{
    double total = 0.0;
    for (size_t kk = terms; kk-- > 0;) {
        double g = 0.0;
        for (size_t i = kk + 1U; i-- > 0;) {
            g = g * u + coefficient[kk][i];
        }
        total = total * rho2 + g;
    }
    return total;
}

static void far_fields(
    const rp_axisymmetric_table *table, double r, double mu, double *a, double *c, double *p
)
{
    const double rho = 1.0 / r;
    const double rho2 = rho * rho;
    const double u = mu * mu;
    *a = far_series(table->far_a, table->far_terms, rho2, u);
    *c = far_series(table->far_c, table->far_terms, rho2, u);
    if (p != NULL) {
        *p = far_series(table->far_p, table->far_terms, rho2, u);
    }
}

/* ------------------------------------------------------------------ */
/* Interpolation                                                      */
/* ------------------------------------------------------------------ */

static inline size_t axis_cell(const rp_axisymmetric_axis *axis, double x)
{
    const uint64_t key = guide_key(x);
    size_t g = key > axis->guide_base ? (size_t)(key - axis->guide_base) : 0U;
    if (g >= axis->guide_count) {
        g = axis->guide_count - 1U;
    }
    size_t i = axis->guide[g];
    const size_t last = axis->count - 2U;
    while (i < last && x >= axis->node[i + 1U]) {
        ++i;
    }
    return i;
}

/* Cubic Hermite weights: values at 0 and 1, slopes (times h) at 0 and 1. */
static inline void hermite_weights(double t, double h, double w[4])
{
    const double s = 1.0 - t;
    const double t2 = t * t;
    w[0] = (1.0 + 2.0 * t) * s * s;
    w[1] = t2 * (3.0 - 2.0 * t);
    w[2] = t * s * s * h;
    w[3] = -t2 * s * h;
}

/* Bicubic Hermite interpolation of A and C inside the grid. */
static inline void grid_ac(
    const rp_axisymmetric_table *table, double r, double mu, double *a, double *c
)
{
    const rp_axisymmetric_axis *ar = &table->radius;
    const rp_axisymmetric_axis *am = &table->cosine;
    const size_t i = axis_cell(ar, r);
    const size_t j = axis_cell(am, mu);
    double wr[4];
    double wm[4];
    hermite_weights((r - ar->node[i]) * ar->inverse_step[i],
                    ar->node[i + 1U] - ar->node[i], wr);
    hermite_weights((mu - am->node[j]) * am->inverse_step[j],
                    am->node[j + 1U] - am->node[j], wm);
    const size_t nr = ar->count;
    const double *n00 = table->field + 8U * (j * nr + i);
    const double *n01 = n00 + 8U * nr;
    const double *node[4] = {n00, n00 + 8, n01, n01 + 8};
    double sa = 0.0;
    double sc = 0.0;
    for (int b = 0; b < 2; ++b) {
        for (int q = 0; q < 2; ++q) {
            const double *v = node[2 * b + q];
            const double w0 = wr[q] * wm[b];
            const double w1 = wr[2 + q] * wm[b];
            const double w2 = wr[q] * wm[2 + b];
            const double w3 = wr[2 + q] * wm[2 + b];
            sa += w0 * v[0] + w1 * v[2] + w2 * v[4] + w3 * v[6];
            sc += w0 * v[1] + w1 * v[3] + w2 * v[5] + w3 * v[7];
        }
    }
    *a = sa;
    *c = sc;
}

static double grid_p(const rp_axisymmetric_table *table, double r, double mu)
{
    const rp_axisymmetric_axis *ar = &table->radius;
    const rp_axisymmetric_axis *am = &table->cosine;
    const size_t i = axis_cell(ar, r);
    const size_t j = axis_cell(am, mu);
    double wr[4];
    double wm[4];
    hermite_weights((r - ar->node[i]) * ar->inverse_step[i],
                    ar->node[i + 1U] - ar->node[i], wr);
    hermite_weights((mu - am->node[j]) * am->inverse_step[j],
                    am->node[j + 1U] - am->node[j], wm);
    const size_t nr = ar->count;
    const double *n00 = table->potential + 4U * (j * nr + i);
    const double *n01 = n00 + 4U * nr;
    const double *node[4] = {n00, n00 + 4, n01, n01 + 4};
    double sp = 0.0;
    for (int b = 0; b < 2; ++b) {
        for (int q = 0; q < 2; ++q) {
            const double *v = node[2 * b + q];
            sp += wr[q] * wm[b] * v[0] + wr[2 + q] * wm[b] * v[1]
                + wr[q] * wm[2 + b] * v[2] + wr[2 + q] * wm[2 + b] * v[3];
        }
    }
    return sp;
}

static inline double inner_scale(const rp_axisymmetric_table *table, double r)
{
    const double x = r / table->radial_min;
    return table->inner_exponent == 1.0 ? x : pow(x, table->inner_exponent);
}

static inline double blend_weight(const rp_axisymmetric_table *table, double r)
{
    const double s = (r - table->radial_match) * table->inverse_blend_width;
    return s * s * s * (10.0 + s * (-15.0 + 6.0 * s));
}

static inline void table_ac(
    const rp_axisymmetric_table *table, double r, double mu, double *a, double *c
)
{
    if (r < table->radial_match) {
        if (r >= table->radial_min) {
            grid_ac(table, r, mu, a, c);
            return;
        }
        grid_ac(table, table->radial_min, mu, a, c);
        const double scale = inner_scale(table, r);
        *a *= scale;
        *c *= scale;
        return;
    }
    if (r >= table->radial_far) {
        far_fields(table, r, mu, a, c, NULL);
        return;
    }
    double ga, gc, fa, fc;
    grid_ac(table, r, mu, &ga, &gc);
    far_fields(table, r, mu, &fa, &fc, NULL);
    const double w = blend_weight(table, r);
    *a = ga + w * (fa - ga);
    *c = gc + w * (fc - gc);
}

void rp_axisymmetric_table_fields(
    const rp_axisymmetric_table *table, double r, double mu, double out[3]
)
{
    double a, c, p;
    if (r < table->radial_min) {
        /* A, C ~ r^k, hence a_r ~ r^(k - 2); integrate a_r inwards for Phi. */
        const double r0 = table->radial_min;
        const double k = table->inner_exponent;
        double a0, c0;
        grid_ac(table, r0, mu, &a0, &c0);
        const double p0 = grid_p(table, r0, mu);
        const double u = mu * mu;
        const double radial0 = -(a0 * (1.0 - u) + c0 * u) / (r0 * r0);
        const double ratio = r / r0;
        const double integral = fabs(k - 1.0) < 1.0e-12
            ? -radial0 * r0 * log(ratio)
            : radial0 * r0 * (1.0 - pow(ratio, k - 1.0)) / (k - 1.0);
        const double scale = inner_scale(table, r);
        a = a0 * scale;
        c = c0 * scale;
        p = -r * (-p0 / r0 + integral);
    } else if (r < table->radial_match) {
        grid_ac(table, r, mu, &a, &c);
        p = grid_p(table, r, mu);
    } else if (r >= table->radial_far) {
        far_fields(table, r, mu, &a, &c, &p);
    } else {
        double ga, gc, fa, fc, fp;
        grid_ac(table, r, mu, &ga, &gc);
        const double gp = grid_p(table, r, mu);
        far_fields(table, r, mu, &fa, &fc, &fp);
        const double w = blend_weight(table, r);
        a = ga + w * (fa - ga);
        c = gc + w * (fc - gc);
        p = gp + w * (fp - gp);
    }
    out[0] = a;
    out[1] = c;
    out[2] = p;
}

/* ------------------------------------------------------------------ */
/* Build                                                              */
/* ------------------------------------------------------------------ */

static int finite_positive(double x)
{
    return isfinite(x) && x > 0.0;
}

static rp_axisymmetric_status validate_spec(const rp_axisymmetric_spec *spec)
{
    if (!finite_positive(spec->radial_min)
        || !(spec->radial_match > spec->radial_min)
        || !(spec->radial_far > spec->radial_match) || !isfinite(spec->radial_far)
        || !(spec->radial_step > 0.0 && spec->radial_step <= 0.5)
        || !(spec->angular_step > 0.0 && spec->angular_step <= 0.5)
        || !finite_positive(spec->floor_step)
        || !finite_positive(spec->plane_width)
        || !finite_positive(spec->inner_exponent)
        || spec->edge_count > RP_AXISYMMETRIC_MAX_EDGES
        || spec->far_order > RP_AXISYMMETRIC_MAX_FAR_ORDER
        || spec->far_order % 2U != 0U) {
        return RP_AXISYMMETRIC_STATUS_INVALID_ARGUMENT;
    }
    for (size_t e = 0; e < spec->edge_count; ++e) {
        if (!finite_positive(spec->edge_radius[e]) || !finite_positive(spec->edge_width[e])) {
            return RP_AXISYMMETRIC_STATUS_INVALID_ARGUMENT;
        }
    }
    for (size_t kk = 0; kk <= spec->far_order / 2U; ++kk) {
        for (size_t i = 0; i <= kk; ++i) {
            if (!isfinite(spec->far_coefficients[kk][i])) {
                return RP_AXISYMMETRIC_STATUS_INVALID_ARGUMENT;
            }
        }
    }
    return RP_AXISYMMETRIC_STATUS_OK;
}

static void set_far_series(rp_axisymmetric_table *table, const rp_axisymmetric_spec *spec)
{
    /* With P = sum r^-k G_k(u): A = sum r^-k [(k + 1) G_k + 2 u G_k'] and
     * C = A - sum r^-k 2 G_k', primes taken in u = mu^2. */
    table->far_terms = spec->far_order / 2U + 1U;
    for (size_t kk = 0; kk < table->far_terms; ++kk) {
        const double k = 2.0 * (double)kk;
        for (size_t i = 0; i <= kk; ++i) {
            const double g = spec->far_coefficients[kk][i];
            const double g_next = i < kk ? spec->far_coefficients[kk][i + 1U] : 0.0;
            table->far_p[kk][i] = g;
            table->far_a[kk][i] = (k + 1.0 + 2.0 * (double)i) * g;
            table->far_c[kk][i] = table->far_a[kk][i] - 2.0 * (double)(i + 1U) * g_next;
        }
    }
}

static void *aligned_buffer(size_t bytes)
{
    const size_t rounded = (bytes + RP_AXISYMMETRIC_ALIGNMENT - 1U)
        / RP_AXISYMMETRIC_ALIGNMENT * RP_AXISYMMETRIC_ALIGNMENT;
    return aligned_alloc(RP_AXISYMMETRIC_ALIGNMENT, rounded);
}

rp_axisymmetric_status rp_axisymmetric_table_build(
    const rp_axisymmetric_spec *spec,
    rp_axisymmetric_field_fn field,
    void *context,
    rp_axisymmetric_table **table_out
)
{
    if (table_out == NULL) {
        return RP_AXISYMMETRIC_STATUS_NULL_POINTER;
    }
    *table_out = NULL;
    if (spec == NULL || field == NULL) {
        return RP_AXISYMMETRIC_STATUS_NULL_POINTER;
    }
    rp_axisymmetric_status status = validate_spec(spec);
    if (status != RP_AXISYMMETRIC_STATUS_OK) {
        return status;
    }
    rp_axisymmetric_table *table = calloc(1U, sizeof *table);
    if (table == NULL) {
        return RP_AXISYMMETRIC_STATUS_ALLOCATION_FAILURE;
    }
    double *scratch = NULL;
    table->bytes = sizeof *table;
    status = build_radial_axis(spec, &table->radius, &table->bytes);
    if (status == RP_AXISYMMETRIC_STATUS_OK) {
        status = build_angular_axis(spec, &table->cosine, &table->bytes);
    }
    if (status != RP_AXISYMMETRIC_STATUS_OK) {
        goto fail;
    }
    const size_t nr = table->radius.count;
    const size_t nm = table->cosine.count;
    const size_t nodes = nr * nm;
    const size_t longest = nr > nm ? nr : nm;
    /* Twelve node arrays (value and three derivatives of A, C, P) plus work. */
    scratch = malloc((12U * nodes + 2U * longest) * sizeof *scratch);
    table->field = aligned_buffer(8U * nodes * sizeof *table->field);
    table->potential = aligned_buffer(4U * nodes * sizeof *table->potential);
    if (scratch == NULL || table->field == NULL || table->potential == NULL) {
        status = RP_AXISYMMETRIC_STATUS_ALLOCATION_FAILURE;
        goto fail;
    }
    table->bytes += 12U * nodes * sizeof(double); /* field (8) + potential (4) */
    double *value[3] = {scratch, scratch + nodes, scratch + 2U * nodes};
    double *deriv[3][3];
    for (int q = 0; q < 3; ++q) {
        for (int k = 0; k < 3; ++k) {
            deriv[q][k] = scratch + (size_t)(3 + 3 * q + k) * nodes;
        }
    }
    double *work = scratch + 12U * nodes;

    for (size_t j = 0; j < nm; ++j) {
        const double mu = table->cosine.node[j];
        const double sine = sqrt((1.0 - mu) * (1.0 + mu));
        for (size_t i = 0; i < nr; ++i) {
            const double r = table->radius.node[i];
            double f[3];
            if (field(r * sine, r * mu, context, f) != 0
                || !isfinite(f[0]) || !isfinite(f[1]) || !isfinite(f[2])) {
                status = RP_AXISYMMETRIC_STATUS_FIELD_FAILURE;
                goto fail;
            }
            const size_t n = j * nr + i;
            value[0][n] = -r * r * r * f[0];
            value[1][n] = -r * r * r * f[1];
            value[2][n] = -r * f[2];
        }
    }
    for (int q = 0; q < 3; ++q) {
        spline_derivatives(nr, nm, table->radius.node, table->cosine.node,
                           value[q], deriv[q][0], deriv[q][1], deriv[q][2], work);
    }
    for (size_t n = 0; n < nodes; ++n) {
        double *f = table->field + 8U * n;
        double *p = table->potential + 4U * n;
        for (int k = 0; k < 4; ++k) {
            f[2 * k] = k == 0 ? value[0][n] : deriv[0][k - 1][n];
            f[2 * k + 1] = k == 0 ? value[1][n] : deriv[1][k - 1][n];
            p[k] = k == 0 ? value[2][n] : deriv[2][k - 1][n];
        }
    }
    free(scratch);

    table->radial_min = spec->radial_min;
    table->radial_match = spec->radial_match;
    table->radial_far = spec->radial_far;
    table->inverse_blend_width = 1.0 / (spec->radial_far - spec->radial_match);
    table->inner_exponent = spec->inner_exponent;
    set_far_series(table, spec);
    *table_out = table;
    return RP_AXISYMMETRIC_STATUS_OK;

fail:
    free(scratch);
    rp_axisymmetric_table_free(table);
    return status;
}

void rp_axisymmetric_table_free(rp_axisymmetric_table *table)
{
    if (table == NULL) {
        return;
    }
    free_axis(&table->radius);
    free_axis(&table->cosine);
    free(table->field);
    free(table->potential);
    free(table);
}

/* ------------------------------------------------------------------ */
/* Placed model                                                       */
/* ------------------------------------------------------------------ */

rp_axisymmetric_status rp_axisymmetric_rotation_from_normal(
    const double normal[3], double rotation[3][3]
)
{
    if (normal == NULL || rotation == NULL) {
        return RP_AXISYMMETRIC_STATUS_NULL_POINTER;
    }
    const double norm = sqrt(normal[0] * normal[0] + normal[1] * normal[1]
                             + normal[2] * normal[2]);
    if (!isfinite(norm) || !(fabs(norm - 1.0) <= 1.0e-10)) {
        return RP_AXISYMMETRIC_STATUS_INVALID_ARGUMENT;
    }
    const double e3[3] = {normal[0] / norm, normal[1] / norm, normal[2] / norm};
    /* Helper axis: the coordinate axis least aligned with the normal. */
    size_t least = 0;
    for (size_t k = 1; k < 3; ++k) {
        if (fabs(e3[k]) < fabs(e3[least])) {
            least = k;
        }
    }
    double e1[3] = {0.0, 0.0, 0.0};
    e1[least] = 1.0;
    const double dot = e3[least];
    for (size_t k = 0; k < 3; ++k) {
        e1[k] -= dot * e3[k];
    }
    const double n1 = sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
    for (size_t k = 0; k < 3; ++k) {
        e1[k] /= n1;
    }
    const double e2[3] = {
        e3[1] * e1[2] - e3[2] * e1[1],
        e3[2] * e1[0] - e3[0] * e1[2],
        e3[0] * e1[1] - e3[1] * e1[0],
    };
    for (size_t k = 0; k < 3; ++k) {
        rotation[0][k] = e1[k];
        rotation[1][k] = e2[k];
        rotation[2][k] = e3[k];
    }
    return RP_AXISYMMETRIC_STATUS_OK;
}

rp_axisymmetric_status rp_axisymmetric_model_init(
    rp_axisymmetric_model *model,
    const rp_axisymmetric_table *table,
    double mass,
    double length_scale,
    const double normal[3]
)
{
    if (model == NULL || table == NULL || normal == NULL) {
        return RP_AXISYMMETRIC_STATUS_NULL_POINTER;
    }
    if (!isfinite(mass) || !finite_positive(length_scale)
        || !finite_positive(1.0 / length_scale)) {
        return RP_AXISYMMETRIC_STATUS_INVALID_ARGUMENT;
    }
    double rotation[3][3];
    const rp_axisymmetric_status status =
        rp_axisymmetric_rotation_from_normal(normal, rotation);
    if (status != RP_AXISYMMETRIC_STATUS_OK) {
        return status;
    }
    model->table = table;
    model->mass = mass;
    model->length_scale = length_scale;
    model->inverse_length_scale = 1.0 / length_scale;
    memcpy(model->rotation, rotation, sizeof rotation);
    return RP_AXISYMMETRIC_STATUS_OK;
}

int rp_axisymmetric_acceleration(
    const double position[3],
    const double velocity[3],
    const void *opaque_model,
    double acceleration[3]
)
{
    (void)velocity;
    const rp_axisymmetric_model *model = opaque_model;
    if (model == NULL || position == NULL || acceleration == NULL) {
        return 1;
    }
    const double x = position[0];
    const double y = position[1];
    const double z = position[2];
    const double r2 = x * x + y * y + z * z;
    if (!(r2 > 0.0 && r2 < INFINITY)) {
        return 1;
    }
    const double *axis = model->rotation[2];
    const double height = axis[0] * x + axis[1] * y + axis[2] * z;
    const double r = sqrt(r2);
    const double inverse_r = 1.0 / r;
    double mu = fabs(height) * inverse_r;
    if (mu > 1.0) {
        mu = 1.0;
    }
    double a, c;
    table_ac(model->table, r * model->inverse_length_scale, mu, &a, &c);
    /* a = -M (A x + (C - A) height n) / r^3 */
    const double factor = -model->mass * inverse_r * inverse_r * inverse_r;
    const double along = (c - a) * height;
    acceleration[0] = factor * (a * x + along * axis[0]);
    acceleration[1] = factor * (a * y + along * axis[1]);
    acceleration[2] = factor * (a * z + along * axis[2]);
    return 0;
}

int rp_axisymmetric_potential(
    const double position[3], const void *opaque_model, double *potential
)
{
    const rp_axisymmetric_model *model = opaque_model;
    if (model == NULL || position == NULL || potential == NULL) {
        return 1;
    }
    const double x = position[0];
    const double y = position[1];
    const double z = position[2];
    const double r2 = x * x + y * y + z * z;
    if (!(r2 > 0.0 && r2 < INFINITY)) {
        return 1;
    }
    const double *axis = model->rotation[2];
    const double height = axis[0] * x + axis[1] * y + axis[2] * z;
    const double r = sqrt(r2);
    const double mu = fmin(fabs(height) / r, 1.0);
    double fields[3];
    rp_axisymmetric_table_fields(model->table, r * model->inverse_length_scale, mu, fields);
    *potential = -model->mass * fields[2] / r;
    return 0;
}
