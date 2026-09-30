/**
 * @file kerr.c
 * @brief Kerr physical evaluators for the native reference.
 *
 * This private translation unit contains the physical formulas used by the
 * orchestration layer in `../kerr.c`.  It does not define the public-facing
 * experimental entry points, own caller buffers, or select integration
 * policy.
 *
 * The metric uses `Delta = r^2 - 2 M r + a^2`.  The connection is evaluated
 * separately as a direct translation of the explicit generated expressions
 * in the frozen legacy `Kerr._get_christoffel_symbols` method.  It is kept as
 * a compatibility reference and is not derived from metric derivatives.  Its
 * 20 independent expressions and 12 explicit symmetry copies populate 32
 * tensor entries; all other entries remain zero.
 *
 * Coordinate-time velocities are converted to a future-directed
 * four-velocity with the optimized expression derived from the frozen legacy
 * metric.  Because that metric has the opposite global signature, the result
 * is normalized to `g_(mu nu) u^mu u^nu = -1` under the metric exposed here.
 *
 * `RP_KERR_STATUS_PHYSICAL_SINGULARITY` distinguishes the Kerr ring from
 * chart singularities; `RP_KERR_STATUS_NUMERICAL_RANGE` reports finite
 * physical inputs whose intermediate or final values cannot be represented
 * as finite `double` values.
 */

#include "kerr.h"
#include "../../utils/numeric.h"
#include "../../utils/tensor.h"

#include <math.h>

rp_kerr_status rp_kerr_physic_prepare_point(
    double mass,
    double spin,
    struct rp_kerr_point *point
)
{
    double radius;
    double radius_squared;
    double spin_squared;
    double sigma_scale;
    double delta_scale;

    radius = point->x[1];
    radius_squared = radius * radius;
    spin_squared = spin * spin;
    if (!isfinite(radius_squared) || !isfinite(spin_squared)
        || !isfinite(mass * radius)) {
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    point->sin_theta = sin(point->x[2]);
    point->cos_theta = cos(point->x[2]);
    point->sigma = radius_squared
        + spin_squared * point->cos_theta * point->cos_theta;
    point->delta = radius_squared - 2.0 * mass * radius + spin_squared;

    if (!isfinite(point->sigma) || !isfinite(point->delta)) {
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }

    sigma_scale = radius_squared + spin_squared;
    delta_scale = radius_squared + 2.0 * mass * fabs(radius) + spin_squared;
    if (!isfinite(sigma_scale) || !isfinite(delta_scale)) {
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    if (rp_effectively_zero(point->sigma, sigma_scale)) {
        return RP_KERR_STATUS_PHYSICAL_SINGULARITY;
    }
    if (rp_bl_polar_axis_singular(point->x[2])
        || rp_effectively_zero(point->delta, delta_scale)) {
        return RP_KERR_STATUS_COORDINATE_SINGULARITY;
    }

    return RP_KERR_STATUS_OK;
}

void rp_kerr_physic_evaluate_metric(
    double mass,
    double spin,
    const struct rp_kerr_point *point,
    double metric[RP_KERR_DIM][RP_KERR_DIM]
)
{
    const double radius = point->x[1];
    const double radius_squared = radius * radius;
    const double spin_squared = spin * spin;
    const double sin_squared = point->sin_theta * point->sin_theta;
    const double common = 2.0 * mass * radius / point->sigma;

    rp_tensor_zero(&metric[0][0], RP_KERR_DIM * RP_KERR_DIM);
    metric[0][0] = -(1.0 - common);
    metric[0][3] = -common * spin * sin_squared;
    metric[3][0] = metric[0][3];
    metric[1][1] = point->sigma / point->delta;
    metric[2][2] = point->sigma;
    metric[3][3] = (radius_squared + spin_squared
        + common * spin_squared * sin_squared) * sin_squared;
}

rp_kerr_status rp_kerr_physic_evaluate_inverse_metric(
    double mass,
    double spin,
    const struct rp_kerr_point *point,
    double inverse_metric[RP_KERR_DIM][RP_KERR_DIM]
)
{
    double metric[RP_KERR_DIM][RP_KERR_DIM];
    double block_determinant;
    double block_scale;

    rp_tensor_zero(&inverse_metric[0][0], RP_KERR_DIM * RP_KERR_DIM);
    rp_kerr_physic_evaluate_metric(mass, spin, point, metric);
    block_determinant = metric[0][0] * metric[3][3]
        - metric[0][3] * metric[0][3];
    block_scale = fabs(metric[0][0] * metric[3][3])
        + metric[0][3] * metric[0][3];
    if (!isfinite(block_determinant) || !isfinite(block_scale)) {
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    if (rp_effectively_zero(block_determinant, block_scale)) {
        return RP_KERR_STATUS_COORDINATE_SINGULARITY;
    }

    inverse_metric[0][0] = metric[3][3] / block_determinant;
    inverse_metric[0][3] = -metric[0][3] / block_determinant;
    inverse_metric[3][0] = inverse_metric[0][3];
    inverse_metric[3][3] = metric[0][0] / block_determinant;
    inverse_metric[1][1] = point->delta / point->sigma;
    inverse_metric[2][2] = 1.0 / point->sigma;

    if (!rp_tensor_is_finite(
        &inverse_metric[0][0], RP_KERR_DIM * RP_KERR_DIM
    )) {
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    return RP_KERR_STATUS_OK;
}

void rp_kerr_physic_evaluate_christoffel(
    double mass,
    double spin,
    const struct rp_kerr_point *point,
    double christoffel[RP_KERR_DIM][RP_KERR_DIM][RP_KERR_DIM]
)
{
    const double r = point->x[1]; /* Boyer--Lindquist radial coordinate. */

    const double r2 = r * r; /* Consecutive powers of the radial coordinate. */
    const double r3 = r2 * r;
    const double r4 = r2 * r2;
    const double r5 = r4 * r;
    const double r6 = r3 * r3;
    const double r7 = r6 * r;

    const double rs = 2.0 * mass; /* Schwarzschild radius. */
    const double rs2 = rs * rs;
    const double rs_r = rs * r;
    const double rs_r3 = rs * r3;

    const double spin2 = spin * spin; /* Even powers of the spin length. */
    const double spin4 = spin2 * spin2;
    const double spin6 = spin4 * spin2;
    const double spin2_r2 = spin2 * r2;

    const double cos_th = point->cos_theta; /* Cosine values and powers. */
    const double cos2_th = cos_th * cos_th;
    const double cos4_th = cos2_th * cos2_th;

    const double sin_th = point->sin_theta; /* Sine values and powers. */
    const double sin2_th = sin_th * sin_th;
    const double sin4_th = sin2_th * sin2_th;

    const double sin_2th = 2.0 * sin_th * cos_th; /* Multiple-angle identities. */
    const double cos_2th = cos2_th - sin2_th;
    const double one_minus_cos_4th = 8.0 * sin2_th * cos2_th;
    const double cot_th = cos_th / sin_th;

    const double spin2_cos2_th = spin2 * cos2_th; /* Reused products. */
    const double spin4_cos2_th = spin4 * cos2_th;
    const double spin2_cos_2th = spin2 * cos_2th;
    const double rs_spin2_cos2_th_r = rs * spin2_cos2_th * r;
    const double two_r2 = 2.0 * r2;

    const double neg_rs_r = -rs_r; /* Reused negative aliases. */
    const double neg_spin2_cos2_th_r2 = -spin2_cos2_th * r2;
    const double neg_spin2_r2 = -spin2_r2;
    const double neg_r4 = -r4;

    const double sigma = spin2_cos2_th + r2; /* Kerr denominator terms. */
    const double two_sigma = two_r2 + spin2 + spin2_cos_2th;
    const double delta = neg_rs_r + spin2 + r2;

    const double neg_spin2_sin_2th_over_two_sigma = -spin2 * sin_2th / two_sigma; /* Reused signed ratio. */

    christoffel[0][0][1] =
        rs
        * (spin2_r2 * sin2_th
            + spin4 * sin2_th
            - spin4
            + r4)
        / (2.0 * sigma * sigma * delta);
    christoffel[0][0][2] =
        2.0 * rs
        * neg_spin2_sin_2th_over_two_sigma * r
        / two_sigma;
    christoffel[0][1][3] =
        rs * spin * sin2_th
        * (neg_spin2_cos2_th_r2
            + neg_spin2_r2
            + 3.0 * neg_r4
            + spin4_cos2_th)
        / (2.0 * sigma * sigma * delta);
    christoffel[0][2][3] =
        rs * spin * spin2 * cos_th
        * sin_th * sin2_th * r
        / (sigma * sigma);

    christoffel[0][1][0] = christoffel[0][0][1];
    christoffel[0][2][0] = christoffel[0][0][2];
    christoffel[0][3][1] = christoffel[0][1][3];
    christoffel[0][3][2] = christoffel[0][2][3];

    christoffel[1][0][0] =
        rs
        * (rs_spin2_cos2_th_r
            - rs_r3
            + neg_spin2_cos2_th_r2
            + neg_spin2_r2
            + spin4_cos2_th
            + r4)
        / (2.0 * sigma * sigma * sigma);
    christoffel[1][0][3] =
        rs * spin * sin2_th
        * (-rs_spin2_cos2_th_r
            + rs_r3
            + neg_r4
            + spin2_cos2_th * r2
            + spin2_r2
            - spin4_cos2_th)
        / (2.0 * sigma * sigma * sigma);
    christoffel[1][1][1] =
        (-rs * spin2_cos2_th / 2.0
            + rs * r2 / 2.0
            + spin2 * r
            + spin2_cos2_th * r)
        / (rs_spin2_cos2_th_r
            + rs_r3
            + neg_spin2_cos2_th_r2
            + neg_r4
            + spin2_r2
            + spin4_cos2_th);
    christoffel[1][1][2] =
        neg_spin2_sin_2th_over_two_sigma;
    christoffel[1][2][2] =
        r * (rs * r
            + spin2 - r2) / sigma;
    christoffel[1][3][3] =
        sin2_th
        * (rs * spin2 * sin2_th
                * r4 / 2.0
            + 2.0 * rs * spin2_cos2_th
                * r4
            + rs * spin4 * cos4_th
                * r2
            - rs * spin4 * sin2_th
                * r2 / 2.0
            - rs * spin4 * r2
                * one_minus_cos_4th / 16.0
            + rs * spin6
                * one_minus_cos_4th / 16.0
            + rs * r6
            - rs2 * spin2 * sin2_th
                * r3 / 2.0
            + rs2 * spin4 * r
                * one_minus_cos_4th / 16.0
            + spin2 * r5
            - 2.0 * spin2_cos2_th * r5
            - spin4 * cos4_th * r3
            + 2.0 * spin4_cos2_th * r3
            + spin6 * cos4_th * r
            - r7)
        / (sigma * sigma * sigma);

    christoffel[1][2][1] = christoffel[1][1][2];
    christoffel[1][3][0] = christoffel[1][0][3];

    christoffel[2][0][0] =
        4.0 * rs
        * neg_spin2_sin_2th_over_two_sigma * r
        / (two_sigma * two_sigma);
    christoffel[2][0][3] =
        4.0 * rs * spin * sin_2th * r
        * (spin2 + r2)
        / (two_sigma * two_sigma * two_sigma);
    christoffel[2][1][1] =
        -spin2 * cos_th * sin_th
        / (rs_spin2_cos2_th_r
            + rs_r3
            + neg_spin2_cos2_th_r2
            - neg_spin2_r2
            + neg_r4
            + spin4_cos2_th);
    christoffel[2][1][2] = r / sigma;
    christoffel[2][2][2] =
        neg_spin2_sin_2th_over_two_sigma;
    christoffel[2][3][3] =
        cos_th * sin_th
        * (-2.0 * rs_r3 * spin2 * sin2_th
            - two_r2 * spin4_cos2_th
            + neg_rs_r * spin4 * sin4_th
            + neg_rs_r * spin4
                * one_minus_cos_4th / 4.0
            + neg_r4 * spin2
            + 2.0 * neg_r4 * spin2_cos2_th
            - spin4 * cos4_th * r2
            - spin6 * cos4_th
            - r6)
        / (sigma * sigma * sigma);

    christoffel[2][2][1] = christoffel[2][1][2];
    christoffel[2][3][0] = christoffel[2][0][3];

    christoffel[3][0][1] =
        rs * spin
        * (-spin2_cos2_th + r2)
        / (2.0 * sigma * sigma * delta);
    christoffel[3][0][2] =
        cot_th * neg_rs_r * spin
        / (sigma * sigma);
    christoffel[3][1][3] =
        (rs
                * neg_spin2_cos2_th_r2
            + rs * neg_spin2_r2
                * sin2_th / 2.0
            + rs * neg_r4
            + rs * spin4
                * one_minus_cos_4th / 16.0
            + 2.0 * spin2_cos2_th * r3
            + spin4 * cos4_th * r
            + r5)
        / (sigma * sigma * delta);
    christoffel[3][2][3] =
        cot_th
        * (rs * spin2 * sin2_th * r
            + two_r2 * spin2
            + 2.0 * neg_spin2_r2 * sin2_th
            - 2.0 * spin4 * sin2_th
            + spin4 * sin4_th
            + spin4
            + r4)
        / (two_r2 * spin2_cos2_th
            + spin4 * cos4_th
            + r4);

    christoffel[3][1][0] = christoffel[3][0][1];
    christoffel[3][2][0] = christoffel[3][0][2];
    christoffel[3][3][1] = christoffel[3][1][3];
    christoffel[3][3][2] = christoffel[3][2][3];
}

rp_kerr_status rp_kerr_physic_evaluate_four_velocity(
    double mass,
    double spin,
    const struct rp_kerr_point *point,
    const double coordinate_velocity[3],
    double four_velocity[RP_KERR_DIM]
)
{
    const double r = point->x[1];  /* Boyer--Lindquist radial coordinate. */

    const double r2 = r * r; /* Radial coordinate squared. */

    const double rs = 2.0 * mass; /* Schwarzschild radius. */
    const double rs_r = rs * r;

    const double spin2 = spin * spin; /* Spin length squared. */

    const double sin_th = point->sin_theta; /* Sine value and its square. */
    const double sin2_th = sin_th * sin_th;

    const double v_r = coordinate_velocity[0]; /* Coordinate-time velocities. */
    const double v_th = coordinate_velocity[1];
    const double v_ph = coordinate_velocity[2];
    const double v_r2 = v_r * v_r;
    const double v_th2 = v_th * v_th;

    const double sigma = point->sigma; /* Kerr denominator terms. */
    const double r2_plus_spin2 = r2 + spin2;
    const double delta = r2_plus_spin2 - rs_r;

    const double rs_r_over_sigma = rs_r / sigma; /* Reused products. */
    const double rs_r_spin_over_sigma = rs_r_over_sigma * spin;
    const double sin2_th_v_ph = sin2_th * v_ph;

    const double normalization =
        -rs_r_over_sigma
        + rs_r_spin_over_sigma * sin2_th_v_ph
        - sigma * v_r2 / delta
        - sigma * v_th2
        + sin2_th_v_ph
            * (rs_r_spin_over_sigma
                - v_ph
                    * (r2_plus_spin2
                        + rs_r_over_sigma * sin2_th * spin2))
        + 1.0;
    double u0;

    if (!isfinite(normalization)) {
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    if (!(normalization > 0.0)) {
        return RP_KERR_STATUS_NON_TIMELIKE_VELOCITY;
    }

    u0 = 1.0 / sqrt(normalization);
    four_velocity[0] = u0;
    four_velocity[1] = u0 * v_r;
    four_velocity[2] = u0 * v_th;
    four_velocity[3] = u0 * v_ph;

    if (!rp_tensor_is_finite(four_velocity, RP_KERR_DIM)) {
        rp_tensor_zero(four_velocity, RP_KERR_DIM);
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    return RP_KERR_STATUS_OK;
}
