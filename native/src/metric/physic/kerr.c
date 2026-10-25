/**
 * @file kerr.c
 * @brief Kerr physical evaluators for the native reference.
 *
 * This private translation unit contains the physical formulas used by the
 * orchestration layer in `../kerr.c`.  It does not define the public-facing
 * experimental entry points, own caller buffers, or select integration
 * policy.
 *
 * The metric uses `Delta = r^2 - 2 M r + a^2`.  The Levi--Civita
 * connection of this metric is evaluated by the geodesic layer
 * (`geodesic/physic/`), not here.
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
