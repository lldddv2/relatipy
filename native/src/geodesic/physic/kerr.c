/**
 * @file kerr.c
 * @brief Physical Kerr geodesic right-hand-side evaluator.
 *
 * This private translation unit evaluates both the semantic C99 candidate
 * derived in `thesis/docs/v1/004-geodesica_optimizada.ipynb` and the
 * corrected null-geodesic path.  The former contracts the frozen-legacy
 * connection; the latter derives the Levi--Civita contraction from analytic
 * derivatives of the corrected Kerr metric.  Neither path materializes a
 * rank-three Christoffel tensor.  Input/output orchestration remains in
 * `../kerr.c`.
 */

#include "kerr.h"
#include "../../utils/numeric.h"
#include "../../utils/tensor.h"

#include <math.h>
#include <stddef.h>

rp_kerr_status rp_kerr_geodesic_physic_prepare_state(
    double mass,
    double spin,
    struct rp_kerr_geodesic_state *state
)
{
    double radius;
    double radius_squared;
    double spin_squared;
    double sigma_scale;
    double delta_scale;

    radius = state->x[1];
    radius_squared = radius * radius;
    spin_squared = spin * spin;
    if (!isfinite(radius_squared) || !isfinite(spin_squared)
        || !isfinite(mass * radius)) {
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }

    state->sin_theta = sin(state->x[2]);
    state->cos_theta = cos(state->x[2]);
    state->sigma = radius_squared
        + spin_squared * state->cos_theta * state->cos_theta;
    state->delta = radius_squared - 2.0 * mass * radius + spin_squared;
    if (!isfinite(state->sigma) || !isfinite(state->delta)) {
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }

    sigma_scale = radius_squared + spin_squared;
    delta_scale = radius_squared + 2.0 * mass * fabs(radius) + spin_squared;
    if (!isfinite(sigma_scale) || !isfinite(delta_scale)) {
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    if (rp_effectively_zero(state->sigma, sigma_scale)) {
        return RP_KERR_STATUS_PHYSICAL_SINGULARITY;
    }
    if (rp_bl_polar_axis_singular(state->x[2])
        || rp_effectively_zero(state->delta, delta_scale)) {
        return RP_KERR_STATUS_COORDINATE_SINGULARITY;
    }
    return RP_KERR_STATUS_OK;
}

void rp_kerr_geodesic_physic_evaluate_rhs(
    double mass,
    double spin,
    const struct rp_kerr_geodesic_state *state,
    double dx[RP_KERR_DIM],
    double du[RP_KERR_DIM]
)
{
    const double r = state->x[1]; /* Boyer--Lindquist radial coordinate. */

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
    const double rs_r4 = rs * r4;
    const double half_rs = 0.5 * rs;
    const double half_rs_r2 = half_rs * r2;

    const double spin2 = spin * spin; /* Powers of the spin length. */
    const double spin3 = spin2 * spin;
    const double spin4 = spin2 * spin2;
    const double spin6 = spin4 * spin2;

    const double cos_th = state->cos_theta; /* Cosine values and powers. */
    const double cos2_th = cos_th * cos_th;
    const double cos4_th = cos2_th * cos2_th;

    const double sin_th = state->sin_theta; /* Sine values and powers. */
    const double sin2_th = sin_th * sin_th;
    const double sin3_th = sin2_th * sin_th;
    const double sin4_th = sin2_th * sin2_th;

    const double sin_2th = 2.0 * sin_th * cos_th; /* Multiple-angle identities. */
    const double cos_2th = cos2_th - sin2_th;
    const double one_minus_cos_4th = 8.0 * sin2_th * cos2_th;
    const double cot_th = cos_th / sin_th;

    const double u0 = state->u[0]; /* Contravariant four-velocity. */
    const double u1 = state->u[1];
    const double u2 = state->u[2];
    const double u3 = state->u[3];
    const double u0_2 = u0 * u0;
    const double u1_2 = u1 * u1;
    const double u2_2 = u2 * u2;
    const double u3_2 = u3 * u3;
    const double two_u1 = 2.0 * u1;

    const double sigma = state->sigma; /* Kerr denominator terms. */
    const double delta = state->delta;
    const double inv_sigma2 = 1.0 / (sigma * sigma);
    const double inv_sigma3 = inv_sigma2 / sigma;
    const double r_over_sigma = r / sigma;
    const double inv_delta_sigma2 = inv_sigma2 / delta;
    const double two_r2 = 2.0 * r2;
    const double two_sigma = spin2 * cos_2th + spin2 + two_r2;
    const double inv_two_sigma = 1.0 / two_sigma;
    const double two_sigma_cubed = two_sigma * two_sigma * two_sigma;
    const double r2_plus_spin2 = r2 + spin2;

    const double spin2_cos2_th = cos2_th * spin2; /* Reused products. */
    const double spin4_cos2_th = cos2_th * spin4;
    const double spin4_cos4_th = cos4_th * spin4;
    const double spin6_cos4_th = cos4_th * spin6;
    const double spin2_cos2_th_r2 = r2 * spin2_cos2_th;
    const double spin4_cos4_th_r2 = r2 * spin4_cos4_th;
    const double two_spin2_cos2_th = 2.0 * spin2_cos2_th;
    const double spin2_r2 = r2 * spin2;
    const double spin2_r4 = r4 * spin2;
    const double two_spin2_r2 = 2.0 * spin2_r2;
    const double spin2_sin2_th = sin2_th * spin2;
    const double spin2_sin_2th = sin_2th * spin2;
    const double spin4_sin2_th = sin2_th * spin4;
    const double spin4_sin4_th = sin4_th * spin4;
    const double spin2_r2_sin2_th = sin2_th * spin2_r2;
    const double cos_th_sin_th = cos_th * sin_th;

    const double rs_u0 = rs * u0; /* Velocity-dependent products. */
    const double spin_u3 = spin * u3;
    const double spin_sin2_th_u3 = sin2_th * spin_u3;
    const double u1_over_delta_sigma2 = inv_delta_sigma2 * u1;
    const double rs_u0_u1_over_delta_sigma2 =
        rs_u0 * u1_over_delta_sigma2;
    const double u3_2_over_sigma3 = inv_sigma3 * u3_2;

    const double spin2_cos2_th_spin2_minus_r2_plus_rs_r =
        rs_r * spin2_cos2_th
        - spin2_cos2_th_r2
        + spin4_cos2_th;
    const double radial_poly_a =
        r4
        - rs_r3
        + spin2_cos2_th_spin2_minus_r2_plus_rs_r
        - spin2_r2;
    const double radial_poly_b =
        -r4
        + rs_r3
        + spin2_cos2_th_spin2_minus_r2_plus_rs_r
        + spin2_r2;
    const double u1_2_over_radial_poly_b = u1_2 / radial_poly_b;

    /* Angular products derived without additional trigonometric calls. */
    const double one_minus_cos_4th_over_16 = one_minus_cos_4th / 16.0;
    const double rs_one_minus_cos_4th_over_16 =
        one_minus_cos_4th_over_16 * rs;
    const double rs_spin4_one_minus_cos_4th_over_16 =
        rs_one_minus_cos_4th_over_16 * spin4;
    const double rs_r_over_two_sigma_cubed = rs_r / two_sigma_cubed;

    dx[0] = u0;
    dx[1] = u1;
    dx[2] = u2;
    dx[3] = u3;

    du[0] =
        -2.0 * cos_th * inv_sigma2 * rs_r * sin3_th * spin3 * u2 * u3
        + 4.0 * r * rs * sin_2th * spin2 * u0 * u2
            / (two_sigma * two_sigma)
        - rs * spin_sin2_th_u3 * u1_over_delta_sigma2
            * (-3.0 * r4
                - spin2_cos2_th_r2
                - spin2_r2
                + spin4_cos2_th)
        - rs_u0_u1_over_delta_sigma2
            * (r4
                + spin2_r2_sin2_th
                - spin4
                + spin4_sin2_th);

    du[1] =
        -half_rs * inv_sigma3 * radial_poly_a * u0_2
        + inv_sigma3 * radial_poly_a * rs_u0 * spin_sin2_th_u3
        + 2.0 * inv_two_sigma * sin_2th * spin2 * u1 * u2
        - r_over_sigma * u2_2 * (-r2 + rs_r + spin2)
        - sin2_th * u3_2_over_sigma3
            * (half_rs * sin2_th * spin2_r4
                - half_rs_r2 * spin4_sin2_th
                + one_minus_cos_4th_over_16 * r * rs2 * spin4
                - r7
                + r * spin6_cos4_th
                - r2 * rs_spin4_one_minus_cos_4th_over_16
                - 0.5 * r3 * rs2 * spin2_sin2_th
                + 2.0 * r3 * spin4_cos2_th
                - r3 * spin4_cos4_th
                + r5 * spin2
                - r5 * two_spin2_cos2_th
                + r6 * rs
                + rs * spin4_cos4_th_r2
                + rs_one_minus_cos_4th_over_16 * spin6
                + rs_r4 * two_spin2_cos2_th)
        - u1_2_over_radial_poly_b
            * (-half_rs * spin2_cos2_th
                + half_rs_r2
                + r * spin2
                + r * spin2_cos2_th);

    du[2] =
        cos_th_sin_th * spin2 * u1_2_over_radial_poly_b
        - cos_th_sin_th * u3_2_over_sigma3
            * (-0.25 * one_minus_cos_4th * rs_r * spin4
                - r4 * two_spin2_cos2_th
                - r6
                - rs_r * spin4_sin4_th
                - 2.0 * rs_r3 * spin2_sin2_th
                - spin2_r4
                - spin4_cos2_th * two_r2
                - spin4_cos4_th_r2
                - spin6_cos4_th)
        + inv_two_sigma * spin2_sin_2th * u2_2
        - 8.0 * r2_plus_spin2 * rs_r_over_two_sigma_cubed
            * sin_2th * spin_u3 * u0
        - r_over_sigma * two_u1 * u2
        + 4.0 * rs_r_over_two_sigma_cubed * spin2_sin_2th * u0_2;

    du[3] =
        2.0 * cot_th * inv_sigma2 * r * rs * spin * u0 * u2
        - 2.0 * cot_th * u2 * u3
            * (r4
                + rs_r * spin2_sin2_th
                - sin2_th * two_spin2_r2
                + spin4
                - 2.0 * spin4_sin2_th
                + spin4_sin4_th
                + two_spin2_r2)
            / (r4
                + spin2_cos2_th * two_r2
                + spin4_cos4_th)
        - inv_delta_sigma2 * two_u1 * u3
            * (-half_rs * spin2_r2_sin2_th
                + r * spin4_cos4_th
                + r3 * two_spin2_cos2_th
                + r5
                - rs * spin2_cos2_th_r2
                - rs_r4
                + rs_spin4_one_minus_cos_4th_over_16)
        - rs_u0_u1_over_delta_sigma2 * spin
            * (r2 - spin2_cos2_th);
}

rp_kerr_status rp_kerr_geodesic_physic_evaluate_null_tangent(
    double spin,
    const struct rp_kerr_geodesic_state *state,
    double energy,
    double axial_angular_momentum,
    double carter_constant,
    int radial_direction,
    int polar_direction,
    double tangent[RP_KERR_DIM]
)
{
    const double radius = state->x[1];
    const double radius_squared = radius * radius;
    const double spin_squared = spin * spin;
    const double sin_squared = state->sin_theta * state->sin_theta;
    const double cos_squared = state->cos_theta * state->cos_theta;
    const double p = energy * (radius_squared + spin_squared)
        - spin * axial_angular_momentum;
    const double shifted_angular_momentum =
        axial_angular_momentum - spin * energy;
    const double radial_second_term = state->delta
        * (shifted_angular_momentum * shifted_angular_momentum
            + carter_constant);
    const double radial_scale = p * p + fabs(radial_second_term);
    const double polar_parenthesis =
        axial_angular_momentum * axial_angular_momentum / sin_squared
        - spin_squared * energy * energy;
    const double polar_second_term = cos_squared * polar_parenthesis;
    /*
     * Theta is commonly supplied as the closest double to pi/2.  Its cosine
     * is then O(DBL_EPSILON), so scaling only by cos(theta)^2 would reject an
     * exactly equatorial Q = 0 ray because of representation roundoff.
     */
    const double polar_scale = fabs(carter_constant)
        + fabs(axial_angular_momentum * axial_angular_momentum / sin_squared)
        + fabs(spin_squared * energy * energy);
    double radial_potential = p * p - radial_second_term;
    double polar_potential = carter_constant - polar_second_term;

    if (!isfinite(p) || !isfinite(radial_second_term)
        || !isfinite(radial_scale) || !isfinite(polar_parenthesis)
        || !isfinite(polar_second_term) || !isfinite(polar_scale)
        || !isfinite(radial_potential) || !isfinite(polar_potential)) {
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    if (radial_potential < 0.0) {
        if (!rp_effectively_zero(radial_potential, radial_scale)) {
            return RP_KERR_STATUS_NO_REAL_NULL_TANGENT;
        }
        radial_potential = 0.0;
    }
    if (polar_potential < 0.0) {
        if (!rp_effectively_zero(polar_potential, polar_scale)) {
            return RP_KERR_STATUS_NO_REAL_NULL_TANGENT;
        }
        polar_potential = 0.0;
    }

    tangent[0] = (
        -spin * (spin * energy * sin_squared - axial_angular_momentum)
        + (radius_squared + spin_squared) * p / state->delta
    ) / state->sigma;
    tangent[1] = (double)radial_direction
        * sqrt(radial_potential) / state->sigma;
    tangent[2] = (double)polar_direction
        * sqrt(polar_potential) / state->sigma;
    tangent[3] = (
        -(spin * energy - axial_angular_momentum / sin_squared)
        + spin * p / state->delta
    ) / state->sigma;

    if (!rp_tensor_is_finite(tangent, RP_KERR_DIM)) {
        rp_tensor_zero(tangent, RP_KERR_DIM);
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    return RP_KERR_STATUS_OK;
}

static void evaluate_corrected_inverse_metric(
    double mass,
    double spin,
    const struct rp_kerr_geodesic_state *state,
    double inverse_metric[RP_KERR_DIM][RP_KERR_DIM]
)
{
    const double radius = state->x[1];
    const double radius_squared = radius * radius;
    const double spin_squared = spin * spin;
    const double sin_squared = state->sin_theta * state->sin_theta;
    const double radius_spin_squared = radius_squared + spin_squared;
    const double common_numerator = 2.0 * mass * radius;
    const double block_numerator = radius_spin_squared * radius_spin_squared
        - spin_squared * state->delta * sin_squared;

    rp_tensor_zero(
        &inverse_metric[0][0], RP_KERR_DIM * RP_KERR_DIM
    );
    inverse_metric[0][0] = -block_numerator
        / (state->sigma * state->delta);
    inverse_metric[0][3] = -common_numerator * spin
        / (state->sigma * state->delta);
    inverse_metric[3][0] = inverse_metric[0][3];
    inverse_metric[1][1] = state->delta / state->sigma;
    inverse_metric[2][2] = 1.0 / state->sigma;
    inverse_metric[3][3] = (state->delta - spin_squared * sin_squared)
        / (state->sigma * state->delta * sin_squared);
}

static void evaluate_corrected_metric_derivatives(
    double mass,
    double spin,
    const struct rp_kerr_geodesic_state *state,
    double derivative[RP_KERR_DIM][RP_KERR_DIM][RP_KERR_DIM]
)
{
    const double radius = state->x[1];
    const double radius_squared = radius * radius;
    const double spin_squared = spin * spin;
    const double sin_theta = state->sin_theta;
    const double cos_theta = state->cos_theta;
    const double sin_squared = sin_theta * sin_theta;
    const double two_sin_cos = 2.0 * sin_theta * cos_theta;
    const double sigma_radial = 2.0 * radius;
    const double sigma_polar = -2.0 * spin_squared * sin_theta * cos_theta;
    const double delta_radial = 2.0 * (radius - mass);
    const double common = 2.0 * mass * radius / state->sigma;
    const double common_radial = 2.0 * mass
        * (state->sigma - radius * sigma_radial)
        / (state->sigma * state->sigma);
    const double common_polar = -2.0 * mass * radius * sigma_polar
        / (state->sigma * state->sigma);
    const double azimuthal_factor = radius_squared + spin_squared
        + spin_squared * sin_squared * common;
    const double azimuthal_factor_radial = 2.0 * radius
        + spin_squared * sin_squared * common_radial;
    const double azimuthal_factor_polar = spin_squared
        * (two_sin_cos * common + sin_squared * common_polar);

    rp_tensor_zero(
        &derivative[0][0][0],
        RP_KERR_DIM * RP_KERR_DIM * RP_KERR_DIM
    );

    derivative[1][0][0] = common_radial;
    derivative[2][0][0] = common_polar;

    derivative[1][0][3] = -spin * sin_squared * common_radial;
    derivative[1][3][0] = derivative[1][0][3];
    derivative[2][0][3] = -spin
        * (two_sin_cos * common + sin_squared * common_polar);
    derivative[2][3][0] = derivative[2][0][3];

    derivative[1][1][1] = (
        sigma_radial * state->delta - state->sigma * delta_radial
    ) / (state->delta * state->delta);
    derivative[2][1][1] = sigma_polar / state->delta;

    derivative[1][2][2] = sigma_radial;
    derivative[2][2][2] = sigma_polar;

    derivative[1][3][3] = azimuthal_factor_radial * sin_squared;
    derivative[2][3][3] = azimuthal_factor_polar * sin_squared
        + azimuthal_factor * two_sin_cos;
}

void rp_kerr_geodesic_physic_evaluate_corrected_rhs(
    double mass,
    double spin,
    const struct rp_kerr_geodesic_state *state,
    double dx[RP_KERR_DIM],
    double du[RP_KERR_DIM]
)
{
    double inverse_metric[RP_KERR_DIM][RP_KERR_DIM];
    double metric_derivative[RP_KERR_DIM][RP_KERR_DIM][RP_KERR_DIM];
    double contracted_connection[RP_KERR_DIM];
    size_t alpha;
    size_t beta;
    size_t lambda;
    size_t sigma;

    evaluate_corrected_inverse_metric(mass, spin, state, inverse_metric);
    evaluate_corrected_metric_derivatives(
        mass, spin, state, metric_derivative
    );
    rp_tensor_zero(contracted_connection, RP_KERR_DIM);

    for (lambda = 0U; lambda < RP_KERR_DIM; ++lambda) {
        dx[lambda] = state->u[lambda];
        for (sigma = 0U; sigma < RP_KERR_DIM; ++sigma) {
            double lower_contraction = 0.0;

            for (alpha = 0U; alpha < RP_KERR_DIM; ++alpha) {
                for (beta = 0U; beta < RP_KERR_DIM; ++beta) {
                    lower_contraction += (
                        metric_derivative[alpha][sigma][beta]
                        - 0.5 * metric_derivative[sigma][alpha][beta]
                    ) * state->u[alpha] * state->u[beta];
                }
            }
            contracted_connection[lambda] +=
                inverse_metric[lambda][sigma] * lower_contraction;
        }
        du[lambda] = -contracted_connection[lambda];
    }
}

/**
 * Second derivatives of the same corrected metric differentiated above.
 * Blocks are (rr, r theta, theta theta); stationarity and axial symmetry
 * make every other coordinate block zero. No metric is reevaluated here.
 */
static void evaluate_corrected_metric_second_derivatives(
    double mass,
    double spin,
    const struct rp_kerr_geodesic_state *state,
    double second[3][RP_KERR_DIM][RP_KERR_DIM]
)
{
    const double radius = state->x[1];
    const double spin_squared = spin * spin;
    const double sine = state->sin_theta;
    const double cosine = state->cos_theta;
    const double sin_squared = sine * sine;
    const double sin_squared_polar = 2.0 * sine * cosine;
    const double sin_squared_polar_polar =
        2.0 * (cosine * cosine - sin_squared);
    const double sigma = state->sigma;
    const double sigma_radial = 2.0 * radius;
    const double sigma_polar = -spin_squared * sin_squared_polar;
    const double sigma_polar_polar =
        -spin_squared * sin_squared_polar_polar;
    const double delta = state->delta;
    const double delta_radial = 2.0 * (radius - mass);
    const double inverse_sigma_squared = 1.0 / (sigma * sigma);
    const double inverse_sigma_cubed = inverse_sigma_squared / sigma;
    const double common = 2.0 * mass * radius / sigma;
    const double common_radial = 2.0 * mass
        * (sigma - radius * sigma_radial) * inverse_sigma_squared;
    const double common_polar =
        -2.0 * mass * radius * sigma_polar * inverse_sigma_squared;
    const double common_radial_radial =
        -4.0 * mass * sigma_radial * inverse_sigma_squared
        -4.0 * mass * radius * inverse_sigma_squared
        +4.0 * mass * radius * sigma_radial * sigma_radial
            * inverse_sigma_cubed;
    const double common_radial_polar =
        -2.0 * mass * sigma_polar * inverse_sigma_squared
        +4.0 * mass * radius * sigma_radial * sigma_polar
            * inverse_sigma_cubed;
    const double common_polar_polar =
        -2.0 * mass * radius * sigma_polar_polar * inverse_sigma_squared
        +4.0 * mass * radius * sigma_polar * sigma_polar
            * inverse_sigma_cubed;
    const double inverse_delta_squared = 1.0 / (delta * delta);
    const double radius_spin_squared = radius * radius + spin_squared;

    rp_tensor_zero(&second[0][0][0], 3U * RP_KERR_DIM * RP_KERR_DIM);
    second[0][0][0] = common_radial_radial;
    second[1][0][0] = common_radial_polar;
    second[2][0][0] = common_polar_polar;

    second[0][0][3] = -spin * sin_squared * common_radial_radial;
    second[1][0][3] = -spin * (sin_squared_polar * common_radial
        + sin_squared * common_radial_polar);
    second[2][0][3] = -spin * (
        sin_squared_polar_polar * common
        + 2.0 * sin_squared_polar * common_polar
        + sin_squared * common_polar_polar
    );
    second[0][3][0] = second[0][0][3];
    second[1][3][0] = second[1][0][3];
    second[2][3][0] = second[2][0][3];

    second[0][1][1] = 2.0 / delta
        - (2.0 * sigma_radial * delta_radial + 2.0 * sigma)
            * inverse_delta_squared
        + 2.0 * sigma * delta_radial * delta_radial
            * inverse_delta_squared / delta;
    second[1][1][1] =
        -sigma_polar * delta_radial * inverse_delta_squared;
    second[2][1][1] = sigma_polar_polar / delta;
    second[0][2][2] = 2.0;
    second[2][2][2] = sigma_polar_polar;

    /* g_phi_phi = sin(theta)^2 (r^2 + a^2)
     *             + a^2 sin(theta)^4 common. */
    second[0][3][3] = 2.0 * sin_squared
        + spin_squared * sin_squared * sin_squared * common_radial_radial;
    second[1][3][3] = 2.0 * radius * sin_squared_polar
        + spin_squared * (
            2.0 * sin_squared * sin_squared_polar * common_radial
            + sin_squared * sin_squared * common_radial_polar
        );
    second[2][3][3] = sin_squared_polar_polar * radius_spin_squared
        + spin_squared * (
            2.0 * (sin_squared_polar * sin_squared_polar
                + sin_squared * sin_squared_polar_polar) * common
            + 4.0 * sin_squared * sin_squared_polar * common_polar
            + sin_squared * sin_squared * common_polar_polar
        );
}

void rp_kerr_geodesic_physic_evaluate_corrected_jacobian(
    double mass,
    double spin,
    const struct rp_kerr_geodesic_state *state,
    double jacobian[2 * RP_KERR_DIM][2 * RP_KERR_DIM]
)
{
    double inverse_metric[RP_KERR_DIM][RP_KERR_DIM];
    double derivative[RP_KERR_DIM][RP_KERR_DIM][RP_KERR_DIM];
    double second[3][RP_KERR_DIM][RP_KERR_DIM];
    double lower_contraction[RP_KERR_DIM] = {0.0};
    double raised_contraction[RP_KERR_DIM] = {0.0};
    double lower_velocity_derivative[RP_KERR_DIM][RP_KERR_DIM] = {{0.0}};
    size_t alpha;
    size_t beta;
    size_t lambda;
    size_t sigma;
    size_t column;

    evaluate_corrected_inverse_metric(mass, spin, state, inverse_metric);
    evaluate_corrected_metric_derivatives(mass, spin, state, derivative);
    evaluate_corrected_metric_second_derivatives(mass, spin, state, second);
    rp_tensor_zero(&jacobian[0][0], 4U * RP_KERR_DIM * RP_KERR_DIM);

    /* C_sigma = (d_alpha g_sigma_beta - d_sigma g_alpha_beta / 2)
     *           u_alpha u_beta, and acceleration = -g_inverse C.
     * Contract before raising indices; never materialize Gamma or d Gamma. */
    for (sigma = 0U; sigma < RP_KERR_DIM; ++sigma) {
        for (alpha = 0U; alpha < RP_KERR_DIM; ++alpha) {
            for (beta = 0U; beta < RP_KERR_DIM; ++beta) {
                lower_contraction[sigma] += (
                    derivative[alpha][sigma][beta]
                    - 0.5 * derivative[sigma][alpha][beta]
                ) * state->u[alpha] * state->u[beta];
            }
            for (column = 0U; column < RP_KERR_DIM; ++column) {
                lower_velocity_derivative[sigma][column] += (
                    derivative[column][sigma][alpha]
                    + derivative[alpha][sigma][column]
                    - derivative[sigma][column][alpha]
                ) * state->u[alpha];
            }
        }
    }
    for (lambda = 0U; lambda < RP_KERR_DIM; ++lambda) {
        jacobian[lambda][lambda + RP_KERR_DIM] = 1.0;
        for (sigma = 0U; sigma < RP_KERR_DIM; ++sigma) {
            raised_contraction[lambda] +=
                inverse_metric[lambda][sigma] * lower_contraction[sigma];
            for (column = 0U; column < RP_KERR_DIM; ++column) {
                jacobian[lambda + RP_KERR_DIM][column + RP_KERR_DIM] -=
                    inverse_metric[lambda][sigma]
                    * lower_velocity_derivative[sigma][column];
            }
        }
    }

    for (column = 1U; column <= 2U; ++column) {
        double coordinate_contraction[RP_KERR_DIM] = {0.0};
        for (sigma = 0U; sigma < RP_KERR_DIM; ++sigma) {
            for (alpha = 0U; alpha < RP_KERR_DIM; ++alpha) {
                for (beta = 0U; beta < RP_KERR_DIM; ++beta) {
                    double differentiated = 0.0;
                    if (alpha == 1U || alpha == 2U) {
                        differentiated +=
                            second[column + alpha - 2U][sigma][beta];
                    }
                    if (sigma == 1U || sigma == 2U) {
                        differentiated -=
                            0.5 * second[column + sigma - 2U][alpha][beta];
                    }
                    coordinate_contraction[sigma] +=
                        differentiated * state->u[alpha] * state->u[beta];
                }
            }
            /* d_k g_inverse = -g_inverse (d_k g) g_inverse. Apply this
             * identity to C, with the remaining inverse applied below. */
            for (alpha = 0U; alpha < RP_KERR_DIM; ++alpha) {
                coordinate_contraction[sigma] -= derivative[column][sigma][alpha]
                    * raised_contraction[alpha];
            }
        }
        for (lambda = 0U; lambda < RP_KERR_DIM; ++lambda) {
            for (sigma = 0U; sigma < RP_KERR_DIM; ++sigma) {
                jacobian[lambda + RP_KERR_DIM][column] -=
                    inverse_metric[lambda][sigma]
                    * coordinate_contraction[sigma];
            }
        }
    }
}
