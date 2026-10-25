/**
 * @file kerr.c
 * @brief ODE callback adapters for the Kerr geodesic evaluators.
 */

#include "kerr.h"

#include "relatipy/kerr_geodesic.h"
#include "relatipy/kerr_geometry.h"
#include "geodesic/jacobian.h"

#include <float.h>
#include <math.h>
#include <stddef.h>

void rp_kerr_timelike_constants_from_momentum(
    double spin,
    double theta,
    const double momentum[4],
    double constants[3]
)
{
    const double sine = sin(theta);
    const double cosine = cos(theta);
    const double energy = -momentum[0];
    const double angular_momentum = momentum[3];

    constants[0] = energy;
    constants[1] = angular_momentum;
    constants[2] = momentum[2] * momentum[2]
        + cosine * cosine * (spin * spin * (1.0 - energy * energy)
            + angular_momentum * angular_momentum / (sine * sine));
}

int rp_kerr_timelike_constants(
    double mass,
    double spin,
    const double state[8],
    double constants[3]
)
{
    double metric[RP_KERR_DIM][RP_KERR_DIM];
    double momentum[RP_KERR_DIM] = {0.0, 0.0, 0.0, 0.0};
    rp_kerr_status status;
    size_t mu;
    size_t nu;

    if (constants == NULL) {
        return (int)RP_KERR_STATUS_NULL_POINTER;
    }
    constants[0] = constants[1] = constants[2] = 0.0;
    if (state == NULL) {
        return (int)RP_KERR_STATUS_NULL_POINTER;
    }
    for (mu = 0U; mu < 2U * RP_KERR_DIM; ++mu) {
        if (!isfinite(state[mu])) {
            return (int)RP_KERR_STATUS_NONFINITE_INPUT;
        }
    }
    status = rp_kerr_metric(mass, spin, state, metric);
    if (status != RP_KERR_STATUS_OK) {
        return (int)status;
    }
    for (mu = 0U; mu < RP_KERR_DIM; ++mu) {
        for (nu = 0U; nu < RP_KERR_DIM; ++nu) {
            momentum[mu] += metric[mu][nu] * state[RP_KERR_DIM + nu];
        }
    }
    rp_kerr_timelike_constants_from_momentum(spin, state[2], momentum, constants);
    if (!isfinite(constants[0]) || !isfinite(constants[1])
        || !isfinite(constants[2])) {
        constants[0] = constants[1] = constants[2] = 0.0;
        return (int)RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    return 0;
}

int rp_kerr_integrator_prepare_projection(
    rp_kerr_integrator_context *context,
    const double state[8]
)
{
    double metric[RP_KERR_DIM][RP_KERR_DIM];
    double momentum[RP_KERR_DIM] = {0.0, 0.0, 0.0, 0.0};
    double norm = 0.0;
    double constants[3];
    double horizon;
    rp_kerr_status status;
    size_t mu;
    size_t nu;

    if (context == NULL || state == NULL) {
        return (int)RP_KERR_STATUS_NULL_POINTER;
    }
    for (mu = 0U; mu < 2U * RP_KERR_DIM; ++mu) {
        if (!isfinite(state[mu])) {
            return (int)RP_KERR_STATUS_NONFINITE_INPUT;
        }
    }
    status = rp_kerr_metric(context->mass, context->spin, state, metric);
    if (status != RP_KERR_STATUS_OK) {
        return (int)status;
    }
    horizon = context->mass
        + sqrt((context->mass - fabs(context->spin))
            * (context->mass + fabs(context->spin)));
    if (!isfinite(horizon)) {
        return (int)RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    if (state[1] <= horizon) {
        return (int)RP_KERR_STATUS_INVALID_PARAMETER;
    }
    for (mu = 0U; mu < RP_KERR_DIM; ++mu) {
        for (nu = 0U; nu < RP_KERR_DIM; ++nu) {
            momentum[mu] += metric[mu][nu] * state[RP_KERR_DIM + nu];
        }
        norm += momentum[mu] * state[RP_KERR_DIM + mu];
    }
    if (!isfinite(norm) || fabs(norm + 1.0) > 1.0e-8) {
        return (int)RP_KERR_STATUS_NON_TIMELIKE_VELOCITY;
    }
    rp_kerr_timelike_constants_from_momentum(
        context->spin, state[2], momentum, constants
    );
    if (!isfinite(constants[0]) || !isfinite(constants[1])
        || !isfinite(constants[2])) {
        return (int)RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    context->E0 = constants[0];
    context->Lz0 = constants[1];
    context->Q0 = constants[2];
    context->projection_ready = 1;
    return 0;
}

/* Projection free variables: (r, theta, u^t, u^r, u^theta, u^phi), packed
 * state indices 1, 2, 4, 5, 6, 7. t and phi do not enter the invariants. */
static const size_t PROJECTION_FREE[6] = {1U, 2U, 4U, 5U, 6U, 7U};
/* J J^T pivots below this fraction of its largest diagonal are dropped. */
#define PROJECTION_PIVOT_CUTOFF 1.0e-12
#define PROJECTION_MAX_ITERATIONS 16

/*
 * Residuals (norm + 1, E - E0, Lz - Lz0, Q - Q0) of y, the summed magnitude
 * of the terms of each (its rounding scale), and their analytic Jacobian
 * with respect to y, in one pass over the Boyer--Lindquist metric and its r
 * and theta derivatives. The caller has validated the point.
 */
static void projection_residuals(
    const rp_kerr_integrator_context *context,
    const double y[6],
    double residual[4],
    double magnitude[4],
    double jacobian[4][6]
)
{
    const double mass = context->mass;
    const double spin = context->spin;
    const double r = y[0];
    const double ut = y[2];
    const double ur = y[3];
    const double uth = y[4];
    const double uph = y[5];
    const double s = sin(y[1]);
    const double c = cos(y[1]);
    const double a2 = spin * spin;
    const double sigma = r * r + a2 * c * c;
    const double sigma_r = 2.0 * r;
    const double sigma_th = -2.0 * a2 * s * c;
    const double delta = r * r - 2.0 * mass * r + a2;
    const double delta_r = 2.0 * r - 2.0 * mass;
    const double inv = 1.0 / sigma;
    const double inv2 = inv * inv;
    const double gtt = -1.0 + 2.0 * mass * r * inv;
    const double gtt_r = 2.0 * mass * (inv - r * sigma_r * inv2);
    const double gtt_th = -2.0 * mass * r * sigma_th * inv2;
    const double gtp = -2.0 * mass * spin * r * s * s * inv;
    const double gtp_r = -2.0 * mass * spin * s * s * (inv - r * sigma_r * inv2);
    const double gtp_th = -2.0 * mass * spin * r
        * (2.0 * s * c * inv - s * s * sigma_th * inv2);
    const double grr = sigma / delta;
    const double grr_r = sigma_r / delta - sigma * delta_r / (delta * delta);
    const double grr_th = sigma_th / delta;
    const double outer = r * r + a2 + 2.0 * mass * a2 * r * s * s * inv;
    const double outer_r = 2.0 * r
        + 2.0 * mass * a2 * s * s * (inv - r * sigma_r * inv2);
    const double outer_th = 2.0 * mass * a2 * r
        * (2.0 * s * c * inv - s * s * sigma_th * inv2);
    const double gpp = outer * s * s;
    const double gpp_r = outer_r * s * s;
    const double gpp_th = outer_th * s * s + 2.0 * outer * s * c;
    const double energy = -(gtt * ut + gtp * uph);
    const double angular_momentum = gtp * ut + gpp * uph;
    const double polar_momentum = sigma * uth;
    const double bracket = a2 * (1.0 - energy * energy)
        + angular_momentum * angular_momentum / (s * s);
    const double dq_denergy = -2.0 * a2 * energy * c * c;
    const double dq_dmomentum = 2.0 * angular_momentum * c * c / (s * s);
    size_t j;

    residual[0] = gtt * ut * ut + 2.0 * gtp * ut * uph + grr * ur * ur
        + sigma * uth * uth + gpp * uph * uph + 1.0;
    residual[1] = energy - context->E0;
    residual[2] = angular_momentum - context->Lz0;
    residual[3] = polar_momentum * polar_momentum + c * c * bracket - context->Q0;
    magnitude[0] = fabs(gtt * ut * ut) + 2.0 * fabs(gtp * ut * uph)
        + fabs(grr * ur * ur) + fabs(sigma * uth * uth) + fabs(gpp * uph * uph)
        + 1.0;
    /* Unit floor: E, Lz and Q are of order one per unit mass with M = 1,
     * and Q vanishes on equatorial orbits. */
    magnitude[1] = fabs(gtt * ut) + fabs(gtp * uph) + 1.0;
    magnitude[2] = fabs(gtp * ut) + fabs(gpp * uph) + 1.0;
    magnitude[3] = polar_momentum * polar_momentum
        + c * c * (a2 * (1.0 + energy * energy)
            + angular_momentum * angular_momentum / (s * s))
        + fabs(context->Q0) + 1.0;

    jacobian[0][0] = gtt_r * ut * ut + 2.0 * gtp_r * ut * uph
        + grr_r * ur * ur + sigma_r * uth * uth + gpp_r * uph * uph;
    jacobian[0][1] = gtt_th * ut * ut + 2.0 * gtp_th * ut * uph
        + grr_th * ur * ur + sigma_th * uth * uth + gpp_th * uph * uph;
    jacobian[0][2] = 2.0 * (gtt * ut + gtp * uph);
    jacobian[0][3] = 2.0 * grr * ur;
    jacobian[0][4] = 2.0 * sigma * uth;
    jacobian[0][5] = 2.0 * (gtp * ut + gpp * uph);
    jacobian[1][0] = -(gtt_r * ut + gtp_r * uph);
    jacobian[1][1] = -(gtt_th * ut + gtp_th * uph);
    jacobian[1][2] = -gtt;
    jacobian[1][3] = 0.0;
    jacobian[1][4] = 0.0;
    jacobian[1][5] = -gtp;
    jacobian[2][0] = gtp_r * ut + gpp_r * uph;
    jacobian[2][1] = gtp_th * ut + gpp_th * uph;
    jacobian[2][2] = gtp;
    jacobian[2][3] = 0.0;
    jacobian[2][4] = 0.0;
    jacobian[2][5] = gpp;
    for (j = 0U; j < 6U; ++j) {
        jacobian[3][j] = dq_denergy * jacobian[1][j]
            + dq_dmomentum * jacobian[2][j];
    }
    jacobian[3][0] += 2.0 * polar_momentum * sigma_r * uth;
    jacobian[3][1] += 2.0 * polar_momentum * sigma_th * uth
        - 2.0 * s * c * bracket
        - 2.0 * c * c * c * angular_momentum * angular_momentum / (s * s * s);
    jacobian[3][4] += 2.0 * polar_momentum * sigma;
}

/* Largest residual relative to the rounding scale of its own terms. */
static double projection_error(
    const double residual[4],
    const double magnitude[4]
)
{
    double error = 0.0;
    size_t i;

    for (i = 0U; i < 4U; ++i) {
        error = fmax(error, fabs(residual[i]) / magnitude[i]);
    }
    return isfinite(error) ? error : INFINITY;
}

/* Every residual at the rounding level of its own terms. */
#define PROJECTION_CONVERGED (8.0 * DBL_EPSILON)
/* A stagnating Newton iteration is accepted below this relative error. */
#define PROJECTION_ACCEPTED 1.0e-10

int rp_kerr_integrator_project(double tau, double state[], void *opaque_context)
{
    const rp_kerr_integrator_context *context = opaque_context;
    double metric[RP_KERR_DIM][RP_KERR_DIM];
    double inverse[RP_KERR_DIM][RP_KERR_DIM];
    double y[6];
    double candidate[6];
    double residual[4];
    double magnitude[4];
    double jacobian[4][6];
    double scratch[4][6];
    double gram[4][4];
    double lower[4][4] = {{0.0}};
    double pivot[4];
    double z[4];
    double increment[4];
    double lambda[4] = {0.0, 0.0, 0.0, 0.0};
    double largest = 0.0;
    double previous_error;
    double horizon;
    rp_kerr_status status;
    size_t i;
    size_t j;
    size_t k;
    int iteration;

    if (context == NULL || state == NULL) {
        return (int)RP_KERR_STATUS_NULL_POINTER;
    }
    if (!isfinite(tau)) {
        return (int)RP_KERR_STATUS_NONFINITE_INPUT;
    }
    if (context->projection_ready != 1 || !(context->mass > 0.0)
        || fabs(context->spin) > context->mass) {
        return (int)RP_KERR_STATUS_INVALID_PARAMETER;
    }
    if (!isfinite(context->mass) || !isfinite(context->spin)
        || !isfinite(context->E0) || !isfinite(context->Lz0)
        || !isfinite(context->Q0)) {
        return (int)RP_KERR_STATUS_NONFINITE_INPUT;
    }
    for (i = 0U; i < 2U * RP_KERR_DIM; ++i) {
        if (!isfinite(state[i])) {
            return (int)RP_KERR_STATUS_NONFINITE_INPUT;
        }
    }

    horizon = context->mass
        + sqrt((context->mass - fabs(context->spin))
            * (context->mass + fabs(context->spin)));
    if (!isfinite(horizon)) {
        return (int)RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    if (state[1] <= horizon) {
        return 0;
    }
    /* Validates the chart point (axis, Delta and ring singularities). */
    status = rp_kerr_metric(context->mass, context->spin, state, metric);
    if (status != RP_KERR_STATUS_OK) {
        return (int)status;
    }

    for (i = 0U; i < 6U; ++i) {
        y[i] = state[PROJECTION_FREE[i]];
        candidate[i] = y[i];
    }
    projection_residuals(context, y, residual, magnitude, jacobian);
    if (projection_error(residual, magnitude) <= PROJECTION_CONVERGED) {
        return 0;
    }

    /* E and Lz are linear in (u^t, u^phi) at fixed position: restore them
     * exactly first, which also resolves large velocity errors at a
     * degenerate (circular) point where the Newton step below cannot. */
    status = rp_kerr_inverse_metric(
        context->mass, context->spin, state, inverse
    );
    if (status != RP_KERR_STATUS_OK) {
        return (int)status;
    }
    y[2] = -inverse[0][0] * context->E0 + inverse[0][3] * context->Lz0;
    y[5] = -inverse[3][0] * context->E0 + inverse[3][3] * context->Lz0;
    if (!isfinite(y[2]) || !isfinite(y[5])) {
        return (int)RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    for (i = 0U; i < 6U; ++i) {
        candidate[i] = y[i];
    }
    projection_residuals(context, y, residual, magnitude, jacobian);
    previous_error = projection_error(residual, magnitude);
    if (previous_error <= PROJECTION_CONVERGED) {
        for (i = 0U; i < 6U; ++i) {
            state[PROJECTION_FREE[i]] = y[i];
        }
        return 0;
    }

    /*
     * Standard orthogonal projection (Hairer, Lubich & Wanner, Geometric
     * Numerical Integration, sec. IV.4): y = y~ + J(y~)^T lambda with
     * simplified Newton iterations on lambda. J J^T is factored once as
     * L D L^T; a negligible pivot drops its constraint direction, so a
     * degenerate invariant set (Q on the equator, circular orbits) does
     * not amplify rounding into the state.
     */
    for (i = 0U; i < 4U; ++i) {
        for (k = 0U; k <= i; ++k) {
            gram[i][k] = 0.0;
            for (j = 0U; j < 6U; ++j) {
                gram[i][k] += jacobian[i][j] * jacobian[k][j];
            }
        }
        largest = fmax(largest, gram[i][i]);
    }
    for (i = 0U; i < 4U; ++i) {
        double value = gram[i][i];

        for (k = 0U; k < i; ++k) {
            value -= lower[i][k] * lower[i][k] * pivot[k];
        }
        pivot[i] = value > PROJECTION_PIVOT_CUTOFF * largest ? value : 0.0;
        lower[i][i] = 1.0;
        for (j = i + 1U; j < 4U; ++j) {
            double entry = gram[j][i];

            for (k = 0U; k < i; ++k) {
                entry -= lower[j][k] * lower[i][k] * pivot[k];
            }
            lower[j][i] = pivot[i] != 0.0 ? entry / pivot[i] : 0.0;
        }
    }

    for (iteration = 0; iteration < PROJECTION_MAX_ITERATIONS; ++iteration) {
        double error;

        for (i = 0U; i < 4U; ++i) {
            z[i] = -residual[i];
            for (k = 0U; k < i; ++k) {
                z[i] -= lower[i][k] * z[k];
            }
        }
        for (i = 0U; i < 4U; ++i) {
            z[i] = pivot[i] != 0.0 ? z[i] / pivot[i] : 0.0;
        }
        for (i = 4U; i-- > 0U;) {
            increment[i] = z[i];
            for (k = i + 1U; k < 4U; ++k) {
                increment[i] -= lower[k][i] * increment[k];
            }
            lambda[i] += increment[i];
        }
        for (i = 0U; i < 6U; ++i) {
            candidate[i] = y[i];
            for (k = 0U; k < 4U; ++k) {
                candidate[i] += jacobian[k][i] * lambda[k];
            }
            if (!isfinite(candidate[i])) {
                return (int)RP_KERR_STATUS_NON_TIMELIKE_VELOCITY;
            }
        }
        if (candidate[0] <= horizon || !(sin(candidate[1]) > 0.0)) {
            return (int)RP_KERR_STATUS_NON_TIMELIKE_VELOCITY;
        }
        projection_residuals(context, candidate, residual, magnitude, scratch);
        error = projection_error(residual, magnitude);
        /* Converged, or stagnating at the rounding floor of a dropped
         * (degenerate) direction with a small error: accept. */
        if (error <= PROJECTION_CONVERGED
            || (error > 0.5 * previous_error && error <= PROJECTION_ACCEPTED)) {
            for (i = 0U; i < 6U; ++i) {
                state[PROJECTION_FREE[i]] = candidate[i];
            }
            return 0;
        }
        previous_error = error;
    }
    if (previous_error <= PROJECTION_ACCEPTED) {
        for (i = 0U; i < 6U; ++i) {
            state[PROJECTION_FREE[i]] = candidate[i];
        }
        return 0;
    }
    return (int)RP_KERR_STATUS_NON_TIMELIKE_VELOCITY;
}

int rp_kerr_integrator_jacobian(
    double affine_parameter,
    const double state[],
    double *jacobian,
    size_t dimension,
    size_t row_stride,
    void *opaque_context
)
{
    const rp_kerr_integrator_context *context = opaque_context;
    double values[2 * RP_KERR_DIM][2 * RP_KERR_DIM];
    rp_kerr_status status;
    size_t row;
    size_t column;

    (void)affine_parameter;
    if (context == NULL || state == NULL || jacobian == NULL) {
        return (int)RP_KERR_STATUS_NULL_POINTER;
    }
    if (dimension != 2U * RP_KERR_DIM || row_stride < dimension) {
        return (int)RP_KERR_STATUS_INVALID_PARAMETER;
    }
    status = rp_kerr_geodesic_jacobian(
        context->mass, context->spin, state, state + RP_KERR_DIM, values
    );
    if (status != RP_KERR_STATUS_OK) {
        return (int)status;
    }
    for (row = 0U; row < dimension; ++row) {
        for (column = 0U; column < dimension; ++column) {
            jacobian[row * row_stride + column] = values[row][column];
        }
    }
    return 0;
}

int rp_kerr_integrator_rhs(
    double affine_parameter,
    const double state[],
    double derivative[],
    void *opaque_context
)
{
    const rp_kerr_integrator_context *context = opaque_context;

    (void)affine_parameter;
    if (context == NULL) {
        return (int)RP_KERR_STATUS_NULL_POINTER;
    }
    /* The corrected metric contraction applies to timelike tangents too. */
    return (int)rp_kerr_null_geodesic_rhs(
        context->mass,
        context->spin,
        state,
        state + RP_KERR_DIM,
        derivative,
        derivative + RP_KERR_DIM
    );
}

int rp_kerr_null_integrator_rhs(
    double affine_parameter,
    const double state[],
    double derivative[],
    void *opaque_context
)
{
    const rp_kerr_integrator_context *context = opaque_context;

    (void)affine_parameter;
    if (context == NULL) {
        return (int)RP_KERR_STATUS_NULL_POINTER;
    }
    return (int)rp_kerr_null_geodesic_rhs(
        context->mass,
        context->spin,
        state,
        state + RP_KERR_DIM,
        derivative,
        derivative + RP_KERR_DIM
    );
}
