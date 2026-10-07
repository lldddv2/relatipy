/**
 * @file initial.c
 * @brief Allocation-free Kerr null states, invariants and coordinate views.
 *
 * G = c = M = 1. All outputs are caller-owned; no pointer is retained.
 */

#include "relatipy/kerr_null.h"
#include "relatipy/kerr_geodesic.h"
#include "../initial/convert.h"
#include "../solution/reconstruct.h"
#include "../../utils/numeric.h"

#include <math.h>
#include <stdint.h>
#include <string.h>

static void clear_values(double *values, size_t count)
{
    size_t index;
    for (index = 0U; index < count; ++index) {
        values[index] = 0.0;
    }
}

static rp_kerr_status validate_spin(double spin)
{
    if (!isfinite(spin)) {
        return RP_KERR_STATUS_NONFINITE_INPUT;
    }
    return spin >= 0.0 && spin <= 1.0
        ? RP_KERR_STATUS_OK : RP_KERR_STATUS_INVALID_PARAMETER;
}

static double outer_horizon(double spin)
{
    return 1.0 + sqrt((1.0 - spin) * (1.0 + spin));
}

/* This chart check does not impose the stronger null-initial horizon margin. */
static rp_kerr_status validate_bl(double spin, const double *values, size_t count)
{
    rp_kerr_status status;
    size_t index;

    if (values == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    status = validate_spin(spin);
    if (status != RP_KERR_STATUS_OK) {
        return status;
    }
    for (index = 0U; index < count; ++index) {
        if (!isfinite(values[index])) {
            return RP_KERR_STATUS_NONFINITE_INPUT;
        }
    }
    if (!(values[1] > outer_horizon(spin))) {
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    if (rp_bl_polar_axis_singular(values[2])) {
        return RP_KERR_STATUS_COORDINATE_SINGULARITY;
    }
    return RP_KERR_STATUS_OK;
}

double rp_kerr_null_horizon_threshold(double spin)
{
    return validate_spin(spin) == RP_KERR_STATUS_OK
        ? outer_horizon(spin) * (1.0 + RP_KERR_NULL_HORIZON_MARGIN) : NAN;
}

rp_kerr_status rp_kerr_null_family_to_bl(
    double spin,
    rp_kerr_null_family family,
    const double row[RP_KERR_NULL_FAMILY_DIM],
    double bl[RP_KERR_NULL_FAMILY_DIM]
)
{
    double cartesian[RP_KERR_NULL_FAMILY_DIM];
    rp_kerr_status status;

    if (bl == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    clear_values(bl, RP_KERR_NULL_FAMILY_DIM);
    if (row == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    status = validate_spin(spin);
    if (status != RP_KERR_STATUS_OK) {
        return status;
    }
    switch (family) {
        case RP_KERR_NULL_FAMILY_CARTESIAN:
            return rp_initial_cartesian_to_bl(spin, row, bl);
        case RP_KERR_NULL_FAMILY_SPHERICAL:
            status = rp_initial_spherical_to_cartesian(row, cartesian);
            return status == RP_KERR_STATUS_OK
                ? rp_initial_cartesian_to_bl(spin, cartesian, bl) : status;
        case RP_KERR_NULL_FAMILY_BOYER_LINDQUIST:
            status = validate_bl(spin, row, RP_KERR_NULL_FAMILY_DIM);
            if (status == RP_KERR_STATUS_OK) {
                memcpy(bl, row, RP_KERR_NULL_FAMILY_DIM * sizeof(*bl));
            }
            return status;
        default:
            return RP_KERR_STATUS_INVALID_PARAMETER;
    }
}

rp_kerr_status rp_kerr_null_state_from_direction(
    double spin,
    rp_kerr_null_family family,
    const double row[RP_KERR_NULL_FAMILY_DIM],
    double state[RP_KERR_NULL_STATE_DIM]
)
{
    double bl[RP_KERR_NULL_FAMILY_DIM];
    double metric[RP_KERR_DIM][RP_KERR_DIM];
    double direction[3];
    double maximum;
    double a = 0.0;
    double b = 0.0;
    double c;
    double scale;
    double discriminant;
    double q;
    double roots[2];
    double root;
    rp_kerr_status status;
    size_t i;
    size_t j;

    if (state == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    clear_values(state, RP_KERR_NULL_STATE_DIM);
    status = rp_kerr_null_family_to_bl(spin, family, row, bl);
    if (status != RP_KERR_STATUS_OK) {
        return status;
    }
    if (!(bl[1] > rp_kerr_null_horizon_threshold(spin))) {
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    maximum = fmax(fabs(bl[4]), fmax(fabs(bl[5]), fabs(bl[6])));
    if (maximum == 0.0) {
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    /* Direction magnitude is irrelevant. Scaling prevents squared overflow. */
    for (i = 0U; i < 3U; ++i) {
        direction[i] = bl[i + RP_KERR_DIM] / maximum;
    }
    status = rp_kerr_metric(1.0, spin, bl, metric);
    if (status != RP_KERR_STATUS_OK) {
        return status;
    }
    for (i = 0U; i < 3U; ++i) {
        b += metric[0][i + 1U] * direction[i];
        for (j = 0U; j < 3U; ++j) {
            a += metric[i + 1U][j + 1U] * direction[i] * direction[j];
        }
    }
    c = metric[0][0];
    if (!(a > 0.0) || !isfinite(a) || !isfinite(b) || !isfinite(c)) {
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    /* Uniform coefficient scaling also protects B^2 - A*C from overflow. */
    scale = fmax(a, fmax(fabs(b), fabs(c)));
    a /= scale;
    b /= scale;
    c /= scale;
    discriminant = fma(-a, c, b * b);
    if (!isfinite(discriminant)) {
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    if (discriminant < 0.0) {
        return RP_KERR_STATUS_NO_REAL_NULL_TANGENT;
    }
    q = -(b + copysign(sqrt(discriminant), b));
    if (q == 0.0) {
        return RP_KERR_STATUS_NO_REAL_NULL_TANGENT;
    }
    roots[0] = q / a;
    roots[1] = discriminant == 0.0 ? roots[0] : c / q;
    if (!isfinite(roots[0]) || !isfinite(roots[1])) {
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    if (roots[0] > 0.0 && roots[1] > 0.0 && roots[0] != roots[1]) {
        /* Selecting between two ergoregion branches is not authorized. */
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    root = roots[0] > 0.0 ? roots[0] : roots[1];
    if (!(root > 0.0)) {
        return RP_KERR_STATUS_NO_REAL_NULL_TANGENT;
    }
    for (i = 0U; i < RP_KERR_DIM; ++i) {
        state[i] = bl[i];
    }
    state[4] = 1.0;
    for (i = 0U; i < 3U; ++i) {
        state[i + 5U] = root * direction[i];
        if (!isfinite(state[i + 5U])) {
            clear_values(state, RP_KERR_NULL_STATE_DIM);
            return RP_KERR_STATUS_NUMERICAL_RANGE;
        }
    }
    return RP_KERR_STATUS_OK;
}

rp_kerr_status rp_kerr_null_state_from_constants(
    double spin,
    const double bl_position[RP_KERR_DIM],
    double impact_parameter,
    double eta,
    int radial_sign,
    int polar_sign,
    double state[RP_KERR_NULL_STATE_DIM]
)
{
    double tangent[RP_KERR_DIM];
    rp_kerr_status status;
    size_t index;

    if (state == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    clear_values(state, RP_KERR_NULL_STATE_DIM);
    status = validate_bl(spin, bl_position, RP_KERR_DIM);
    if (status != RP_KERR_STATUS_OK) {
        return status;
    }
    if (!isfinite(impact_parameter) || !isfinite(eta)) {
        return RP_KERR_STATUS_NONFINITE_INPUT;
    }
    if (!(bl_position[1] > rp_kerr_null_horizon_threshold(spin))
        || (radial_sign != -1 && radial_sign != 1)
        || (polar_sign != -1 && polar_sign != 1)) {
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    status = rp_kerr_null_tangent(1.0, spin, bl_position, 1.0,
        impact_parameter, eta, radial_sign, polar_sign, tangent);
    if (status != RP_KERR_STATUS_OK) {
        return status;
    }
    if (!(tangent[0] > 0.0) || !isfinite(tangent[0])) {
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    for (index = 0U; index < RP_KERR_DIM; ++index) {
        state[index] = bl_position[index];
        state[index + RP_KERR_DIM] = tangent[index] / tangent[0];
        if (!isfinite(state[index + RP_KERR_DIM])) {
            clear_values(state, RP_KERR_NULL_STATE_DIM);
            return RP_KERR_STATUS_NUMERICAL_RANGE;
        }
    }
    state[4] = 1.0;
    return RP_KERR_STATUS_OK;
}

rp_kerr_status rp_kerr_null_invariants_evaluate(
    double spin,
    const double state[RP_KERR_NULL_STATE_DIM],
    rp_kerr_null_invariants *invariants
)
{
    double metric[RP_KERR_DIM][RP_KERR_DIM];
    long double covariant[RP_KERR_DIM] = {0.0L};
    long double norm = 0.0L;
    long double denominator = 0.0L;
    long double energy;
    long double angular_momentum;
    long double carter;
    const double *tangent;
    double sine;
    double cosine;
    rp_kerr_null_invariants result;
    rp_kerr_status status;
    size_t mu;
    size_t nu;

    if (invariants == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    memset(invariants, 0, sizeof(*invariants));
    status = validate_bl(spin, state, RP_KERR_NULL_STATE_DIM);
    if (status != RP_KERR_STATUS_OK) {
        return status;
    }
    status = rp_kerr_metric(1.0, spin, state, metric);
    if (status != RP_KERR_STATUS_OK) {
        return status;
    }
    tangent = state + RP_KERR_DIM;
    /* Wider intermediates retain cancellation and avoid artificial E^2
     * overflow/underflow. Every required final double is range-checked. */
    for (mu = 0U; mu < RP_KERR_DIM; ++mu) {
        for (nu = 0U; nu < RP_KERR_DIM; ++nu) {
            const long double term = (long double) metric[mu][nu]
                * tangent[mu] * tangent[nu];
            covariant[mu] += (long double) metric[mu][nu] * tangent[nu];
            norm += term;
            denominator += fabsl(term);
        }
    }
    if (denominator == 0.0L) {
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    energy = -covariant[0];
    angular_momentum = covariant[3];
    sine = sin(state[2]);
    cosine = cos(state[2]);
    carter = covariant[2] * covariant[2]
        + (long double) cosine * cosine
            * (angular_momentum * angular_momentum / ((long double) sine * sine)
                - (long double) spin * spin * energy * energy);
    result.energy = (double) energy;
    result.axial_angular_momentum = (double) angular_momentum;
    result.carter_constant = (double) carter;
    result.norm = (double) norm;
    result.relative_norm = (double) (fabsl(norm) / denominator);
    result.impact_parameter = energy == 0.0L ? NAN
        : (double) (angular_momentum / energy);
    result.eta = energy == 0.0L ? NAN : (double) (carter / (energy * energy));
    if (!isfinite(result.energy) || !isfinite(result.axial_angular_momentum)
        || !isfinite(result.carter_constant) || !isfinite(result.norm)
        || !isfinite(result.relative_norm)
        || (energy != 0.0L
            && (!isfinite(result.impact_parameter) || !isfinite(result.eta)))) {
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    *invariants = result;
    return RP_KERR_STATUS_OK;
}

rp_kerr_status rp_kerr_null_views_batch(
    double spin,
    rp_kerr_null_family family,
    const double *states,
    size_t count,
    double *views,
    rp_kerr_status *row_status
)
{
    rp_kerr_status global_status = RP_KERR_STATUS_OK;
    rp_kerr_status first_status = RP_KERR_STATUS_OK;
    size_t row;
    size_t column;

    if (count == 0U) {
        return RP_KERR_STATUS_OK;
    }
    if (views == NULL || row_status == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    if (count > SIZE_MAX / (RP_KERR_NULL_STATE_DIM * sizeof(*states))
        || count > SIZE_MAX / (RP_KERR_NULL_FAMILY_DIM * sizeof(*views))
        || count > SIZE_MAX / sizeof(*row_status)) {
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    if (states == NULL) {
        global_status = RP_KERR_STATUS_NULL_POINTER;
    } else {
        global_status = validate_spin(spin);
        if (global_status == RP_KERR_STATUS_OK
            && family != RP_KERR_NULL_FAMILY_CARTESIAN
            && family != RP_KERR_NULL_FAMILY_SPHERICAL
            && family != RP_KERR_NULL_FAMILY_BOYER_LINDQUIST) {
            global_status = RP_KERR_STATUS_INVALID_PARAMETER;
        }
    }
    for (row = 0U; row < count; ++row) {
        double *view = views + row * RP_KERR_NULL_FAMILY_DIM;
        rp_kerr_status status = global_status;
        clear_values(view, RP_KERR_NULL_FAMILY_DIM);
        if (status == RP_KERR_STATUS_OK) {
            const double *state = states + row * RP_KERR_NULL_STATE_DIM;
            status = validate_bl(spin, state, RP_KERR_NULL_STATE_DIM);
            if (status == RP_KERR_STATUS_OK && !(state[4] > 0.0)) {
                status = RP_KERR_STATUS_INVALID_PARAMETER;
            }
            if (status == RP_KERR_STATUS_OK) {
                if (family == RP_KERR_NULL_FAMILY_BOYER_LINDQUIST) {
                    for (column = 0U; column < RP_KERR_DIM; ++column) {
                        view[column] = state[column];
                    }
                    for (column = 0U; column < 3U; ++column) {
                        view[column + RP_KERR_DIM] = state[column + 5U] / state[4];
                        if (!isfinite(view[column + RP_KERR_DIM])) {
                            status = RP_KERR_STATUS_NUMERICAL_RANGE;
                        }
                    }
                } else {
                    const rp_solution_family selected =
                        family == RP_KERR_NULL_FAMILY_CARTESIAN
                            ? RP_SOLUTION_FAMILY_CARTESIAN
                            : RP_SOLUTION_FAMILY_SPHERICAL;
                    /* This path preserves the supplied tangent and never
                     * solves or checks timelike normalization. */
                    status = rp_solution_reconstruct_canonical_family_batch(
                        spin, state, 1U, selected, view, row_status + row);
                }
            }
        }
        if (status != RP_KERR_STATUS_OK) {
            clear_values(view, RP_KERR_NULL_FAMILY_DIM);
            if (first_status == RP_KERR_STATUS_OK) {
                first_status = status;
            }
        }
        row_status[row] = status;
    }
    return first_status;
}
