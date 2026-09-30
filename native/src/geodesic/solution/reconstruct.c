/**
 * @file reconstruct.c
 * @brief Allocation-free Kerr state reconstruction after Cartesian interpolation.
 */

#include "reconstruct.h"
#include "../initial/convert.h"

#include <float.h>
#include <math.h>
#include <stdint.h>

static void clear_row(double *row)
{
    size_t column;

    for (column = 0U; column < RP_SOLUTION_RECONSTRUCTED_DIM; ++column) {
        row[column] = 0.0;
    }
}

static double positive_angle(double angle)
{
    const double full_turn = 2.0 * acos(-1.0);

    if (fabs(angle) <= 64.0 * DBL_EPSILON) {
        return 0.0;
    }
    if (angle < 0.0) {
        angle += full_turn;
    }
    return full_turn - angle <= 64.0 * DBL_EPSILON ? 0.0 : angle;
}

/* Osculating Newtonian conic in the spin-aligned Cartesian frame, mu = 1. */
static rp_kerr_status reconstruct_elements(
    double x,
    double y,
    double z,
    double vx,
    double vy,
    double vz,
    double distance,
    double *output
)
{
    const double angular_momentum[3] = {
        y * vz - z * vy,
        z * vx - x * vz,
        x * vy - y * vx
    };
    double h;
    double h_scale;
    double node;
    double eccentricity_vector[3];
    double eccentricity;
    double energy;
    double semimajor;
    double inclination;
    double ascending_node;
    double node_p[3];
    double node_q[3];
    double peri_p[3];
    double peri_q[3];
    double periapsis_argument;
    double anomaly;
    const double circular_tolerance = 128.0 * DBL_EPSILON;
    size_t column;

    h = hypot(hypot(angular_momentum[0], angular_momentum[1]),
        angular_momentum[2]);
    h_scale = distance * hypot(hypot(vx, vy), vz);
    if (!isfinite(h) || !isfinite(h_scale)) {
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    energy = 0.5 * (vx * vx + vy * vy + vz * vz) - 1.0 / distance;
    if (!isfinite(energy)) {
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    semimajor = energy == 0.0 ? INFINITY : -0.5 / energy;
    eccentricity_vector[0] = vy * angular_momentum[2]
        - vz * angular_momentum[1] - x / distance;
    eccentricity_vector[1] = vz * angular_momentum[0]
        - vx * angular_momentum[2] - y / distance;
    eccentricity_vector[2] = vx * angular_momentum[1]
        - vy * angular_momentum[0] - z / distance;
    eccentricity = hypot(hypot(eccentricity_vector[0],
        eccentricity_vector[1]), eccentricity_vector[2]);
    if ((energy != 0.0 && (!isfinite(semimajor) || semimajor == 0.0))
        || !isfinite(eccentricity)) {
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    output[RP_SOL_SEMIMAJOR] = semimajor;
    output[RP_SOL_ECCENTRICITY] = eccentricity;
    if (!(h > 128.0 * DBL_EPSILON * h_scale)) {
        /* An unresolved orbital plane leaves all four angles undefined. */
        output[RP_SOL_INCLINATION] = NAN;
        output[RP_SOL_ASCENDING_NODE] = NAN;
        output[RP_SOL_PERIAPSIS_ARGUMENT] = NAN;
        output[RP_SOL_TRUE_ANOMALY] = NAN;
        return RP_KERR_STATUS_OK;
    }
    node = hypot(angular_momentum[0], angular_momentum[1]);
    inclination = atan2(node, angular_momentum[2]);
    if (node <= circular_tolerance * h) {
        ascending_node = 0.0;
    } else {
        ascending_node = positive_angle(atan2(
            angular_momentum[0], -angular_momentum[1]
        ));
    }
    node_p[0] = cos(ascending_node);
    node_p[1] = sin(ascending_node);
    node_p[2] = 0.0;
    node_q[0] = -sin(ascending_node) * cos(inclination);
    node_q[1] = cos(ascending_node) * cos(inclination);
    node_q[2] = sin(inclination);

    if (eccentricity <= circular_tolerance) {
        eccentricity = 0.0;
        periapsis_argument = 0.0;
        for (column = 0U; column < 3U; ++column) {
            peri_p[column] = node_p[column];
            peri_q[column] = node_q[column];
        }
    } else {
        periapsis_argument = positive_angle(atan2(
            eccentricity_vector[0] * node_q[0]
                + eccentricity_vector[1] * node_q[1]
                + eccentricity_vector[2] * node_q[2],
            eccentricity_vector[0] * node_p[0]
                + eccentricity_vector[1] * node_p[1]
                + eccentricity_vector[2] * node_p[2]
        ));
        for (column = 0U; column < 3U; ++column) {
            peri_p[column] = eccentricity_vector[column] / eccentricity;
        }
        peri_q[0] = angular_momentum[1] / h * peri_p[2]
            - angular_momentum[2] / h * peri_p[1];
        peri_q[1] = angular_momentum[2] / h * peri_p[0]
            - angular_momentum[0] / h * peri_p[2];
        peri_q[2] = angular_momentum[0] / h * peri_p[1]
            - angular_momentum[1] / h * peri_p[0];
    }
    anomaly = positive_angle(atan2(
        x * peri_q[0] + y * peri_q[1] + z * peri_q[2],
        x * peri_p[0] + y * peri_p[1] + z * peri_p[2]
    ));

    if (!isfinite(inclination) || !isfinite(ascending_node)
        || !isfinite(periapsis_argument) || !isfinite(anomaly)) {
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    output[RP_SOL_ECCENTRICITY] = eccentricity;
    output[RP_SOL_INCLINATION] = inclination;
    output[RP_SOL_ASCENDING_NODE] = ascending_node;
    output[RP_SOL_PERIAPSIS_ARGUMENT] = periapsis_argument;
    output[RP_SOL_TRUE_ANOMALY] = anomaly;
    return RP_KERR_STATUS_OK;
}

static rp_kerr_status assemble_row(
    const double *input,
    const double canonical[RP_INITIAL_CANONICAL_DIM],
    double *output
)
{
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
    double transverse;
    double transverse_velocity;
    double spherical_radius;
    double spherical_velocity[3];
    double velocity[3];
    rp_kerr_status status;
    size_t column;

    for (column = 0U; column < RP_SOLUTION_CARTESIAN_DIM; ++column) {
        if (!isfinite(input[column])) {
            return RP_KERR_STATUS_NONFINITE_INPUT;
        }
    }

    x = input[1];
    y = input[2];
    z = input[3];
    vx = input[4];
    vy = input[5];
    vz = input[6];
    transverse = hypot(x, y);
    if (!isfinite(transverse)) {
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    if (transverse == 0.0) {
        return RP_KERR_STATUS_COORDINATE_SINGULARITY;
    }

    transverse_velocity = (x / transverse) * vx + (y / transverse) * vy;
    velocity[0] = canonical[5] / canonical[4];
    velocity[1] = canonical[6] / canonical[4];
    velocity[2] = canonical[7] / canonical[4];
    for (column = 0U; column < 3U; ++column) {
        if (!isfinite(velocity[column])) {
            return RP_KERR_STATUS_NUMERICAL_RANGE;
        }
    }

    spherical_radius = hypot(transverse, z);
    if (!(spherical_radius > 0.0) || !isfinite(spherical_radius)) {
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    spherical_velocity[0] = (transverse * transverse_velocity + z * vz)
        / spherical_radius;
    spherical_velocity[1] = (z * transverse_velocity - transverse * vz)
        / (spherical_radius * spherical_radius);
    spherical_velocity[2] = velocity[2];
    for (column = 0U; column < 3U; ++column) {
        if (!isfinite(spherical_velocity[column])) {
            return RP_KERR_STATUS_NUMERICAL_RANGE;
        }
    }

    status = reconstruct_elements(
        x, y, z, vx, vy, vz, spherical_radius, output
    );
    if (status != RP_KERR_STATUS_OK) {
        clear_row(output);
        return status;
    }

    for (column = 0U; column < RP_KERR_DIM; ++column) {
        output[column] = canonical[column];
        output[column + 7U] = canonical[column + RP_KERR_DIM];
    }
    for (column = 0U; column < 3U; ++column) {
        output[column + 4U] = velocity[column];
        output[RP_SOL_SPH_VR + column] = spherical_velocity[column];
        output[RP_SOL_SPH_UR + column] = canonical[4]
            * spherical_velocity[column];
    }
    output[RP_SOL_SPH_R] = spherical_radius;
    output[RP_SOL_SPH_THETA] = atan2(transverse, z);
    output[RP_SOL_SPH_PHI] = canonical[3];
    output[RP_SOL_UX] = canonical[4] * vx;
    output[RP_SOL_UY] = canonical[4] * vy;
    output[RP_SOL_UZ] = canonical[4] * vz;
    /* Elements were validated above, with only their declared sentinels allowed. */
    for (column = 0U; column < RP_SOL_SEMIMAJOR; ++column) {
        if (!isfinite(output[column])) {
            clear_row(output);
            return RP_KERR_STATUS_NUMERICAL_RANGE;
        }
    }
    return RP_KERR_STATUS_OK;
}

static rp_kerr_status reconstruct_one(
    double spin,
    const double *input,
    double *output
)
{
    double canonical[RP_INITIAL_CANONICAL_DIM];
    rp_kerr_status status = rp_initial_cartesian_to_canonical(
        spin, input, canonical
    );
    if (status != RP_KERR_STATUS_OK) {
        return status;
    }
    return assemble_row(input, canonical, output);
}

rp_kerr_status rp_solution_reconstruct_batch(
    double spin,
    const double *cartesian,
    size_t count,
    double *reconstructed,
    rp_kerr_status *row_status
)
{
    rp_kerr_status global_status = RP_KERR_STATUS_OK;
    rp_kerr_status first_status = RP_KERR_STATUS_OK;
    rp_kerr_status status;
    size_t row;
    double *output;

    if (count == 0U) {
        return RP_KERR_STATUS_OK;
    }
    if (reconstructed == NULL || row_status == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    if (count > SIZE_MAX / RP_SOLUTION_RECONSTRUCTED_DIM
        || count > SIZE_MAX / RP_SOLUTION_CARTESIAN_DIM) {
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    if (cartesian == NULL) {
        global_status = RP_KERR_STATUS_NULL_POINTER;
    } else if (!isfinite(spin)) {
        global_status = RP_KERR_STATUS_NONFINITE_INPUT;
    } else if (spin < 0.0 || spin > 1.0) {
        global_status = RP_KERR_STATUS_INVALID_PARAMETER;
    }

    for (row = 0U; row < count; ++row) {
        output = reconstructed + row * RP_SOLUTION_RECONSTRUCTED_DIM;
        clear_row(output);
        status = global_status == RP_KERR_STATUS_OK
            ? reconstruct_one(
                spin,
                cartesian + row * RP_SOLUTION_CARTESIAN_DIM,
                output
            )
            : global_status;
        row_status[row] = status;
        if (first_status == RP_KERR_STATUS_OK && status != RP_KERR_STATUS_OK) {
            first_status = status;
        }
    }
    return first_status;
}

rp_kerr_status rp_solution_reconstruct_canonical_batch(
    double spin,
    const double *canonical,
    size_t count,
    double *reconstructed,
    double *cartesian,
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
    if (reconstructed == NULL || cartesian == NULL || row_status == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    if (count > SIZE_MAX / RP_SOLUTION_RECONSTRUCTED_DIM
        || count > SIZE_MAX / RP_SOLUTION_CARTESIAN_DIM
        || count > SIZE_MAX / RP_INITIAL_CANONICAL_DIM) {
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    if (canonical == NULL) {
        global_status = RP_KERR_STATUS_NULL_POINTER;
    } else if (!isfinite(spin)) {
        global_status = RP_KERR_STATUS_NONFINITE_INPUT;
    } else if (spin < 0.0 || spin > 1.0) {
        global_status = RP_KERR_STATUS_INVALID_PARAMETER;
    }
    for (row = 0U; row < count; ++row) {
        double *output = reconstructed + row * RP_SOLUTION_RECONSTRUCTED_DIM;
        double *cartesian_row = cartesian + row * RP_SOLUTION_CARTESIAN_DIM;
        rp_kerr_status status = global_status;
        clear_row(output);
        for (column = 0U; column < RP_SOLUTION_CARTESIAN_DIM; ++column) {
            cartesian_row[column] = 0.0;
        }
        if (status == RP_KERR_STATUS_OK) {
            const double *state = canonical + row * RP_INITIAL_CANONICAL_DIM;
            status = rp_initial_canonical_to_cartesian(spin, state, cartesian_row);
            if (status == RP_KERR_STATUS_OK) {
                /* Integrated u is preserved, including its numerical drift.
                 * Re-normalizing coordinate velocities here hides errors and
                 * loses valid partial output near the BL horizon. */
                status = assemble_row(cartesian_row, state, output);
            }
        }
        if (status != RP_KERR_STATUS_OK) {
            clear_row(output);
            for (column = 0U; column < RP_SOLUTION_CARTESIAN_DIM; ++column) {
                cartesian_row[column] = 0.0;
            }
            if (first_status == RP_KERR_STATUS_OK) {
                first_status = status;
            }
        }
        row_status[row] = status;
    }
    return first_status;
}

rp_kerr_status rp_solution_reconstruct_canonical_family_batch(
    double spin,
    const double *canonical,
    size_t count,
    rp_solution_family family,
    double *output,
    rp_kerr_status *row_status
)
{
    rp_kerr_status global_status = RP_KERR_STATUS_OK;
    rp_kerr_status first_status = RP_KERR_STATUS_OK;
    const size_t width = family == RP_SOLUTION_FAMILY_ELEMENTS
        ? RP_SOLUTION_ELEMENTS_DIM : RP_SOLUTION_CARTESIAN_DIM;
    size_t row;
    size_t column;

    if (count == 0U) {
        return RP_KERR_STATUS_OK;
    }
    if (output == NULL || row_status == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    if (count > SIZE_MAX / width || count > SIZE_MAX / RP_INITIAL_CANONICAL_DIM) {
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    if (canonical == NULL) {
        global_status = RP_KERR_STATUS_NULL_POINTER;
    } else if (!isfinite(spin)) {
        global_status = RP_KERR_STATUS_NONFINITE_INPUT;
    } else if (spin < 0.0 || spin > 1.0
        || (family != RP_SOLUTION_FAMILY_CARTESIAN
            && family != RP_SOLUTION_FAMILY_SPHERICAL
            && family != RP_SOLUTION_FAMILY_ELEMENTS)) {
        global_status = RP_KERR_STATUS_INVALID_PARAMETER;
    }
    for (row = 0U; row < count; ++row) {
        double *selected = output + row * width;
        rp_kerr_status status = global_status;
        for (column = 0U; column < width; ++column) {
            selected[column] = 0.0;
        }
        if (status == RP_KERR_STATUS_OK) {
            const double *state = canonical + row * RP_INITIAL_CANONICAL_DIM;
            double cartesian[RP_SOLUTION_CARTESIAN_DIM];
            status = rp_initial_canonical_to_cartesian(spin, state, cartesian);
            if (status == RP_KERR_STATUS_OK) {
                if (family == RP_SOLUTION_FAMILY_CARTESIAN) {
                    for (column = 0U; column < width; ++column) {
                        selected[column] = cartesian[column];
                    }
                } else {
                    const double x = cartesian[1];
                    const double y = cartesian[2];
                    const double z = cartesian[3];
                    const double transverse = hypot(x, y);
                    const double distance = hypot(transverse, z);
                    if (!isfinite(transverse) || !isfinite(distance)) {
                        status = RP_KERR_STATUS_NUMERICAL_RANGE;
                    } else if (transverse == 0.0) {
                        status = RP_KERR_STATUS_COORDINATE_SINGULARITY;
                    } else if (!(distance > 0.0)) {
                        status = RP_KERR_STATUS_NUMERICAL_RANGE;
                    } else if (family == RP_SOLUTION_FAMILY_SPHERICAL) {
                        const double radial_velocity =
                            (x / transverse) * cartesian[4]
                            + (y / transverse) * cartesian[5];
                        const double vr = (transverse * radial_velocity
                            + z * cartesian[6]) / distance;
                        const double vtheta = (z * radial_velocity
                            - transverse * cartesian[6]) / (distance * distance);
                        const double vphi = state[7] / state[4];
                        if (!isfinite(vr) || !isfinite(vtheta) || !isfinite(vphi)) {
                            status = RP_KERR_STATUS_NUMERICAL_RANGE;
                        } else {
                            selected[0] = state[0];
                            selected[1] = distance;
                            selected[2] = atan2(transverse, z);
                            selected[3] = state[3];
                            selected[4] = vr;
                            selected[5] = vtheta;
                            selected[6] = vphi;
                        }
                    } else {
                        double elements[RP_SOLUTION_RECONSTRUCTED_DIM] = {0.0};
                        status = reconstruct_elements(x, y, z, cartesian[4],
                            cartesian[5], cartesian[6], distance, elements);
                        if (status == RP_KERR_STATUS_OK) {
                            for (column = 0U; column < width; ++column) {
                                selected[column] = elements[RP_SOL_SEMIMAJOR + column];
                            }
                        }
                    }
                }
            }
        }
        if (status != RP_KERR_STATUS_OK) {
            for (column = 0U; column < width; ++column) {
                selected[column] = 0.0;
            }
            if (first_status == RP_KERR_STATUS_OK) {
                first_status = status;
            }
        }
        row_status[row] = status;
    }
    return first_status;
}
