/**
 * @file convert.c
 * @brief Kerr oblate Cartesian/Boyer--Lindquist state conversion.
 */

#include "convert.h"
#include "../../utils/numeric.h"

#include <math.h>
#include <stddef.h>

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

static rp_kerr_status cartesian_to_bl_kinematics(
    double spin,
    const double cartesian[RP_INITIAL_CARTESIAN_DIM],
    double coordinates[RP_KERR_DIM],
    double velocity[3]
)
{
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
    double transverse;
    double distance_squared;
    double discriminant_base;
    double discriminant;
    double radius_squared;
    double radius;
    double radial_scale;
    double sin_theta;
    double cos_theta;
    double sigma;
    double transverse_velocity;
    rp_kerr_status status;
    size_t index;

    if (cartesian == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    status = validate_spin(spin);
    if (status != RP_KERR_STATUS_OK) {
        return status;
    }
    for (index = 0U; index < RP_INITIAL_CARTESIAN_DIM; ++index) {
        if (!isfinite(cartesian[index])) {
            return RP_KERR_STATUS_NONFINITE_INPUT;
        }
    }

    x = cartesian[1];
    y = cartesian[2];
    z = cartesian[3];
    vx = cartesian[4];
    vy = cartesian[5];
    vz = cartesian[6];
    transverse = hypot(x, y);
    if (!isfinite(transverse)) {
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    if (transverse == 0.0) {
        return RP_KERR_STATUS_COORDINATE_SINGULARITY;
    }
    distance_squared = transverse * transverse + z * z;
    if (!isfinite(distance_squared)) {
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    discriminant_base = distance_squared - spin * spin;
    discriminant = hypot(discriminant_base, 2.0 * spin * z);
    if (!isfinite(discriminant)) {
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    if (discriminant_base >= 0.0) {
        radius_squared = 0.5 * (discriminant_base + discriminant);
    } else {
        radius_squared = 2.0 * spin * spin * z * z
            / (discriminant - discriminant_base);
    }
    radius = sqrt(radius_squared);
    radial_scale = hypot(radius, spin);
    if (!(radius > 0.0) || !isfinite(radial_scale)) {
        return RP_KERR_STATUS_PHYSICAL_SINGULARITY;
    }
    if (radius <= outer_horizon(spin)) {
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    sin_theta = transverse / radial_scale;
    cos_theta = z / radius;
    sigma = radius_squared + spin * spin * cos_theta * cos_theta;
    if (!(sigma > 0.0) || !isfinite(sigma)) {
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }

    coordinates[0] = cartesian[0];
    coordinates[1] = radius;
    coordinates[2] = atan2(sin_theta, cos_theta);
    coordinates[3] = atan2(y, x);
    if (rp_bl_polar_axis_singular(coordinates[2])) {
        return RP_KERR_STATUS_COORDINATE_SINGULARITY;
    }
    transverse_velocity = (x / transverse) * vx + (y / transverse) * vy;
    velocity[0] = radial_scale / sigma
        * (radius * sin_theta * transverse_velocity
            + radial_scale * cos_theta * vz);
    velocity[1] = (radial_scale * cos_theta * transverse_velocity
        - radius * sin_theta * vz) / sigma;
    velocity[2] = ((x / transverse) * vy - (y / transverse) * vx)
        / transverse;
    for (index = 0U; index < 3U; ++index) {
        if (!isfinite(velocity[index])) {
            return RP_KERR_STATUS_NUMERICAL_RANGE;
        }
    }

    return RP_KERR_STATUS_OK;
}

/** Pure chart kinematics shared by timelike and null initial-state builders. */
rp_kerr_status rp_initial_cartesian_to_bl(
    double spin,
    const double cartesian[RP_INITIAL_CARTESIAN_DIM],
    double bl[RP_INITIAL_CARTESIAN_DIM]
)
{
    double coordinates[RP_KERR_DIM];
    double velocity[3];
    rp_kerr_status status;
    size_t index;

    if (bl == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    clear_values(bl, RP_INITIAL_CARTESIAN_DIM);
    status = cartesian_to_bl_kinematics(spin, cartesian, coordinates, velocity);
    if (status != RP_KERR_STATUS_OK) {
        return status;
    }
    for (index = 0U; index < RP_KERR_DIM; ++index) {
        bl[index] = coordinates[index];
    }
    for (index = 0U; index < 3U; ++index) {
        bl[index + RP_KERR_DIM] = velocity[index];
    }
    return RP_KERR_STATUS_OK;
}

rp_kerr_status rp_initial_cartesian_to_canonical(
    double spin,
    const double cartesian[RP_INITIAL_CARTESIAN_DIM],
    double canonical[RP_INITIAL_CANONICAL_DIM]
)
{
    double coordinates[RP_KERR_DIM];
    double velocity[3];
    double four_velocity[RP_KERR_DIM];
    rp_kerr_status status;
    size_t index;

    if (canonical == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    clear_values(canonical, RP_INITIAL_CANONICAL_DIM);
    status = cartesian_to_bl_kinematics(spin, cartesian, coordinates, velocity);
    if (status != RP_KERR_STATUS_OK) {
        return status;
    }
    status = rp_kerr_four_velocity(
        1.0, spin, coordinates, velocity, four_velocity
    );
    if (status != RP_KERR_STATUS_OK) {
        return status;
    }
    for (index = 0U; index < RP_KERR_DIM; ++index) {
        canonical[index] = coordinates[index];
        canonical[index + RP_KERR_DIM] = four_velocity[index];
    }
    return RP_KERR_STATUS_OK;
}

rp_kerr_status rp_initial_canonical_to_cartesian(
    double spin,
    const double canonical[RP_INITIAL_CANONICAL_DIM],
    double cartesian[RP_INITIAL_CARTESIAN_DIM]
)
{
    double radius;
    double theta;
    double phi;
    double radial_scale;
    double sin_theta;
    double cos_theta;
    double sin_phi;
    double cos_phi;
    double vr;
    double vtheta;
    double vphi;
    double radial_scale_dot;
    rp_kerr_status status;
    size_t index;

    if (cartesian == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    clear_values(cartesian, RP_INITIAL_CARTESIAN_DIM);
    if (canonical == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    status = validate_spin(spin);
    if (status != RP_KERR_STATUS_OK) {
        return status;
    }
    for (index = 0U; index < RP_INITIAL_CANONICAL_DIM; ++index) {
        if (!isfinite(canonical[index])) {
            return RP_KERR_STATUS_NONFINITE_INPUT;
        }
    }

    radius = canonical[1];
    theta = canonical[2];
    phi = canonical[3];
    if (radius <= outer_horizon(spin)) {
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    if (rp_bl_polar_axis_singular(theta)) {
        return RP_KERR_STATUS_COORDINATE_SINGULARITY;
    }
    if (!(canonical[4] > 0.0)) {
        return RP_KERR_STATUS_NON_TIMELIKE_VELOCITY;
    }
    radial_scale = hypot(radius, spin);
    sin_theta = sin(theta);
    cos_theta = cos(theta);
    sin_phi = sin(phi);
    cos_phi = cos(phi);
    vr = canonical[5] / canonical[4];
    vtheta = canonical[6] / canonical[4];
    vphi = canonical[7] / canonical[4];
    radial_scale_dot = radius / radial_scale * vr;

    cartesian[0] = canonical[0];
    cartesian[1] = radial_scale * sin_theta * cos_phi;
    cartesian[2] = radial_scale * sin_theta * sin_phi;
    cartesian[3] = radius * cos_theta;
    cartesian[4] = radial_scale_dot * sin_theta * cos_phi
        + radial_scale * cos_theta * cos_phi * vtheta
        - radial_scale * sin_theta * sin_phi * vphi;
    cartesian[5] = radial_scale_dot * sin_theta * sin_phi
        + radial_scale * cos_theta * sin_phi * vtheta
        + radial_scale * sin_theta * cos_phi * vphi;
    cartesian[6] = vr * cos_theta - radius * sin_theta * vtheta;
    for (index = 0U; index < RP_INITIAL_CARTESIAN_DIM; ++index) {
        if (!isfinite(cartesian[index])) {
            clear_values(cartesian, RP_INITIAL_CARTESIAN_DIM);
            return RP_KERR_STATUS_NUMERICAL_RANGE;
        }
    }
    return RP_KERR_STATUS_OK;
}

rp_kerr_status rp_initial_from_bl(
    double spin,
    const double bl[RP_INITIAL_CARTESIAN_DIM],
    double canonical[RP_INITIAL_CANONICAL_DIM]
)
{
    double four_velocity[RP_KERR_DIM];
    rp_kerr_status status;
    size_t index;

    if (canonical == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    clear_values(canonical, RP_INITIAL_CANONICAL_DIM);
    if (bl == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    status = validate_spin(spin);
    if (status != RP_KERR_STATUS_OK) {
        return status;
    }
    for (index = 0U; index < RP_INITIAL_CARTESIAN_DIM; ++index) {
        if (!isfinite(bl[index])) {
            return RP_KERR_STATUS_NONFINITE_INPUT;
        }
    }
    if (bl[1] <= outer_horizon(spin)) {
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    if (rp_bl_polar_axis_singular(bl[2])) {
        return RP_KERR_STATUS_COORDINATE_SINGULARITY;
    }
    status = rp_kerr_four_velocity(
        1.0, spin, bl, bl + RP_KERR_DIM, four_velocity
    );
    if (status != RP_KERR_STATUS_OK) {
        return status;
    }
    for (index = 0U; index < RP_KERR_DIM; ++index) {
        canonical[index] = bl[index];
        canonical[index + RP_KERR_DIM] = four_velocity[index];
    }
    return RP_KERR_STATUS_OK;
}

rp_kerr_status rp_initial_spherical_to_cartesian(
    const double spherical[RP_INITIAL_CARTESIAN_DIM],
    double cartesian[RP_INITIAL_CARTESIAN_DIM]
)
{
    double radius;
    double theta;
    double phi;
    double vr;
    double vtheta;
    double vphi;
    double sin_theta;
    double cos_theta;
    double sin_phi;
    double cos_phi;
    size_t index;

    if (cartesian == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    clear_values(cartesian, RP_INITIAL_CARTESIAN_DIM);
    if (spherical == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    for (index = 0U; index < RP_INITIAL_CARTESIAN_DIM; ++index) {
        if (!isfinite(spherical[index])) {
            return RP_KERR_STATUS_NONFINITE_INPUT;
        }
    }
    radius = spherical[1];
    theta = spherical[2];
    phi = spherical[3];
    vr = spherical[4];
    vtheta = spherical[5];
    vphi = spherical[6];
    if (!(radius > 0.0)) {
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    if (rp_bl_polar_axis_singular(theta)) {
        return RP_KERR_STATUS_COORDINATE_SINGULARITY;
    }
    sin_theta = sin(theta);
    cos_theta = cos(theta);
    sin_phi = sin(phi);
    cos_phi = cos(phi);
    cartesian[0] = spherical[0];
    cartesian[1] = radius * sin_theta * cos_phi;
    cartesian[2] = radius * sin_theta * sin_phi;
    cartesian[3] = radius * cos_theta;
    cartesian[4] = vr * sin_theta * cos_phi
        + radius * vtheta * cos_theta * cos_phi
        - radius * vphi * sin_theta * sin_phi;
    cartesian[5] = vr * sin_theta * sin_phi
        + radius * vtheta * cos_theta * sin_phi
        + radius * vphi * sin_theta * cos_phi;
    cartesian[6] = vr * cos_theta - radius * vtheta * sin_theta;
    for (index = 0U; index < RP_INITIAL_CARTESIAN_DIM; ++index) {
        if (!isfinite(cartesian[index])) {
            clear_values(cartesian, RP_INITIAL_CARTESIAN_DIM);
            return RP_KERR_STATUS_NUMERICAL_RANGE;
        }
    }
    return RP_KERR_STATUS_OK;
}

rp_kerr_status rp_initial_from_spherical(
    double spin,
    const double spherical[RP_INITIAL_CARTESIAN_DIM],
    double canonical[RP_INITIAL_CANONICAL_DIM]
)
{
    double cartesian[RP_INITIAL_CARTESIAN_DIM];
    rp_kerr_status status;

    if (canonical == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    clear_values(canonical, RP_INITIAL_CANONICAL_DIM);
    status = validate_spin(spin);
    if (status != RP_KERR_STATUS_OK) {
        return status;
    }
    status = rp_initial_spherical_to_cartesian(spherical, cartesian);
    if (status != RP_KERR_STATUS_OK) {
        return status;
    }
    return rp_initial_cartesian_to_canonical(spin, cartesian, canonical);
}

rp_kerr_status rp_initial_elements_to_cartesian(
    const double elements[RP_INITIAL_CARTESIAN_DIM],
    double cartesian[RP_INITIAL_CARTESIAN_DIM]
)
{
    double semimajor;
    double eccentricity;
    double inclination;
    double ascending_node;
    double periapsis_argument;
    double anomaly;
    double semilatus;
    double denominator;
    double radius;
    double speed_scale;
    double cos_node;
    double sin_node;
    double cos_inc;
    double sin_inc;
    double cos_peri;
    double sin_peri;
    double cos_anomaly;
    double sin_anomaly;
    double peri_p[3];
    double peri_q[3];
    size_t index;

    if (cartesian == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    clear_values(cartesian, RP_INITIAL_CARTESIAN_DIM);
    if (elements == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    for (index = 0U; index < RP_INITIAL_CARTESIAN_DIM; ++index) {
        if (!isfinite(elements[index])) {
            return RP_KERR_STATUS_NONFINITE_INPUT;
        }
    }
    semimajor = elements[1];
    eccentricity = elements[2];
    inclination = elements[3];
    ascending_node = elements[4];
    periapsis_argument = elements[5];
    anomaly = elements[6];
    if (!((semimajor > 0.0 && eccentricity >= 0.0
                && eccentricity < 1.0)
            || (semimajor < 0.0 && eccentricity > 1.0))) {
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    if (!(inclination >= 0.0 && inclination <= acos(-1.0))) {
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    semilatus = semimajor * (1.0 - eccentricity * eccentricity);
    cos_anomaly = cos(anomaly);
    sin_anomaly = sin(anomaly);
    denominator = 1.0 + eccentricity * cos_anomaly;
    if (!(semilatus > 0.0 && denominator > 0.0)) {
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    radius = semilatus / denominator;
    speed_scale = 1.0 / sqrt(semilatus);
    if (!isfinite(radius) || !isfinite(speed_scale)) {
        return RP_KERR_STATUS_NUMERICAL_RANGE;
    }
    cos_node = cos(ascending_node);
    sin_node = sin(ascending_node);
    cos_inc = cos(inclination);
    sin_inc = sin(inclination);
    cos_peri = cos(periapsis_argument);
    sin_peri = sin(periapsis_argument);
    peri_p[0] = cos_node * cos_peri - sin_node * sin_peri * cos_inc;
    peri_p[1] = sin_node * cos_peri + cos_node * sin_peri * cos_inc;
    peri_p[2] = sin_peri * sin_inc;
    peri_q[0] = -cos_node * sin_peri - sin_node * cos_peri * cos_inc;
    peri_q[1] = -sin_node * sin_peri + cos_node * cos_peri * cos_inc;
    peri_q[2] = cos_peri * sin_inc;
    cartesian[0] = elements[0];
    for (index = 0U; index < 3U; ++index) {
        cartesian[1U + index] = radius
            * (cos_anomaly * peri_p[index]
                + sin_anomaly * peri_q[index]);
        cartesian[4U + index] = speed_scale
            * (-sin_anomaly * peri_p[index]
                + (eccentricity + cos_anomaly) * peri_q[index]);
        if (!isfinite(cartesian[1U + index])
            || !isfinite(cartesian[4U + index])) {
            clear_values(cartesian, RP_INITIAL_CARTESIAN_DIM);
            return RP_KERR_STATUS_NUMERICAL_RANGE;
        }
    }
    return RP_KERR_STATUS_OK;
}

rp_kerr_status rp_initial_from_elements(
    double spin,
    const double elements[RP_INITIAL_CARTESIAN_DIM],
    double canonical[RP_INITIAL_CANONICAL_DIM]
)
{
    double cartesian[RP_INITIAL_CARTESIAN_DIM];
    rp_kerr_status status;

    if (canonical == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    clear_values(canonical, RP_INITIAL_CANONICAL_DIM);
    status = validate_spin(spin);
    if (status != RP_KERR_STATUS_OK) {
        return status;
    }
    status = rp_initial_elements_to_cartesian(elements, cartesian);
    if (status != RP_KERR_STATUS_OK) {
        return status;
    }
    return rp_initial_cartesian_to_canonical(spin, cartesian, canonical);
}

rp_kerr_status rp_initial_observer_cartesian_to_canonical(
    double spin,
    const double rotation[3][3],
    const double observer[RP_INITIAL_CARTESIAN_DIM],
    double canonical[RP_INITIAL_CANONICAL_DIM]
)
{
    double body[RP_INITIAL_CARTESIAN_DIM];
    size_t index;
    size_t component;

    if (canonical == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    clear_values(canonical, RP_INITIAL_CANONICAL_DIM);
    if (rotation == NULL || observer == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    body[0] = observer[0];
    for (index = 0U; index < 3U; ++index) {
        body[1U + index] = 0.0;
        body[4U + index] = 0.0;
        for (component = 0U; component < 3U; ++component) {
            body[1U + index] += rotation[component][index]
                * observer[1U + component];
            body[4U + index] += rotation[component][index]
                * observer[4U + component];
        }
    }
    return rp_initial_cartesian_to_canonical(spin, body, canonical);
}
