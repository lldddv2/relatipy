#include "kerr_properties.h"

#include <math.h>

static double clamp_unit(double value)
{
    return value < -1.0 ? -1.0 : (value > 1.0 ? 1.0 : value);
}

rp_kerr_status rp_kerr_characteristic_radii(double spin, double radii[6])
{
    double z1;
    double z2;
    double radical;
    double shift;
    size_t i;

    if (radii == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    for (i = 0U; i < 6U; ++i) {
        radii[i] = 0.0;
    }
    if (!isfinite(spin) || spin < 0.0 || spin > 1.0) {
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    radical = sqrt((1.0 - spin) * (1.0 + spin));
    z1 = 1.0 + cbrt((1.0 - spin) * (1.0 + spin))
        * (cbrt(1.0 + spin) + cbrt(1.0 - spin));
    z2 = sqrt(3.0 * spin * spin + z1 * z1);
    shift = sqrt(fmax(0.0, (3.0 - z1) * (3.0 + z1 + 2.0 * z2)));
    radii[0] = 1.0 + radical;
    radii[1] = 1.0 - radical;
    radii[2] = 3.0 + z2 - shift;
    radii[3] = 3.0 + z2 + shift;
    radii[4] = 2.0 * (1.0 + cos((2.0 / 3.0) * acos(clamp_unit(-spin))));
    radii[5] = 2.0 * (1.0 + cos((2.0 / 3.0) * acos(clamp_unit(spin))));
    return RP_KERR_STATUS_OK;
}

rp_kerr_status rp_kerr_ergosurface_radii(
    double spin, const double *theta, size_t count, double *radii
)
{
    size_t i;

    if (count != 0U && (theta == NULL || radii == NULL)) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    if (!isfinite(spin) || spin < 0.0 || spin > 1.0) {
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    for (i = 0U; i < count; ++i) {
        double sine;
        double radicand;
        if (!isfinite(theta[i]) || theta[i] < 0.0 || theta[i] > acos(-1.0)) {
            return RP_KERR_STATUS_INVALID_PARAMETER;
        }
        sine = sin(theta[i]);
        radicand = (1.0 - spin) * (1.0 + spin) + spin * spin * sine * sine;
        radii[i] = 1.0 + sqrt(radicand);
    }
    return RP_KERR_STATUS_OK;
}

rp_kerr_status rp_kerr_surface_profile(
    double spin, rp_kerr_surface surface, const double *theta, size_t count,
    double *rho, double *z
)
{
    double characteristic[6];
    double cartesian_radius = 0.0;
    double radius = 0.0;
    rp_kerr_status status;
    size_t i;

    if (count != 0U && (theta == NULL || rho == NULL || z == NULL)) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    if (!isfinite(spin) || spin < 0.0 || spin > 1.0) {
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    if (surface != RP_KERR_SURFACE_OUTER_HORIZON
        && surface != RP_KERR_SURFACE_ERGOSURFACE) {
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    for (i = 0U; i < count; ++i) {
        if (!isfinite(theta[i]) || theta[i] < 0.0 || theta[i] > acos(-1.0)) {
            return RP_KERR_STATUS_INVALID_PARAMETER;
        }
    }
    if (surface == RP_KERR_SURFACE_OUTER_HORIZON) {
        status = rp_kerr_characteristic_radii(spin, characteristic);
        if (status != RP_KERR_STATUS_OK) {
            return status;
        }
        radius = characteristic[0];
        cartesian_radius = hypot(radius, spin);
    }
    for (i = 0U; i < count; ++i) {
        const double angle = theta[i];
        if (surface == RP_KERR_SURFACE_ERGOSURFACE) {
            status = rp_kerr_ergosurface_radii(spin, &angle, 1U, &radius);
            if (status != RP_KERR_STATUS_OK) {
                return status;
            }
            cartesian_radius = hypot(radius, spin);
        }
        rho[i] = cartesian_radius * sin(angle);
        z[i] = radius * cos(angle);
    }
    return RP_KERR_STATUS_OK;
}
