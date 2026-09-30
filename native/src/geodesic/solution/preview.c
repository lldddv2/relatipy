#include "preview.h"

#include "geodesic/initial/convert.h"
#include "geodesic/solution/reconstruct.h"
#include "metric/kerr_properties.h"

#include <math.h>
#include <stdint.h>

rp_kerr_status rp_solution_osculating_preview(
    double spin, const double canonical[8], size_t count,
    double *xyz, double references[3]
)
{
    double cartesian[7];
    double row[RP_SOLUTION_RECONSTRUCTED_DIM];
    double radii[6];
    double a, e, inc, node, omega, p, limit;
    double cp, sp, ci, si, cn, sn;
    double basis_p[3], basis_q[3];
    rp_kerr_status status, row_status;
    size_t i, j;

    if (canonical == NULL || xyz == NULL || references == NULL) {
        return RP_KERR_STATUS_NULL_POINTER;
    }
    if (count < 2U || count > SIZE_MAX / (3U * sizeof(double))) {
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    status = rp_initial_canonical_to_cartesian(spin, canonical, cartesian);
    if (status != RP_KERR_STATUS_OK) {
        return status;
    }
    status = rp_solution_reconstruct_batch(spin, cartesian, 1U, row, &row_status);
    if (status != RP_KERR_STATUS_OK) {
        return status;
    }
    status = rp_kerr_characteristic_radii(spin, radii);
    if (status != RP_KERR_STATUS_OK) {
        return status;
    }
    for (j = 0U; j < 3U; ++j) {
        size_t index = j == 0U ? 0U : (j == 1U ? 2U : 3U);
        references[j] = hypot(radii[index], spin);
    }
    a = row[RP_SOL_SEMIMAJOR];
    e = row[RP_SOL_ECCENTRICITY];
    inc = row[RP_SOL_INCLINATION];
    node = row[RP_SOL_ASCENDING_NODE];
    omega = row[RP_SOL_PERIAPSIS_ARGUMENT];
    p = a * (1.0 - e * e);
    if (!(p > 0.0) || !isfinite(p)) {
        return RP_KERR_STATUS_INVALID_PARAMETER;
    }
    cp = cos(omega); sp = sin(omega);
    ci = cos(inc); si = sin(inc);
    cn = cos(node); sn = sin(node);
    basis_p[0] = cn * cp - sn * sp * ci;
    basis_p[1] = sn * cp + cn * sp * ci;
    basis_p[2] = sp * si;
    basis_q[0] = -cn * sp - sn * cp * ci;
    basis_q[1] = -sn * sp + cn * cp * ci;
    basis_q[2] = cp * si;
    limit = e < 1.0 ? 2.0 * acos(-1.0) : 0.9 * acos(-1.0 / e);
    for (i = 0U; i < count; ++i) {
        double anomaly = e < 1.0
            ? limit * (double)i / (double)(count - 1U)
            : -limit + 2.0 * limit * (double)i / (double)(count - 1U);
        double cosine = cos(anomaly);
        double sine = sin(anomaly);
        double radius = p / (1.0 + e * cosine);
        if (!(radius > 0.0) || !isfinite(radius)) {
            return RP_KERR_STATUS_NUMERICAL_RANGE;
        }
        for (j = 0U; j < 3U; ++j) {
            xyz[3U * i + j] = radius
                * (cosine * basis_p[j] + sine * basis_q[j]);
            if (!isfinite(xyz[3U * i + j])) {
                return RP_KERR_STATUS_NUMERICAL_RANGE;
            }
        }
    }
    return RP_KERR_STATUS_OK;
}
