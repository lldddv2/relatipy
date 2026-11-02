#include "force/models/thin_disk.h"

#include "utils/integrators/axisymmetric/ring_sum.h"

#include <math.h>

void rp_thin_disk_default_params(rp_thin_disk_params *params)
{
    if (params == NULL) {
        return;
    }
    params->alpha = -1.0;
    params->inner_ratio = 0.0;
    params->softening = 0.01;
    params->quadrature_points = 16U;
    params->radial_min = 1.0e-5;
    params->radial_match = 4.0;
    params->radial_far = 5.0;
    params->radial_step = 0.04;
    params->angular_step = 0.02;
    params->floor_step = 0.08;
    params->far_order = 24U;
}

/* The disk's mass: Sigma ~ R^alpha on [R_in / R_out, 1], softening H = eps R'.
 * A NULL density is the power law of ring_sum.h, Sigma ~ R^(beta - 2). */
static rp_thin_mass disk_mass(const rp_thin_disk_params *params)
{
    const rp_thin_mass mass = {
        .sigma = NULL,
        .context = NULL,
        .inner_ratio = params->inner_ratio,
        .softening = params->softening,
        .beta = params->alpha + 2.0,
        .inner_exponent = params->alpha + 2.0,
        .quadrature_points = params->quadrature_points,
    };
    return mass;
}

static rp_axisymmetric_grid disk_grid(const rp_thin_disk_params *params)
{
    const rp_axisymmetric_grid grid = {
        .radial_min = params->radial_min,
        .radial_match = params->radial_match,
        .radial_far = params->radial_far,
        .radial_step = params->radial_step,
        .angular_step = params->angular_step,
        .floor_step = params->floor_step,
        .far_order = params->far_order,
    };
    return grid;
}

rp_axisymmetric_status rp_thin_disk_direct_field(
    const rp_thin_disk_params *params, double d, double z, double field[3]
)
{
    if (params == NULL || field == NULL) {
        return RP_AXISYMMETRIC_STATUS_NULL_POINTER;
    }
    if (!(d >= 0.0) || !isfinite(d + z) || !isfinite(params->alpha)) {
        return RP_AXISYMMETRIC_STATUS_INVALID_ARGUMENT;
    }
    const rp_thin_mass mass = disk_mass(params);
    rp_thin_ring_sum sum;
    const rp_axisymmetric_status status = rp_thin_ring_sum_init(&mass, &sum);
    if (status != RP_AXISYMMETRIC_STATUS_OK) {
        return status;
    }
    return rp_thin_ring_sum_field(d, fabs(z), &sum, field) == 0
        ? RP_AXISYMMETRIC_STATUS_OK : RP_AXISYMMETRIC_STATUS_FIELD_FAILURE;
}

rp_axisymmetric_status rp_thin_disk_table_spec(
    const rp_thin_disk_params *params, rp_axisymmetric_spec *spec
)
{
    if (params == NULL || spec == NULL) {
        return RP_AXISYMMETRIC_STATUS_NULL_POINTER;
    }
    const rp_thin_mass mass = disk_mass(params);
    const rp_axisymmetric_grid grid = disk_grid(params);
    rp_thin_ring_sum sum;
    const rp_axisymmetric_status status = rp_thin_ring_sum_init(&mass, &sum);
    if (status != RP_AXISYMMETRIC_STATUS_OK) {
        return status;
    }
    return rp_thin_ring_sum_table_spec(&sum, &grid, spec);
}

rp_axisymmetric_status rp_thin_disk_table_build(
    const rp_thin_disk_params *params, rp_axisymmetric_table **table
)
{
    if (params == NULL || table == NULL) {
        if (table != NULL) {
            *table = NULL;
        }
        return RP_AXISYMMETRIC_STATUS_NULL_POINTER;
    }
    const rp_thin_mass mass = disk_mass(params);
    const rp_axisymmetric_grid grid = disk_grid(params);
    return rp_thin_ring_sum_table_build(&mass, &grid, table);
}

void rp_thin_disk_table_free(rp_axisymmetric_table *table)
{
    rp_axisymmetric_table_free(table);
}

rp_axisymmetric_status rp_thin_disk_rotation_from_normal(
    const double normal[3], double rotation[3][3]
)
{
    return rp_axisymmetric_rotation_from_normal(normal, rotation);
}

rp_axisymmetric_status rp_thin_disk_model_init(
    rp_thin_disk_model *model,
    const rp_axisymmetric_table *table,
    double mass,
    double outer_radius,
    const double normal[3]
)
{
    return rp_axisymmetric_model_init(model, table, mass, outer_radius, normal);
}

int rp_thin_disk_acceleration(
    const double position[3],
    const double velocity[3],
    const void *model,
    double acceleration[3]
)
{
    return rp_axisymmetric_acceleration(position, velocity, model, acceleration);
}

int rp_thin_disk_potential(
    const double position[3], const void *model, double *potential
)
{
    return rp_axisymmetric_potential(position, model, potential);
}
