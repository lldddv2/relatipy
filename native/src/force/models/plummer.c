/**
 * @file plummer.c
 * @brief Plummer sphere potential and acceleration (see plummer.h).
 */

#include "force/models/plummer.h"

#include <math.h>
#include <stddef.h>

/* Validate (m, b) and return |x|^2 + b^2, or a negative value on failure. */
static double softened_square(const double position[3], const double *parameters)
{
    const double mass = parameters[0];
    const double scale = parameters[1];

    if (!isfinite(mass) || !isfinite(scale) || !(scale > 0.0)) {
        return -1.0;
    }
    return position[0] * position[0] + position[1] * position[1]
        + position[2] * position[2] + scale * scale;
}

int rp_force_plummer_acceleration(
    const double position[3],
    const double velocity[3],
    const void *model,
    double acceleration[3]
)
{
    const double *parameters = model;
    double d2;
    double factor;

    (void)velocity;
    acceleration[0] = 0.0;
    acceleration[1] = 0.0;
    acceleration[2] = 0.0;
    if (parameters == NULL) {
        return 1;
    }
    d2 = softened_square(position, parameters);
    if (!(d2 > 0.0) || !isfinite(d2)) {
        return 1;
    }
    factor = -parameters[0] / (d2 * sqrt(d2));
    if (!isfinite(factor)) {
        return 1;
    }
    acceleration[0] = factor * position[0];
    acceleration[1] = factor * position[1];
    acceleration[2] = factor * position[2];
    return 0;
}

int rp_force_plummer_potential(
    const double position[3],
    const void *model,
    double *potential
)
{
    const double *parameters = model;
    double d2;
    double value;

    *potential = 0.0;
    if (parameters == NULL) {
        return 1;
    }
    d2 = softened_square(position, parameters);
    if (!(d2 > 0.0) || !isfinite(d2)) {
        return 1;
    }
    value = -parameters[0] / sqrt(d2);
    if (!isfinite(value)) {
        return 1;
    }
    *potential = value;
    return 0;
}
