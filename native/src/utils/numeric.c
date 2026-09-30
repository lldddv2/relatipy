/**
 * @file numeric.c
 * @brief Internal scalar numerical predicates.
 */

#include "numeric.h"

#include <float.h>
#include <math.h>

int rp_effectively_zero(double value, double scale)
{
    return fabs(value) <= 64.0 * DBL_EPSILON * scale;
}

int rp_bl_polar_axis_singular(double theta)
{
    const double tolerance = 64.0 * DBL_EPSILON;

    return theta <= tolerance || acos(-1.0) - theta <= tolerance;
}
