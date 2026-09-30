/**
 * @file tensor.c
 * @brief Internal operations on contiguous tensor storage.
 */

#include "tensor.h"

#include <math.h>

void rp_tensor_zero(double *components, size_t component_count)
{
    size_t index;

    for (index = 0; index < component_count; ++index) {
        components[index] = 0.0;
    }
}

int rp_tensor_is_finite(
    const double *components,
    size_t component_count
)
{
    size_t index;

    for (index = 0; index < component_count; ++index) {
        if (!isfinite(components[index])) {
            return 0;
        }
    }
    return 1;
}
