/**
 * @file tensor.h
 * @brief Internal helpers for contiguous tensor storage.
 *
 * This header belongs to the native implementation and is not an installed or
 * stable C API.  The caller owns the storage; these helpers allocate nothing
 * and retain no pointers.
 */

#ifndef RELATIPY_NATIVE_TENSOR_H
#define RELATIPY_NATIVE_TENSOR_H

#include <stddef.h>

/**
 * Set every component in a contiguous tensor buffer to floating-point zero.
 *
 * @param components First component of caller-owned contiguous storage.
 * @param component_count Number of `double` components to clear.
 *
 * @pre `components` is non-null when `component_count` is greater than zero.
 * @note A zero component count performs no memory access.
 */
void rp_tensor_zero(double *components, size_t component_count);

/**
 * Report whether every component in a contiguous tensor buffer is finite.
 *
 * The tensor rank and shape are intentionally external to this helper so the
 * same operation applies to metrics, inverse metrics, connections, and future
 * metric-specific tensors.
 *
 * @param components First component of caller-owned contiguous storage.
 * @param component_count Number of `double` components to inspect.
 * @return Nonzero if every requested component is finite; zero otherwise.
 *
 * @pre `components` is non-null when `component_count` is greater than zero.
 * @note A zero component count returns nonzero without accessing memory.
 */
int rp_tensor_is_finite(
    const double *components,
    size_t component_count
);

#endif /* RELATIPY_NATIVE_TENSOR_H */
