/**
 * @file numeric.h
 * @brief Internal scalar numerical predicates.
 */

#ifndef RELATIPY_NATIVE_NUMERIC_H
#define RELATIPY_NATIVE_NUMERIC_H

/**
 * Report whether a value is negligible relative to a nonnegative scale.
 *
 * The threshold is `64 * DBL_EPSILON * scale`, preserving the tolerance used
 * by the native geometry reference for singularity checks.
 *
 * @param value Scalar to test.
 * @param scale Nonnegative reference magnitude.
 * @return Nonzero when `abs(value)` does not exceed the scaled threshold.
 */
int rp_effectively_zero(double value, double scale);

/**
 * Report a Boyer--Lindquist polar angle on or within the axis guard.
 *
 * Callers validate finiteness separately. The guard uses angular distance
 * `64 * DBL_EPSILON` from either axis and rejects angles outside `(0, pi)`.
 */
int rp_bl_polar_axis_singular(double theta);

#endif /* RELATIPY_NATIVE_NUMERIC_H */
