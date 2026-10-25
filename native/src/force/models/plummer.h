/**
 * @file plummer.h
 * @brief Plummer sphere centred on the black hole as an external force.
 *
 * Potential and acceleration, in the black-hole Cartesian frame with the
 * black-hole mass set to one,
 *
 *     Phi = -m / sqrt(|x|^2 + b^2),
 *     a   = -m x / (|x|^2 + b^2)^(3/2),
 *
 * where `m` is the sphere mass in black-hole masses and `b` its scale radius
 * in `G M / c^2`.  The model data are `RP_FORCE_PLUMMER_PARAMETERS` doubles,
 * `model[0] = m` and `model[1] = b`, the `const double *` layout shared with
 * forces compiled at run time, so a caller can change both values between
 * integrations without rebuilding anything.  `m` may be any finite value
 * (zero gives an identically zero acceleration); `b` must be positive.
 *
 * Nothing here allocates, retains a pointer or uses mutable global state.
 * The header is private to the native layer and is not a stable ABI.
 */

#ifndef RELATIPY_NATIVE_FORCE_MODELS_PLUMMER_H
#define RELATIPY_NATIVE_FORCE_MODELS_PLUMMER_H

#include "force/external_force.h"

/** Number of doubles in the Plummer model data. */
#define RP_FORCE_PLUMMER_PARAMETERS 2

/**
 * `rp_external_acceleration_fn` of a Plummer sphere.
 *
 * `model` points to `(m, b)`.  Returns zero on success and one, writing a
 * zero acceleration, for a NULL model, non-finite parameters, `b <= 0` or a
 * non-finite result.  `velocity` is not read.
 */
int rp_force_plummer_acceleration(
    const double position[3],
    const double velocity[3],
    const void *model,
    double acceleration[3]
);

/**
 * Potential of the same model at one position.
 *
 * Returns zero and writes `*potential` on success; returns one and writes
 * zero under the conditions of `rp_force_plummer_acceleration`.
 */
int rp_force_plummer_potential(
    const double position[3],
    const void *model,
    double *potential
);

#endif /* RELATIPY_NATIVE_FORCE_MODELS_PLUMMER_H */
