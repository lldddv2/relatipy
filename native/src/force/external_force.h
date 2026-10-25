/**
 * @file external_force.h
 * @brief Native contract of an external Newtonian force on a Kerr geodesic.
 *
 * Positions, velocities and accelerations use the black-hole Cartesian frame
 * defined in `docs/development/external-forces.rst`, in geometric units with
 * the black-hole mass set to one.  Every force is one callback with
 * read-only model data; no Python callback runs inside an integrator.
 */

#ifndef RELATIPY_NATIVE_FORCE_EXTERNAL_FORCE_H
#define RELATIPY_NATIVE_FORCE_EXTERNAL_FORCE_H

/**
 * Evaluate a Newtonian acceleration.
 * @param position Cartesian position `(x, y, z)` in the black-hole frame.
 * @param velocity Coordinate-time velocity `dx/dt` in the same frame.
 * @param model Read-only model data owned by the caller; NULL for models
 *        without data.  It must not change during an integration.
 * @param acceleration Caller-owned output of three doubles.
 * @return Zero on success; any nonzero value is reported as an RHS failure.
 */
typedef int (*rp_external_acceleration_fn)(
    const double position[3],
    const double velocity[3],
    const void *model,
    double acceleration[3]
);

/** One external force: its callback and the model data passed to it. */
typedef struct {
    rp_external_acceleration_fn evaluate;
    const void *model;
} rp_external_force;

#endif /* RELATIPY_NATIVE_FORCE_EXTERNAL_FORCE_H */
