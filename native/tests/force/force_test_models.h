/**
 * @file force_test_models.h
 * @brief Analytic external force models used only by the native force tests.
 *
 * Each model is a static potential with its exact gradient, or a tilted and
 * offset combination with a velocity-dependent drag, all in the black-hole
 * Cartesian frame with mass one.  Callbacks follow `rp_external_acceleration_fn`
 * and read immutable model data.
 */

#ifndef RELATIPY_NATIVE_TESTS_FORCE_TEST_MODELS_H
#define RELATIPY_NATIVE_TESTS_FORCE_TEST_MODELS_H

#include "force/external_force.h"

#include <math.h>

/** Plummer sphere `Phi = -m / sqrt(|x - c|^2 + b^2)`. */
typedef struct {
    double mass;
    double scale;
    double center[3];
} test_plummer;

/** Miyamoto--Nagai disk with symmetry axis z, centered at the origin. */
typedef struct {
    double mass;
    double radial_scale;
    double vertical_scale;
} test_miyamoto_nagai;

/**
 * Tilted, offset, non-axisymmetric model: a Miyamoto--Nagai disk whose
 * frame is `rotation` (disk coordinates `q = R^T x`), an offset Plummer
 * sphere, and linear drag `-drag * dx/dt`.
 */
typedef struct {
    test_miyamoto_nagai disk;
    double rotation[3][3];
    test_plummer clump;
    double drag;
} test_tilted;

static inline double test_plummer_potential(
    const test_plummer *model,
    const double position[3]
)
{
    const double dx = position[0] - model->center[0];
    const double dy = position[1] - model->center[1];
    const double dz = position[2] - model->center[2];
    return -model->mass
        / sqrt(dx * dx + dy * dy + dz * dz + model->scale * model->scale);
}

static inline int test_plummer_acceleration(
    const double position[3],
    const double velocity[3],
    const void *opaque_model,
    double acceleration[3]
)
{
    const test_plummer *model = opaque_model;
    const double dx = position[0] - model->center[0];
    const double dy = position[1] - model->center[1];
    const double dz = position[2] - model->center[2];
    const double d2 = dx * dx + dy * dy + dz * dz
        + model->scale * model->scale;
    const double factor = -model->mass / (d2 * sqrt(d2));

    (void)velocity;
    acceleration[0] = factor * dx;
    acceleration[1] = factor * dy;
    acceleration[2] = factor * dz;
    return 0;
}

static inline double test_miyamoto_nagai_potential(
    const test_miyamoto_nagai *model,
    const double position[3]
)
{
    const double zeta = sqrt(position[2] * position[2]
        + model->vertical_scale * model->vertical_scale);
    const double s = model->radial_scale + zeta;
    return -model->mass / sqrt(position[0] * position[0]
        + position[1] * position[1] + s * s);
}

static inline int test_miyamoto_nagai_acceleration(
    const double position[3],
    const double velocity[3],
    const void *opaque_model,
    double acceleration[3]
)
{
    const test_miyamoto_nagai *model = opaque_model;
    const double zeta = sqrt(position[2] * position[2]
        + model->vertical_scale * model->vertical_scale);
    const double s = model->radial_scale + zeta;
    const double d2 = position[0] * position[0] + position[1] * position[1]
        + s * s;
    const double factor = -model->mass / (d2 * sqrt(d2));

    (void)velocity;
    acceleration[0] = factor * position[0];
    acceleration[1] = factor * position[1];
    acceleration[2] = factor * position[2] * s / zeta;
    return 0;
}

static inline int test_tilted_acceleration(
    const double position[3],
    const double velocity[3],
    const void *opaque_model,
    double acceleration[3]
)
{
    const test_tilted *model = opaque_model;
    double local[3];
    double local_acceleration[3];
    double clump[3];
    int i;
    int j;

    for (i = 0; i < 3; ++i) {
        local[i] = 0.0;
        for (j = 0; j < 3; ++j) {
            local[i] += model->rotation[j][i] * position[j];
        }
    }
    test_miyamoto_nagai_acceleration(
        local, velocity, &model->disk, local_acceleration
    );
    test_plummer_acceleration(position, velocity, &model->clump, clump);
    for (i = 0; i < 3; ++i) {
        acceleration[i] = clump[i] - model->drag * velocity[i];
        for (j = 0; j < 3; ++j) {
            acceleration[i] += model->rotation[i][j] * local_acceleration[j];
        }
    }
    return 0;
}

/** Callback that always fails, for status propagation tests. */
static inline int test_failing_acceleration(
    const double position[3],
    const double velocity[3],
    const void *model,
    double acceleration[3]
)
{
    (void)position;
    (void)velocity;
    (void)model;
    acceleration[0] = acceleration[1] = acceleration[2] = 0.0;
    return 7;
}

/** Rotation by `angle` about the unit vector `axis`. */
static inline void test_axis_rotation(
    const double axis[3],
    double angle,
    double rotation[3][3]
)
{
    const double c = cos(angle);
    const double s = sin(angle);
    const double t = 1.0 - c;
    const double x = axis[0];
    const double y = axis[1];
    const double z = axis[2];

    rotation[0][0] = t * x * x + c;
    rotation[0][1] = t * x * y - s * z;
    rotation[0][2] = t * x * z + s * y;
    rotation[1][0] = t * x * y + s * z;
    rotation[1][1] = t * y * y + c;
    rotation[1][2] = t * y * z - s * x;
    rotation[2][0] = t * x * z - s * y;
    rotation[2][1] = t * y * z + s * x;
    rotation[2][2] = t * z * z + c;
}

#endif /* RELATIPY_NATIVE_TESTS_FORCE_TEST_MODELS_H */
