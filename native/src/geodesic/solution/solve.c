#include "solve.h"

#include "geodesic/integrators/kerr.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    rp_kerr_trajectory *trajectory;
    double horizon;
    int store_steps;
} rp_collection_context;

static int reserve(rp_kerr_trajectory *trajectory, size_t capacity)
{
    double *taus;
    double *states;
    if (capacity > SIZE_MAX / (8U * sizeof(double))) {
        return -1;
    }
    taus = realloc(trajectory->taus, capacity * sizeof(double));
    if (taus == NULL) {
        return -1;
    }
    trajectory->taus = taus;
    states = realloc(trajectory->states, capacity * 8U * sizeof(double));
    if (states == NULL) {
        return -1;
    }
    trajectory->states = states;
    trajectory->capacity = capacity;
    return 0;
}

static int collect_step(double tau, const double state[], void *opaque_context)
{
    rp_collection_context *context = opaque_context;
    rp_kerr_trajectory *trajectory = context->trajectory;
    size_t index;

    /* The event is observed only when an accepted step actually crosses. */
    if (state[1] <= context->horizon) {
        trajectory->crossed_outer_horizon = 1;
        return 1;
    }
    index = context->store_steps ? trajectory->count : 1U;
    if (index >= trajectory->capacity) {
        size_t capacity = trajectory->capacity < 2U ? 2U : trajectory->capacity * 2U;
        if (capacity <= trajectory->capacity || reserve(trajectory, capacity) != 0) {
            trajectory->allocation_failed = 1;
            return -1;
        }
    }
    trajectory->taus[index] = tau;
    memcpy(trajectory->states + 8U * index, state, 8U * sizeof(double));
    trajectory->count = index + 1U;
    return 0;
}

rp_integrator_status rp_kerr_trajectory_solve(
    double spin,
    const double initial_state[8],
    double tau_initial,
    double tau_final,
    rp_integrator_method method,
    double rtol,
    double scalar_atol,
    const double *vector_atol,
    double first_step,
    double max_step,
    int store_steps,
    rp_kerr_trajectory *trajectory
)
{
    rp_integrator_config config;
    rp_kerr_integrator_context kerr_context;
    rp_collection_context collection;
    rp_integrator_status result;
    double current[8];
    double horizon;
    size_t i;

    if (trajectory == NULL || initial_state == NULL) {
        return RP_INTEGRATOR_STATUS_NULL_POINTER;
    }
    memset(trajectory, 0, sizeof(*trajectory));
    if (!isfinite(spin) || spin < 0.0 || spin > 1.0
        || !isfinite(tau_initial) || !isfinite(tau_final)
        || tau_final < tau_initial) {
        return RP_INTEGRATOR_STATUS_INVALID_ARGUMENT;
    }
    horizon = 1.0 + sqrt((1.0 - spin) * (1.0 + spin));
    for (i = 0U; i < 8U; ++i) {
        if (!isfinite(initial_state[i])) {
            return RP_INTEGRATOR_STATUS_NONFINITE_VALUE;
        }
        current[i] = initial_state[i];
    }
    if (current[1] <= horizon) {
        return RP_INTEGRATOR_STATUS_INVALID_ARGUMENT;
    }
    memset(&kerr_context, 0, sizeof(kerr_context));
    kerr_context.mass = 1.0;
    kerr_context.spin = spin;
    if (method == RP_INTEGRATOR_METHOD_PROJECTION_RADAU
        && rp_kerr_integrator_prepare_projection(&kerr_context, current) != 0) {
        trajectory->integrator_status = RP_INTEGRATOR_STATUS_PROJECTION_FAILURE;
        return RP_INTEGRATOR_STATUS_PROJECTION_FAILURE;
    }
    if (reserve(trajectory, store_steps ? 32U : 2U) != 0) {
        trajectory->allocation_failed = 1;
        return RP_INTEGRATOR_STATUS_OBSERVER_FAILURE;
    }
    trajectory->taus[0] = tau_initial;
    memcpy(trajectory->states, current, sizeof(current));
    trajectory->count = 1U;
    config.method = method;
    config.dimension = 8U;
    config.relative_tolerance = rtol;
    config.absolute_tolerance = scalar_atol;
    config.initial_step = first_step;
    config.maximum_step = max_step;
    config.maximum_steps = 1000000U;
    config.absolute_tolerances = vector_atol;
    config.step_observer = collect_step;
    collection.trajectory = trajectory;
    collection.horizon = horizon;
    collection.store_steps = store_steps != 0;
    config.step_observer_context = &collection;
    config.jacobian = rp_kerr_integrator_jacobian;
    config.projector = method == RP_INTEGRATOR_METHOD_PROJECTION_RADAU
        ? rp_kerr_integrator_project : NULL;
    result = rp_integrator_integrate(
        &config, rp_kerr_integrator_rhs, &kerr_context,
        tau_initial, tau_final, current, &trajectory->stats
    );
    trajectory->integrator_status = result;
    return result;
}

void rp_kerr_trajectory_free(rp_kerr_trajectory *trajectory)
{
    if (trajectory != NULL) {
        free(trajectory->taus);
        free(trajectory->states);
        memset(trajectory, 0, sizeof(*trajectory));
    }
}
