/**
 * @file integrator.c
 * @brief Validation and registry dispatch for native geodesic integrators.
 */

#include "integrator.h"
#include "dop853.h"
#include "dp45.h"
#include "radau.h"

#include <math.h>
#include <stddef.h>

typedef rp_integrator_status (*rp_integrator_implementation)(
    const rp_integrator_config *,
    rp_integrator_rhs,
    void *,
    double,
    double,
    double[],
    rp_integrator_stats *
);

typedef struct {
    const char *name;
    rp_integrator_implementation integrate;
} rp_integrator_registry_entry;

static const rp_integrator_registry_entry RP_INTEGRATOR_REGISTRY[] = {
    {"dop853", rp_dop853_integrate},
    {"radau", rp_radau_integrate},
    {"dp45", rp_dp45_integrate},
    {"projection_radau", rp_radau_integrate}
};

static int state_is_finite(const double state[], size_t dimension)
{
    size_t index;

    for (index = 0U; index < dimension; ++index) {
        if (!isfinite(state[index])) {
            return 0;
        }
    }
    return 1;
}

static void clear_stats(rp_integrator_stats *stats, double initial_value)
{
    if (stats == NULL) {
        return;
    }
    stats->accepted_steps = 0U;
    stats->rejected_steps = 0U;
    stats->rhs_evaluations = 0U;
    stats->jacobian_evaluations = 0U;
    stats->linear_solves = 0U;
    stats->final_independent_variable = initial_value;
}

rp_integrator_status rp_integrator_report_step(
    const rp_integrator_config *config,
    double independent_variable,
    const double state[]
)
{
    int result;

    if (config->step_observer == NULL) {
        return RP_INTEGRATOR_STATUS_OK;
    }
    result = config->step_observer(
        independent_variable, state, config->step_observer_context
    );
    if (result > 0) {
        return RP_INTEGRATOR_STATUS_OBSERVER_STOPPED;
    }
    if (result < 0) {
        return RP_INTEGRATOR_STATUS_OBSERVER_FAILURE;
    }
    return RP_INTEGRATOR_STATUS_OK;
}

const char *rp_integrator_method_name(rp_integrator_method method)
{
    if ((int)method < 0 || method >= RP_INTEGRATOR_METHOD_COUNT) {
        return NULL;
    }
    return RP_INTEGRATOR_REGISTRY[(size_t)method].name;
}

rp_integrator_status rp_integrator_integrate(
    const rp_integrator_config *config,
    rp_integrator_rhs rhs,
    void *context,
    double initial_independent_variable,
    double final_independent_variable,
    double state[],
    rp_integrator_stats *stats
)
{
    clear_stats(stats, initial_independent_variable);

    if (config == NULL || rhs == NULL || state == NULL) {
        return RP_INTEGRATOR_STATUS_NULL_POINTER;
    }
    if ((int)config->method < 0
        || config->method >= RP_INTEGRATOR_METHOD_COUNT
        || config->dimension == 0U
        || config->dimension > RP_INTEGRATOR_MAX_DIMENSION
        || config->maximum_steps == 0U
        || (config->method == RP_INTEGRATOR_METHOD_PROJECTION_RADAU
            && config->projector == NULL)) {
        return RP_INTEGRATOR_STATUS_INVALID_ARGUMENT;
    }
    if (!isfinite(initial_independent_variable)
        || !isfinite(final_independent_variable)
        || !isfinite(config->relative_tolerance)
        || !isfinite(config->initial_step)
        || !isfinite(config->maximum_step)
        || !state_is_finite(state, config->dimension)) {
        return RP_INTEGRATOR_STATUS_NONFINITE_VALUE;
    }
    if (!(config->relative_tolerance > 0.0)
        || config->initial_step < 0.0
        || config->maximum_step < 0.0) {
        return RP_INTEGRATOR_STATUS_INVALID_ARGUMENT;
    }
    if (config->absolute_tolerances == NULL) {
        if (!isfinite(config->absolute_tolerance)) {
            return RP_INTEGRATOR_STATUS_NONFINITE_VALUE;
        }
        if (config->absolute_tolerance < 0.0) {
            return RP_INTEGRATOR_STATUS_INVALID_ARGUMENT;
        }
    } else {
        size_t component;

        for (component = 0U; component < config->dimension; ++component) {
            const double tolerance = config->absolute_tolerances[component];

            if (!isfinite(tolerance)) {
                return RP_INTEGRATOR_STATUS_NONFINITE_VALUE;
            }
            if (tolerance < 0.0) {
                return RP_INTEGRATOR_STATUS_INVALID_ARGUMENT;
            }
        }
    }
    if (initial_independent_variable == final_independent_variable) {
        return RP_INTEGRATOR_STATUS_OK;
    }

    return RP_INTEGRATOR_REGISTRY[(size_t)config->method].integrate(
        config,
        rhs,
        context,
        initial_independent_variable,
        final_independent_variable,
        state,
        stats
    );
}
