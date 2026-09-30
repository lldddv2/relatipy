/**
 * @file integrator.h
 * @brief Internal, extensible endpoint-integration contract.
 *
 * This header is private to the native backend.  It deliberately does not
 * choose a public default method, a canonical geodesic state, event semantics,
 * or dense output.
 */

#ifndef RELATIPY_GEODESIC_INTEGRATORS_INTEGRATOR_H
#define RELATIPY_GEODESIC_INTEGRATORS_INTEGRATOR_H

#include <stddef.h>
#include <float.h>
#include <math.h>

#define RP_INTEGRATOR_MAX_DIMENSION 16U

typedef enum {
    RP_INTEGRATOR_STATUS_OK = 0,
    RP_INTEGRATOR_STATUS_NULL_POINTER,
    RP_INTEGRATOR_STATUS_INVALID_ARGUMENT,
    RP_INTEGRATOR_STATUS_NONFINITE_VALUE,
    RP_INTEGRATOR_STATUS_RHS_FAILURE,
    RP_INTEGRATOR_STATUS_STEP_UNDERFLOW,
    RP_INTEGRATOR_STATUS_MAXIMUM_STEPS,
    RP_INTEGRATOR_STATUS_CONVERGENCE_FAILURE,
    RP_INTEGRATOR_STATUS_OBSERVER_STOPPED,
    RP_INTEGRATOR_STATUS_OBSERVER_FAILURE,
    RP_INTEGRATOR_STATUS_PROJECTION_FAILURE
} rp_integrator_status;

typedef enum {
    RP_INTEGRATOR_METHOD_DOP853 = 0,
    RP_INTEGRATOR_METHOD_RADAU = 1,
    RP_INTEGRATOR_METHOD_DP45 = 2,
    RP_INTEGRATOR_METHOD_PROJECTION_RADAU = 3,
    RP_INTEGRATOR_METHOD_COUNT
} rp_integrator_method;

/**
 * Evaluate an ODE right-hand side entirely in native code.
 *
 * @param independent_variable Current independent variable.
 * @param state Read-only state with the configured dimension.
 * @param derivative Caller-owned output with the configured dimension.
 * @param context Opaque native context retained by the caller.
 * @return Zero on success; any nonzero value is reported as an RHS failure.
 */
typedef int (*rp_integrator_rhs)(
    double independent_variable,
    const double state[],
    double derivative[],
    void *context
);

/**
 * Fill an optional analytic Jacobian entirely in native code.
 * Output is caller-owned row-major storage; write only dimension columns
 * per row, separated by row_stride doubles. Context is the RHS context.
 * Return zero on success; nonzero is an RHS failure, without an FD retry.
 */
typedef int (*rp_integrator_jacobian)(
    double independent_variable,
    const double state[],
    double *jacobian,
    size_t dimension,
    size_t row_stride,
    void *context
);

/**
 * Inspect a newly accepted state without crossing the Python boundary.
 * The initial state is not reported: the caller already owns it.  Return 0
 * to continue, a positive value to stop successfully at this accepted state,
 * or a negative value to report an observer failure.  No pointer is retained.
 */
typedef int (*rp_integrator_step_observer)(
    double independent_variable,
    const double state[],
    void *context
);

/**
 * Project a candidate accepted state before it replaces the current state.
 * The callback receives the same caller-owned context as the RHS. Return zero
 * on success. A nonzero return leaves the last accepted state intact and is
 * reported as PROJECTION_FAILURE. Only projection_radau invokes this callback.
 */
typedef int (*rp_integrator_projector)(
    double independent_variable,
    double state[],
    void *context
);

typedef struct {
    rp_integrator_method method;
    size_t dimension;
    double relative_tolerance;
    double absolute_tolerance;
    double initial_step;
    double maximum_step;
    size_t maximum_steps;
    /** NULL selects absolute_tolerance for every component. */
    const double *absolute_tolerances;
    rp_integrator_step_observer step_observer;
    void *step_observer_context;
    /** NULL selects finite differences; used only by Radau. */
    rp_integrator_jacobian jacobian;
    /** Required by projection_radau; ignored by the other methods. */
    rp_integrator_projector projector;
} rp_integrator_config;

typedef struct {
    size_t accepted_steps;
    size_t rejected_steps;
    size_t rhs_evaluations;
    size_t jacobian_evaluations;
    size_t linear_solves;
    double final_independent_variable;
} rp_integrator_stats;

/**
 * Integrate a caller-owned state to one endpoint.
 *
 * The operation performs no allocation, retains no pointer, has no mutable
 * global state, and invokes only the supplied C callback.  The state is
 * changed only after each accepted step; on failure it contains the last
 * accepted state.  absolute_tolerances, when non-NULL, points to dimension
 * caller-owned nonnegative finite values and is read only for this call.
 * step_observer receives each newly accepted state, including the last one,
 * after state/stats are updated.  On a positive observer return, this function
 * returns OBSERVER_STOPPED with that accepted state as the endpoint; on a
 * negative return it returns OBSERVER_FAILURE, also with that last accepted
 * state. projection_radau invokes projector before committing and observing
 * each accepted candidate; projection failure keeps the prior accepted state.
 * The caller owns and retains all pointers.
 */
rp_integrator_status rp_integrator_integrate(
    const rp_integrator_config *config,
    rp_integrator_rhs rhs,
    void *context,
    double initial_independent_variable,
    double final_independent_variable,
    double state[],
    rp_integrator_stats *stats
);

/** Return the stable internal spelling of a valid method, or NULL. */
const char *rp_integrator_method_name(rp_integrator_method method);

/** Internal helpers used by the registered method implementations. */
static inline double rp_integrator_atol(
    const rp_integrator_config *config,
    size_t component
)
{
    return config->absolute_tolerances == NULL
        ? config->absolute_tolerance
        : config->absolute_tolerances[component];
}
static inline double rp_integrator_scale(
    const rp_integrator_config *config,
    size_t component,
    double magnitude
)
{
    const double scale = rp_integrator_atol(config, component)
        + config->relative_tolerance * magnitude;
    return fmax(scale, DBL_MIN);
}
rp_integrator_status rp_integrator_report_step(
    const rp_integrator_config *config,
    double independent_variable,
    const double state[]
);

#endif /* RELATIPY_GEODESIC_INTEGRATORS_INTEGRATOR_H */
