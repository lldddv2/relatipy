/**
 * @file solve.c
 * @brief Reentrant affine integration with coordinate-time Hermite sampling.
 */

#include "solve.h"

#include "geodesic/integrators/kerr.h"
#include "relatipy/kerr_null.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define RP_NULL_MAXIMUM_STEPS 1000000U
/* Reject clearly non-null inputs, not the numerical drift of continuation. */
#define RP_KERR_NULL_SOLVE_NORM_SANITY 1.0e-3

typedef struct {
    rp_kerr_null_trajectory *trajectory;
    rp_kerr_integrator_context *kerr;
    const double *t_eval;
    size_t n_eval;
    size_t next_eval;
    double horizon;
    double t_origin;
    double absolute_t_final;
    double t_final;
    double r_escape;
    double previous[8];
    double previous_lambda;
    double previous_rhs[8];
    int previous_rhs_valid;
    int store_steps;
    int reached_final;
} rp_null_collection;

/* A partial realloc failure still leaves both pointers owned and freeable. */
static int reserve(rp_kerr_null_trajectory *trajectory, size_t capacity)
{
    double *values;
    if (capacity > SIZE_MAX / (8U * sizeof(double))) {
        return -1;
    }
    values = realloc(trajectory->lambdas, capacity * sizeof(double));
    if (values == NULL) {
        return -1;
    }
    trajectory->lambdas = values;
    values = realloc(trajectory->states, capacity * 8U * sizeof(double));
    if (values == NULL) {
        return -1;
    }
    trajectory->states = values;
    trajectory->capacity = capacity;
    return 0;
}

static int store(rp_null_collection *context, double lambda, const double state[8],
    double absolute_time)
{
    rp_kerr_null_trajectory *trajectory = context->trajectory;
    const size_t index = context->t_eval != NULL || context->store_steps
        || trajectory->count == 0U ? trajectory->count : 1U;
    if (index >= trajectory->capacity) {
        const size_t capacity = trajectory->capacity < 2U
            ? 2U : trajectory->capacity * 2U;
        if (capacity <= trajectory->capacity || reserve(trajectory, capacity) != 0) {
            trajectory->allocation_failed = 1;
            return -1;
        }
    }
    trajectory->lambdas[index] = lambda;
    memcpy(trajectory->states + 8U * index, state, 8U * sizeof(double));
    trajectory->states[8U * index] = absolute_time;
    trajectory->count = index + 1U;
    return 0;
}

static void set_final(rp_kerr_null_trajectory *trajectory,
    double lambda, const double state[8], double absolute_time)
{
    memcpy(trajectory->final_state, state, 8U * sizeof(double));
    trajectory->final_state[0] = absolute_time;
    trajectory->final_lambda = lambda;
}

/* Difference form avoids cancellation of large, nearly equal endpoints. */
static double hermite(double left, double right, double dleft, double dright,
    double width, double fraction)
{
    const double complement = 1.0 - fraction;
    return left + fraction * fraction * (3.0 - 2.0 * fraction) * (right - left)
        + width * fraction * complement
            * (complement * dleft - fraction * dright);
}

/* Quintic Hermite in difference form, with affine derivatives at both ends.
 * Factored bases preserve the endpoint zeros. Optional slope is d/ds, not
 * d/dlambda, for Newton's iteration in the dimensionless fraction s.
 */
static double hermite_quintic(double left, double right,
    double dleft, double dright, double ddleft, double ddright,
    double width, double fraction, double *slope)
{
    const double s = fraction;
    const double c = 1.0 - s;
    const double s2 = s * s;
    const double s3 = s2 * s;
    const double c2 = c * c;
    const double c3 = c2 * c;
    const double difference = right - left;
    const double vleft = width * dleft;
    const double vright = width * dright;
    const double aleft = width * (width * ddleft);
    const double aright = width * (width * ddright);

    if (slope != NULL) {
        *slope = 30.0 * s2 * c2 * difference
            + c2 * (1.0 + 2.0 * s - 15.0 * s2) * vleft
            + s2 * (-12.0 + 28.0 * s - 15.0 * s2) * vright
            + 0.5 * s * c2 * (2.0 - 5.0 * s) * aleft
            + 0.5 * s2 * c * (3.0 - 5.0 * s) * aright;
    }
    return left + s3 * (10.0 + s * (-15.0 + 6.0 * s)) * difference
        + s * c3 * (1.0 + 3.0 * s) * vleft
        - s3 * c * (4.0 - 3.0 * s) * vright
        + 0.5 * s2 * c3 * aleft + 0.5 * s3 * c2 * aright;
}

/* Locate the time root in the accepted affine bracket, without extrapolation.
 * Endpoint RHS values are computed lazily, once per sampled accepted step.
 * A cached previous right endpoint supplies the next left endpoint, if valid.
 * Their evaluations contribute to the reported RHS counter too.
 */
static int interpolate(rp_null_collection *context,
    double lambda, const double state[8], double target,
    double derivatives[2][8], int *derivatives_ready,
    double *sample_lambda, double sample[8])
{
    const double *previous = context->previous;
    const double width = lambda - context->previous_lambda;
    double lower = 0.0;
    double upper = 1.0;
    double fraction;
    double ulp;
    double tolerance;
    size_t iteration;
    size_t component;

    if (target == previous[0] || target == state[0]) {
        const int at_left = target == previous[0];
        memcpy(sample, at_left ? previous : state, 8U * sizeof(double));
        *sample_lambda = at_left ? context->previous_lambda : lambda;
        sample[0] = target;
        return 0;
    }
    if (!(width > 0.0) || !isfinite(width)
        || !(previous[0] < target && target < state[0])) {
        return -1;
    }
    if (!*derivatives_ready) {
        if (context->previous_rhs_valid) {
            memcpy(derivatives[0], context->previous_rhs, sizeof(derivatives[0]));
        } else {
            ++context->trajectory->stats.rhs_evaluations;
            if (rp_kerr_null_integrator_rhs(context->previous_lambda, previous,
                    derivatives[0], context->kerr) != 0) {
                return -1;
            }
        }
        ++context->trajectory->stats.rhs_evaluations;
        if (rp_kerr_null_integrator_rhs(lambda, state,
                derivatives[1], context->kerr) != 0) {
            return -1;
        }
        *derivatives_ready = 1;
    }
    ulp = fabs(nextafter(target, target >= 0.0 ? -INFINITY : INFINITY) - target);
    tolerance = 4.0 * ulp;
    fraction = (target - previous[0]) / (state[0] - previous[0]);
    if (!(fraction > lower && fraction < upper)) {
        fraction = 0.5;
    }
    /* Safeguarded Newton, with enough bisections even for subnormal roots. */
    for (iteration = 0U; iteration < 1100U; ++iteration) {
        double slope;
        const double residual = hermite_quintic(previous[0] - target,
            state[0] - target, previous[4], state[4],
            derivatives[0][4], derivatives[1][4], width, fraction, &slope);
        double candidate;
        if (!isfinite(residual)) {
            return -1;
        }
        if (fabs(residual) <= tolerance) {
            break;
        }
        if (residual < 0.0) {
            lower = fraction;
        } else {
            upper = fraction;
        }
        candidate = fraction - residual / slope;
        if (!isfinite(candidate) || !(candidate > lower && candidate < upper)) {
            candidate = lower + 0.5 * (upper - lower);
        }
        if (candidate == fraction || candidate == lower || candidate == upper) {
            /* The bracket has reached floating-point resolution. */
            break;
        }
        fraction = candidate;
    }
    if (iteration == 1100U) {
        return -1;
    }
    for (component = 0U; component < 8U; ++component) {
        sample[component] = component < 4U
            ? hermite_quintic(previous[component], state[component],
                previous[component + 4U], state[component + 4U],
                derivatives[0][component + 4U], derivatives[1][component + 4U],
                width, fraction, NULL)
            : hermite(previous[component], state[component],
                derivatives[0][component], derivatives[1][component], width, fraction);
        if (!isfinite(sample[component])) {
            return -1;
        }
    }
    sample[0] = target;
    *sample_lambda = context->previous_lambda + fraction * width;
    if (!isfinite(*sample_lambda) || !(sample[4] > 0.0)) {
        return -1;
    }
    return 0;
}

static int collect_step(double lambda, const double state[], void *opaque)
{
    rp_null_collection *context = opaque;
    rp_kerr_null_trajectory *trajectory = context->trajectory;
    double derivatives[2][8];
    int derivatives_ready = 0;
    const int reached_final = state[0] >= context->t_final;
    const double sample_limit = reached_final ? context->t_final : state[0];
    double endpoint[8];
    double endpoint_lambda = lambda;
    size_t component;

    /* This event takes precedence even over a time target in the same step. */
    if (state[1] <= context->horizon) {
        trajectory->termination = RP_KERR_NULL_TERMINATION_HORIZON;
        return 1;
    }
    for (component = 0U; component < 8U; ++component) {
        if (!isfinite(state[component])) {
            return -1;
        }
    }
    if (!isfinite(lambda) || !(lambda > context->previous_lambda)
        || !(state[4] > 0.0) || !(state[0] > context->previous[0])) {
        return -1;
    }
    /* Relative time must advance strictly to provide a nonzero time bracket.
     * An equal time is numerical stagnation, even with positive k^t. */
    /* Storage/interpolation failures still report this valid accepted state. */
    set_final(trajectory, lambda, state, context->t_origin + state[0]);
    if (reached_final) {
        if (interpolate(context, lambda, state, context->t_final,
                derivatives, &derivatives_ready, &endpoint_lambda, endpoint) != 0) {
            return -1;
        }
    } else {
        memcpy(endpoint, state, sizeof(endpoint));
    }
    while (context->next_eval < context->n_eval
        && context->t_eval[context->next_eval] - context->t_origin <= sample_limit) {
        double sample[8];
        double sample_lambda;
        const double absolute_time = context->t_eval[context->next_eval];
        if (interpolate(context, lambda, state, absolute_time - context->t_origin,
                derivatives, &derivatives_ready, &sample_lambda, sample) != 0
            || store(context, sample_lambda, sample, absolute_time) != 0) {
            return -1;
        }
        ++context->next_eval;
    }
    if (context->t_eval == NULL && store(context, endpoint_lambda, endpoint,
            reached_final ? context->absolute_t_final
                : context->t_origin + endpoint[0]) != 0) {
        return -1;
    }
    if (reached_final) {
        set_final(trajectory, endpoint_lambda, endpoint, context->absolute_t_final);
        /* A step beyond t_final can cross escape only after the requested
         * interval. Decide escape on the cut state, never on that step end. */
        if (endpoint[1] >= context->r_escape) {
            trajectory->termination = RP_KERR_NULL_TERMINATION_ESCAPE;
            return 1;
        }
        context->reached_final = 1;
        return 1;
    }
    if (state[1] >= context->r_escape) {
        trajectory->termination = RP_KERR_NULL_TERMINATION_ESCAPE;
        return 1;
    }
    memcpy(context->previous, state, sizeof(context->previous));
    context->previous_lambda = lambda;
    /* Never reuse a derivative belonging to an older accepted state. */
    context->previous_rhs_valid = derivatives_ready;
    if (derivatives_ready) {
        memcpy(context->previous_rhs, derivatives[1], sizeof(context->previous_rhs));
    }
    return 0;
}

static void accumulate(rp_integrator_stats *total, const rp_integrator_stats *part)
{
    total->accepted_steps += part->accepted_steps;
    total->rejected_steps += part->rejected_steps;
    total->rhs_evaluations += part->rhs_evaluations;
    total->jacobian_evaluations += part->jacobian_evaluations;
    total->linear_solves += part->linear_solves;
    total->final_independent_variable = part->final_independent_variable;
}

static rp_integrator_status validate(double spin, const double initial[8],
    double t_final, const double *t_eval, size_t n_eval, double r_escape,
    rp_integrator_method method, double rtol, double scalar_atol,
    const double *vector_atol)
{
    rp_kerr_null_invariants invariants;
    rp_kerr_status status;
    size_t index;
    if (!isfinite(spin) || !isfinite(t_final) || !isfinite(rtol) || isnan(r_escape)) {
        return RP_INTEGRATOR_STATUS_NONFINITE_VALUE;
    }
    for (index = 0U; index < 8U; ++index) {
        if (!isfinite(initial[index])) {
            return RP_INTEGRATOR_STATUS_NONFINITE_VALUE;
        }
        if (!isfinite(vector_atol == NULL ? scalar_atol : vector_atol[index])) {
            return RP_INTEGRATOR_STATUS_NONFINITE_VALUE;
        }
        if ((vector_atol == NULL ? scalar_atol : vector_atol[index]) < 0.0) {
            return RP_INTEGRATOR_STATUS_INVALID_ARGUMENT;
        }
    }
    if (spin < 0.0 || spin > 1.0 || !(t_final > initial[0]) || !(rtol > 0.0)
        || !(initial[4] > 0.0) || initial[1] <= rp_kerr_null_horizon_threshold(spin)
        || (method != RP_INTEGRATOR_METHOD_RADAU
            && method != RP_INTEGRATOR_METHOD_DOP853 && method != RP_INTEGRATOR_METHOD_DP45)
        || (r_escape > 0.0 && isfinite(r_escape) && r_escape <= initial[1])
        || (t_eval == NULL ? n_eval != 0U : n_eval == 0U)
        || n_eval > SIZE_MAX / sizeof(double)) {
        return RP_INTEGRATOR_STATUS_INVALID_ARGUMENT;
    }
    for (index = 0U; index < n_eval; ++index) {
        if (!isfinite(t_eval[index])) {
            return RP_INTEGRATOR_STATUS_NONFINITE_VALUE;
        }
        if (t_eval[index] < initial[0] || t_eval[index] > t_final
            || (index > 0U && t_eval[index] <= t_eval[index - 1U])) {
            return RP_INTEGRATOR_STATUS_INVALID_ARGUMENT;
        }
    }
    if (!isfinite(t_final - initial[0])) {
        return RP_INTEGRATOR_STATUS_NONFINITE_VALUE;
    }
    status = rp_kerr_null_invariants_evaluate(spin, initial, &invariants);
    if (status == RP_KERR_STATUS_NUMERICAL_RANGE) {
        return RP_INTEGRATOR_STATUS_NONFINITE_VALUE;
    }
    if (status != RP_KERR_STATUS_OK
        || invariants.relative_norm > RP_KERR_NULL_SOLVE_NORM_SANITY) {
        return RP_INTEGRATOR_STATUS_INVALID_ARGUMENT;
    }
    return RP_INTEGRATOR_STATUS_OK;
}

rp_integrator_status rp_kerr_null_trajectory_solve(
    double spin, const double initial_state[8], double t_final,
    const double *t_eval, size_t n_eval, double r_escape,
    rp_integrator_method method, double rtol, double scalar_atol,
    const double *vector_atol, int store_steps, rp_kerr_null_trajectory *trajectory)
{
    rp_kerr_integrator_context kerr = {0};
    rp_null_collection collection = {0};
    rp_integrator_config config = {0};
    rp_integrator_status result;
    double current[8];
    double lambda = 0.0;

    if (trajectory == NULL) {
        return RP_INTEGRATOR_STATUS_NULL_POINTER;
    }
    memset(trajectory, 0, sizeof(*trajectory));
    if (initial_state == NULL) {
        trajectory->integrator_status = RP_INTEGRATOR_STATUS_NULL_POINTER;
        return RP_INTEGRATOR_STATUS_NULL_POINTER;
    }
    result = validate(spin, initial_state, t_final, t_eval, n_eval, r_escape,
        method, rtol, scalar_atol, vector_atol);
    if (result != RP_INTEGRATOR_STATUS_OK) {
        trajectory->integrator_status = result;
        return result;
    }
    memcpy(current, initial_state, sizeof(current));
    /* Kerr is stationary: subtracting the time origin leaves the RHS and
     * Jacobian unchanged, while preserving small increments at large t0. */
    current[0] = 0.0;
    set_final(trajectory, 0.0, current, initial_state[0]);
    kerr.mass = 1.0;
    kerr.spin = spin;
    collection.trajectory = trajectory;
    collection.kerr = &kerr;
    collection.t_eval = t_eval;
    collection.n_eval = n_eval;
    collection.horizon = rp_kerr_null_horizon_threshold(spin);
    collection.t_origin = initial_state[0];
    collection.absolute_t_final = t_final;
    collection.t_final = t_final - initial_state[0];
    collection.r_escape = r_escape <= 0.0 ? INFINITY : r_escape;
    collection.store_steps = store_steps != 0;
    memcpy(collection.previous, current, sizeof(current));
    config.method = method;
    config.dimension = 8U;
    config.relative_tolerance = rtol;
    config.absolute_tolerance = scalar_atol;
    config.absolute_tolerances = vector_atol;
    config.step_observer = collect_step;
    config.step_observer_context = &collection;
    config.jacobian = rp_kerr_integrator_jacobian;
    /* initial_step and maximum_step remain zero (automatic); no projector. */
    if (t_eval == NULL || t_eval[0] == initial_state[0]) {
        if (store(&collection, 0.0, current, initial_state[0]) != 0) {
            trajectory->integrator_status = RP_INTEGRATOR_STATUS_OBSERVER_FAILURE;
            return RP_INTEGRATOR_STATUS_OBSERVER_FAILURE;
        }
        collection.next_eval = t_eval == NULL ? 0U : 1U;
    }
    for (;;) {
        rp_integrator_stats part;
        const size_t used = trajectory->stats.accepted_steps + trajectory->stats.rejected_steps;
        double endpoint;
        long double span;
        if (used >= RP_NULL_MAXIMUM_STEPS) {
            result = RP_INTEGRATOR_STATUS_MAXIMUM_STEPS;
            break;
        }
        /* No global upper bound on lambda(t) exists. Restart finite segments,
         * estimating their extent from the current positive k^t. Clipping in
         * long double avoids overflow for large but finite input times/scales.
         * Every segment shares the accepted-step observer and total step cap;
         * an integrator OK alone never means the time target was reached.
         */
        span = fmaxl(((long double)collection.t_final - current[0]) / current[4], 1.0L);
        span = fminl(span, ((long double)DBL_MAX - lambda) / 2.0L);
        endpoint = (double)((long double)lambda + span);
        if (!isfinite(endpoint) || !(endpoint > lambda)) {
            result = RP_INTEGRATOR_STATUS_STEP_UNDERFLOW;
            break;
        }
        config.maximum_steps = RP_NULL_MAXIMUM_STEPS - used;
        result = rp_integrator_integrate(&config, rp_kerr_null_integrator_rhs, &kerr,
            lambda, endpoint, current, &part);
        accumulate(&trajectory->stats, &part);
        if (result == RP_INTEGRATOR_STATUS_OBSERVER_STOPPED && collection.reached_final) {
            result = RP_INTEGRATOR_STATUS_OK;
            break;
        }
        if (result != RP_INTEGRATOR_STATUS_OK) {
            break;
        }
        lambda = endpoint;
    }
    trajectory->integrator_status = result;
    return result;
}

void rp_kerr_null_trajectory_free(rp_kerr_null_trajectory *trajectory)
{
    if (trajectory != NULL) {
        free(trajectory->lambdas);
        free(trajectory->states);
        memset(trajectory, 0, sizeof(*trajectory));
    }
}
