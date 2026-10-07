/**
 * @file solve.h
 * @brief Private Kerr null-trajectory solver in coordinate time (section 6.13).
 *
 * Internal to RelatiPy; not a stable ABI.  Units: G = c = M = 1.
 */

#ifndef RELATIPY_GEODESIC_NULL_SOLVE_H
#define RELATIPY_GEODESIC_NULL_SOLVE_H

#include <stddef.h>

#include "geodesic/integrators/integrator.h"

/** Terminal reason of a null trajectory (null-specific; not the Orbit enum). */
typedef enum rp_kerr_null_termination {
    /** No terminal event: t_final reached or a failure occurred. */
    RP_KERR_NULL_TERMINATION_NONE = 0,
    /** Accepted step with r <= r_+ (1 + RP_KERR_NULL_HORIZON_MARGIN). */
    RP_KERR_NULL_TERMINATION_HORIZON = 1,
    /** Accepted escape state, or the t_final cut state with r >= r_escape. */
    RP_KERR_NULL_TERMINATION_ESCAPE = 2
} rp_kerr_null_termination;

/**
 * Heap-allocated result of rp_kerr_null_trajectory_solve.
 *
 * `states` holds `count` rows of eight doubles
 * `(t, r, theta, phi, k^t, k^r, k^theta, k^phi)`, row-major; `lambdas` holds
 * the matching internal affine parameters (lambda = 0 at the initial state).
 * `lambdas` is internal and must not become public API.  Both arrays are
 * allocated with realloc by the solver, have `capacity >= count` rows, are
 * owned by this structure and are released only by
 * rp_kerr_null_trajectory_free.
 *
 * `final_state`/`final_lambda` always hold the last valid state reached:
 * the state interpolated exactly at t_final, the last exterior accepted
 * state on a horizon termination, the accepted escape state (cut at t_final
 * if that step also reaches the time target), or the last
 * accepted state on a failure (the initial state if no step was accepted).
 * `termination` is the terminal reason; `integrator_status` copies the
 * returned status; `allocation_failed` is nonzero when storage could not
 * grow.
 */
typedef struct rp_kerr_null_trajectory {
    double *lambdas;
    double *states;
    size_t count;
    size_t capacity;
    double final_state[8];
    double final_lambda;
    rp_integrator_stats stats;
    rp_integrator_status integrator_status;
    rp_kerr_null_termination termination;
    int allocation_failed;
} rp_kerr_null_trajectory;

/**
 * Integrate one exterior Kerr null geodesic forward in coordinate time.
 *
 * `spin` is a/M in [0, 1].  `initial_state` is a null state (see
 * relatipy/kerr_null.h) with finite components, `k^t > 0`,
 * `r > r_+ (1 + RP_KERR_NULL_HORIZON_MARGIN)`, the polar guard satisfied and
 * relative norm (robust definition of rp_kerr_null_invariants) at most the
 * internal RP_KERR_NULL_SOLVE_NORM_SANITY = 1e-3. This rejects clearly
 * non-null states; it is not an accuracy criterion. Construction enforces
 * the tighter null constraint; continuation may retain numerical drift.
 * No projection or rescaling is applied to the supplied tangent.
 * Integration runs in lambda from 0 with the corrected null RHS; the
 * analytic Jacobian is supplied to Radau.  Only the future direction is
 * supported: `t_final > initial_state[0]`.
 * The internal time component starts at zero and evolves as t - t0;
 * t_final and t_eval are shifted by the same origin. The finite time span
 * must be representable as a double. Accepted relative times must increase
 * strictly; stagnation fails the observer rather than forming a zero-width
 * interpolation bracket. Output times restore t0 + relative time; sampled
 * rows and the cut endpoint copy t_eval[i] and t_final exactly instead.
 *
 * Termination, checked after each accepted step in this order:
 * 1. horizon: r <= r_+ (1 + eps_h); the interior state is not stored and the
 *    last exterior state is kept (no crossing localization);
 * 2. t_final: the first accepted step with t >= t_final is cut; t_final is
 *    located inside the step by quintic Hermite interpolation of positions
 *    in lambda (x, k and dk/dlambda at both ends); tangents use cubic Hermite
 *    (k and dk/dlambda at both ends). The time root uses the quintic, and
 *    the interpolated state, with t set exactly to t_final, is the endpoint.
 *    If escape is enabled and this cut state has r >= r_escape, termination
 *    is ESCAPE; otherwise reaching t_final returns OK, even when the accepted
 *    step end has r >= r_escape;
 * 3. escape before t_final: when enabled, r >= r_escape stores that accepted
 *    state without localizing the radial crossing.
 *
 * `r_escape <= 0` or `+INFINITY` disables escape; otherwise it must be
 * finite and greater than the initial r.  `method` is radau, dop853 or dp45;
 * projection_radau is rejected with INVALID_ARGUMENT.  `rtol`,
 * `scalar_atol` and the optional eight-element `vector_atol` (NULL selects
 * the scalar) follow the integrator contract, on the null state vector.
 *
 * Output rows: with `t_eval == NULL` (then `n_eval` must be 0), row 0 is the
 * initial state; with `store_steps` nonzero every accepted step is appended
 * (the step that crosses t_final is replaced by the interpolated t_final
 * state); with zero only the initial and the latest state are kept
 * (count <= 2).  With `t_eval != NULL` (`n_eval > 0`), `store_steps` is
 * ignored and only requested samples are stored: `t_eval` must be finite,
 * strictly increasing and inside `[t0, t_final]`; a sample equal to t0 is the
 * initial state, others use the same quintic-position/cubic-tangent Hermite
 * interpolation inside the accepted step
 * containing them.  Samples beyond a terminal event are not produced.
 *
 * Returns OK when t_final is reached (termination NONE);
 * OBSERVER_STOPPED with termination HORIZON or ESCAPE; NULL_POINTER;
 * INVALID_ARGUMENT or NONFINITE_VALUE for rejected inputs (nothing
 * allocated); OBSERVER_FAILURE when storage cannot grow (allocation_failed)
 * or interpolation fails; any other integrator failure as returned.
 *
 * `*trajectory` is zeroed on entry, so it must not own storage from an
 * earlier call.  After any return other than NULL_POINTER the caller must
 * call rp_kerr_null_trajectory_free.  Reentrant: all working storage lives
 * in the call; a trajectory must not be shared between concurrent calls.
 */
rp_integrator_status rp_kerr_null_trajectory_solve(
    double spin,
    const double initial_state[8],
    double t_final,
    const double *t_eval,
    size_t n_eval,
    double r_escape,
    rp_integrator_method method,
    double rtol,
    double scalar_atol,
    const double *vector_atol,
    int store_steps,
    rp_kerr_null_trajectory *trajectory
);

/** Free both arrays and zero the structure. Accepts NULL; idempotent. */
void rp_kerr_null_trajectory_free(rp_kerr_null_trajectory *trajectory);

#endif /* RELATIPY_GEODESIC_NULL_SOLVE_H */
