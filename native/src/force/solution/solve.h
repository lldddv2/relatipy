/**
 * @file solve.h
 * @brief Forced Kerr timelike trajectory solve.
 *
 * Parallel to `geodesic/solution/solve.h`: same inputs, the same
 * `rp_kerr_trajectory` output and the same ownership, plus one external
 * force.  The Python layer can therefore build the same Solution objects.
 */

#ifndef RELATIPY_NATIVE_FORCE_SOLUTION_SOLVE_H
#define RELATIPY_NATIVE_FORCE_SOLUTION_SOLVE_H

#include "force/external_force.h"
#include "geodesic/solution/solve.h"

/*
 * Integrate one forced timelike Kerr trajectory from tau_initial to tau_final.
 *
 * Every argument except `force` has the meaning, units and validation of
 * rp_kerr_trajectory_solve: `spin` is a/M in [0, 1], `initial_state` is
 * (t, r, theta, phi, u^t, u^r, u^theta, u^phi) with r > r_+, the observer
 * stops at the first accepted step with r <= r_+ (not stored), and
 * `store_steps` selects all accepted steps or only the latest. `*trajectory`
 * is overwritten with zeros on entry and, after any return other than
 * NULL_POINTER, must be released with rp_kerr_trajectory_free.
 *
 * `force` is caller owned, read only and may be NULL; it and its model must
 * stay valid and unchanged for the call. A NULL force, or one whose
 * `evaluate` is NULL, gives a trajectory bitwise identical to
 * rp_kerr_trajectory_solve for dop853, dp45 and radau (radau keeps the
 * analytic Kerr Jacobian). With an active force radau uses the integrator's
 * finite-difference Jacobian. `projection_radau` is always rejected with
 * INVALID_ARGUMENT, because Carter's constant is not conserved. A failing
 * force callback is reported by the integrator as RHS_FAILURE.
 */
rp_integrator_status rp_kerr_forced_trajectory_solve(
    double spin,
    const rp_external_force *force,
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
);

#endif /* RELATIPY_NATIVE_FORCE_SOLUTION_SOLVE_H */
