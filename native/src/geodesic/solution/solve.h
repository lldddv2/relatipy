#ifndef RELATIPY_GEODESIC_SOLUTION_SOLVE_H
#define RELATIPY_GEODESIC_SOLUTION_SOLVE_H

#include <stddef.h>
#include "geodesic/integrators/integrator.h"

/*
 * Heap-allocated Kerr trajectory produced by rp_kerr_trajectory_solve.
 *
 * Units are normalized geometric units with M = 1 (proper time in GM/c^3).
 * `taus` holds `count` proper times and `states` holds `count` rows of eight
 * doubles in Boyer--Lindquist order (t, r, theta, phi, u^t, u^r, u^theta,
 * u^phi), row-major. `capacity` is the allocated number of rows (at least
 * `count`). Both arrays are allocated with realloc by the solver and are
 * owned by the structure; release them only with rp_kerr_trajectory_free.
 * `stats` and `integrator_status` copy the integrator result;
 * `crossed_outer_horizon` is nonzero when an accepted step reached
 * r <= r_+; `allocation_failed` is nonzero when storage could not grow.
 */
typedef struct {
    double *taus;
    double *states;
    size_t count;
    size_t capacity;
    rp_integrator_stats stats;
    rp_integrator_status integrator_status;
    int crossed_outer_horizon;
    int allocation_failed;
} rp_kerr_trajectory;

/*
 * Integrate one timelike Kerr geodesic from tau_initial to tau_final.
 *
 * `spin` is a/M in [0, 1]; `initial_state` is the eight-component
 * Boyer--Lindquist state above, with r outside the outer horizon. `rtol`,
 * `scalar_atol` and the optional eight-element `vector_atol` (NULL for the
 * scalar value) are forwarded to the integrator, as are `first_step` and
 * `max_step`. With `store_steps` nonzero, row 0 is the initial state and
 * every accepted step is appended; with zero, only the initial state and the
 * latest accepted state are kept (count <= 2). The observer stops at the
 * first accepted step with r <= r_+ and that step is not stored, so the last
 * row is the last valid exterior state.
 *
 * `*trajectory` is overwritten with zeros on entry. It must therefore not
 * own storage from an earlier call: call rp_kerr_trajectory_free first or
 * the earlier arrays leak. After any return other than NULL_POINTER the
 * caller must call rp_kerr_trajectory_free, whatever the status, because
 * storage may have been allocated before a failure. Returns the integrator
 * status: INVALID_ARGUMENT for a spin outside [0, 1], non-finite or
 * decreasing times, or an initial radius at or inside r_+;
 * NONFINITE_VALUE for a non-finite initial state; PROJECTION_FAILURE when
 * projection_radau cannot be prepared; OBSERVER_FAILURE when storage cannot
 * be allocated. A `trajectory` must not be shared between concurrent calls.
 */
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
);

/* Free both arrays and zero the structure. Accepts NULL; safe to repeat. */
void rp_kerr_trajectory_free(rp_kerr_trajectory *trajectory);

#endif
