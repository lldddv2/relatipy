/**
 * @file kerr_mcmc.h
 * @brief Forced one-call Kerr orbit to geometric observables evaluator.
 *
 * The single implementation of the observable evaluators: the homogeneous
 * `geodesic/kerr_mcmc.h` calls it with no force.  Same status codes, time
 * kinds and observation terms: orbital elements give the osculating state
 * at true anomaly zero, rotated from the observer frame into the
 * spin-aligned Kerr frame; light follows the same straight-line Rømer
 * model; redshift, offsets, drifts and statistics are unchanged.  The orbit
 * additionally feels one external force through the forced RHS of
 * `force/kerr_forced.h`, in coordinate time (arrival and coordinate-time
 * samples) or proper time (proper-time samples).
 */

#ifndef RELATIPY_NATIVE_FORCE_KERR_MCMC_H
#define RELATIPY_NATIVE_FORCE_KERR_MCMC_H

#include "force/external_force.h"
#include "geodesic/kerr_mcmc.h"

#include <stddef.h>

/**
 * Compute alpha, delta and v_los at sorted arrival times under a force.
 *
 * Every argument except `force` has the meaning and validation of
 * `rp_kerr_mcmc_evaluate`.  `force` is caller owned, read only and may be
 * NULL; it and its model must stay valid and unchanged for the call.
 * With an active force, `solver->method` must be accepted by
 * `rp_kerr_forced_check_method` (`dop853`, `dp45` or `radau`); the
 * projection_* methods, which restore E, Lz and Q, are accepted only with an
 * inactive force and otherwise return `RP_KERR_MCMC_INVALID_INPUT`.  Their
 * invariants are those of the initial state, in either clock.
 * `solver->jacobian` and `solver->projector` are ignored (Radau uses finite
 * differences; the projector is `rp_kerr_integrator_project`) because the
 * RHS context differs from the homogeneous one.  A failing force callback is
 * reported as
 * `RP_KERR_MCMC_INTEGRATION_FAILURE`.  The caller owns all buffers; no
 * allocation or pointer retention occurs; any failure zeroes the output.
 */
rp_kerr_mcmc_status rp_kerr_forced_mcmc_evaluate(
    const rp_integrator_config *solver,
    double spin,
    const rp_external_force *force,
    const double rotation[3][3],
    double semi_major_axis,
    double eccentricity,
    double inclination,
    double ascending_node,
    double periapsis_argument,
    const double *arrival_times,
    size_t sample_count,
    const rp_kerr_mcmc_observation *observation,
    double *output,
    rp_integrator_stats *statistics
);

/**
 * Evaluate a caller-supplied normalized state at endpoints under a force.
 *
 * Every argument except `force` has the meaning and validation of
 * `rp_kerr_mcmc_evaluate_state`; `force` and the solver rules are those of
 * `rp_kerr_forced_mcmc_evaluate`.
 */
rp_kerr_mcmc_status rp_kerr_forced_mcmc_evaluate_state(
    const rp_integrator_config *solver,
    double spin,
    const rp_external_force *force,
    const double rotation[3][3],
    const double initial[8],
    double tau_initial,
    rp_kerr_mcmc_time_kind time_kind,
    const double *sorted_times,
    size_t sample_count,
    double angular_scale,
    double *output,
    rp_integrator_stats *statistics
);

#endif /* RELATIPY_NATIVE_FORCE_KERR_MCMC_H */
