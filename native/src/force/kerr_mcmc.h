/**
 * @file kerr_mcmc.h
 * @brief Forced one-call Kerr orbit to geometric observables evaluator.
 *
 * Parallel to `geodesic/kerr_mcmc.h`, with the same status codes, time
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
 * `solver->method` must be `dop853`, `dp45` or `radau`; `projection_radau`
 * returns `RP_KERR_MCMC_INVALID_INPUT`.  `solver->jacobian` and
 * `solver->projector` are ignored (Radau uses finite differences) because
 * the RHS context differs from the homogeneous one.  With an inactive force
 * and a NULL solver Jacobian, the output and statistics are bitwise those of
 * `rp_kerr_mcmc_evaluate`.  A failing force callback is reported as
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
 * `rp_kerr_forced_mcmc_evaluate`.  With an inactive force and a NULL solver
 * Jacobian, output and statistics are bitwise those of
 * `rp_kerr_mcmc_evaluate_state`.
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
