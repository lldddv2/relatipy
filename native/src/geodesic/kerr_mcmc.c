/**
 * @file kerr_mcmc.c
 * @brief Evaluate an osculating Kerr orbit at observed arrival times.
 *
 * The homogeneous evaluator is the forced one of `force/kerr_mcmc.c` with
 * no force: an inactive force reproduces the homogeneous RHS bit for bit, so
 * the initial state, Rømer arrival solve and observation terms have a single
 * implementation.
 */

#include "kerr_mcmc.h"
#include "force/kerr_mcmc.h"

#include <stddef.h>

rp_kerr_mcmc_status rp_kerr_mcmc_evaluate(
    const rp_integrator_config *solver,
    double spin,
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
)
{
    return rp_kerr_forced_mcmc_evaluate(
        solver, spin, NULL, rotation, semi_major_axis, eccentricity,
        inclination, ascending_node, periapsis_argument, arrival_times,
        sample_count, observation, output, statistics
    );
}

rp_kerr_mcmc_status rp_kerr_mcmc_evaluate_state(
    const rp_integrator_config *solver,
    double spin,
    const double rotation[3][3],
    const double initial[8],
    double tau_initial,
    rp_kerr_mcmc_time_kind time_kind,
    const double *sorted_times,
    size_t sample_count,
    double angular_scale,
    double *output,
    rp_integrator_stats *statistics
)
{
    return rp_kerr_forced_mcmc_evaluate_state(
        solver, spin, NULL, rotation, initial, tau_initial, time_kind,
        sorted_times, sample_count, angular_scale, output, statistics
    );
}
