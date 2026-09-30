/**
 * @file kerr_mcmc.h
 * @brief Private one-call Kerr orbit to geometric observables evaluator.
 *
 * Orbital elements determine the same osculating Kepler state as Kerr.orbit
 * at true anomaly zero. Its observer-frame position and coordinate velocity
 * are rotated into the spin-aligned Kerr frame before evolution. Light
 * follows a straight-line Rømer model.
 */

#ifndef RELATIPY_NATIVE_GEODESIC_KERR_MCMC_H
#define RELATIPY_NATIVE_GEODESIC_KERR_MCMC_H

#include "integrators/integrator.h"

#include <stddef.h>

typedef enum {
    RP_KERR_MCMC_OK = 0,
    RP_KERR_MCMC_INVALID_INPUT,
    RP_KERR_MCMC_INITIAL_STATE_FAILURE,
    RP_KERR_MCMC_INTEGRATION_FAILURE,
    RP_KERR_MCMC_ARRIVAL_FAILURE
} rp_kerr_mcmc_status;

typedef enum {
    RP_KERR_MCMC_TIME_COORDINATE = 0,
    RP_KERR_MCMC_TIME_PROPER = 1,
    RP_KERR_MCMC_TIME_ARRIVAL = 2
} rp_kerr_mcmc_time_kind;

typedef struct {
    double angular_scale;   /* gravitational radius / source distance */
    double alpha_offset;    /* radians */
    double delta_offset;    /* radians */
    double alpha_drift;     /* radians per geometric time */
    double delta_drift;     /* radians per geometric time */
    double velocity_offset; /* v_LSR / c, subtracted from v_los/c */
    double reference_time;  /* dimensionless time from periapsis */
} rp_kerr_mcmc_observation;

/**
 * Compute alpha, delta and v_los at sorted arrival times.
 *
 * `rotation` is a proper body-to-observer rotation matrix. `spin` is the
 * nonnegative Kerr spin length in units of M=1. Orbital `semi_major_axis`,
 * `eccentricity`, `inclination`, `ascending_node`, and `periapsis_argument`
 * specify the osculating Cartesian position and coordinate velocity at
 * true anomaly zero in the observer frame. Angles are radians. The initial
 * four-velocity is normalized with the corrected Kerr metric.
 *
 * `arrival_times` are dimensionless `(t_obs - t_p)/(GM/c^3)`, sorted in
 * ascending order. `observation` supplies dimensionless scales, angular
 * offsets, geometric-time drifts and redshift offset. `output` has three
 * doubles per sample: `(alpha_rad, delta_rad, v_los_over_c)`.
 * The caller owns all buffers. No allocation or pointer retention occurs.
 * On any failure, all output elements are zeroed. The integrator
 * statistics aggregate all endpoint integrations, including arrival solves.
 */
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
);

/**
 * Evaluate a caller-supplied normalized Boyer--Lindquist state at endpoints.
 *
 * `initial` is `(t,r,theta,phi,u^t,u^r,u^theta,u^phi)` at `tau_initial`,
 * with `G=c=M=1`, `r>r_+`, future-directed timelike velocity and norm -1.
 * `sorted_times` are ascending coordinate, proper or arrival times according
 * to `time_kind`; arrival time is `t+Z` in the observer frame. Coordinate and
 * arrival times use the absolute `initial[0]` origin; proper times use the
 * supplied `tau_initial` origin. `angular_scale` is gravitational radius over
 * source distance. Each output row contains `(alpha_rad,delta_rad,v_los/c)`.
 * Input and output buffers must be disjoint. The caller owns all buffers;
 * no allocation or pointer retention occurs.
 * A non-null output is cleared on failure. Statistics aggregate integrations.
 */
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
);

#endif /* RELATIPY_NATIVE_GEODESIC_KERR_MCMC_H */
