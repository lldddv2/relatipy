/**
 * @file kerr_mcmc.c
 * @brief Evaluate an osculating Kerr orbit at observed arrival times.
 */

#include "kerr_mcmc.h"
#include "physic/kerr_observables.h"
#include "initial/convert.h"
#include "relatipy/kerr_geodesic.h"
#include "relatipy/kerr_geometry.h"

#include <math.h>
#include <stddef.h>
#include <stdint.h>

#define STATE_DIM 8U
#define OUTPUT_DIM 3U
#define GEOMETRIC_DIM 4U
#define ARRIVAL_ITERATIONS 12U

typedef struct {
    double spin;
} rhs_context;

static void zero_output(double *output, size_t sample_count)
{
    size_t index;
    for (index = 0U; index < sample_count * OUTPUT_DIM; ++index) {
        output[index] = 0.0;
    }
}

static int coordinate_time_rhs(
    double coordinate_time,
    const double state[],
    double derivative[],
    void *opaque_context
)
{
    const rhs_context *context = (const rhs_context *)opaque_context;
    double dx[4];
    double du[4];
    double inverse_ut;
    rp_kerr_status status;
    size_t index;

    (void)coordinate_time;
    if (context == NULL || !(state[4] > 0.0) || !isfinite(state[4])) {
        return -1;
    }
    /* The corrected Levi--Civita contraction applies to timelike as well as
       null tangents; this operation does not enforce a null constraint. */
    status = rp_kerr_null_geodesic_rhs(
        1.0, context->spin, state, state + 4U, dx, du
    );
    if (status != RP_KERR_STATUS_OK) {
        return -1;
    }
    inverse_ut = 1.0 / state[4];
    for (index = 0U; index < 4U; ++index) {
        derivative[index] = dx[index] * inverse_ut;
        derivative[index + 4U] = du[index] * inverse_ut;
        if (!isfinite(derivative[index])
            || !isfinite(derivative[index + 4U])) {
            return -1;
        }
    }
    return 0;
}

static int proper_time_rhs(
    double proper_time,
    const double state[],
    double derivative[],
    void *opaque_context
)
{
    const rhs_context *context = (const rhs_context *)opaque_context;
    rp_kerr_status status;
    size_t index;

    (void)proper_time;
    if (context == NULL || !(state[4] > 0.0) || !isfinite(state[4])) {
        return -1;
    }
    status = rp_kerr_null_geodesic_rhs(
        1.0, context->spin, state, state + 4U,
        derivative, derivative + 4U
    );
    if (status != RP_KERR_STATUS_OK) {
        return -1;
    }
    for (index = 0U; index < STATE_DIM; ++index) {
        if (!isfinite(derivative[index])) {
            return -1;
        }
    }
    return 0;
}

static void transpose_multiply(
    const double rotation[3][3],
    const double observer_vector[3],
    double body_vector[3]
)
{
    size_t index;
    size_t component;
    for (index = 0U; index < 3U; ++index) {
        body_vector[index] = 0.0;
        for (component = 0U; component < 3U; ++component) {
            body_vector[index] += rotation[component][index]
                * observer_vector[component];
        }
    }
}

static rp_kerr_mcmc_status initial_state(
    double spin,
    const double rotation[3][3],
    double semi_major_axis,
    double eccentricity,
    double inclination,
    double ascending_node,
    double periapsis_argument,
    double state[STATE_DIM]
)
{
    const double elements[RP_INITIAL_CARTESIAN_DIM] = {
        0.0, semi_major_axis, eccentricity, inclination,
        ascending_node, periapsis_argument, 0.0
    };
    double observer[RP_INITIAL_CARTESIAN_DIM];
    double body[RP_INITIAL_CARTESIAN_DIM];
    rp_kerr_status status;

    if (!(semi_major_axis > 0.0) || !(eccentricity >= 0.0)
        || !(eccentricity < 1.0) || !(spin >= 0.0) || !(spin <= 1.0)) {
        return RP_KERR_MCMC_INVALID_INPUT;
    }
    /* Keep orbital angles in the observer frame; Kerr.orbit uses the same
       osculating Cartesian state in the spin-aligned frame. */
    status = rp_initial_elements_to_cartesian(elements, observer);
    if (status != RP_KERR_STATUS_OK) {
        return RP_KERR_MCMC_INVALID_INPUT;
    }
    body[0] = observer[0];
    transpose_multiply(rotation, observer + 1U, body + 1U);
    transpose_multiply(rotation, observer + 4U, body + 4U);
    status = rp_initial_cartesian_to_canonical(spin, body, state);
    return status == RP_KERR_STATUS_OK
        ? RP_KERR_MCMC_OK : RP_KERR_MCMC_INITIAL_STATE_FAILURE;
}

static void add_statistics(
    rp_integrator_stats *total,
    const rp_integrator_stats *step
)
{
    total->accepted_steps += step->accepted_steps;
    total->rejected_steps += step->rejected_steps;
    total->rhs_evaluations += step->rhs_evaluations;
    total->jacobian_evaluations += step->jacobian_evaluations;
    total->linear_solves += step->linear_solves;
    total->final_independent_variable = step->final_independent_variable;
}

static rp_kerr_mcmc_status solve_arrival(
    const rp_integrator_config *solver,
    rhs_context *context,
    const double rotation[3][3],
    double target,
    double state[STATE_DIM],
    double result[GEOMETRIC_DIM],
    rp_integrator_stats *statistics
)
{
    double current_time = state[0];
    double next_time;
    double derivative;
    double residual;
    double tolerance;
    rp_integrator_stats step_statistics;
    rp_integrator_status solver_status;
    rp_kerr_status projection_status;
    size_t iteration;

    projection_status = rp_kerr_observables_project(
        state, context->spin, rotation, 1.0, result
    );
    if (projection_status != RP_KERR_STATUS_OK) {
        return RP_KERR_MCMC_INITIAL_STATE_FAILURE;
    }
    tolerance = 1e-9 * (1.0 + fabs(target));
    for (iteration = 0U; iteration < ARRIVAL_ITERATIONS; ++iteration) {
        residual = result[0] - target;
        if (fabs(residual) <= tolerance) {
            return RP_KERR_MCMC_OK;
        }
        derivative = (result[3] + 1.0) / state[4];
        if (!(derivative > 0.0) || !isfinite(derivative)) {
            return RP_KERR_MCMC_ARRIVAL_FAILURE;
        }
        next_time = current_time - residual / derivative;
        if (!isfinite(next_time)) {
            return RP_KERR_MCMC_ARRIVAL_FAILURE;
        }
        solver_status = rp_integrator_integrate(
            solver,
            coordinate_time_rhs,
            context,
            current_time,
            next_time,
            state,
            &step_statistics
        );
        add_statistics(statistics, &step_statistics);
        if (solver_status != RP_INTEGRATOR_STATUS_OK) {
            return RP_KERR_MCMC_INTEGRATION_FAILURE;
        }
        current_time = next_time;
        state[0] = current_time;
        projection_status = rp_kerr_observables_project(
            state, context->spin, rotation, 1.0, result
        );
        if (projection_status != RP_KERR_STATUS_OK) {
            return RP_KERR_MCMC_INTEGRATION_FAILURE;
        }
    }
    return RP_KERR_MCMC_ARRIVAL_FAILURE;
}

static int write_observation(
    const double geometric[GEOMETRIC_DIM],
    double arrival_time,
    const rp_kerr_mcmc_observation *observation,
    double output[OUTPUT_DIM]
)
{
    double elapsed = arrival_time - observation->reference_time;
    output[0] = geometric[1] * observation->angular_scale
        + observation->alpha_offset + observation->alpha_drift * elapsed;
    output[1] = geometric[2] * observation->angular_scale
        + observation->delta_offset + observation->delta_drift * elapsed;
    output[2] = geometric[3] - observation->velocity_offset;
    return isfinite(output[0]) && isfinite(output[1]) && isfinite(output[2]);
}

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
    rhs_context context;
    double initial[STATE_DIM];
    double state[STATE_DIM];
    double initial_observable[GEOMETRIC_DIM];
    double observable[GEOMETRIC_DIM];
    rp_kerr_mcmc_status status;
    rp_kerr_status projection_status;
    size_t split;
    size_t index;
    size_t component;

    if (statistics != NULL) {
        statistics->accepted_steps = 0U;
        statistics->rejected_steps = 0U;
        statistics->rhs_evaluations = 0U;
        statistics->jacobian_evaluations = 0U;
        statistics->linear_solves = 0U;
        statistics->final_independent_variable = 0.0;
    }
    if (output != NULL && sample_count <= SIZE_MAX / OUTPUT_DIM) {
        zero_output(output, sample_count);
    }
    if (solver == NULL || rotation == NULL || arrival_times == NULL
        || observation == NULL
        || output == NULL || statistics == NULL
        || sample_count > SIZE_MAX / OUTPUT_DIM
        || solver->dimension != STATE_DIM) {
        return RP_KERR_MCMC_INVALID_INPUT;
    }
    if (!isfinite(spin) || !isfinite(semi_major_axis)
        || !isfinite(eccentricity) || !isfinite(inclination)
        || !isfinite(ascending_node) || !isfinite(periapsis_argument)
        || !isfinite(observation->angular_scale)
        || !(observation->angular_scale > 0.0)
        || !isfinite(observation->alpha_offset)
        || !isfinite(observation->delta_offset)
        || !isfinite(observation->alpha_drift)
        || !isfinite(observation->delta_drift)
        || !isfinite(observation->velocity_offset)
        || !isfinite(observation->reference_time)) {
        return RP_KERR_MCMC_INVALID_INPUT;
    }
    for (index = 0U; index < sample_count; ++index) {
        if (!isfinite(arrival_times[index])
            || (index > 0U && arrival_times[index] < arrival_times[index - 1U])) {
            return RP_KERR_MCMC_INVALID_INPUT;
        }
    }

    status = initial_state(
        spin, rotation, semi_major_axis, eccentricity,
        inclination, ascending_node, periapsis_argument, initial
    );
    if (status != RP_KERR_MCMC_OK) {
        return status;
    }
    context.spin = spin;
    projection_status = rp_kerr_observables_project(
        initial, spin, rotation, 1.0, initial_observable
    );
    if (projection_status != RP_KERR_STATUS_OK) {
        return RP_KERR_MCMC_INITIAL_STATE_FAILURE;
    }

    split = 0U;
    while (split < sample_count
        && arrival_times[split] < initial_observable[0]) {
        ++split;
    }
    for (component = 0U; component < STATE_DIM; ++component) {
        state[component] = initial[component];
    }
    for (index = split; index > 0U; --index) {
        status = solve_arrival(
            solver, &context, rotation, arrival_times[index - 1U],
            state, observable, statistics
        );
        if (status != RP_KERR_MCMC_OK) {
            zero_output(output, sample_count);
            return status;
        }
        if (!write_observation(
            observable, arrival_times[index - 1U], observation,
            output + (index - 1U) * OUTPUT_DIM
        )) {
            zero_output(output, sample_count);
            return RP_KERR_MCMC_INTEGRATION_FAILURE;
        }
    }

    for (component = 0U; component < STATE_DIM; ++component) {
        state[component] = initial[component];
    }
    for (index = split; index < sample_count; ++index) {
        status = solve_arrival(
            solver, &context, rotation, arrival_times[index],
            state, observable, statistics
        );
        if (status != RP_KERR_MCMC_OK) {
            zero_output(output, sample_count);
            return status;
        }
        if (!write_observation(
            observable, arrival_times[index], observation,
            output + index * OUTPUT_DIM
        )) {
            zero_output(output, sample_count);
            return RP_KERR_MCMC_INTEGRATION_FAILURE;
        }
    }
    return RP_KERR_MCMC_OK;
}

static int valid_initial_state(double spin, const double initial[STATE_DIM])
{
    double metric[RP_KERR_DIM][RP_KERR_DIM];
    double norm = 0.0;
    double horizon;
    size_t mu;
    size_t nu;

    if (!(spin >= 0.0) || !(spin <= 1.0)) {
        return 0;
    }
    for (mu = 0U; mu < STATE_DIM; ++mu) {
        if (!isfinite(initial[mu])) {
            return 0;
        }
    }
    horizon = 1.0 + sqrt((1.0 - spin) * (1.0 + spin));
    if (!(initial[1] > horizon) || !(initial[4] > 0.0)
        || rp_kerr_metric(1.0, spin, initial, metric) != RP_KERR_STATUS_OK) {
        return 0;
    }
    for (mu = 0U; mu < RP_KERR_DIM; ++mu) {
        for (nu = 0U; nu < RP_KERR_DIM; ++nu) {
            norm += initial[4U + mu] * metric[mu][nu]
                * initial[4U + nu];
        }
    }
    return isfinite(norm) && fabs(norm + 1.0) <= 1e-8;
}

static rp_kerr_mcmc_status evaluate_state_endpoint(
    const rp_integrator_config *solver,
    rhs_context *context,
    const double rotation[3][3],
    rp_kerr_mcmc_time_kind time_kind,
    double target,
    double *current_time,
    double state[STATE_DIM],
    double angular_scale,
    double output[OUTPUT_DIM],
    rp_integrator_stats *statistics
)
{
    double observable[GEOMETRIC_DIM];
    rp_integrator_stats step_statistics;
    rp_integrator_status integration_status;
    rp_kerr_status projection_status;
    rp_kerr_mcmc_status status;

    if (time_kind == RP_KERR_MCMC_TIME_ARRIVAL) {
        status = solve_arrival(
            solver, context, rotation, target, state, observable, statistics
        );
        if (status != RP_KERR_MCMC_OK) {
            return status;
        }
        output[0] = observable[1] * angular_scale;
        output[1] = observable[2] * angular_scale;
        output[2] = observable[3];
    } else {
        if (target != *current_time) {
            integration_status = rp_integrator_integrate(
                solver,
                time_kind == RP_KERR_MCMC_TIME_PROPER
                    ? proper_time_rhs : coordinate_time_rhs,
                context, *current_time, target, state, &step_statistics
            );
            add_statistics(statistics, &step_statistics);
            if (integration_status != RP_INTEGRATOR_STATUS_OK) {
                return RP_KERR_MCMC_INTEGRATION_FAILURE;
            }
        }
        *current_time = target;
        if (time_kind == RP_KERR_MCMC_TIME_COORDINATE) {
            state[0] = target;
        }
        projection_status = rp_kerr_observables_project(
            state, context->spin, rotation, angular_scale, observable
        );
        if (projection_status != RP_KERR_STATUS_OK) {
            return RP_KERR_MCMC_INTEGRATION_FAILURE;
        }
        output[0] = observable[1];
        output[1] = observable[2];
        output[2] = observable[3];
    }
    return isfinite(output[0]) && isfinite(output[1]) && isfinite(output[2])
        ? RP_KERR_MCMC_OK : RP_KERR_MCMC_INTEGRATION_FAILURE;
}

rp_kerr_mcmc_status rp_kerr_mcmc_evaluate_state(
    const rp_integrator_config *solver,
    double spin,
    const double rotation[3][3],
    const double initial[STATE_DIM],
    double tau_initial,
    rp_kerr_mcmc_time_kind time_kind,
    const double *sorted_times,
    size_t sample_count,
    double angular_scale,
    double *output,
    rp_integrator_stats *statistics
)
{
    rhs_context context;
    double state[STATE_DIM];
    double initial_observable[GEOMETRIC_DIM];
    double initial_time;
    double current_time;
    size_t split;
    size_t index;
    size_t component;
    rp_kerr_mcmc_status status;

    if (statistics != NULL) {
        statistics->accepted_steps = 0U;
        statistics->rejected_steps = 0U;
        statistics->rhs_evaluations = 0U;
        statistics->jacobian_evaluations = 0U;
        statistics->linear_solves = 0U;
        statistics->final_independent_variable = 0.0;
    }
    if (output != NULL && sample_count <= SIZE_MAX / OUTPUT_DIM) {
        zero_output(output, sample_count);
    }
    if (solver == NULL || rotation == NULL || initial == NULL
        || sorted_times == NULL || output == NULL || statistics == NULL
        || sample_count > SIZE_MAX / OUTPUT_DIM
        || solver->dimension != STATE_DIM
        || (time_kind != RP_KERR_MCMC_TIME_COORDINATE
            && time_kind != RP_KERR_MCMC_TIME_PROPER
            && time_kind != RP_KERR_MCMC_TIME_ARRIVAL)
        || !isfinite(tau_initial) || !(angular_scale > 0.0)
        || !isfinite(angular_scale) || !valid_initial_state(spin, initial)) {
        return RP_KERR_MCMC_INVALID_INPUT;
    }
    for (index = 0U; index < sample_count; ++index) {
        if (!isfinite(sorted_times[index])
            || (index > 0U && sorted_times[index] < sorted_times[index - 1U])) {
            return RP_KERR_MCMC_INVALID_INPUT;
        }
    }
    if (rp_kerr_observables_project(
        initial, spin, rotation, angular_scale, initial_observable
    ) != RP_KERR_STATUS_OK) {
        return RP_KERR_MCMC_INVALID_INPUT;
    }
    context.spin = spin;
    initial_time = time_kind == RP_KERR_MCMC_TIME_COORDINATE
        ? initial[0]
        : time_kind == RP_KERR_MCMC_TIME_PROPER
            ? tau_initial : initial_observable[0];
    split = 0U;
    while (split < sample_count && sorted_times[split] < initial_time) {
        ++split;
    }
    for (component = 0U; component < STATE_DIM; ++component) {
        state[component] = initial[component];
    }
    current_time = initial_time;
    for (index = split; index > 0U; --index) {
        status = evaluate_state_endpoint(
            solver, &context, rotation, time_kind,
            sorted_times[index - 1U], &current_time, state, angular_scale,
            output + (index - 1U) * OUTPUT_DIM, statistics
        );
        if (status != RP_KERR_MCMC_OK) {
            zero_output(output, sample_count);
            return status;
        }
    }
    for (component = 0U; component < STATE_DIM; ++component) {
        state[component] = initial[component];
    }
    current_time = initial_time;
    for (index = split; index < sample_count; ++index) {
        status = evaluate_state_endpoint(
            solver, &context, rotation, time_kind,
            sorted_times[index], &current_time, state, angular_scale,
            output + index * OUTPUT_DIM, statistics
        );
        if (status != RP_KERR_MCMC_OK) {
            zero_output(output, sample_count);
            return status;
        }
    }
    return RP_KERR_MCMC_OK;
}
