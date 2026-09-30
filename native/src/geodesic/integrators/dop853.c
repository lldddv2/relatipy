/**
 * @file dop853.c
 * @brief Allocation-free DOP853 endpoint integrator.
 *
 * The tableau is isolated under native/vendor/scipy-dop853.  This controller
 * follows the published 8(5,3) error estimator and does not implement dense
 * output, which remains outside the current internal endpoint contract.
 */

#include "dop853.h"
#include "../../../vendor/scipy-dop853/dop853_coefficients.h"

#include <float.h>
#include <math.h>
#include <stddef.h>

#define RP_DOP853_STAGES 12U
#define RP_DOP853_ERROR_STAGES 13U
#define RP_DOP853_SAFETY 0.9
#define RP_DOP853_MINIMUM_FACTOR 0.2
#define RP_DOP853_MAXIMUM_FACTOR 10.0

static int vector_is_finite(const double values[], size_t dimension)
{
    size_t index;

    for (index = 0U; index < dimension; ++index) {
        if (!isfinite(values[index])) {
            return 0;
        }
    }
    return 1;
}

static rp_integrator_status evaluate_rhs(
    rp_integrator_rhs rhs,
    void *context,
    double independent_variable,
    const double state[],
    double derivative[],
    size_t dimension,
    rp_integrator_stats *stats
)
{
    ++stats->rhs_evaluations;
    if (rhs(independent_variable, state, derivative, context) != 0) {
        return RP_INTEGRATOR_STATUS_RHS_FAILURE;
    }
    if (!vector_is_finite(derivative, dimension)) {
        return RP_INTEGRATOR_STATUS_NONFINITE_VALUE;
    }
    return RP_INTEGRATOR_STATUS_OK;
}

static double scaled_root_mean_square(
    const double values[],
    const double scales[],
    size_t dimension
)
{
    double sum = 0.0;
    size_t index;

    for (index = 0U; index < dimension; ++index) {
        const double scaled_value = values[index] / scales[index];
        sum += scaled_value * scaled_value;
    }
    return sqrt(sum / (double)dimension);
}

static rp_integrator_status select_initial_step(
    const rp_integrator_config *config,
    rp_integrator_rhs rhs,
    void *context,
    double initial_value,
    double final_value,
    const double state[],
    const double derivative[],
    double direction,
    rp_integrator_stats *stats,
    double *step
)
{
    double scales[RP_INTEGRATOR_MAX_DIMENSION];
    double trial_state[RP_INTEGRATOR_MAX_DIMENSION];
    double trial_derivative[RP_INTEGRATOR_MAX_DIMENSION];
    double difference[RP_INTEGRATOR_MAX_DIMENSION];
    const double interval = fabs(final_value - initial_value);
    double state_norm;
    double derivative_norm;
    double curvature_norm;
    double trial_step;
    double asymptotic_step;
    size_t index;
    rp_integrator_status status;

    if (config->initial_step > 0.0) {
        *step = fmin(config->initial_step, interval);
        if (config->maximum_step > 0.0) {
            *step = fmin(*step, config->maximum_step);
        }
        return RP_INTEGRATOR_STATUS_OK;
    }

    for (index = 0U; index < config->dimension; ++index) {
        scales[index] = rp_integrator_scale(
            config, index, fabs(state[index])
        );
    }
    state_norm = scaled_root_mean_square(state, scales, config->dimension);
    derivative_norm = scaled_root_mean_square(
        derivative, scales, config->dimension
    );
    trial_step = (state_norm < 1.0e-5 || derivative_norm < 1.0e-5)
        ? 1.0e-6
        : 0.01 * state_norm / derivative_norm;
    trial_step = fmin(trial_step, interval);
    if (config->maximum_step > 0.0) {
        trial_step = fmin(trial_step, config->maximum_step);
    }

    for (index = 0U; index < config->dimension; ++index) {
        trial_state[index] = state[index]
            + direction * trial_step * derivative[index];
    }
    status = evaluate_rhs(
        rhs,
        context,
        initial_value + direction * trial_step,
        trial_state,
        trial_derivative,
        config->dimension,
        stats
    );
    if (status != RP_INTEGRATOR_STATUS_OK) {
        return status;
    }
    for (index = 0U; index < config->dimension; ++index) {
        difference[index] = (trial_derivative[index] - derivative[index])
            / trial_step;
    }
    curvature_norm = scaled_root_mean_square(
        difference, scales, config->dimension
    );
    if (fmax(derivative_norm, curvature_norm) <= 1.0e-15) {
        asymptotic_step = fmax(1.0e-6, trial_step * 1.0e-3);
    } else {
        asymptotic_step = pow(
            0.01 / fmax(derivative_norm, curvature_norm), 1.0 / 8.0
        );
    }
    *step = fmin(100.0 * trial_step, asymptotic_step);
    *step = fmin(*step, interval);
    if (config->maximum_step > 0.0) {
        *step = fmin(*step, config->maximum_step);
    }
    return RP_INTEGRATOR_STATUS_OK;
}

static rp_integrator_status dop853_step(
    const rp_integrator_config *config,
    rp_integrator_rhs rhs,
    void *context,
    double independent_variable,
    const double state[],
    const double derivative[],
    double step,
    double candidate[],
    double candidate_derivative[],
    double *error_norm,
    rp_integrator_stats *stats
)
{
    double stages[RP_DOP853_ERROR_STAGES][RP_INTEGRATOR_MAX_DIMENSION];
    double stage_state[RP_INTEGRATOR_MAX_DIMENSION];
    double error_five_squared = 0.0;
    double error_three_squared = 0.0;
    size_t stage;
    size_t component;
    size_t previous_stage;
    rp_integrator_status status;

    for (component = 0U; component < config->dimension; ++component) {
        stages[0][component] = derivative[component];
    }
    for (stage = 1U; stage < RP_DOP853_STAGES; ++stage) {
        for (component = 0U; component < config->dimension; ++component) {
            double increment = 0.0;

            for (previous_stage = 0U; previous_stage < stage;
                 ++previous_stage) {
                increment += RP_DOP853_A[stage][previous_stage]
                    * stages[previous_stage][component];
            }
            stage_state[component] = state[component] + step * increment;
        }
        status = evaluate_rhs(
            rhs,
            context,
            independent_variable + RP_DOP853_C[stage] * step,
            stage_state,
            stages[stage],
            config->dimension,
            stats
        );
        if (status != RP_INTEGRATOR_STATUS_OK) {
            return status;
        }
    }
    for (component = 0U; component < config->dimension; ++component) {
        double increment = 0.0;

        for (stage = 0U; stage < RP_DOP853_STAGES; ++stage) {
            increment += RP_DOP853_B[stage] * stages[stage][component];
        }
        candidate[component] = state[component] + step * increment;
    }
    if (!vector_is_finite(candidate, config->dimension)) {
        return RP_INTEGRATOR_STATUS_NONFINITE_VALUE;
    }
    status = evaluate_rhs(
        rhs,
        context,
        independent_variable + step,
        candidate,
        candidate_derivative,
        config->dimension,
        stats
    );
    if (status != RP_INTEGRATOR_STATUS_OK) {
        return status;
    }
    for (component = 0U; component < config->dimension; ++component) {
        double error_five = 0.0;
        double error_three = 0.0;
        const double scale = rp_integrator_scale(
            config, component,
            fmax(fabs(state[component]), fabs(candidate[component]))
        );

        stages[12][component] = candidate_derivative[component];
        for (stage = 0U; stage < RP_DOP853_ERROR_STAGES; ++stage) {
            error_five += RP_DOP853_E5[stage] * stages[stage][component];
            error_three += RP_DOP853_E3[stage] * stages[stage][component];
        }
        error_five /= scale;
        error_three /= scale;
        error_five_squared += error_five * error_five;
        error_three_squared += error_three * error_three;
    }
    if (error_five_squared == 0.0 && error_three_squared == 0.0) {
        *error_norm = 0.0;
    } else {
        const double denominator = error_five_squared
            + 0.01 * error_three_squared;
        *error_norm = fabs(step) * error_five_squared
            / sqrt(denominator * (double)config->dimension);
    }
    return isfinite(*error_norm)
        ? RP_INTEGRATOR_STATUS_OK
        : RP_INTEGRATOR_STATUS_NONFINITE_VALUE;
}

rp_integrator_status rp_dop853_integrate(
    const rp_integrator_config *config,
    rp_integrator_rhs rhs,
    void *context,
    double initial_independent_variable,
    double final_independent_variable,
    double state[],
    rp_integrator_stats *stats
)
{
    rp_integrator_stats local_stats = {0};
    double derivative[RP_INTEGRATOR_MAX_DIMENSION];
    double candidate[RP_INTEGRATOR_MAX_DIMENSION];
    double candidate_derivative[RP_INTEGRATOR_MAX_DIMENSION];
    const double direction = final_independent_variable
        > initial_independent_variable ? 1.0 : -1.0;
    double independent_variable = initial_independent_variable;
    double step_magnitude;
    int previous_step_rejected = 0;
    size_t component;
    rp_integrator_status status;

    if (stats == NULL) {
        stats = &local_stats;
    }
    status = evaluate_rhs(
        rhs,
        context,
        independent_variable,
        state,
        derivative,
        config->dimension,
        stats
    );
    if (status != RP_INTEGRATOR_STATUS_OK) {
        return status;
    }
    status = select_initial_step(
        config,
        rhs,
        context,
        initial_independent_variable,
        final_independent_variable,
        state,
        derivative,
        direction,
        stats,
        &step_magnitude
    );
    if (status != RP_INTEGRATOR_STATUS_OK) {
        return status;
    }

    while (direction * (final_independent_variable - independent_variable) > 0.0) {
        const double remaining = fabs(
            final_independent_variable - independent_variable
        );
        const double minimum_step = 10.0 * fabs(
            nextafter(independent_variable, direction * INFINITY)
                - independent_variable
        );
        double step;
        double error_norm;
        double factor;

        if (stats->accepted_steps + stats->rejected_steps
            >= config->maximum_steps) {
            return RP_INTEGRATOR_STATUS_MAXIMUM_STEPS;
        }
        step_magnitude = fmin(step_magnitude, remaining);
        if (config->maximum_step > 0.0) {
            step_magnitude = fmin(step_magnitude, config->maximum_step);
        }
        if (!(step_magnitude >= minimum_step) || step_magnitude == 0.0) {
            return RP_INTEGRATOR_STATUS_STEP_UNDERFLOW;
        }
        step = direction * step_magnitude;
        status = dop853_step(
            config,
            rhs,
            context,
            independent_variable,
            state,
            derivative,
            step,
            candidate,
            candidate_derivative,
            &error_norm,
            stats
        );
        if (status != RP_INTEGRATOR_STATUS_OK) {
            return status;
        }
        if (error_norm < 1.0) {
            factor = error_norm == 0.0
                ? RP_DOP853_MAXIMUM_FACTOR
                : fmin(
                    RP_DOP853_MAXIMUM_FACTOR,
                    RP_DOP853_SAFETY * pow(error_norm, -1.0 / 8.0)
                );
            if (previous_step_rejected) {
                factor = fmin(1.0, factor);
            }
            independent_variable += step;
            for (component = 0U; component < config->dimension; ++component) {
                state[component] = candidate[component];
                derivative[component] = candidate_derivative[component];
            }
            ++stats->accepted_steps;
            stats->final_independent_variable = independent_variable;
            status = rp_integrator_report_step(
                config, independent_variable, state
            );
            if (status != RP_INTEGRATOR_STATUS_OK) {
                return status;
            }
            step_magnitude *= factor;
            previous_step_rejected = 0;
        } else {
            factor = fmax(
                RP_DOP853_MINIMUM_FACTOR,
                RP_DOP853_SAFETY * pow(error_norm, -1.0 / 8.0)
            );
            step_magnitude *= factor;
            ++stats->rejected_steps;
            previous_step_rejected = 1;
        }
    }
    stats->final_independent_variable = final_independent_variable;
    return RP_INTEGRATOR_STATUS_OK;
}
