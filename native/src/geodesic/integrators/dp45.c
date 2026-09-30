/**
 * @file dp45.c
 * @brief Allocation-free Dormand--Prince embedded 5(4) endpoint integrator.
 *
 * Method: J. R. Dormand and P. J. Prince, "A family of embedded
 * Runge-Kutta formulae", J. Comput. Appl. Math. 6 (1980), pp. 19--26;
 * https://doi.org/10.1016/0771-050X(80)90013-3 .  Fractions were checked
 * against official SciPy v1.16.2, scipy/integrate/_ivp/rk.py, RK45.C/A/B/E:
 * https://github.com/scipy/scipy/blob/v1.16.2/scipy/integrate/_ivp/rk.py
 * SciPy's BSD-3-Clause notice ships in
 * native/vendor/scipy-dop853/LICENSE-SCIPY.
 * The implementation and controller below are original to RelatiPy.  It
 * advances with order 5 and estimates local error with the embedded order 4.
 * Dense output is outside the internal endpoint contract.
 */

#include "dp45.h"

#include <math.h>
#include <stdint.h>

#define RP_DP45_STAGES 6U
#define RP_DP45_SAFETY 0.9
#define RP_DP45_MINIMUM_FACTOR 0.2
#define RP_DP45_MAXIMUM_FACTOR 10.0

static const double RP_DP45_C[RP_DP45_STAGES] = {
    0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0
};

static const double RP_DP45_A[RP_DP45_STAGES][RP_DP45_STAGES] = {
    {0.0},
    {1.0 / 5.0},
    {3.0 / 40.0, 9.0 / 40.0},
    {44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0},
    {19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0,
     -212.0 / 729.0},
    {9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0,
     49.0 / 176.0, -5103.0 / 18656.0}
};

static const double RP_DP45_B[RP_DP45_STAGES] = {
    35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0,
    -2187.0 / 6784.0, 11.0 / 84.0
};

/* Fifth-order minus fourth-order weights; SciPy's RK45.E has opposite sign. */
static const double RP_DP45_E[RP_DP45_STAGES + 1U] = {
    71.0 / 57600.0, 0.0, -71.0 / 16695.0, 71.0 / 1920.0,
    -17253.0 / 339200.0, 22.0 / 525.0, -1.0 / 40.0
};

static int vector_is_finite(const double values[], size_t dimension)
{
    size_t component;

    for (component = 0U; component < dimension; ++component) {
        if (!isfinite(values[component])) {
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
    if (stats->rhs_evaluations == SIZE_MAX) {
        return RP_INTEGRATOR_STATUS_MAXIMUM_STEPS;
    }
    ++stats->rhs_evaluations;
    if (rhs(independent_variable, state, derivative, context) != 0) {
        return RP_INTEGRATOR_STATUS_RHS_FAILURE;
    }
    return vector_is_finite(derivative, dimension)
        ? RP_INTEGRATOR_STATUS_OK
        : RP_INTEGRATOR_STATUS_NONFINITE_VALUE;
}

static double scaled_root_mean_square(
    const double values[],
    const double scales[],
    size_t dimension
)
{
    double sum = 0.0;
    size_t component;

    for (component = 0U; component < dimension; ++component) {
        const double scaled_value = values[component] / scales[component];
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
    const double interval = fabs(final_value - initial_value);
    double scales[RP_INTEGRATOR_MAX_DIMENSION];
    double trial_state[RP_INTEGRATOR_MAX_DIMENSION];
    double trial_derivative[RP_INTEGRATOR_MAX_DIMENSION];
    double difference[RP_INTEGRATOR_MAX_DIMENSION];
    double state_norm;
    double derivative_norm;
    double curvature_norm;
    double trial_step;
    double asymptotic_step;
    size_t component;
    rp_integrator_status status;

    if (config->initial_step > 0.0) {
        *step = fmin(config->initial_step, interval);
        if (config->maximum_step > 0.0) {
            *step = fmin(*step, config->maximum_step);
        }
        return RP_INTEGRATOR_STATUS_OK;
    }
    for (component = 0U; component < config->dimension; ++component) {
        scales[component] = rp_integrator_scale(
            config, component, fabs(state[component])
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
    if (!(trial_step > 0.0) || !isfinite(trial_step)) {
        return RP_INTEGRATOR_STATUS_STEP_UNDERFLOW;
    }
    for (component = 0U; component < config->dimension; ++component) {
        trial_state[component] = state[component]
            + direction * trial_step * derivative[component];
    }
    if (!vector_is_finite(trial_state, config->dimension)) {
        return RP_INTEGRATOR_STATUS_NONFINITE_VALUE;
    }
    status = evaluate_rhs(
        rhs, context, initial_value + direction * trial_step,
        trial_state, trial_derivative, config->dimension, stats
    );
    if (status != RP_INTEGRATOR_STATUS_OK) {
        return status;
    }
    for (component = 0U; component < config->dimension; ++component) {
        difference[component] =
            (trial_derivative[component] - derivative[component]) / trial_step;
    }
    curvature_norm = scaled_root_mean_square(
        difference, scales, config->dimension
    );
    if (fmax(derivative_norm, curvature_norm) <= 1.0e-15) {
        asymptotic_step = fmax(1.0e-6, trial_step * 1.0e-3);
    } else {
        asymptotic_step = pow(
            0.01 / fmax(derivative_norm, curvature_norm), 1.0 / 5.0
        );
    }
    *step = fmin(100.0 * trial_step, asymptotic_step);
    *step = fmin(*step, interval);
    if (config->maximum_step > 0.0) {
        *step = fmin(*step, config->maximum_step);
    }
    return RP_INTEGRATOR_STATUS_OK;
}

static rp_integrator_status dp45_step(
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
    double stages[RP_DP45_STAGES + 1U][RP_INTEGRATOR_MAX_DIMENSION];
    double stage_state[RP_INTEGRATOR_MAX_DIMENSION];
    double error_squared = 0.0;
    size_t stage;
    size_t previous_stage;
    size_t component;
    rp_integrator_status status;

    for (component = 0U; component < config->dimension; ++component) {
        stages[0][component] = derivative[component];
    }
    for (stage = 1U; stage < RP_DP45_STAGES; ++stage) {
        for (component = 0U; component < config->dimension; ++component) {
            double increment = 0.0;

            for (previous_stage = 0U; previous_stage < stage;
                 ++previous_stage) {
                increment += RP_DP45_A[stage][previous_stage]
                    * stages[previous_stage][component];
            }
            stage_state[component] = state[component] + step * increment;
        }
        if (!vector_is_finite(stage_state, config->dimension)) {
            return RP_INTEGRATOR_STATUS_NONFINITE_VALUE;
        }
        status = evaluate_rhs(
            rhs, context, independent_variable + RP_DP45_C[stage] * step,
            stage_state, stages[stage], config->dimension, stats
        );
        if (status != RP_INTEGRATOR_STATUS_OK) {
            return status;
        }
    }
    for (component = 0U; component < config->dimension; ++component) {
        double increment = 0.0;

        for (stage = 0U; stage < RP_DP45_STAGES; ++stage) {
            increment += RP_DP45_B[stage] * stages[stage][component];
        }
        candidate[component] = state[component] + step * increment;
    }
    if (!vector_is_finite(candidate, config->dimension)) {
        return RP_INTEGRATOR_STATUS_NONFINITE_VALUE;
    }
    status = evaluate_rhs(
        rhs, context, independent_variable + step, candidate,
        candidate_derivative, config->dimension, stats
    );
    if (status != RP_INTEGRATOR_STATUS_OK) {
        return status;
    }
    for (component = 0U; component < config->dimension; ++component) {
        double error = 0.0;
        const double scale = rp_integrator_scale(
            config, component,
            fmax(fabs(state[component]), fabs(candidate[component]))
        );

        stages[RP_DP45_STAGES][component] = candidate_derivative[component];
        for (stage = 0U; stage <= RP_DP45_STAGES; ++stage) {
            error += RP_DP45_E[stage] * stages[stage][component];
        }
        error *= step / scale;
        error_squared += error * error;
    }
    *error_norm = sqrt(error_squared / (double)config->dimension);
    return isfinite(*error_norm)
        ? RP_INTEGRATOR_STATUS_OK
        : RP_INTEGRATOR_STATUS_NONFINITE_VALUE;
}

rp_integrator_status rp_dp45_integrate(
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
        rhs, context, independent_variable, state, derivative,
        config->dimension, stats
    );
    if (status != RP_INTEGRATOR_STATUS_OK) {
        return status;
    }
    status = select_initial_step(
        config, rhs, context, initial_independent_variable,
        final_independent_variable, state, derivative, direction,
        stats, &step_magnitude
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

        if (stats->rejected_steps >= config->maximum_steps
            || stats->accepted_steps >=
                config->maximum_steps - stats->rejected_steps) {
            return RP_INTEGRATOR_STATUS_MAXIMUM_STEPS;
        }
        step_magnitude = fmin(step_magnitude, remaining);
        if (config->maximum_step > 0.0) {
            step_magnitude = fmin(step_magnitude, config->maximum_step);
        }
        if (!(step_magnitude > 0.0) || !isfinite(step_magnitude)
            || (step_magnitude < minimum_step
                && !(step_magnitude == remaining
                    && independent_variable + direction * step_magnitude
                        == final_independent_variable))) {
            return RP_INTEGRATOR_STATUS_STEP_UNDERFLOW;
        }
        step = direction * step_magnitude;
        status = dp45_step(
            config, rhs, context, independent_variable, state, derivative,
            step, candidate, candidate_derivative, &error_norm, stats
        );
        if (status != RP_INTEGRATOR_STATUS_OK) {
            return status;
        }
        if (error_norm < 1.0) {
            factor = error_norm == 0.0
                ? RP_DP45_MAXIMUM_FACTOR
                : fmin(RP_DP45_MAXIMUM_FACTOR,
                    RP_DP45_SAFETY * pow(error_norm, -1.0 / 5.0));
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
            factor = fmax(RP_DP45_MINIMUM_FACTOR,
                RP_DP45_SAFETY * pow(error_norm, -1.0 / 5.0));
            step_magnitude *= factor;
            ++stats->rejected_steps;
            previous_step_rejected = 1;
        }
    }
    stats->final_independent_variable = final_independent_variable;
    return RP_INTEGRATOR_STATUS_OK;
}
