/**
 * @file radau.c
 * @brief Allocation-free three-stage Radau IIA endpoint integrator.
 *
 * The method is the fifth-order, L-stable Radau IIA collocation formula.  A
 * caller-supplied analytic Jacobian, or finite differences when omitted, and
 * simplified Newton iteration solve each stage system. Adaptive error control
 * uses step doubling; dense output and events are intentionally outside this
 * private endpoint prototype.
 */

#include "radau.h"

#include <float.h>
#include <math.h>
#include <stddef.h>

#define RP_RADAU_STAGES 3U
#define RP_RADAU_ORDER 5.0
#define RP_RADAU_SYSTEM_DIMENSION (RP_RADAU_STAGES * RP_INTEGRATOR_MAX_DIMENSION)
#define RP_RADAU_MAXIMUM_NEWTON_ITERATIONS 8U
#define RP_RADAU_SAFETY 0.9
#define RP_RADAU_MINIMUM_FACTOR 0.2
#define RP_RADAU_MAXIMUM_FACTOR 5.0

static const double RP_RADAU_C[RP_RADAU_STAGES] = {
    0.15505102572168219,
    0.64494897427831781,
    1.0
};

static const double RP_RADAU_A[RP_RADAU_STAGES][RP_RADAU_STAGES] = {
    {0.19681547722366044, -0.06553542585019838, 0.02377097434822015},
    {0.39442431473908729, 0.29207341166522843, -0.04154875212599793},
    {0.37640306270046728, 0.51248582618842161, 0.11111111111111111}
};

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

static int lu_decompose(
    double matrix[RP_RADAU_SYSTEM_DIMENSION][RP_RADAU_SYSTEM_DIMENSION],
    size_t pivots[RP_RADAU_SYSTEM_DIMENSION],
    size_t dimension
)
{
    size_t column;

    for (column = 0U; column < dimension; ++column) {
        size_t pivot = column;
        double pivot_magnitude = fabs(matrix[column][column]);
        size_t row;
        size_t next_column;

        for (row = column + 1U; row < dimension; ++row) {
            const double magnitude = fabs(matrix[row][column]);
            if (magnitude > pivot_magnitude) {
                pivot = row;
                pivot_magnitude = magnitude;
            }
        }
        if (pivot_magnitude <= 64.0 * DBL_MIN) {
            return 0;
        }
        pivots[column] = pivot;
        if (pivot != column) {
            for (next_column = 0U; next_column < dimension; ++next_column) {
                const double temporary = matrix[column][next_column];
                matrix[column][next_column] = matrix[pivot][next_column];
                matrix[pivot][next_column] = temporary;
            }
        }
        for (row = column + 1U; row < dimension; ++row) {
            matrix[row][column] /= matrix[column][column];
            for (next_column = column + 1U;
                 next_column < dimension;
                 ++next_column) {
                matrix[row][next_column] -= matrix[row][column]
                    * matrix[column][next_column];
            }
        }
    }
    return 1;
}

static void lu_solve(
    double matrix[RP_RADAU_SYSTEM_DIMENSION][RP_RADAU_SYSTEM_DIMENSION],
    const size_t pivots[RP_RADAU_SYSTEM_DIMENSION],
    double right_hand_side[RP_RADAU_SYSTEM_DIMENSION],
    size_t dimension
)
{
    size_t row;

    for (row = 0U; row < dimension; ++row) {
        const size_t pivot = pivots[row];
        size_t column;

        if (pivot != row) {
            const double temporary = right_hand_side[row];
            right_hand_side[row] = right_hand_side[pivot];
            right_hand_side[pivot] = temporary;
        }
        for (column = 0U; column < row; ++column) {
            right_hand_side[row] -= matrix[row][column]
                * right_hand_side[column];
        }
    }
    for (row = dimension; row-- > 0U;) {
        size_t column;

        for (column = row + 1U; column < dimension; ++column) {
            right_hand_side[row] -= matrix[row][column]
                * right_hand_side[column];
        }
        right_hand_side[row] /= matrix[row][row];
    }
}

static rp_integrator_status finite_difference_jacobian(
    const rp_integrator_config *config,
    rp_integrator_rhs rhs,
    void *context,
    double independent_variable,
    const double state[],
    const double derivative[],
    double jacobian[RP_INTEGRATOR_MAX_DIMENSION][RP_INTEGRATOR_MAX_DIMENSION],
    rp_integrator_stats *stats
)
{
    double perturbed_state[RP_INTEGRATOR_MAX_DIMENSION];
    double perturbed_derivative[RP_INTEGRATOR_MAX_DIMENSION];
    size_t row;
    size_t column;

    for (row = 0U; row < config->dimension; ++row) {
        perturbed_state[row] = state[row];
    }
    for (column = 0U; column < config->dimension; ++column) {
        const double perturbation = sqrt(DBL_EPSILON)
            * fmax(1.0, fabs(state[column]));
        rp_integrator_status status;

        perturbed_state[column] += perturbation;
        status = evaluate_rhs(
            rhs,
            context,
            independent_variable,
            perturbed_state,
            perturbed_derivative,
            config->dimension,
            stats
        );
        perturbed_state[column] = state[column];
        if (status != RP_INTEGRATOR_STATUS_OK) {
            return status;
        }
        for (row = 0U; row < config->dimension; ++row) {
            jacobian[row][column] = (
                perturbed_derivative[row] - derivative[row]
            ) / perturbation;
        }
    }
    ++stats->jacobian_evaluations;
    return RP_INTEGRATOR_STATUS_OK;
}

/** Evaluate the active Jacobian without changing RHS evaluation accounting. */
static rp_integrator_status evaluate_jacobian(
    const rp_integrator_config *config,
    rp_integrator_rhs rhs,
    void *context,
    double independent_variable,
    const double state[],
    const double derivative[],
    double jacobian[RP_INTEGRATOR_MAX_DIMENSION][RP_INTEGRATOR_MAX_DIMENSION],
    rp_integrator_stats *stats
)
{
    size_t row;
    size_t column;

    if (config->jacobian == NULL) {
        return finite_difference_jacobian(
            config, rhs, context, independent_variable, state, derivative,
            jacobian, stats
        );
    }
    /* Missing active entries must fail rather than use indeterminate values. */
    for (row = 0U; row < config->dimension; ++row) {
        for (column = 0U; column < config->dimension; ++column) {
            jacobian[row][column] = NAN;
        }
    }
    if (config->jacobian(
            independent_variable, state, &jacobian[0][0], config->dimension,
            RP_INTEGRATOR_MAX_DIMENSION, context
        ) != 0) {
        return RP_INTEGRATOR_STATUS_RHS_FAILURE;
    }
    for (row = 0U; row < config->dimension; ++row) {
        if (!vector_is_finite(jacobian[row], config->dimension)) {
            return RP_INTEGRATOR_STATUS_NONFINITE_VALUE;
        }
    }
    ++stats->jacobian_evaluations;
    return RP_INTEGRATOR_STATUS_OK;
}

static rp_integrator_status radau_collocation_step(
    const rp_integrator_config *config,
    rp_integrator_rhs rhs,
    void *context,
    double independent_variable,
    const double state[],
    double step,
    double candidate[],
    rp_integrator_stats *stats
)
{
    double derivative[RP_INTEGRATOR_MAX_DIMENSION];
    double jacobian[RP_INTEGRATOR_MAX_DIMENSION][RP_INTEGRATOR_MAX_DIMENSION];
    double stages[RP_RADAU_STAGES][RP_INTEGRATOR_MAX_DIMENSION];
    double stage_derivatives[RP_RADAU_STAGES][RP_INTEGRATOR_MAX_DIMENSION];
    double matrix[RP_RADAU_SYSTEM_DIMENSION][RP_RADAU_SYSTEM_DIMENSION];
    double correction[RP_RADAU_SYSTEM_DIMENSION];
    size_t pivots[RP_RADAU_SYSTEM_DIMENSION];
    const size_t system_dimension = RP_RADAU_STAGES * config->dimension;
    const double newton_tolerance = fmax(
        10.0 * DBL_EPSILON / config->relative_tolerance,
        fmin(0.03, sqrt(config->relative_tolerance))
    );
    size_t stage;
    size_t component;
    size_t coupled_stage;
    size_t coupled_component;
    size_t iteration;
    rp_integrator_status status;

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
    status = evaluate_jacobian(
        config,
        rhs,
        context,
        independent_variable,
        state,
        derivative,
        jacobian,
        stats
    );
    if (status != RP_INTEGRATOR_STATUS_OK) {
        return status;
    }

    for (stage = 0U; stage < RP_RADAU_STAGES; ++stage) {
        for (component = 0U; component < config->dimension; ++component) {
            const size_t row = stage * config->dimension + component;

            stages[stage][component] = state[component]
                + RP_RADAU_C[stage] * step * derivative[component];
            for (coupled_stage = 0U;
                 coupled_stage < RP_RADAU_STAGES;
                 ++coupled_stage) {
                for (coupled_component = 0U;
                     coupled_component < config->dimension;
                     ++coupled_component) {
                    const size_t column = coupled_stage * config->dimension
                        + coupled_component;
                    const double identity = stage == coupled_stage
                            && component == coupled_component
                        ? 1.0
                        : 0.0;
                    matrix[row][column] = identity
                        - step * RP_RADAU_A[stage][coupled_stage]
                            * jacobian[component][coupled_component];
                }
            }
        }
    }
    if (!lu_decompose(matrix, pivots, system_dimension)) {
        return RP_INTEGRATOR_STATUS_CONVERGENCE_FAILURE;
    }

    for (iteration = 0U;
         iteration < RP_RADAU_MAXIMUM_NEWTON_ITERATIONS;
         ++iteration) {
        double correction_norm_squared = 0.0;

        for (stage = 0U; stage < RP_RADAU_STAGES; ++stage) {
            status = evaluate_rhs(
                rhs,
                context,
                independent_variable + RP_RADAU_C[stage] * step,
                stages[stage],
                stage_derivatives[stage],
                config->dimension,
                stats
            );
            if (status != RP_INTEGRATOR_STATUS_OK) {
                return status;
            }
        }
        for (stage = 0U; stage < RP_RADAU_STAGES; ++stage) {
            for (component = 0U; component < config->dimension; ++component) {
                const size_t row = stage * config->dimension + component;
                double residual = stages[stage][component] - state[component];

                for (coupled_stage = 0U;
                     coupled_stage < RP_RADAU_STAGES;
                     ++coupled_stage) {
                    residual -= step * RP_RADAU_A[stage][coupled_stage]
                        * stage_derivatives[coupled_stage][component];
                }
                correction[row] = -residual;
            }
        }
        lu_solve(matrix, pivots, correction, system_dimension);
        ++stats->linear_solves;
        for (stage = 0U; stage < RP_RADAU_STAGES; ++stage) {
            for (component = 0U; component < config->dimension; ++component) {
                const size_t row = stage * config->dimension + component;
                const double scale = rp_integrator_scale(
                    config, component,
                    fmax(fabs(state[component]),
                         fabs(stages[stage][component]))
                );
                const double scaled_correction = correction[row] / scale;

                stages[stage][component] += correction[row];
                correction_norm_squared += scaled_correction
                    * scaled_correction;
            }
        }
        if (sqrt(correction_norm_squared / (double)system_dimension)
            <= newton_tolerance) {
            for (component = 0U; component < config->dimension; ++component) {
                candidate[component] = stages[RP_RADAU_STAGES - 1U][component];
            }
            return vector_is_finite(candidate, config->dimension)
                ? RP_INTEGRATOR_STATUS_OK
                : RP_INTEGRATOR_STATUS_NONFINITE_VALUE;
        }
    }
    return RP_INTEGRATOR_STATUS_CONVERGENCE_FAILURE;
}

static double choose_initial_step(
    const rp_integrator_config *config,
    const double state[],
    const double derivative[],
    double interval
)
{
    double state_norm_squared = 0.0;
    double derivative_norm_squared = 0.0;
    double step;
    size_t component;

    if (config->initial_step > 0.0) {
        step = config->initial_step;
    } else {
        for (component = 0U; component < config->dimension; ++component) {
            const double scale = rp_integrator_scale(
                config, component, fabs(state[component])
            );
            const double scaled_state = state[component] / scale;
            const double scaled_derivative = derivative[component] / scale;

            state_norm_squared += scaled_state * scaled_state;
            derivative_norm_squared += scaled_derivative * scaled_derivative;
        }
        if (state_norm_squared < 1.0e-10
            || derivative_norm_squared < 1.0e-10) {
            step = 1.0e-6;
        } else {
            step = 0.01 * sqrt(state_norm_squared / derivative_norm_squared);
        }
    }
    step = fmin(step, interval);
    if (config->maximum_step > 0.0) {
        step = fmin(step, config->maximum_step);
    }
    return step;
}

rp_integrator_status rp_radau_integrate(
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
    double full_step_state[RP_INTEGRATOR_MAX_DIMENSION];
    double first_half_state[RP_INTEGRATOR_MAX_DIMENSION];
    double second_half_state[RP_INTEGRATOR_MAX_DIMENSION];
    const double direction = final_independent_variable
        > initial_independent_variable ? 1.0 : -1.0;
    double independent_variable = initial_independent_variable;
    double step_magnitude;
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
    step_magnitude = choose_initial_step(
        config,
        state,
        derivative,
        fabs(final_independent_variable - initial_independent_variable)
    );

    while (direction * (final_independent_variable - independent_variable) > 0.0) {
        const double remaining = fabs(
            final_independent_variable - independent_variable
        );
        const double minimum_step = 10.0 * fabs(
            nextafter(independent_variable, direction * INFINITY)
                - independent_variable
        );
        double step;
        double error_norm_squared = 0.0;
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

        status = radau_collocation_step(
            config,
            rhs,
            context,
            independent_variable,
            state,
            step,
            full_step_state,
            stats
        );
        if (status == RP_INTEGRATOR_STATUS_OK) {
            status = radau_collocation_step(
                config,
                rhs,
                context,
                independent_variable,
                state,
                0.5 * step,
                first_half_state,
                stats
            );
        }
        if (status == RP_INTEGRATOR_STATUS_OK) {
            status = radau_collocation_step(
                config,
                rhs,
                context,
                independent_variable + 0.5 * step,
                first_half_state,
                0.5 * step,
                second_half_state,
                stats
            );
        }
        if (status == RP_INTEGRATOR_STATUS_CONVERGENCE_FAILURE) {
            step_magnitude *= 0.25;
            ++stats->rejected_steps;
            continue;
        }
        if (status != RP_INTEGRATOR_STATUS_OK) {
            return status;
        }

        for (component = 0U; component < config->dimension; ++component) {
            const double scale = rp_integrator_scale(
                config, component,
                fmax(fabs(state[component]),
                     fabs(second_half_state[component]))
            );
            const double scaled_error = (
                second_half_state[component] - full_step_state[component]
            ) / (31.0 * scale);
            error_norm_squared += scaled_error * scaled_error;
        }
        error_norm = sqrt(
            error_norm_squared / (double)config->dimension
        );
        if (!isfinite(error_norm)) {
            return RP_INTEGRATOR_STATUS_NONFINITE_VALUE;
        }
        factor = error_norm == 0.0
            ? RP_RADAU_MAXIMUM_FACTOR
            : RP_RADAU_SAFETY * pow(error_norm, -1.0 / (RP_RADAU_ORDER + 1.0));
        factor = fmax(
            RP_RADAU_MINIMUM_FACTOR,
            fmin(RP_RADAU_MAXIMUM_FACTOR, factor)
        );

        if (error_norm <= 1.0) {
            if (config->method == RP_INTEGRATOR_METHOD_PROJECTION_RADAU) {
                if (config->projector(
                        independent_variable + step,
                        second_half_state,
                        context
                    ) != 0
                    || !vector_is_finite(
                        second_half_state, config->dimension
                    )) {
                    return RP_INTEGRATOR_STATUS_PROJECTION_FAILURE;
                }
            }
            independent_variable += step;
            for (component = 0U; component < config->dimension; ++component) {
                state[component] = second_half_state[component];
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
        } else {
            step_magnitude *= factor;
            ++stats->rejected_steps;
        }
    }
    stats->final_independent_variable = final_independent_variable;
    return RP_INTEGRATOR_STATUS_OK;
}
