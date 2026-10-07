/**
 * @file test_kerr_null_solve.c
 * @brief Contract tests for the private coordinate-time Kerr photon solver.
 *
 * Units: G = c = M = 1. Tests use only the fixed native interfaces and own
 * every trajectory until its matching free call. No external test framework
 * is needed. Numerical tolerances below distinguish integration error from
 * the cubic-Hermite interpolation error required by the solver contract.
 */

#include "geodesic/null/solve.h"
#include "relatipy/kerr_null.h"

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

static unsigned int failures;

#define CHECK(condition) do { \
    if (!(condition)) { \
        fprintf(stderr, "%s:%d: CHECK failed: %s\n", \
            __FILE__, __LINE__, #condition); \
        ++failures; \
    } \
} while (0)

static const rp_integrator_method methods[] = {
    RP_INTEGRATOR_METHOD_RADAU,
    RP_INTEGRATOR_METHOD_DOP853,
    RP_INTEGRATOR_METHOD_DP45
};

static void report_case(const char *name, rp_integrator_method method,
    unsigned int before)
{
    printf("%s %s: %s\n", name, rp_integrator_method_name(method),
        failures == before ? "PASS" : "FAIL");
}

static rp_integrator_status solve(double spin, const double state[8],
    double t_final, const double *t_eval, size_t n_eval, double r_escape,
    rp_integrator_method method, int store_steps,
    rp_kerr_null_trajectory *trajectory)
{
    return rp_kerr_null_trajectory_solve(spin, state, t_final, t_eval,
        n_eval, r_escape, method, 1.0e-10, 1.0e-12, NULL,
        store_steps, trajectory);
}

static int state_from_constants(double spin, double t0, double radius,
    double theta, double b, double eta, int radial_sign, double state[8])
{
    const double position[4] = {t0, radius, theta, 0.25};
    const rp_kerr_status status = rp_kerr_null_state_from_constants(
        spin, position, b, eta, radial_sign, 1, state);
    CHECK(status == RP_KERR_STATUS_OK);
    return status == RP_KERR_STATUS_OK;
}

static void radial_state(double t0, double radius, double state[8])
{
    const double initial[8] = {
        t0, radius, 0.5 * acos(-1.0), 0.0,
        1.0, -(1.0 - 2.0 / radius), 0.0, 0.0
    };
    memcpy(state, initial, sizeof(initial));
}

static int usable_rows(const rp_kerr_null_trajectory *trajectory)
{
    const int usable = trajectory->count > 0U
        && trajectory->states != NULL && trajectory->lambdas != NULL
        && trajectory->capacity >= trajectory->count;
    CHECK(usable);
    return usable;
}

static void check_endpoint(const rp_kerr_null_trajectory *trajectory)
{
    if (usable_rows(trajectory)) {
        const size_t last = trajectory->count - 1U;
        CHECK(memcmp(trajectory->final_state,
            trajectory->states + 8U * last, 8U * sizeof(double)) == 0);
        CHECK(trajectory->final_lambda == trajectory->lambdas[last]);
        CHECK(trajectory->lambdas[0] == 0.0);
    }
}

static void test_radial(rp_integrator_method method)
{
    const unsigned int before = failures;
    const double t0 = 2.5;
    const double r0 = 10.0;
    const double threshold = 2.0 * (1.0 + 1.0e-6);
    rp_kerr_null_trajectory trajectory = {0};
    double state[8];
    double largest_error = 0.0;
    double worst_radius = r0;
    int exterior = 1;
    int finite_rows = 1;
    size_t row;
    rp_integrator_status status;

    radial_state(t0, r0, state);
    status = solve(0.0, state, t0 + 200.0, NULL, 0U, 0.0,
        method, 1, &trajectory);
    CHECK(status == RP_INTEGRATOR_STATUS_OBSERVER_STOPPED);
    CHECK(trajectory.integrator_status == status);
    CHECK(trajectory.termination == RP_KERR_NULL_TERMINATION_HORIZON);
    CHECK((int)trajectory.termination == 1);
    CHECK(trajectory.allocation_failed == 0);
    if (usable_rows(&trajectory)) {
        CHECK(trajectory.count > 2U);
        CHECK(memcmp(trajectory.states, state, sizeof(state)) == 0);
        for (row = 0U; row < trajectory.count; ++row) {
            const double *sample = trajectory.states + 8U * row;
            double exact_elapsed;
            double error;
            size_t component;
            if (!(sample[1] > threshold)) {
                exterior = 0;
                continue;
            }
            /* Integrating dt/dr = -1/(1-2/r) gives this exact relation.
             * max(1, elapsed) avoids division by zero at the initial row.
             * The requested 1e-7 bound allows accumulated error near the BL
             * horizon while still testing every stored exterior sample. */
            exact_elapsed = (r0 - sample[1])
                + 2.0 * log((r0 - 2.0) / (sample[1] - 2.0));
            error = fabs((sample[0] - t0) - exact_elapsed)
                / fmax(1.0, fabs(exact_elapsed));
            if (!isfinite(error)) {
                finite_rows = 0;
            } else if (error > largest_error) {
                largest_error = error;
                worst_radius = sample[1];
            }
            for (component = 0U; component < 8U; ++component) {
                if (!isfinite(sample[component])) {
                    finite_rows = 0;
                }
            }
            if (row > 0U) {
                CHECK(trajectory.lambdas[row] > trajectory.lambdas[row - 1U]);
                CHECK(sample[0] > trajectory.states[8U * (row - 1U)]);
            }
        }
        CHECK(exterior);
        CHECK(finite_rows);
        CHECK(largest_error < 1.0e-7);
        CHECK(trajectory.final_state[1] > threshold);
        check_endpoint(&trajectory);
    }
    printf("radial %s: status=%d rows=%zu max_relative_error=%.17g r=%.17g\n",
        rp_integrator_method_name(method), (int)status, trajectory.count,
        largest_error, worst_radius);
    rp_kerr_null_trajectory_free(&trajectory);
    report_case("radial/HORIZON", method, before);
}

/* Independent analytic inversion for the radial r0=10 Schwarzschild case.
 * t(r) is strictly decreasing on (2,10]; long double keeps reference rounding
 * below the measured solver errors. Newton always stays inside its bracket.
 */
static long double radial_radius_at_time(double time)
{
    long double lower = 2.0L;
    long double upper = 10.0L;
    long double radius = 10.0L - 0.5L * (long double)time;
    size_t iteration;

    if (time == 0.0) {
        return 10.0L;
    }
    for (iteration = 0U; iteration < 256U; ++iteration) {
        const long double residual = (10.0L - radius)
            + 2.0L * logl(8.0L / (radius - 2.0L)) - (long double)time;
        long double candidate;
        if (residual > 0.0L) {
            lower = radius;
        } else {
            upper = radius;
        }
        candidate = radius + residual * (radius - 2.0L) / radius;
        if (!(candidate > lower && candidate < upper)) {
            candidate = lower + 0.5L * (upper - lower);
        }
        if (candidate == radius
            || upper - lower <= 4.0L * LDBL_EPSILON * radius) {
            break;
        }
        radius = candidate;
    }
    CHECK(iteration < 256U);
    return radius;
}

static void test_radial_interpolation(rp_integrator_method method)
{
    const unsigned int before = failures;
    const double tolerances[2] = {1.0e-10, 1.0e-12};
    /* Measured on 400 uniform t_eval in (0,12], G=c=M=1, atol=rtol*1e-2:
     * method     max |r-r_exact| (rtol=1e-10, 1e-12)     final errors
     * radau      7.08155e-9, 5.55309e-11                5.81767e-10, 1.45730e-11
     * dop853     8.90490e-8, 5.47280e-9                 6.08143e-9,  1.25504e-10
     * dp45       1.12822e-10,9.53839e-13                8.19579e-13, 7.92821e-14
     * Bounds below round approximately ten times each measured maximum up.
     * The cubic-position implementation exceeds every sample bound.
     * These are measured regression bounds, not guarantees equal to rtol.
     */
    const double radau_sample[2] = {7.1e-8, 5.6e-10};
    const double radau_final[2] = {5.9e-9, 1.5e-10};
    const double dop853_sample[2] = {9.0e-7, 5.5e-8};
    const double dop853_final[2] = {6.1e-8, 1.3e-9};
    const double dp45_sample[2] = {1.13e-9, 9.54e-12};
    const double dp45_final[2] = {8.2e-12, 8.0e-13};
    const double *sample_bounds = method == RP_INTEGRATOR_METHOD_RADAU
        ? radau_sample : method == RP_INTEGRATOR_METHOD_DOP853
            ? dop853_sample : dp45_sample;
    const double *final_bounds = method == RP_INTEGRATOR_METHOD_RADAU
        ? radau_final : method == RP_INTEGRATOR_METHOD_DOP853
            ? dop853_final : dp45_final;
    double state[8];
    double times[400];
    size_t row;
    size_t tolerance_index;

    radial_state(0.0, 10.0, state);
    for (row = 0U; row < 400U; ++row) {
        times[row] = 12.0 * (double)(row + 1U) / 400.0;
    }
    for (tolerance_index = 0U; tolerance_index < 2U; ++tolerance_index) {
        rp_kerr_null_trajectory samples = {0};
        rp_kerr_null_trajectory endpoint = {0};
        const double rtol = tolerances[tolerance_index];
        double largest_error = 0.0;
        double final_error;
        const rp_integrator_status status = rp_kerr_null_trajectory_solve(
            0.0, state, 12.0, times, 400U, 0.0, method,
            rtol, rtol * 1.0e-2, NULL, 0, &samples);

        CHECK(status == RP_INTEGRATOR_STATUS_OK);
        CHECK(samples.integrator_status == status);
        CHECK(samples.termination == RP_KERR_NULL_TERMINATION_NONE);
        CHECK(samples.count == 400U);
        if (usable_rows(&samples) && samples.count == 400U) {
            for (row = 0U; row < 400U; ++row) {
                const double error = (double)fabsl(
                    (long double)samples.states[8U * row + 1U]
                    - radial_radius_at_time(times[row]));
                CHECK(samples.states[8U * row] == times[row]);
                CHECK(isfinite(error));
                largest_error = fmax(largest_error, error);
            }
            CHECK(largest_error < sample_bounds[tolerance_index]);
            CHECK(samples.final_state[0] == 12.0);
            CHECK(memcmp(samples.final_state, samples.states + 8U * 399U,
                sizeof(samples.final_state)) == 0);
            CHECK(samples.final_lambda == samples.lambdas[399U]);
            final_error = (double)fabsl((long double)samples.final_state[1]
                - radial_radius_at_time(12.0));
            CHECK(isfinite(final_error));
            CHECK(final_error < final_bounds[tolerance_index]);
        }
        /* Without intermediate samples, the cache must be invalidated as
         * accepted steps advance; the same final-time bound still applies. */
        CHECK(rp_kerr_null_trajectory_solve(0.0, state, 12.0, NULL, 0U,
            0.0, method, rtol, rtol * 1.0e-2, NULL, 0, &endpoint)
            == RP_INTEGRATOR_STATUS_OK);
        CHECK(endpoint.termination == RP_KERR_NULL_TERMINATION_NONE);
        CHECK(endpoint.count == 2U);
        CHECK(endpoint.final_state[0] == 12.0);
        final_error = (double)fabsl((long double)endpoint.final_state[1]
            - radial_radius_at_time(12.0));
        CHECK(isfinite(final_error));
        CHECK(final_error < final_bounds[tolerance_index]);
        check_endpoint(&endpoint);
        printf("radial-interpolation %s rtol=%.1e: max_radius_error=%.17g "
            "final_radius_error=%.17g\n", rp_integrator_method_name(method),
            rtol, largest_error, final_error);
        rp_kerr_null_trajectory_free(&endpoint);
        rp_kerr_null_trajectory_free(&samples);
    }
    report_case("radial-interpolation", method, before);
}

static void test_photon_orbit(rp_integrator_method method)
{
    const unsigned int before = failures;
    rp_kerr_null_trajectory trajectory = {0};
    double state[8];
    double largest_error = 0.0;
    size_t row;
    rp_integrator_status status;

    if (!state_from_constants(0.0, 0.0, 3.0, 0.5 * acos(-1.0),
        3.0 * sqrt(3.0), 0.0, 1, state)) {
        report_case("photon-orbit", method, before);
        return;
    }
    status = solve(0.0, state, 20.0, NULL, 0U, 0.0, method, 1, &trajectory);
    CHECK(status == RP_INTEGRATOR_STATUS_OK);
    CHECK(trajectory.termination == RP_KERR_NULL_TERMINATION_NONE);
    if (usable_rows(&trajectory)) {
        for (row = 0U; row < trajectory.count; ++row) {
            const double error = fabs(trajectory.states[8U * row + 1U] - 3.0);
            CHECK(isfinite(error));
            largest_error = fmax(largest_error, error);
        }
        /* The circular photon orbit is unstable. This short t=20 test and
         * absolute 1e-3 radius bound test the known r=3 solution without
         * claiming long-time stability of the numerical orbit. */
        CHECK(largest_error < 1.0e-3);
        CHECK(trajectory.final_state[0] == 20.0);
        check_endpoint(&trajectory);
    }
    printf("photon-orbit %s: max_radius_error=%.17g\n",
        rp_integrator_method_name(method), largest_error);
    rp_kerr_null_trajectory_free(&trajectory);
    report_case("photon-orbit", method, before);
}

static void test_capture_escape(rp_integrator_method method)
{
    const double b_values[2] = {5.0, 6.0};
    const rp_kerr_null_termination expected[2] = {
        RP_KERR_NULL_TERMINATION_HORIZON, RP_KERR_NULL_TERMINATION_ESCAPE
    };
    const unsigned int before = failures;
    size_t case_index;

    for (case_index = 0U; case_index < 2U; ++case_index) {
        rp_kerr_null_trajectory trajectory = {0};
        double state[8];
        rp_integrator_status status;
        if (!state_from_constants(0.0, 0.0, 30.0, 0.5 * acos(-1.0),
            b_values[case_index], 0.0, -1, state)) {
            continue;
        }
        /* Both begin inward. For b=6 the trajectory turns before reaching
         * the horizon and subsequently crosses the accepted-step escape
         * radius; t=500 leaves time for either terminal event. */
        status = solve(0.0, state, 500.0, NULL, 0U, 40.0,
            method, 1, &trajectory);
        CHECK(status == RP_INTEGRATOR_STATUS_OBSERVER_STOPPED);
        CHECK(trajectory.integrator_status == status);
        CHECK(trajectory.termination == expected[case_index]);
        check_endpoint(&trajectory);
        if (case_index == 0U) {
            CHECK((int)trajectory.termination == 1);
            CHECK(trajectory.final_state[1] > 2.0 * (1.0 + 1.0e-6));
        } else {
            CHECK((int)trajectory.termination == 2);
            CHECK(trajectory.final_state[1] >= 40.0);
            CHECK(trajectory.final_state[5] > 0.0);
        }
        printf("capture/escape %s b=%.1f: status=%d termination=%d r=%.17g\n",
            rp_integrator_method_name(method), b_values[case_index],
            (int)status, (int)trajectory.termination, trajectory.final_state[1]);
        rp_kerr_null_trajectory_free(&trajectory);
    }
    report_case("capture/escape", method, before);
}

/* Independent contraction of the metric and null Carter expression.
 * Long-double accumulation reduces cancellation in the diagnostic oracle;
 * production still evolves and stores only double precision. */
static int invariant_oracle(double spin, const double state[8],
    rp_kerr_null_invariants *invariants)
{
    double metric[RP_KERR_DIM][RP_KERR_DIM];
    long double momentum[RP_KERR_DIM] = {0.0L, 0.0L, 0.0L, 0.0L};
    long double norm = 0.0L;
    long double denominator = 0.0L;
    const double sine = sin(state[2]);
    const double cosine = cos(state[2]);
    size_t mu;
    size_t nu;
    const rp_kerr_status status = rp_kerr_metric(1.0, spin, state, metric);
    CHECK(status == RP_KERR_STATUS_OK);
    if (status != RP_KERR_STATUS_OK) {
        return 0;
    }
    for (mu = 0U; mu < RP_KERR_DIM; ++mu) {
        for (nu = 0U; nu < RP_KERR_DIM; ++nu) {
            const long double term = (long double)metric[mu][nu]
                * state[mu + RP_KERR_DIM] * state[nu + RP_KERR_DIM];
            momentum[mu] += (long double)metric[mu][nu]
                * state[nu + RP_KERR_DIM];
            norm += term;
            denominator += fabsl(term);
        }
    }
    invariants->energy = (double)-momentum[0];
    invariants->axial_angular_momentum = (double)momentum[3];
    invariants->carter_constant = (double)(momentum[2] * momentum[2])
        + cosine * cosine * (invariants->axial_angular_momentum
            * invariants->axial_angular_momentum / (sine * sine)
            - spin * spin * invariants->energy * invariants->energy);
    invariants->impact_parameter = invariants->axial_angular_momentum
        / invariants->energy;
    invariants->eta = invariants->carter_constant
        / (invariants->energy * invariants->energy);
    invariants->norm = (double)norm;
    invariants->relative_norm = (double)(fabsl(norm) / denominator);
    return 1;
}

static void test_kerr_invariants(rp_integrator_method method)
{
    const unsigned int before = failures;
    const double spin = 0.9;
    rp_kerr_null_trajectory trajectory = {0};
    rp_kerr_null_invariants initial;
    double state[8];
    double largest[4] = {0.0, 0.0, 0.0, 0.0};
    double accepted_largest[4] = {0.0, 0.0, 0.0, 0.0};
    double final_error[4] = {0.0, 0.0, 0.0, 0.0};
    size_t worst_rows[4] = {0U, 0U, 0U, 0U};
    int finite_diagnostics = 1;
    size_t row;
    rp_integrator_status status;

    /* The requested b=2, eta=20 tangent exists at r=8, theta=1:
     * R/E^2 = (r^2+a^2-a*b)^2 - Delta*((b-a)^2+eta) > 0,
     * Theta/E^2 = eta-cos(theta)^2*(b^2/sin(theta)^2-a^2) > 0.
     * Choose the outward radial sign to exercise the full t=200 interval
     * without the BL horizon's conditioning dominating this drift test. */
    {
        const double delta = 8.0 * 8.0 - 2.0 * 8.0 + spin * spin;
        const double p = 8.0 * 8.0 + spin * spin - spin * 2.0;
        const double radial_potential = p * p
            - delta * ((2.0 - spin) * (2.0 - spin) + 20.0);
        const double polar_potential = 20.0 - cos(1.0) * cos(1.0)
            * (4.0 / (sin(1.0) * sin(1.0)) - spin * spin);
        CHECK(radial_potential > 0.0);
        CHECK(polar_potential > 0.0);
    }
    if (!state_from_constants(spin, 0.0, 8.0, 1.0,
        2.0, 20.0, 1, state) || !invariant_oracle(spin, state, &initial)) {
        report_case("Kerr-invariants", method, before);
        return;
    }
    status = solve(spin, state, 200.0, NULL, 0U, 0.0,
        method, 1, &trajectory);
    CHECK(status == RP_INTEGRATOR_STATUS_OK);
    CHECK(trajectory.termination == RP_KERR_NULL_TERMINATION_NONE);
    if (usable_rows(&trajectory)) {
        CHECK(trajectory.final_state[0] == 200.0);
        for (row = 0U; row < trajectory.count; ++row) {
            rp_kerr_null_invariants sample;
            rp_kerr_null_invariants reported = {0};
            const double *saved = trajectory.states + 8U * row;
            double error[4];
            size_t quantity;
            CHECK(rp_kerr_null_invariants_evaluate(spin, saved, &reported)
                == RP_KERR_STATUS_OK);
            if (!invariant_oracle(spin, saved, &sample)) {
                continue;
            }
            error[0] = fabs(sample.energy - initial.energy) / fabs(initial.energy);
            error[1] = fabs(sample.axial_angular_momentum
                - initial.axial_angular_momentum)
                / fabs(initial.axial_angular_momentum);
            error[2] = fabs(sample.carter_constant - initial.carter_constant)
                / fabs(initial.carter_constant);
            error[3] = fmax(sample.relative_norm, reported.relative_norm);
            for (quantity = 0U; quantity < 4U; ++quantity) {
                if (!isfinite(error[quantity])) {
                    finite_diagnostics = 0;
                }
                if (error[quantity] > largest[quantity]) {
                    largest[quantity] = error[quantity];
                    worst_rows[quantity] = row;
                }
                if (row + 1U < trajectory.count) {
                    accepted_largest[quantity] = fmax(accepted_largest[quantity],
                        error[quantity]);
                } else {
                    final_error[quantity] = error[quantity];
                }
            }
        }
        /* 1e-6 is the requested long-interval relative drift limit, measured
         * against each nonzero initial invariant. It includes the final
         * Hermite sample as well as every accepted step; no projection or
         * renormalization is applied to any saved tangent. */
        CHECK(finite_diagnostics);
        CHECK(largest[0] < 1.0e-6);
        CHECK(largest[1] < 1.0e-6);
        CHECK(largest[2] < 1.0e-6);
        CHECK(largest[3] < 1.0e-6);
        check_endpoint(&trajectory);
    }
    printf("Kerr-invariants %s: drift E=%.17g Lz=%.17g Q=%.17g norm=%.17g\n",
        rp_integrator_method_name(method), largest[0], largest[1],
        largest[2], largest[3]);
    printf("Kerr-invariants %s accepted: E=%.17g Lz=%.17g Q=%.17g norm=%.17g\n",
        rp_integrator_method_name(method), accepted_largest[0],
        accepted_largest[1], accepted_largest[2], accepted_largest[3]);
    printf("Kerr-invariants %s final: E=%.17g Lz=%.17g Q=%.17g norm=%.17g\n",
        rp_integrator_method_name(method), final_error[0], final_error[1],
        final_error[2], final_error[3]);
    if (trajectory.count > 0U && trajectory.states != NULL) {
        size_t quantity;
        for (quantity = 0U; quantity < 4U; ++quantity) {
            printf("Kerr-invariants %s quantity=%zu worst_row=%zu t=%.17g\n",
                rp_integrator_method_name(method), quantity, worst_rows[quantity],
                trajectory.states[8U * worst_rows[quantity]]);
        }
    }
    rp_kerr_null_trajectory_free(&trajectory);
    report_case("Kerr-invariants", method, before);
}

static void test_t_final(rp_integrator_method method)
{
    const unsigned int before = failures;
    const double t_final = 19.123456789;
    double state[8];
    int store_steps;

    if (!state_from_constants(0.9, 3.5, 20.0, 1.0, 3.0, 5.0, 1, state)) {
        report_case("t_final", method, before);
        return;
    }
    for (store_steps = 0; store_steps <= 1; ++store_steps) {
        rp_kerr_null_trajectory trajectory = {0};
        const rp_integrator_status status = solve(0.9, state, t_final,
            NULL, 0U, 0.0, method, store_steps, &trajectory);
        CHECK(status == RP_INTEGRATOR_STATUS_OK);
        CHECK(trajectory.integrator_status == status);
        CHECK(trajectory.termination == RP_KERR_NULL_TERMINATION_NONE);
        CHECK((int)trajectory.termination == 0);
        if (usable_rows(&trajectory)) {
            /* The public target is written exactly after Hermite location.
             * The 1e-12 scaled check also documents the requested bound. */
            CHECK(fabs(trajectory.final_state[0] - t_final)
                <= 1.0e-12 * fmax(1.0, fabs(t_final)));
            CHECK(trajectory.final_state[0] == t_final);
            CHECK(memcmp(trajectory.states, state, sizeof(state)) == 0);
            CHECK(trajectory.count >= 2U);
            if (store_steps == 0) {
                CHECK(trajectory.count == 2U);
            }
            check_endpoint(&trajectory);
        }
        rp_kerr_null_trajectory_free(&trajectory);
    }
    report_case("t_final", method, before);
}

static void test_t_eval(rp_integrator_method method)
{
    const unsigned int before = failures;
    const size_t comparison_rows[3] = {13U, 26U, 41U};
    const double t0 = 3.5;
    const double t_final = 38.5;
    rp_kerr_null_trajectory trajectory = {0};
    double state[8];
    double times[50];
    double largest_comparison_error = 0.0;
    size_t row;
    size_t comparison;
    rp_integrator_status status;

    if (!state_from_constants(0.9, t0, 20.0, 1.0, 3.0, 5.0, 1, state)) {
        report_case("t_eval", method, before);
        return;
    }
    for (row = 0U; row < 50U; ++row) {
        times[row] = t0 + (t_final - t0) * (double)row / 49.0;
    }
    times[0] = t0;
    times[49] = t_final;
    status = solve(0.9, state, t_final, times, 50U, 0.0,
        method, 1, &trajectory);
    CHECK(status == RP_INTEGRATOR_STATUS_OK);
    CHECK(trajectory.termination == RP_KERR_NULL_TERMINATION_NONE);
    CHECK(trajectory.count == 50U);
    if (usable_rows(&trajectory) && trajectory.count == 50U) {
        CHECK(memcmp(trajectory.states, state, sizeof(state)) == 0);
        for (row = 0U; row < 50U; ++row) {
            CHECK(trajectory.states[8U * row] == times[row]);
        }
        check_endpoint(&trajectory);
        for (comparison = 0U; comparison < 3U; ++comparison) {
            rp_kerr_null_trajectory independent = {0};
            const size_t selected = comparison_rows[comparison];
            /* A tighter independent integration changes accepted-step sizes
             * and avoids comparing identical Hermite brackets by accident. */
            const rp_integrator_status reference_status =
                rp_kerr_null_trajectory_solve(0.9, state, times[selected],
                    NULL, 0U, 0.0, method, 1.0e-12, 1.0e-14,
                    NULL, 0, &independent);
            size_t component;
            CHECK(reference_status == RP_INTEGRATOR_STATUS_OK);
            CHECK(independent.count == 2U);
            if (usable_rows(&independent)) {
                for (component = 0U; component < 8U; ++component) {
                    const double reference = independent.final_state[component];
                    const double error = fabs(trajectory.states[8U * selected
                        + component] - reference) / fmax(1.0, fabs(reference));
                    /* Cubic Hermite is lower order than DOP853 and samples
                     * use different containing steps. A 1e-6 componentwise
                     * scaled bound tests consistency, including small tangent
                     * components, without confusing it with rtol=1e-10. */
                    CHECK(isfinite(error));
                    largest_comparison_error = fmax(largest_comparison_error,
                        error);
                }
            }
            rp_kerr_null_trajectory_free(&independent);
        }
        CHECK(largest_comparison_error < 1.0e-6);
    }
    printf("t_eval %s: max_scaled_independent_error=%.17g\n",
        rp_integrator_method_name(method), largest_comparison_error);
    rp_kerr_null_trajectory_free(&trajectory);
    report_case("t_eval", method, before);
}

static void test_continuation(rp_integrator_method method)
{
    const unsigned int before = failures;
    const double position[4] = {0.0, 10.0, 1.1, 0.0};
    rp_kerr_null_trajectory first = {0};
    rp_kerr_null_trajectory continued = {0};
    rp_kerr_null_invariants diagnostics = {0};
    double state[8];
    rp_integrator_status status;

    CHECK(rp_kerr_null_state_from_constants(0.9, position,
        3.0, 5.0, -1, 1, state) == RP_KERR_STATUS_OK);
    status = rp_kerr_null_trajectory_solve(0.9, state, 10.0,
        NULL, 0U, 0.0, method, 1.0e-6, 1.0e-8, NULL, 0, &first);
    CHECK(status == RP_INTEGRATOR_STATUS_OK);
    CHECK(first.termination == RP_KERR_NULL_TERMINATION_NONE);
    CHECK(first.final_state[0] == 10.0);
    CHECK(rp_kerr_null_invariants_evaluate(0.9, first.final_state,
        &diagnostics) == RP_KERR_STATUS_OK);
    /* The first endpoint must actually reproduce the old strict-validation
     * rejection. No tangent projection or reconstruction is permitted. */
    CHECK(diagnostics.relative_norm > 1.0e-8);
    CHECK(diagnostics.relative_norm < 1.0e-3);
    status = rp_kerr_null_trajectory_solve(0.9, first.final_state, 10.01,
        NULL, 0U, 0.0, method, 1.0e-6, 1.0e-8, NULL, 0, &continued);
    CHECK(status == RP_INTEGRATOR_STATUS_OK);
    CHECK(continued.termination == RP_KERR_NULL_TERMINATION_NONE);
    CHECK(continued.final_state[0] == 10.01);
    if (usable_rows(&continued)) {
        CHECK(memcmp(continued.states, first.final_state,
            sizeof(first.final_state)) == 0);
        check_endpoint(&continued);
    }
    printf("continuation %s: first_relative_norm=%.17g status=%d\n",
        rp_integrator_method_name(method), diagnostics.relative_norm,
        (int)status);
    rp_kerr_null_trajectory_free(&continued);
    rp_kerr_null_trajectory_free(&first);
    report_case("continuation", method, before);
}

static void test_large_time_origin(rp_integrator_method method)
{
    const unsigned int before = failures;
    const double t0 = 1.0e16;
    const double t_final = t0 + 10.0;
    const double relative_times[6] = {0.0, 2.0, 4.0, 6.0, 8.0, 10.0};
    double absolute_times[6];
    double reference_state[8];
    double shifted_state[8];
    int sampled;
    size_t row;

    radial_state(0.0, 10.0, reference_state);
    radial_state(t0, 10.0, shifted_state);
    for (row = 0U; row < 6U; ++row) {
        absolute_times[row] = t0 + relative_times[row];
    }
    for (sampled = 0; sampled <= 1; ++sampled) {
        rp_kerr_null_trajectory reference = {0};
        rp_kerr_null_trajectory shifted = {0};
        const size_t count = sampled ? 6U : 0U;
        double radius_error;

        CHECK(solve(0.0, reference_state, 10.0,
            sampled ? relative_times : NULL, count, 0.0,
            method, 1, &reference) == RP_INTEGRATOR_STATUS_OK);
        CHECK(solve(0.0, shifted_state, t_final,
            sampled ? absolute_times : NULL, count, 0.0,
            method, 1, &shifted) == RP_INTEGRATOR_STATUS_OK);
        CHECK(shifted.integrator_status == RP_INTEGRATOR_STATUS_OK);
        CHECK(shifted.termination == RP_KERR_NULL_TERMINATION_NONE);
        CHECK(shifted.final_state[0] == t_final);
        radius_error = fabs(shifted.final_state[1] - reference.final_state[1])
            / fabs(reference.final_state[1]);
        /* Time translation must leave the physical trajectory unchanged.
         * This compares two integrations at identical relative times;
         * 1e-9 allows roundoff, not absolute-time loss of significance. */
        CHECK(isfinite(radius_error));
        CHECK(radius_error < 1.0e-9);
        CHECK(fabsl((long double)reference.final_state[1]
            - radial_radius_at_time(10.0)) < 1.0e-6L);
        if (sampled) {
            CHECK(shifted.count == 6U);
            CHECK(reference.count == 6U);
            if (usable_rows(&shifted) && usable_rows(&reference)
                && shifted.count == 6U && reference.count == 6U) {
                for (row = 0U; row < 6U; ++row) {
                    const double expected_radius = reference.states[8U * row + 1U];
                    CHECK(shifted.states[8U * row] == absolute_times[row]);
                    CHECK(fabs(shifted.states[8U * row + 1U] - expected_radius)
                        / fabs(expected_radius) < 1.0e-9);
                }
            }
        }
        check_endpoint(&shifted);
        printf("large-time-origin %s sampled=%d: relative_radius_error=%.17g "
            "r=%.17g t=%.17g\n", rp_integrator_method_name(method),
            sampled, radius_error, shifted.final_state[1], shifted.final_state[0]);
        rp_kerr_null_trajectory_free(&shifted);
        rp_kerr_null_trajectory_free(&reference);
    }
    report_case("large-time-origin", method, before);
}

static void test_escape_at_final_time(rp_integrator_method method)
{
    const unsigned int before = failures;
    const double t_final = 0.1;
    rp_kerr_null_trajectory baseline = {0};
    double state[8];
    double accepted_radius;
    double escape_after_target;
    size_t direction;

    radial_state(0.0, 10.0, state);
    state[5] = -state[5];
    CHECK(solve(0.0, state, t_final, NULL, 0U, 0.0,
        method, 1, &baseline) == RP_INTEGRATOR_STATUS_OK);
    /* For this radial Schwarzschild null tangent, dk^r/dlambda = 0 and
     * r(lambda) = r0 + k^r lambda. The stats retain the accepted endpoint
     * lambda even though final_lambda is clipped at t_final. Thus this
     * independently proves that the actual last accepted step overshoots
     * the later escape radius; checking only the clipped row cannot do so. */
    accepted_radius = state[1] + state[5]
        * baseline.stats.final_independent_variable;
    CHECK(baseline.stats.final_independent_variable > baseline.final_lambda);
    CHECK(accepted_radius - baseline.final_state[1] > 1.0e-8);
    escape_after_target = baseline.final_state[1]
        + 0.5 * (accepted_radius - baseline.final_state[1]);
    for (direction = 0U; direction < 2U; ++direction) {
        const double r_escape = direction == 0U ? 10.079 : escape_after_target;
        const rp_integrator_status expected_status = direction == 0U
            ? RP_INTEGRATOR_STATUS_OBSERVER_STOPPED : RP_INTEGRATOR_STATUS_OK;
        const rp_kerr_null_termination expected_termination = direction == 0U
            ? RP_KERR_NULL_TERMINATION_ESCAPE : RP_KERR_NULL_TERMINATION_NONE;
        rp_kerr_null_trajectory trajectory = {0};
        const rp_integrator_status status = solve(0.0, state, t_final,
            NULL, 0U, r_escape, method, 1, &trajectory);
        const double actual_accepted_radius = state[1] + state[5]
            * trajectory.stats.final_independent_variable;

        CHECK(status == expected_status);
        CHECK(trajectory.integrator_status == status);
        CHECK(trajectory.termination == expected_termination);
        CHECK(trajectory.final_state[0] == t_final);
        CHECK(trajectory.stats.final_independent_variable
            == baseline.stats.final_independent_variable);
        CHECK(actual_accepted_radius > r_escape);
        if (direction == 0U) {
            CHECK(trajectory.final_state[1] >= r_escape);
        } else {
            CHECK(trajectory.final_state[1] < r_escape);
        }
        check_endpoint(&trajectory);
        printf("escape-at-t_final %s later=%zu: status=%d r=%.17g "
            "r_escape=%.17g accepted_r=%.17g\n",
            rp_integrator_method_name(method), direction, (int)status,
            trajectory.final_state[1], r_escape, actual_accepted_radius);
        rp_kerr_null_trajectory_free(&trajectory);
    }
    rp_kerr_null_trajectory_free(&baseline);
    report_case("escape-at-t_final", method, before);
}

static void expect_rejection(const char *name, double spin,
    const double *state, double t_final, const double *times, size_t count,
    double r_escape, rp_integrator_method method,
    rp_integrator_status expected)
{
    rp_kerr_null_trajectory trajectory = {0};
    rp_integrator_status status;
    /* There is no owned storage. A sentinel verifies entry clearing even
     * when the initial-state pointer is NULL and the result pointer is valid. */
    trajectory.count = 9U;
    trajectory.capacity = 9U;
    status = solve(spin, state, t_final, times, count, r_escape,
        method, 1, &trajectory);
    if (status != expected) {
        fprintf(stderr, "rejection %s: status=%d expected=%d\n",
            name, (int)status, (int)expected);
    }
    CHECK(status == expected);
    CHECK(trajectory.count == 0U);
    CHECK(trajectory.states == NULL);
    CHECK(trajectory.lambdas == NULL);
    rp_kerr_null_trajectory_free(&trajectory);
    rp_kerr_null_trajectory_free(&trajectory);
    CHECK(trajectory.count == 0U);
    CHECK(trajectory.capacity == 0U);
}

static void test_rejections(void)
{
    const unsigned int before = failures;
    const double duplicate[3] = {2.0, 4.0, 4.0};
    const double decreasing[3] = {2.0, 5.0, 4.0};
    const double before_start[2] = {1.0, 4.0};
    const double after_end[2] = {2.0, 11.0};
    const double nonfinite_times[2] = {2.0, NAN};
    const double one_time[1] = {2.0};
    double initial[8];
    double invalid[8];
    rp_kerr_null_trajectory trajectory = {0};
    size_t component;

    radial_state(2.0, 10.0, initial);
    expect_rejection("projection_radau", 0.0, initial, 10.0, NULL, 0U,
        0.0, RP_INTEGRATOR_METHOD_PROJECTION_RADAU,
        RP_INTEGRATOR_STATUS_INVALID_ARGUMENT);
    memcpy(invalid, initial, sizeof(initial));
    invalid[1] = 2.0 * (1.0 + 1.0e-6);
    expect_rejection("threshold equality", 0.0, invalid, 10.0, NULL, 0U,
        0.0, RP_INTEGRATOR_METHOD_RADAU, RP_INTEGRATOR_STATUS_INVALID_ARGUMENT);
    invalid[1] = 2.0;
    expect_rejection("inside threshold", 0.0, invalid, 10.0, NULL, 0U,
        0.0, RP_INTEGRATOR_METHOD_RADAU, RP_INTEGRATOR_STATUS_INVALID_ARGUMENT);
    expect_rejection("t_final equality", 0.0, initial, 2.0, NULL, 0U,
        0.0, RP_INTEGRATOR_METHOD_RADAU, RP_INTEGRATOR_STATUS_INVALID_ARGUMENT);
    expect_rejection("past t_final", 0.0, initial, 1.0, NULL, 0U,
        0.0, RP_INTEGRATOR_METHOD_RADAU, RP_INTEGRATOR_STATUS_INVALID_ARGUMENT);
    expect_rejection("duplicate t_eval", 0.0, initial, 10.0, duplicate, 3U,
        0.0, RP_INTEGRATOR_METHOD_RADAU, RP_INTEGRATOR_STATUS_INVALID_ARGUMENT);
    expect_rejection("decreasing t_eval", 0.0, initial, 10.0, decreasing, 3U,
        0.0, RP_INTEGRATOR_METHOD_RADAU, RP_INTEGRATOR_STATUS_INVALID_ARGUMENT);
    expect_rejection("t_eval before start", 0.0, initial, 10.0,
        before_start, 2U, 0.0, RP_INTEGRATOR_METHOD_RADAU,
        RP_INTEGRATOR_STATUS_INVALID_ARGUMENT);
    expect_rejection("t_eval after end", 0.0, initial, 10.0,
        after_end, 2U, 0.0, RP_INTEGRATOR_METHOD_RADAU,
        RP_INTEGRATOR_STATUS_INVALID_ARGUMENT);
    expect_rejection("NULL t_eval/nonzero n_eval", 0.0, initial, 10.0,
        NULL, 1U, 0.0, RP_INTEGRATOR_METHOD_RADAU,
        RP_INTEGRATOR_STATUS_INVALID_ARGUMENT);
    expect_rejection("non-NULL empty t_eval", 0.0, initial, 10.0,
        one_time, 0U, 0.0, RP_INTEGRATOR_METHOD_RADAU,
        RP_INTEGRATOR_STATUS_INVALID_ARGUMENT);
    expect_rejection("r_escape equality", 0.0, initial, 10.0, NULL, 0U,
        10.0, RP_INTEGRATOR_METHOD_RADAU, RP_INTEGRATOR_STATUS_INVALID_ARGUMENT);
    expect_rejection("r_escape below r0", 0.0, initial, 10.0, NULL, 0U,
        5.0, RP_INTEGRATOR_METHOD_RADAU, RP_INTEGRATOR_STATUS_INVALID_ARGUMENT);
    memcpy(invalid, initial, sizeof(initial));
    invalid[5] *= 0.5;
    {
        rp_kerr_null_invariants diagnostics;
        CHECK(invariant_oracle(0.0, invalid, &diagnostics));
        CHECK(diagnostics.relative_norm > 1.0e-3);
    }
    expect_rejection("non-null state", 0.0, invalid, 10.0, NULL, 0U,
        0.0, RP_INTEGRATOR_METHOD_RADAU, RP_INTEGRATOR_STATUS_INVALID_ARGUMENT);
    memcpy(invalid, initial, sizeof(initial));
    invalid[4] = 1.0 / sqrt(1.0 - 2.0 / invalid[1]);
    invalid[5] = 0.0;
    {
        rp_kerr_null_invariants diagnostics;
        CHECK(invariant_oracle(0.0, invalid, &diagnostics));
        CHECK(fabs(diagnostics.norm + 1.0) < 1.0e-14);
        CHECK(diagnostics.relative_norm > 1.0e-3);
    }
    expect_rejection("timelike norm -1", 0.0, invalid, 10.0, NULL, 0U,
        0.0, RP_INTEGRATOR_METHOD_RADAU, RP_INTEGRATOR_STATUS_INVALID_ARGUMENT);
    memcpy(invalid, initial, sizeof(initial));
    invalid[4] = 0.0;
    expect_rejection("zero k^t", 0.0, invalid, 10.0, NULL, 0U,
        0.0, RP_INTEGRATOR_METHOD_RADAU, RP_INTEGRATOR_STATUS_INVALID_ARGUMENT);
    for (component = 4U; component < 8U; ++component) {
        invalid[component] = -initial[component];
    }
    expect_rejection("past-directed null state", 0.0, invalid, 10.0, NULL, 0U,
        0.0, RP_INTEGRATOR_METHOD_RADAU, RP_INTEGRATOR_STATUS_INVALID_ARGUMENT);
    expect_rejection("negative spin", -0.1, initial, 10.0, NULL, 0U,
        0.0, RP_INTEGRATOR_METHOD_RADAU, RP_INTEGRATOR_STATUS_INVALID_ARGUMENT);
    expect_rejection("spin above one", 1.1, initial, 10.0, NULL, 0U,
        0.0, RP_INTEGRATOR_METHOD_RADAU, RP_INTEGRATOR_STATUS_INVALID_ARGUMENT);
    memcpy(invalid, initial, sizeof(initial));
    invalid[2] = NAN;
    expect_rejection("nonfinite state", 0.0, invalid, 10.0, NULL, 0U,
        0.0, RP_INTEGRATOR_METHOD_RADAU, RP_INTEGRATOR_STATUS_NONFINITE_VALUE);
    expect_rejection("nonfinite t_final", 0.0, initial, NAN, NULL, 0U,
        0.0, RP_INTEGRATOR_METHOD_RADAU, RP_INTEGRATOR_STATUS_NONFINITE_VALUE);
    expect_rejection("nonfinite t_eval", 0.0, initial, 10.0,
        nonfinite_times, 2U, 0.0, RP_INTEGRATOR_METHOD_RADAU,
        RP_INTEGRATOR_STATUS_NONFINITE_VALUE);
    expect_rejection("nonfinite escape", 0.0, initial, 10.0, NULL, 0U,
        NAN, RP_INTEGRATOR_METHOD_RADAU, RP_INTEGRATOR_STATUS_NONFINITE_VALUE);
    expect_rejection("NULL initial state", 0.0, NULL, 10.0, NULL, 0U,
        0.0, RP_INTEGRATOR_METHOD_RADAU, RP_INTEGRATOR_STATUS_NULL_POINTER);
    CHECK(solve(0.0, initial, 10.0, NULL, 0U, 0.0,
        RP_INTEGRATOR_METHOD_RADAU, 1, NULL)
        == RP_INTEGRATOR_STATUS_NULL_POINTER);
    rp_kerr_null_trajectory_free(&trajectory);
    printf("rejections: %s\n", failures == before ? "PASS" : "FAIL");
}

static void test_free(void)
{
    const unsigned int before = failures;
    rp_kerr_null_trajectory trajectory = {0};
    double state[8];
    size_t component;
    radial_state(0.0, 10.0, state);
    CHECK(solve(0.0, state, 1.0, NULL, 0U, 0.0,
        RP_INTEGRATOR_METHOD_DP45, 1, &trajectory) == RP_INTEGRATOR_STATUS_OK);
    CHECK(trajectory.states != NULL);
    CHECK(trajectory.lambdas != NULL);
    rp_kerr_null_trajectory_free(&trajectory);
    rp_kerr_null_trajectory_free(&trajectory);
    rp_kerr_null_trajectory_free(NULL);
    CHECK(trajectory.states == NULL);
    CHECK(trajectory.lambdas == NULL);
    CHECK(trajectory.count == 0U);
    CHECK(trajectory.capacity == 0U);
    CHECK(trajectory.final_lambda == 0.0);
    CHECK(trajectory.allocation_failed == 0);
    CHECK(trajectory.termination == RP_KERR_NULL_TERMINATION_NONE);
    for (component = 0U; component < 8U; ++component) {
        CHECK(trajectory.final_state[component] == 0.0);
    }
    printf("free/idempotent/NULL: %s\n", failures == before ? "PASS" : "FAIL");
}

static void compare_trajectories(const rp_kerr_null_trajectory *left,
    const rp_kerr_null_trajectory *right)
{
    CHECK(left->integrator_status == right->integrator_status);
    CHECK(left->termination == right->termination);
    CHECK(left->allocation_failed == right->allocation_failed);
    CHECK(left->count == right->count);
    CHECK(left->capacity == right->capacity);
    CHECK(memcmp(left->final_state, right->final_state,
        sizeof(left->final_state)) == 0);
    CHECK(memcmp(&left->final_lambda, &right->final_lambda,
        sizeof(left->final_lambda)) == 0);
    CHECK(left->stats.accepted_steps == right->stats.accepted_steps);
    CHECK(left->stats.rejected_steps == right->stats.rejected_steps);
    CHECK(left->stats.rhs_evaluations == right->stats.rhs_evaluations);
    CHECK(left->stats.jacobian_evaluations == right->stats.jacobian_evaluations);
    CHECK(left->stats.linear_solves == right->stats.linear_solves);
    CHECK(memcmp(&left->stats.final_independent_variable,
        &right->stats.final_independent_variable, sizeof(double)) == 0);
    if (usable_rows(left) && usable_rows(right) && left->count == right->count) {
        CHECK(memcmp(left->states, right->states,
            left->count * 8U * sizeof(double)) == 0);
        CHECK(memcmp(left->lambdas, right->lambdas,
            left->count * sizeof(double)) == 0);
    }
}

static void test_reentrancy(rp_integrator_method method)
{
    const unsigned int before = failures;
    rp_kerr_null_trajectory isolated_a = {0};
    rp_kerr_null_trajectory isolated_b = {0};
    rp_kerr_null_trajectory interleaved_a = {0};
    rp_kerr_null_trajectory interleaved_b = {0};
    double state_a[8];
    double state_b[8];
    radial_state(1.0, 10.0, state_a);
    if (!state_from_constants(0.9, 3.5, 20.0, 1.0,
        3.0, 5.0, 1, state_b)) {
        report_case("reentrancy", method, before);
        return;
    }
    /* The fixed interface is an endpoint operation, not a resumable stepper.
     * Interleave two retained result lifetimes in A,B,B,A call order and
     * compare every active double bit-for-bit against isolated baselines.
     * This checks same-thread determinism and independent owned storage;
     * it does not claim a concurrent or step-level interleaving test. */
    CHECK(solve(0.0, state_a, 5.0, NULL, 0U, 0.0, method, 1, &isolated_a)
        == RP_INTEGRATOR_STATUS_OK);
    CHECK(solve(0.9, state_b, 18.0, NULL, 0U, 0.0, method, 1, &isolated_b)
        == RP_INTEGRATOR_STATUS_OK);
    CHECK(solve(0.9, state_b, 18.0, NULL, 0U, 0.0, method, 1, &interleaved_b)
        == RP_INTEGRATOR_STATUS_OK);
    CHECK(solve(0.0, state_a, 5.0, NULL, 0U, 0.0, method, 1, &interleaved_a)
        == RP_INTEGRATOR_STATUS_OK);
    compare_trajectories(&isolated_a, &interleaved_a);
    compare_trajectories(&isolated_b, &interleaved_b);
    rp_kerr_null_trajectory_free(&interleaved_a);
    /* Releasing one result must not corrupt the independently owned B. */
    compare_trajectories(&isolated_b, &interleaved_b);
    rp_kerr_null_trajectory_free(&interleaved_b);
    rp_kerr_null_trajectory_free(&isolated_a);
    rp_kerr_null_trajectory_free(&isolated_b);
    report_case("reentrancy", method, before);
}

int main(void)
{
    size_t method_index;
    for (method_index = 0U;
        method_index < sizeof(methods) / sizeof(methods[0]); ++method_index) {
        const rp_integrator_method method = methods[method_index];
        test_radial(method);
        test_radial_interpolation(method);
        test_photon_orbit(method);
        test_capture_escape(method);
        test_kerr_invariants(method);
        test_t_final(method);
        test_t_eval(method);
        if (method != RP_INTEGRATOR_METHOD_DP45) {
            test_continuation(method);
        }
        test_large_time_origin(method);
        test_escape_at_final_time(method);
        test_reentrancy(method);
    }
    test_rejections();
    test_free();
    printf("null-solve: %s (%u failed checks)\n",
        failures == 0U ? "PASS" : "FAIL", failures);
    return failures == 0U ? 0 : 1;
}
