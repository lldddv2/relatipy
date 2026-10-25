/*
 * rp_kerr_null_geodesic_rhs evaluates the generated corrected Kerr kernel,
 * rp_kerr_geodesic_physic_evaluate_generated_rhs, bit for bit.  Compare it
 * with the generic contraction of analytic metric derivatives,
 * rp_kerr_geodesic_physic_evaluate_corrected_rhs, kept as the reference
 * oracle, on a wide grid: spins -0.99..0.99, radii from just outside the
 * outer horizon to 1e6 M, polar angles from 1e-6 off either pole to the
 * equator, and timelike and null tangents built in the zero-angular-momentum
 * frame, so that they are valid inside the ergoregion too.
 *
 * Each derivative is measured in that frame's orthonormal scales,
 * w = (alpha, sqrt(Sigma/Delta), sqrt(Sigma), sqrt(g_phiphi)).  The largest
 * weighted component difference is compared with the largest weighted size
 * of the terms -Gamma^mu_(alpha beta) u^alpha u^beta summed by the
 * contraction, which bounds its rounding error even where the terms cancel
 * (a circular photon orbit has du = 0).  The difference must stay within
 * 4 eps (r/M + (r^2 + a^2)/Delta) of that size; the measured maximum is
 * at most about 0.16 of this bound at every radius of the grid.  Near the horizon
 * both kernels inherit the rounding of the prepared Delta.  At large r the
 * difference is the error of the oracle: against a 50-digit evaluation at
 * r = 1e6 M the generated timelike derivative is accurate to 5e-15 per
 * component, the contraction only to 3e-9.
 */

#include "relatipy/kerr_geodesic.h"
#include "relatipy/kerr_geometry.h"
#include "geodesic/physic/kerr.h"
#include "geodesic/physic/kerr_rhs_generated.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define TOLERANCE_FACTOR 4.0

static const double PI = 3.14159265358979323846;

static double worst_error;
static double worst_point[4];
static size_t evaluated;

/* Future-directed tangent with frame velocity `speed * direction`. */
static void zamo_tangent(
    double mass,
    double spin,
    double r,
    double theta,
    const double direction[3],
    double speed,
    int null,
    double tangent[4],
    double weights[4]
)
{
    const double s = sin(theta);
    const double c = cos(theta);
    const double sigma = r * r + spin * spin * c * c;
    const double delta = r * r - 2.0 * mass * r + spin * spin;
    const double big_a = (r * r + spin * spin) * (r * r + spin * spin)
        - spin * spin * delta * s * s;
    const double alpha = sqrt(sigma * delta / big_a);
    const double omega = 2.0 * mass * spin * r / big_a;
    const double root_g_phiphi = sqrt(big_a / sigma) * s;
    double norm;
    double v[3];
    double gamma;
    size_t i;

    norm = sqrt(direction[0] * direction[0] + direction[1] * direction[1]
        + direction[2] * direction[2]);
    for (i = 0U; i < 3U; ++i) {
        v[i] = (null ? 1.0 : speed) * direction[i] / norm;
    }
    gamma = null ? 1.0 : 1.0 / sqrt(1.0 - speed * speed);
    tangent[0] = gamma / alpha;
    tangent[1] = gamma * v[0] * sqrt(delta / sigma);
    tangent[2] = gamma * v[1] / sqrt(sigma);
    tangent[3] = gamma * (omega / alpha + v[2] / root_g_phiphi);
    weights[0] = alpha;
    weights[1] = sqrt(sigma / delta);
    weights[2] = sqrt(sigma);
    weights[3] = root_g_phiphi;
}

/* Normalization residual relative to the largest metric term. */
static double norm_residual(
    double mass,
    double spin,
    const double x[4],
    const double u[4],
    double expected
)
{
    double metric[4][4];
    double sum = 0.0;
    double scale = 0.0;
    size_t i;
    size_t j;

    assert(rp_kerr_metric(mass, spin, x, metric) == RP_KERR_STATUS_OK);
    for (i = 0U; i < 4U; ++i) {
        for (j = 0U; j < 4U; ++j) {
            const double term = metric[i][j] * u[i] * u[j];
            sum += term;
            scale = fmax(scale, fabs(term));
        }
    }
    return fabs(sum - expected) / scale;
}

/*
 * Sum over alpha <= beta of |Gamma^k_(alpha beta) u^alpha u^beta| (times
 * the symmetry factor), from the oracle on one- and two-component tangents
 * by polarization.  `state` is prepared; its tangent is overwritten.
 */
static void term_sizes(
    double mass,
    double spin,
    struct rp_kerr_geodesic_state *state,
    const double u[4],
    double sizes[4]
)
{
    double diagonal[4][4];
    double pair[4];
    double dx[4];
    size_t a;
    size_t b;
    size_t k;

    for (k = 0U; k < 4U; ++k) {
        sizes[k] = 0.0;
    }
    for (a = 0U; a < 4U; ++a) {
        memset(state->u, 0, sizeof(state->u));
        state->u[a] = u[a];
        rp_kerr_geodesic_physic_evaluate_corrected_rhs(
            mass, spin, state, dx, diagonal[a]
        );
        for (k = 0U; k < 4U; ++k) {
            sizes[k] += fabs(diagonal[a][k]);
        }
    }
    for (a = 0U; a < 4U; ++a) {
        for (b = a + 1U; b < 4U; ++b) {
            memset(state->u, 0, sizeof(state->u));
            state->u[a] = u[a];
            state->u[b] = u[b];
            rp_kerr_geodesic_physic_evaluate_corrected_rhs(
                mass, spin, state, dx, pair
            );
            for (k = 0U; k < 4U; ++k) {
                sizes[k] += fabs(pair[k] - diagonal[a][k] - diagonal[b][k]);
            }
        }
    }
}

static void check_point(
    double mass,
    double spin,
    double r,
    double theta,
    const double direction[3],
    double speed,
    int null
)
{
    const double x[4] = {12.5, r, theta, 0.7};
    struct rp_kerr_geodesic_state state;
    double u[4];
    double w[4];
    double dx[4];
    double du[4];
    double dx_ref[4];
    double du_ref[4];
    double dx_gen[4];
    double du_gen[4];
    double sizes[4];
    double difference = 0.0;
    double scale = 0.0;
    double bound;
    double error;
    size_t k;

    zamo_tangent(mass, spin, r, theta, direction, speed, null, u, w);
    assert(norm_residual(mass, spin, x, u, null ? 0.0 : -1.0) <= 1e-9);

    assert(rp_kerr_null_geodesic_rhs(mass, spin, x, u, dx, du)
        == RP_KERR_STATUS_OK);

    memcpy(state.x, x, sizeof(state.x));
    memcpy(state.u, u, sizeof(state.u));
    assert(rp_kerr_geodesic_physic_prepare_state(mass, spin, &state)
        == RP_KERR_STATUS_OK);
    rp_kerr_geodesic_physic_evaluate_generated_rhs(
        mass, spin, &state, dx_gen, du_gen
    );
    assert(memcmp(dx, dx_gen, sizeof(dx)) == 0);
    assert(memcmp(du, du_gen, sizeof(du)) == 0);
    rp_kerr_geodesic_physic_evaluate_corrected_rhs(
        mass, spin, &state, dx_ref, du_ref
    );
    term_sizes(mass, spin, &state, u, sizes);

    for (k = 0U; k < 4U; ++k) {
        assert(dx[k] == u[k] && dx_ref[k] == u[k]);
        assert(isfinite(du[k]));
        difference = fmax(difference, w[k] * fabs(du[k] - du_ref[k]));
        scale = fmax(scale, w[k] * sizes[k]);
    }
    assert(scale > 0.0);
    bound = TOLERANCE_FACTOR * DBL_EPSILON * (r / mass
        + (r * r + spin * spin) / (r * r - 2.0 * mass * r + spin * spin));
    error = difference / scale / bound;
    if (error > worst_error) {
        worst_error = error;
        worst_point[0] = spin / mass;
        worst_point[1] = r / mass;
        worst_point[2] = theta;
        worst_point[3] = null ? -1.0 : speed;
    }
    assert(error <= 1.0);
    ++evaluated;
}

static void check_grid(double mass)
{
    const double spins[10] = {
        -0.99, -0.9, -0.6, -0.3, 0.0, 0.1, 0.5, 0.8, 0.95, 0.99
    };
    const double horizon_factors[4] = {1.0 + 1e-6, 1.0 + 1e-3, 1.1, 1.5};
    const double radii[9] = {3.0, 6.0, 10.0, 30.0, 1e2, 1e3, 1e4, 1e5, 1e6};
    const double thetas[12] = {
        1e-6, 1e-4, 1e-2, 0.2, 0.6, 1.0, 1.4,
        0.5 * PI - 1e-6, 0.5 * PI, PI - 0.6, PI - 1e-2, PI - 1e-6
    };
    const double directions[7][3] = {
        {1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}, {0.0, 0.0, -1.0}, {0.6, 0.48, 0.64},
        {-0.3, -0.5, 0.81}
    };
    const double speeds[2] = {0.3, 0.95};
    const double at_rest[3] = {1.0, 0.0, 0.0};
    size_t s;
    size_t i;
    size_t j;
    size_t d;
    size_t v;

    for (s = 0U; s < 10U; ++s) {
        const double spin = spins[s] * mass;
        const double outer = mass + sqrt(mass * mass - spin * spin);
        for (i = 0U; i < 13U; ++i) {
            const double r = i < 4U
                ? horizon_factors[i] * outer : radii[i - 4U] * mass;
            for (j = 0U; j < 12U; ++j) {
                check_point(mass, spin, r, thetas[j], at_rest, 0.0, 0);
                for (d = 0U; d < 7U; ++d) {
                    for (v = 0U; v < 2U; ++v) {
                        check_point(
                            mass, spin, r, thetas[j], directions[d],
                            speeds[v], 0
                        );
                    }
                    check_point(
                        mass, spin, r, thetas[j], directions[d], 1.0, 1
                    );
                }
            }
        }
    }
}

int main(void)
{
    check_grid(1.0);
    check_grid(3.7);
    assert(evaluated == 2U * 10U * 13U * 12U * (1U + 7U * 3U));
    printf(
        "%zu tangents, worst difference %.2f of the bound "
        "(a/M %.2f, r/M %.6g, theta %.3g, speed %g, -1 if null)\n",
        evaluated, worst_error, worst_point[0], worst_point[1],
        worst_point[2], worst_point[3]
    );
    printf("test_kerr_generated_rhs: ok\n");
    return 0;
}
