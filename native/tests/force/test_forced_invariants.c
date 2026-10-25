/*
 * Invariants of forced Kerr orbits.
 *
 * For a static potential Phi axisymmetric about the spin axis,
 * f_t = u_t (u.F) and f_phi = u_phi (u.F) with u.F = -dPhi/dtau, so
 * E' = -u_t exp(Phi) and L' = u_phi exp(Phi) are exact invariants while
 * the plain -u_t and u_phi change.  Every forced orbit keeps u.u = -1 to
 * integrator tolerance, including a tilted, offset, non-axisymmetric model
 * with velocity-dependent drag.  Plummer and Miyamoto--Nagai models are
 * checked at spins 0, 0.5 and 0.9.  The tabulated thin disk holds E' and
 * L' to the accuracy of its separate potential and acceleration
 * interpolants (about 1e-5 apart in the softened plane).
 */

#include "force/kerr_forced.h"
#include "force/models/thin_disk.h"
#include "force_test_models.h"
#include "geodesic/initial/convert.h"
#include "geodesic/integrators/kerr.h"
#include "relatipy/kerr_geometry.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

static const double ELEMENTS[7] = {0.0, 50.0, 0.4, 0.6, 0.3, 1.1, 0.0};

typedef double (*potential_fn)(const void *model, const double position[3]);

typedef struct {
    double norm;
    double energy;
    double angular_momentum;
    double plain_energy;
    double plain_angular_momentum;
    double carter;
} drift;

static double plummer_potential(const void *model, const double position[3])
{
    return test_plummer_potential(model, position);
}

static double disk_potential(const void *model, const double position[3])
{
    return test_miyamoto_nagai_potential(model, position);
}

static double thin_disk_potential(const void *model, const double position[3])
{
    double potential = 0.0;

    assert(rp_thin_disk_potential(position, model, &potential) == 0);
    return potential;
}

static double relative_change(double value, double reference)
{
    return fabs(value - reference) / fabs(reference);
}

/* Observer state: reference values and the largest drift so far. */
typedef struct {
    double spin;
    const rp_external_force *force;
    potential_fn potential;
    size_t count;
    double reference[5];
    drift result;
} drift_observer;

/* Black-hole Cartesian position of a Boyer--Lindquist point. */
static void cartesian_position(
    double spin,
    const double coordinates[4],
    double position[3]
)
{
    struct rp_kerr_geodesic_state prepared;
    rp_force_frame frame;
    size_t k;

    for (k = 0U; k < RP_KERR_DIM; ++k) {
        prepared.x[k] = coordinates[k];
        prepared.u[k] = 0.0;
    }
    assert(rp_kerr_geodesic_physic_prepare_state(1.0, spin, &prepared)
        == RP_KERR_STATUS_OK);
    rp_force_frame_prepare(spin, &prepared, &frame);
    for (k = 0U; k < 3U; ++k) {
        position[k] = frame.position[k];
    }
}

/* Record the drift of each quantity at one accepted state. */
static int observe_drift(double tau, const double state[], void *opaque)
{
    drift_observer *observer = opaque;
    double metric[RP_KERR_DIM][RP_KERR_DIM];
    double lower[RP_KERR_DIM] = {0.0, 0.0, 0.0, 0.0};
    double constants[3];
    double position[3];
    double norm = 0.0;
    double weight = 1.0;
    double values[5];
    drift *result = &observer->result;
    size_t mu;
    size_t nu;

    (void)tau;
    assert(rp_kerr_metric(1.0, observer->spin, state, metric)
        == RP_KERR_STATUS_OK);
    for (mu = 0U; mu < 4U; ++mu) {
        for (nu = 0U; nu < 4U; ++nu) {
            lower[mu] += metric[mu][nu] * state[4U + nu];
        }
        norm += lower[mu] * state[4U + mu];
    }
    if (observer->potential != NULL) {
        cartesian_position(observer->spin, state, position);
        weight = exp(observer->potential(observer->force->model, position));
    }
    rp_kerr_timelike_constants_from_momentum(
        observer->spin, state[2], lower, constants
    );
    values[0] = -lower[0] * weight;
    values[1] = lower[3] * weight;
    values[2] = -lower[0];
    values[3] = lower[3];
    values[4] = constants[2];
    if (observer->count == 0U) {
        for (mu = 0U; mu < 5U; ++mu) {
            observer->reference[mu] = values[mu];
        }
    }
    ++observer->count;
    result->norm = fmax(result->norm, fabs(norm + 1.0));
    result->energy = fmax(
        result->energy, relative_change(values[0], observer->reference[0])
    );
    result->angular_momentum = fmax(
        result->angular_momentum,
        relative_change(values[1], observer->reference[1])
    );
    result->plain_energy = fmax(
        result->plain_energy, relative_change(values[2], observer->reference[2])
    );
    result->plain_angular_momentum = fmax(
        result->plain_angular_momentum,
        relative_change(values[3], observer->reference[3])
    );
    result->carter = fmax(
        result->carter, relative_change(values[4], observer->reference[4])
    );
    return 0;
}

/* Integrate one orbit and record the largest drift of each quantity. */
static drift run_orbit(
    const double elements[7],
    double spin,
    const rp_external_force *force,
    potential_fn potential,
    rp_integrator_method method,
    double rtol,
    double atol,
    double tau_final
)
{
    rp_integrator_config config = {
        RP_INTEGRATOR_METHOD_DOP853, 8U, 0.0, 0.0, 0.0, 0.0, 1000000U,
        NULL, NULL, NULL, NULL, NULL
    };
    rp_kerr_forced_context context;
    drift_observer observer = {
        0.0, NULL, NULL, 0U, {0.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
    };
    double state[8];

    assert(rp_initial_from_elements(spin, elements, state)
        == RP_KERR_STATUS_OK);
    observer.spin = spin;
    observer.force = force;
    observer.potential = potential;
    observe_drift(0.0, state, &observer);
    context.spin = spin;
    context.force = force;
    config.method = method;
    config.relative_tolerance = rtol;
    config.absolute_tolerance = atol;
    config.step_observer = observe_drift;
    config.step_observer_context = &observer;
    assert(rp_integrator_integrate(
        &config, rp_kerr_forced_integrator_rhs, &context,
        0.0, tau_final, state, NULL
    ) == RP_INTEGRATOR_STATUS_OK);
    assert(observer.count > 50U);
    return observer.result;
}

static void check_axisymmetric(void)
{
    const double spins[3] = {0.0, 0.5, 0.9};
    const test_plummer plummer = {0.05, 30.0, {0.0, 0.0, 0.0}};
    const test_miyamoto_nagai disk = {0.05, 20.0, 3.0};
    const rp_external_force forces[2] = {
        {test_plummer_acceleration, &plummer},
        {test_miyamoto_nagai_acceleration, &disk}
    };
    const potential_fn potentials[2] = {plummer_potential, disk_potential};
    const char *names[2] = {"plummer", "miyamoto-nagai"};
    size_t s;
    size_t m;

    for (m = 0U; m < 2U; ++m) {
        for (s = 0U; s < 3U; ++s) {
            const drift d = run_orbit(
                ELEMENTS, spins[s], &forces[m], potentials[m],
                RP_INTEGRATOR_METHOD_DOP853, 1e-12, 1e-14, 6000.0
            );
            printf("%-15s a=%.1f dop853 |u.u+1|=%.1e dE'=%.1e dL'=%.1e"
                " (plain dE=%.1e dL=%.1e dQ=%.1e)\n",
                names[m], spins[s], d.norm, d.energy, d.angular_momentum,
                d.plain_energy, d.plain_angular_momentum, d.carter);
            assert(d.norm < 1e-10);
            assert(d.energy < 1e-10);
            assert(d.angular_momentum < 1e-10);
            /* The invariants are not trivially the homogeneous ones. */
            assert(d.plain_energy > 1e3 * d.energy);
            assert(d.plain_angular_momentum > 1e3 * d.angular_momentum);
        }
    }

    /* The other supported methods, radau with finite differences. */
    {
        const drift d45 = run_orbit(
            ELEMENTS, 0.9, &forces[1], potentials[1], RP_INTEGRATOR_METHOD_DP45,
            1e-11, 1e-13, 6000.0
        );
        const drift dradau = run_orbit(
            ELEMENTS, 0.9, &forces[1], potentials[1], RP_INTEGRATOR_METHOD_RADAU,
            1e-10, 1e-12, 3000.0
        );
        printf("miyamoto-nagai  a=0.9 dp45   |u.u+1|=%.1e dE'=%.1e dL'=%.1e\n",
            d45.norm, d45.energy, d45.angular_momentum);
        printf("miyamoto-nagai  a=0.9 radau  |u.u+1|=%.1e dE'=%.1e dL'=%.1e\n",
            dradau.norm, dradau.energy, dradau.angular_momentum);
        assert(d45.norm < 1e-8 && d45.energy < 1e-8
            && d45.angular_momentum < 1e-8);
        assert(dradau.norm < 1e-7 && dradau.energy < 1e-7
            && dradau.angular_momentum < 1e-7);
    }
}

static void check_tilted(void)
{
    const double axis[3] = {0.6, -0.48, 0.64};
    test_tilted model = {
        {0.05, 20.0, 3.0},
        {{0.0}},
        {0.02, 8.0, {12.0, -7.0, 5.0}},
        2e-6
    };
    const rp_external_force force = {test_tilted_acceleration, &model};
    const double spins[2] = {0.0, 0.9};
    size_t s;

    test_axis_rotation(axis, 0.7, model.rotation);
    for (s = 0U; s < 2U; ++s) {
        const drift d = run_orbit(
            ELEMENTS, spins[s], &force, NULL, RP_INTEGRATOR_METHOD_DOP853,
            1e-12, 1e-14, 6000.0
        );
        printf("tilted+drag     a=%.1f dop853 |u.u+1|=%.1e"
            " (dE=%.1e dL=%.1e dQ=%.1e)\n",
            spins[s], d.norm, d.plain_energy, d.plain_angular_momentum,
            d.carter);
        assert(d.norm < 1e-10);
        assert(d.plain_energy > 1e-6 && d.plain_angular_momentum > 1e-6);
    }
}

static void check_thin_disk(void)
{
    /* Disk of 1e-2 black-hole masses out to 1000 M; orbit a = 300 M. */
    const double elements[7] = {0.0, 300.0, 0.5, 0.6, 0.4, 0.0, 0.0};
    const double normal[3] = {0.0, 0.0, 1.0};
    rp_thin_disk_params params;
    rp_axisymmetric_table *table = NULL;
    rp_thin_disk_model model;
    rp_external_force force;
    drift d;

    rp_thin_disk_default_params(&params);
    assert(rp_thin_disk_table_build(&params, &table)
        == RP_AXISYMMETRIC_STATUS_OK);
    assert(rp_thin_disk_model_init(&model, table, 1e-2, 1000.0, normal)
        == RP_AXISYMMETRIC_STATUS_OK);
    force.evaluate = rp_thin_disk_acceleration;
    force.model = &model;
    d = run_orbit(
        elements, 0.7, &force, thin_disk_potential,
        RP_INTEGRATOR_METHOD_DOP853, 1e-12, 1e-14, 2e5
    );
    printf("thin disk       a=0.7 dop853 |u.u+1|=%.1e dE'=%.1e dL'=%.1e"
        " (plain dE=%.1e dL=%.1e)\n",
        d.norm, d.energy, d.angular_momentum, d.plain_energy,
        d.plain_angular_momentum);
    assert(d.norm < 1e-10);
    assert(d.energy < 1e-6 && d.angular_momentum < 1e-6);
    assert(d.energy < 1e-2 * d.plain_energy);
    assert(d.angular_momentum < 1e-2 * d.plain_energy);
    rp_thin_disk_table_free(table);
}

static void check_callback_failure(void)
{
    const rp_external_force failing = {test_failing_acceleration, NULL};
    const double elements[7] = {0.0, 50.0, 0.4, 0.6, 0.3, 1.1, 0.0};
    const rp_integrator_config config = {
        RP_INTEGRATOR_METHOD_DOP853, 8U, 1e-10, 1e-12, 0.0, 0.0, 1000000U,
        NULL, NULL, NULL, NULL, NULL
    };
    rp_kerr_forced_context context;
    rp_integrator_stats stats;
    double initial[8];
    double state[8];
    size_t k;

    assert(rp_initial_from_elements(0.5, elements, initial)
        == RP_KERR_STATUS_OK);
    for (k = 0U; k < 8U; ++k) {
        state[k] = initial[k];
    }
    context.spin = 0.5;
    context.force = &failing;
    assert(rp_integrator_integrate(
        &config, rp_kerr_forced_integrator_rhs, &context,
        0.0, 100.0, state, &stats
    ) == RP_INTEGRATOR_STATUS_RHS_FAILURE);
    assert(stats.accepted_steps == 0U);
    for (k = 0U; k < 8U; ++k) {
        assert(state[k] == initial[k]);
    }
}

int main(void)
{
    check_axisymmetric();
    check_tilted();
    check_thin_disk();
    check_callback_failure();
    puts("test_forced_invariants: ok");
    return 0;
}
