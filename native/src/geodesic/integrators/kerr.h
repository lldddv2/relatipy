/**
 * @file kerr.h
 * @brief Private adapter from the experimental Kerr RHS to an ODE callback.
 *
 * The packed (x, u) layout is local glue for tests and native experiments.
 * Production has selected normalized (x, u), but this adapter still accepts
 * G = c = 1 inputs with explicit mass and does not implement that conversion.
 */

#ifndef RELATIPY_GEODESIC_INTEGRATORS_KERR_H
#define RELATIPY_GEODESIC_INTEGRATORS_KERR_H

#include <stddef.h>

typedef struct {
    double mass;
    double spin;
    double E0;
    double Lz0;
    double Q0;
    int projection_ready;
} rp_kerr_integrator_context;

/**
 * Evaluate specific E, Lz and Carter Q from covariant momentum at theta.
 *
 * `momentum` is `u_mu = g_(mu nu) u^nu` in geometric units with explicit spin
 * length. Q uses the fixed rest mass `mu^2 = 1`, not the measured norm, so
 * normalization drift is not hidden. It is not `K = Q + (Lz - a E)^2`.
 * Outputs are written to `constants` as `(E, Lz, Q)`; nothing is validated.
 */
void rp_kerr_timelike_constants_from_momentum(
    double spin,
    double theta,
    const double momentum[4],
    double constants[3]
);

/**
 * Evaluate specific E, Lz and Carter Q of one packed `(x, u)` state.
 *
 * Requires finite input at a point accepted by `rp_kerr_metric`. The
 * four-velocity is not renormalized. `constants` is
 * caller owned and zeroed on failure; nothing is allocated or retained.
 */
int rp_kerr_timelike_constants(
    double mass,
    double spin,
    const double state[8],
    double constants[3]
);

/**
 * Capture the four timelike invariants from a unit-normalized initial state.
 *
 * The caller owns context and state. On failure context is unchanged. The
 * position and contravariant velocity use the packed (x, u) layout.
 */
int rp_kerr_integrator_prepare_projection(
    rp_kerr_integrator_context *context,
    const double state[8]
);

/**
 * Restore norm, energy, axial momentum and Carter Q by orthogonal projection.
 *
 * Off the invariant set, first sets (u^t, u^phi) from E and Lz at fixed
 * position, then moves (r, theta, u^t, u^r, u^theta, u^phi) along J^T, J
 * the Jacobian of the four invariants, with simplified Newton steps and an
 * L D L^T solve of J J^T whose negligible pivots are dropped (Hairer,
 * Lubich & Wanner, Geometric Numerical Integration, sec. IV.4). t and phi
 * never change. A state already on the invariant set to rounding is
 * returned unchanged; a Newton iteration that stagnates below a relative
 * invariant error of 1e-10 is accepted.
 * The caller owns state and context; all eight components remain unchanged
 * on failure, reported as RP_KERR_STATUS_NON_TIMELIKE_VELOCITY when the
 * invariant set is not reached. A candidate at or inside the outer horizon
 * is returned unchanged so the integration observer can handle the
 * terminal event.
 */
int rp_kerr_integrator_project(double tau, double state[], void *context);

int rp_kerr_integrator_rhs(
    double affine_parameter,
    const double state[],
    double derivative[],
    void *context
);

/** Analytic Jacobian adapter for the corrected eight-component affine RHS. */
int rp_kerr_integrator_jacobian(
    double affine_parameter,
    const double state[],
    double *jacobian,
    size_t dimension,
    size_t row_stride,
    void *context
);

/** Adapter for the corrected Kerr null-geodesic right-hand side. */
int rp_kerr_null_integrator_rhs(
    double affine_parameter,
    const double state[],
    double derivative[],
    void *context
);

#endif /* RELATIPY_GEODESIC_INTEGRATORS_KERR_H */
