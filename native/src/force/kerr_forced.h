/**
 * @file kerr_forced.h
 * @brief Forced Kerr timelike right-hand sides and their ODE adapters.
 *
 * The forced proper-time derivatives are
 *
 *     dx^mu/dtau = u^mu,
 *     du^mu/dtau = -Gamma^mu_(alpha beta) u^alpha u^beta + f^mu,
 *
 * where the first term is the corrected Kerr Levi--Civita contraction used
 * by every homogeneous production path (`rp_kerr_null_geodesic_rhs`) and
 * `f^mu` is the projected external four-acceleration of `external.h`.  The
 * mass is one and `spin` is the spin length.  An inactive force (NULL, or a
 * NULL callback) reproduces the homogeneous derivatives bit for bit.
 *
 * Nothing here allocates, retains a pointer or uses mutable global state.
 * The header is private to the native layer and is not a stable ABI.
 */

#ifndef RELATIPY_NATIVE_FORCE_KERR_FORCED_H
#define RELATIPY_NATIVE_FORCE_KERR_FORCED_H

#include "force/external.h"
#include "geodesic/integrators/integrator.h"

#include <stddef.h>

/**
 * ODE context of the forced adapters: spin length (mass one) and a
 * caller-owned, read-only force that may be NULL.  The caller keeps both
 * alive and unchanged during an integration.
 */
typedef struct {
    double spin;
    const rp_external_force *force;
} rp_kerr_forced_context;

/**
 * Evaluate the forced proper-time derivatives at one state.
 *
 * Validation, aliasing and output rules are those of
 * `rp_kerr_null_geodesic_rhs` with mass one: `|spin| <= 1`, finite inputs,
 * the two four-component caller-owned outputs must not overlap each other
 * but may overlap the inputs, and every validated output is zeroed on
 * failure.  Returns zero on success, an `rp_kerr_status` value, or
 * `RP_FORCE_STATUS_CALLBACK_FAILURE` when the force callback fails.
 */
int rp_kerr_forced_rhs(
    double spin,
    const rp_external_force *force,
    const double coordinates[RP_KERR_DIM],
    const double four_velocity[RP_KERR_DIM],
    double coordinate_derivative[RP_KERR_DIM],
    double four_velocity_derivative[RP_KERR_DIM]
);

/**
 * `rp_integrator_rhs` adapter in proper time for the packed eight-component
 * state `(x^mu, u^mu)`.  `context` points to `rp_kerr_forced_context`.
 * Returns zero or the nonzero status of `rp_kerr_forced_rhs`; a NULL
 * context returns `RP_KERR_STATUS_NULL_POINTER`.
 */
int rp_kerr_forced_integrator_rhs(
    double proper_time,
    const double state[],
    double derivative[],
    void *context
);

/**
 * `rp_integrator_rhs` adapter in coordinate time for the same packed state:
 * every proper-time derivative, force term included, is divided by `u^t`,
 * exactly as the coordinate-time RHS of `geodesic/kerr_mcmc.c`.  Requires a
 * finite `u^t > 0`, otherwise returns `RP_KERR_STATUS_NON_TIMELIKE_VELOCITY`;
 * non-finite results return `RP_KERR_STATUS_NUMERICAL_RANGE`.
 */
int rp_kerr_forced_coordinate_time_rhs(
    double coordinate_time,
    const double state[],
    double derivative[],
    void *context
);

/**
 * Analytic `rp_integrator_jacobian` of the homogeneous term only, for an
 * `rp_kerr_forced_context`.  Valid for an inactive force (it equals
 * `rp_kerr_integrator_jacobian` with mass one); with an active force it
 * returns `RP_KERR_STATUS_INVALID_PARAMETER`, and callers leave the config
 * Jacobian NULL so that Radau uses finite differences.
 */
int rp_kerr_forced_integrator_jacobian(
    double proper_time,
    const double state[],
    double *jacobian,
    size_t dimension,
    size_t row_stride,
    void *context
);

/**
 * Report whether the forced layer accepts an integration method.
 *
 * `dop853`, `dp45` and `radau` return `RP_INTEGRATOR_STATUS_OK`.
 * `projection_radau` and unknown methods return
 * `RP_INTEGRATOR_STATUS_INVALID_ARGUMENT`: its projector restores Carter's
 * constant, which an external force does not conserve.
 */
rp_integrator_status rp_kerr_forced_check_method(rp_integrator_method method);

#endif /* RELATIPY_NATIVE_FORCE_KERR_FORCED_H */
