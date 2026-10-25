/**
 * @file external.h
 * @brief Black-hole Cartesian frame and projected external four-acceleration.
 *
 * Private boundary of the external-force layer described in
 * `docs/development/external-forces.rst`.  Boyer--Lindquist points
 * `(t, r, theta, phi)` map to the black-hole Cartesian frame
 *
 *     x = sqrt(r^2 + a^2) sin(theta) cos(phi),
 *     y = sqrt(r^2 + a^2) sin(theta) sin(phi),
 *     z = r cos(theta),
 *
 * with `a` the spin length and the black-hole mass set to one.  A Newtonian
 * acceleration `A` evaluated in that frame gives the covector `F_t = 0`,
 * `F_i = A . dx/dx^i` (i = r, theta, phi), and the four-acceleration
 * `f^mu = (g^(mu nu) + u^mu u^nu) F_nu` added to the Kerr geodesic equation.
 *
 * The functions act on a state already prepared and validated by the forced
 * right-hand side (`force/kerr_forced.h`).  Every function allocates
 * nothing, retains no pointer and uses no mutable global state.  This header
 * is not installed and is not a stable ABI.
 */

#ifndef RELATIPY_NATIVE_FORCE_EXTERNAL_H
#define RELATIPY_NATIVE_FORCE_EXTERNAL_H

#include "force/external_force.h"
#include "geodesic/physic/kerr.h"

#include <stddef.h>

/**
 * Status returned when a force callback reports a nonzero value.
 *
 * Forced operations return either an `rp_kerr_status` value or this code,
 * which lies outside the `rp_kerr_status` range.
 */
#define RP_FORCE_STATUS_CALLBACK_FAILURE 32

/** Cartesian position and frame Jacobian at one Boyer--Lindquist point. */
typedef struct {
    /** Black-hole Cartesian position `(x, y, z)`. */
    double position[3];
    /** `jacobian[k][j] = d x^k / d x^(j + 1)`, columns `(r, theta, phi)`. */
    double jacobian[3][3];
} rp_force_frame;

/** Report whether a force is set: non-NULL with a non-NULL callback. */
static inline int rp_force_is_active(const rp_external_force *force)
{
    return force != NULL && force->evaluate != NULL;
}

/**
 * Fill the frame map at a prepared geodesic state.
 *
 * `state` must come from `rp_kerr_geodesic_physic_prepare_state` (its
 * `sin_theta` and `cos_theta` are reused); only `phi` costs two further
 * trigonometric calls.  `spin` is the spin length with the mass set to one.
 * Caller-owned `frame` is fully overwritten; nothing is validated.
 */
void rp_force_frame_prepare(
    double spin,
    const struct rp_kerr_geodesic_state *state,
    rp_force_frame *frame
);

/**
 * Coordinate-time velocity `dx/dt = J (u^r, u^theta, u^phi) / u^t`.
 *
 * `four_velocity` is `(u^t, u^r, u^theta, u^phi)`; the caller owns the three
 * output doubles and checks finiteness (`u^t = 0` gives non-finite values).
 */
void rp_force_frame_velocity(
    const rp_force_frame *frame,
    const double four_velocity[RP_KERR_DIM],
    double velocity[3]
);

/**
 * Evaluate the projected four-acceleration at a prepared geodesic state.
 *
 * Builds the frame, the coordinate-time velocity, calls `force->evaluate`,
 * forms `F_i` and returns `f^mu = g^(mu i) F_i + u^mu (u^i F_i)` with the
 * closed-form Boyer--Lindquist inverse metric (mass one) of the prepared
 * state.  `force` must be active (`rp_force_is_active`).  On success the
 * four caller-owned outputs are finite and the return value is zero.  A
 * nonzero callback value returns `RP_FORCE_STATUS_CALLBACK_FAILURE`; a
 * non-finite velocity, acceleration or result returns
 * `RP_KERR_STATUS_NUMERICAL_RANGE`.  Every failure zeroes the output.
 */
int rp_force_evaluate_prepared(
    double spin,
    const struct rp_kerr_geodesic_state *state,
    const rp_external_force *force,
    double four_acceleration[RP_KERR_DIM]
);

#endif /* RELATIPY_NATIVE_FORCE_EXTERNAL_H */
