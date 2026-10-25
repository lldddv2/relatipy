"""Specification of ``native/src/force/generated/projection.{c,h}``.

Projects a Newtonian acceleration ``a`` given in the black-hole Cartesian
frame

    x = sqrt(r^2 + a_*^2) sin(theta) cos(phi),
    y = sqrt(r^2 + a_*^2) sin(theta) sin(phi),
    z = r cos(theta)

onto the Kerr four-force ``f^mu = (g^(mu nu) + u^mu u^nu) F_nu`` with
``F_t = 0`` and ``F_i = a . dx/dx^i`` for ``i = r, theta, phi``
(``docs/development/external-forces.rst``).  The inverse metric uses the
prepared ``Sigma`` and ``Delta`` of the Boyer--Lindquist state.
"""

from __future__ import annotations

import sympy as sp

from relatipy._codegen import optimize
from relatipy._codegen.kerr import inverse_metric, kerr_symbols, metric

GUARD = "RELATIPY_NATIVE_FORCE_GENERATED_PROJECTION_H"
PARAMETERS = [
    "double mass",
    "double spin",
    "const struct rp_kerr_geodesic_state *state",
    "double sin_phi",
    "double cos_phi",
    "const double acceleration[3]",
    "double force[4]",
]


def symbols():
    s = kerr_symbols()
    phi = sp.Symbol("phi", real=True)
    acceleration = sp.symbols("a_x a_y a_z", real=True)
    return s, phi, acceleration


def frame_jacobian(s, phi):
    """``dx_j / dq_i`` of the Cartesian frame map, ``q = (r, theta, phi)``."""
    spheroidal = sp.sqrt(s.r**2 + s.spin**2)
    position = sp.Matrix([
        spheroidal * sp.sin(s.theta) * sp.cos(phi),
        spheroidal * sp.sin(s.theta) * sp.sin(phi),
        s.r * sp.cos(s.theta),
    ])
    return position.jacobian([s.r, s.theta, phi])


def covector(s, phi, acceleration):
    """``(F_t, F_r, F_theta, F_phi)`` with ``F_t = 0``."""
    jacobian = frame_jacobian(s, phi)
    spatial = [sum(acceleration[j] * jacobian[j, i] for j in range(3)) for i in range(3)]
    return [sp.S.Zero, *spatial]


def expressions():
    """``f^mu`` in terms of the prepared state, ``sin/cos phi`` and ``a``."""
    s, phi, acceleration = symbols()
    force_covector = covector(s, phi, acceleration)
    inverse = inverse_metric(metric(s, s.sigma, s.delta), s, s.delta)
    u = s.u
    contraction = sum(u[nu] * force_covector[nu] for nu in range(4))
    return [
        sum(inverse[mu, nu] * force_covector[nu] for nu in range(4)) + u[mu] * contraction
        for mu in range(4)
    ]


def build() -> dict[str, str]:
    s, phi, acceleration = symbols()
    inputs = {
        s.mass: None,
        s.spin: None,
        s.r: "state->x[1]",
        s.sigma: "state->sigma",
        s.delta: "state->delta",
        **{component: f"state->u[{index}]" for index, component in enumerate(s.u)},
        **{component: f"acceleration[{index}]" for index, component in enumerate(acceleration)},
    }
    code = optimize(
        sp.Matrix(expressions()),
        inputs=inputs,
        outputs="force",
        angles={
            s.theta: ("state->sin_theta", "state->cos_theta"),
            phi: ("sin_phi", "cos_phi"),
        },
        transforms=("none", "factor_terms", "cancel", "factor"),
    )
    function = code.c_function("rp_force_generated_projection", PARAMETERS)
    header = f"""/**
 * @file projection.h
 * @brief Generated projection of a Cartesian acceleration onto a Kerr force.
 *
 * Evaluates `f^mu = (g^(mu nu) + u^mu u^nu) F_nu` with `F_t = 0` and
 * `F_i = a . dx/dx^i` for the black-hole Cartesian frame of
 * `docs/development/external-forces.rst`.  Boyer--Lindquist coordinates,
 * signature (-,+,+,+), geometric units with explicit mass and spin length.
 * The function allocates nothing, retains no pointer, keeps no global state
 * and does not validate its inputs.
 */

#ifndef {GUARD}
#define {GUARD}

#include "../../geodesic/physic/kerr.h"

/**
 * Four-force of a Newtonian acceleration at a prepared state.
 * @param mass Black-hole mass (one in the normalized production units).
 * @param spin Spin length `a_*`.
 * @param state State prepared by `rp_kerr_geodesic_physic_prepare_state`:
 *        `x[1]`, `u`, `sin_theta`, `cos_theta`, `sigma` and `delta` are read.
 * @param sin_phi Sine of the azimuth `state->x[3]`.
 * @param cos_phi Cosine of the same azimuth.
 * @param acceleration Cartesian acceleration `(a_x, a_y, a_z)`, three doubles.
 * @param force Caller-owned output of four doubles, `f^mu`; it must not
 *        alias `acceleration` or the state.
 */
{code.prototype("rp_force_generated_projection", PARAMETERS)};

#endif /* {GUARD} */
"""
    source = f"""/**
 * @file projection.c
 * @brief Generated Kerr projection of a Cartesian acceleration; see projection.h.
 */

#include "projection.h"

#include <math.h>

{function}"""
    return {
        "native/src/force/generated/projection.h": header,
        "native/src/force/generated/projection.c": source,
    }
