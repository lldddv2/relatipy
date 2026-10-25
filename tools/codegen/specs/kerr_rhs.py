"""Specification of ``native/src/geodesic/physic/kerr_rhs_generated.{c,h}``.

Optimized affine right-hand side ``dx^mu = u^mu``, ``du^mu = -Gamma^mu_(alpha
beta) u^alpha u^beta`` of the corrected (Levi--Civita) Kerr connection,
derived symbolically from the Boyer--Lindquist metric.  Each connection
coefficient is simplified with explicit ``Sigma`` and ``Delta`` (so its
numerator has no large-radius cancellation) and then written with the
prepared ``Sigma`` and ``Delta`` factors of the state, so the kernel has the
inputs and the result of ``rp_kerr_geodesic_physic_evaluate_corrected_rhs``.
"""

from __future__ import annotations

import sympy as sp

from relatipy._codegen import optimize
from relatipy._codegen.kerr import kerr_symbols, simplified_geodesic_acceleration

GUARD = "RELATIPY_NATIVE_GEODESIC_PHYSIC_KERR_RHS_GENERATED_H"
FUNCTION = "rp_kerr_geodesic_physic_evaluate_generated_rhs"
PARAMETERS = [
    "double mass",
    "double spin",
    "const struct rp_kerr_geodesic_state *state",
    "double dx[RP_KERR_DIM]",
    "double du[RP_KERR_DIM]",
]


def optimized():
    """The :class:`OptimizedCode` of the generated kernel."""
    s = kerr_symbols()
    acceleration = simplified_geodesic_acceleration(s)
    inputs = {
        s.mass: None,
        s.spin: None,
        s.r: "state->x[1]",
        s.sigma: "state->sigma",
        s.delta: "state->delta",
        **{component: f"state->u[{index}]" for index, component in enumerate(s.u)},
    }
    return optimize(
        {"dx": sp.Matrix(s.u), "du": sp.Matrix(acceleration)},
        inputs=inputs,
        angles={s.theta: ("state->sin_theta", "state->cos_theta")},
        polynomial_in=s.u,
        transforms=("none", "factor_terms", "cancel", "factor"),
    )


def build() -> dict[str, str]:
    code = optimized()
    function = code.c_function(FUNCTION, PARAMETERS)
    header = f"""/**
 * @file kerr_rhs_generated.h
 * @brief Generated corrected Kerr geodesic right-hand side.
 *
 * Same contract as `rp_kerr_geodesic_physic_evaluate_corrected_rhs`: the
 * affine derivatives `dx^mu/dlambda = u^mu` and
 * `du^mu/dlambda = -Gamma^mu_(alpha beta) u^alpha u^beta` of the
 * Levi--Civita connection of the Boyer--Lindquist Kerr metric, signature
 * (-,+,+,+), geometric units with explicit mass and spin length.  The
 * expression is norm-agnostic (timelike and null tangents).
 * `rp_kerr_null_geodesic_rhs`, and the Kerr integration paths built on it,
 * evaluate this kernel; the generic contraction stays as its test
 * reference.  The function allocates nothing, retains no pointer,
 * keeps no global state and does not validate its inputs.
 */

#ifndef {GUARD}
#define {GUARD}

#include "kerr.h"

/**
 * Affine derivatives at a prepared state.
 * @param state State prepared by `rp_kerr_geodesic_physic_prepare_state`
 *        with status OK: `x[1]`, `u`, `sin_theta`, `cos_theta`, `sigma` and
 *        `delta` are read.
 * @param dx Caller-owned output of four doubles, `dx^mu/dlambda`.
 * @param du Caller-owned output of four doubles, `du^mu/dlambda`.  Neither
 *        output may alias the state.
 */
{code.prototype(FUNCTION, PARAMETERS)};

#endif /* {GUARD} */
"""
    source = f"""/**
 * @file kerr_rhs_generated.c
 * @brief Generated corrected Kerr geodesic RHS; see kerr_rhs_generated.h.
 */

#include "kerr_rhs_generated.h"

#include <math.h>

{function}"""
    return {
        "native/src/geodesic/physic/kerr_rhs_generated.h": header,
        "native/src/geodesic/physic/kerr_rhs_generated.c": source,
    }
