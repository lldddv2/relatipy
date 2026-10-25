"""Specification of ``native/src/force/generated/ring_field.{c,h}``.

Field of a thin ring of mass ``dM`` and radius ``R`` in its own cylindrical
frame, at cylindrical radius ``d > 0`` and height ``z``, with the softening
rule ``z^2 -> z^2 + H^2``:

    p = sqrt((d + R)^2 + z^2 + H^2),   q^2 = (d - R)^2 + z^2 + H^2,
    m = 4 R d / p^2,
    a_d = -(dM / (pi p d)) [K + E (d^2 - R^2 - z^2 - H^2) / q^2],
    a_z = -(2 dM / (pi p q^2)) E z,
    Phi = -2 dM K / (pi p),

with ``K = K(m)`` and ``E = E(m)`` the complete elliptic integrals of
parameter ``m``, supplied by the caller.
"""

from __future__ import annotations

import sympy as sp

from relatipy._codegen import optimize

GUARD = "RELATIPY_NATIVE_FORCE_GENERATED_RING_FIELD_H"


def symbols():
    d, ring_radius, softening, ring_mass = sp.symbols(
        "d ring_radius softening ring_mass", positive=True
    )
    z = sp.Symbol("z", real=True)
    elliptic_k, elliptic_e = sp.symbols("elliptic_k elliptic_e", positive=True)
    return d, z, ring_radius, softening, ring_mass, elliptic_k, elliptic_e


def expressions():
    """``(m, a_d, a_z, Phi)`` as SymPy expressions of :func:`symbols`."""
    d, z, R, H, dM, K, E = symbols()
    z2 = z**2 + H**2
    p = sp.sqrt((d + R) ** 2 + z2)
    q2 = (d - R) ** 2 + z2
    parameter = 4 * R * d / p**2
    radial = -(dM / (sp.pi * p * d)) * (K + E * (d**2 - R**2 - z2) / q2)
    vertical = -(2 * dM / (sp.pi * p * q2)) * E * z
    potential = -2 * dM * K / (sp.pi * p)
    return parameter, radial, vertical, potential


def build() -> dict[str, str]:
    d, z, R, H, dM, K, E = symbols()
    parameter, radial, vertical, potential = expressions()
    geometry = (d, z, R, H)
    parameter_code = optimize(
        parameter, inputs=geometry, outputs="parameter", transforms=("none", "factor_terms"),
    )
    field_code = optimize(
        [radial, vertical, potential], inputs=(*geometry, dM, K, E), outputs="field",
    )
    geometry_parameters = [
        "double d", "double z", "double ring_radius", "double softening",
    ]
    parameter_function = parameter_code.c_function(
        "rp_force_generated_ring_parameter",
        geometry_parameters,
        return_type="double",
        prologue=("double parameter;",),
        epilogue=("return parameter;",),
    )
    field_function = field_code.c_function(
        "rp_force_generated_ring_field",
        [*geometry_parameters, "double ring_mass", "double elliptic_k", "double elliptic_e",
         "double field[3]"],
    )
    header = f"""/**
 * @file ring_field.h
 * @brief Softened field of one thin ring from precomputed elliptic integrals.
 *
 * The ring has mass `ring_mass` and radius `ring_radius` and lies in the
 * plane `z = 0` of its own frame.  `d > 0` is the cylindrical radius of the
 * field point and `z` its height; the softening length `softening = H`
 * replaces `z^2` by `z^2 + H^2` in every distance.  Units are geometric with
 * the black-hole mass set to one.  Functions allocate nothing, retain no
 * pointer and keep no global state; they do not validate their inputs.
 */

#ifndef {GUARD}
#define {GUARD}

/**
 * Parameter `m = 4 R d / ((d + R)^2 + z^2 + H^2)` of the complete elliptic
 * integrals `K(m)` and `E(m)` (parameter, not modulus: `m = k^2`).
 * @return `m`, in `[0, 1)` for `d > 0` and `H > 0`.
 */
{parameter_code.prototype("rp_force_generated_ring_parameter", geometry_parameters, return_type="double")};

/**
 * Field of one ring at `(d, z)`.
 * @param elliptic_k `K(m)` with `m` from `rp_force_generated_ring_parameter`.
 * @param elliptic_e `E(m)` at the same `m`.
 * @param field Caller-owned output of three doubles: `field[0]` is the
 *        cylindrical-radial acceleration `a_d`, `field[1]` the vertical
 *        acceleration `a_z` and `field[2]` the potential `Phi`.
 */
{field_code.prototype("rp_force_generated_ring_field", [*geometry_parameters, "double ring_mass", "double elliptic_k", "double elliptic_e", "double field[3]"])};

#endif /* {GUARD} */
"""
    source = f"""/**
 * @file ring_field.c
 * @brief Generated softened thin-ring field; see ring_field.h.
 */

#include "ring_field.h"

#include <math.h>

{parameter_function}
{field_function}"""
    return {
        "native/src/force/generated/ring_field.h": header,
        "native/src/force/generated/ring_field.c": source,
    }
