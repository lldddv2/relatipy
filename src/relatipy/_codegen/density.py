"""C source of a surface density ``Sigma(R)`` for the native ring sum.

:func:`build_surface_density_source` accepts ``Sigma`` as a string in the
coordinates of the disk plane and returns one C translation unit that
implements the ``rp_surface_density_fn`` of
``native/src/utils/integrators/axisymmetric/ring_sum.h``::

    double name(double radius, const void *context);

``radius`` is the cylindrical radius in units of the outer radius of the
disk; declared parameters are read from ``context`` as ``const double *``.
The density must be axisymmetric: written with ``R`` (or ``r``, equal to
``R`` in the plane), or with ``x, y`` or ``phi`` combined so that the
azimuth cancels.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from ._sympy import require_sympy

sp = require_sympy()

from ._naming import C_RESERVED, IDENTIFIER  # noqa: E402
from .forces import _CONSTANTS, _FUNCTIONS, parse_expression  # noqa: E402
from .optimizer import optimize  # noqa: E402

__all__ = ["DensitySource", "build_surface_density_source"]

PLANE_NAMES = ("R", "r", "phi", "x", "y")
_RESERVED = frozenset({"radius", "context", "parameters", "value", "z", "theta"})


@dataclass(frozen=True)
class DensitySource:
    """Generated surface density.

    Attributes
    ----------
    source : str
        Complete C11 translation unit; it includes only ``<math.h>``.
    function_name : str
        Name of the exported ``rp_surface_density_fn``.
    parameters : tuple of str
        Parameter names; ``context`` must point to ``len(parameters)``
        doubles in this order (it may be NULL when there are none).
    expression : sympy.Expr
        ``Sigma`` in terms of the positive symbol ``R`` and the parameters.
    """

    source: str
    function_name: str
    parameters: tuple[str, ...]
    expression: object


def _validate_parameters(parameters: Sequence[str]) -> tuple[str, ...]:
    names = tuple(parameters)
    for name in names:
        if not isinstance(name, str) or not IDENTIFIER.match(name):
            raise ValueError(f"parameter {name!r} is not a valid identifier")
        if name in C_RESERVED or name in _RESERVED or name in PLANE_NAMES:
            raise ValueError(f"parameter name {name!r} is reserved")
        if name in _FUNCTIONS or name in _CONSTANTS:
            raise ValueError(f"parameter name {name!r} shadows a function or constant")
        if name.startswith("__"):
            raise ValueError(f"parameter name {name!r} must not start with '__'")
    if len(set(names)) != len(names):
        raise ValueError("parameter names must be unique")
    return names


def _azimuth_free(expression, radius, azimuth, values: Mapping) -> bool:
    """True when ``expression`` does not depend on ``azimuth``."""
    if azimuth not in expression.free_symbols:
        return True
    derivative = sp.simplify(sp.diff(expression, azimuth))
    if derivative == 0:
        return True
    # SymPy may fail to prove zero: sample the derivative instead.
    numeric = sp.lambdify((radius, azimuth), derivative.xreplace(values), "math")
    scale = abs(complex(sp.lambdify((radius, azimuth), expression.xreplace(values), "math")(0.5, 0.3)))
    for radius_value in (0.013, 0.21, 0.5, 0.87):
        for angle in (0.1, 1.3, 2.9, 4.4):
            try:
                value = complex(numeric(radius_value, angle))
            except (ValueError, ZeroDivisionError, OverflowError):
                return False
            if not abs(value) <= 1e-12 * max(scale, 1.0):
                return False
    return True


def build_surface_density_source(
    text: str,
    *,
    parameters: Sequence[str] = (),
    values: Mapping[str, float] | None = None,
    name: str = "rp_jit_surface_density",
) -> DensitySource:
    """Build the C source of an axisymmetric surface density.

    Parameters
    ----------
    text
        ``Sigma`` as text parsed by :func:`~relatipy._codegen.forces.parse_expression`
        with the names ``R, r, phi, x, y`` (lengths in units of the outer
        radius) and the declared parameters.
    parameters
        Names of the shape parameters, read in this order from ``context``.
    values
        Parameter values used only to check axisymmetry when SymPy cannot
        prove it.
    name
        C name of the generated function.

    Raises
    ------
    ValueError
        For an invalid name, undeclared names, or a density that depends on
        the azimuth.
    """
    if not IDENTIFIER.match(name) or name in C_RESERVED:
        raise ValueError(f"{name!r} is not a usable C function name")
    parameter_names = _validate_parameters(parameters)
    radius = sp.Symbol("R", positive=True)
    azimuth = sp.Symbol("phi", real=True)
    symbols = {
        "R": radius, "r": radius, "phi": azimuth,
        "x": radius * sp.cos(azimuth), "y": radius * sp.sin(azimuth),
    }
    parameter_symbols = tuple(sp.Symbol(item, real=True) for item in parameter_names)
    namespace = symbols | {symbol.name: symbol for symbol in parameter_symbols}
    expression = sp.sympify(parse_expression(text, namespace))
    unknown = expression.free_symbols - {radius, azimuth, *parameter_symbols}
    if unknown:
        raise ValueError("surface density uses undeclared names: "
                         + ", ".join(sorted(item.name for item in unknown)))
    numeric_values = {
        symbol: float(values[symbol.name]) if values and symbol.name in values else 0.37 + 0.1 * k
        for k, symbol in enumerate(parameter_symbols)
    }
    if not _azimuth_free(expression, radius, azimuth, numeric_values):
        raise ValueError("surface density must be axisymmetric (independent of phi)")
    expression = sp.simplify(expression.subs(azimuth, 0)) if azimuth in expression.free_symbols else expression

    inputs = {radius: "radius"} | {symbol: f"parameters[{i}]" for i, symbol in enumerate(parameter_symbols)}
    code = optimize(expression, inputs=inputs, outputs="value",
                    reserved=("radius", "context", "parameters", "value", name))
    prologue = ["double value;"]
    if parameter_symbols:
        prologue.insert(0, "const double *parameters = (const double *)context;")
    else:
        prologue.insert(0, "(void)context;")
    if radius not in expression.free_symbols:
        prologue.append("(void)radius;")
    layout = ", ".join(f"context[{i}] = {item}" for i, item in enumerate(parameter_names)) or "none"
    function = code.c_function(
        name,
        ["double radius", "const void *context"],
        return_type="double",
        doc=("Run-time surface density Sigma(R), R in units of the outer radius.\n"
             f"Parameters read from context as const double: {layout}."),
        prologue=prologue,
        epilogue=("return value;",),
    )
    header = (
        "/* Generated at run time by relatipy._codegen.density; do not edit. */\n"
        "#include <math.h>\n\n"
    )
    return DensitySource(source=header + function, function_name=name,
                         parameters=parameter_names, expression=expression)


def central_slope(expression, values: Mapping[str, float]) -> float:
    """Logarithmic slope ``d ln Sigma / d ln R`` at ``R -> 0``.

    Returns ``inf`` for a density that vanishes faster than any power, and
    ``nan`` when the slope cannot be determined.
    """
    radius = sp.Symbol("R", positive=True)
    substituted = expression.xreplace({sp.Symbol(k, real=True): v for k, v in values.items()})
    try:
        # Floats mislead SymPy's limit (R**-1.7 gives -oo): use exact rationals.
        exact = sp.nsimplify(substituted, rational=True)
        slope = sp.limit(radius * sp.diff(exact, radius) / exact, radius, 0, "+")
        if slope.is_real and slope.is_finite:
            return float(slope)
        if slope is sp.oo or slope is -sp.oo:
            return math.copysign(math.inf, float(slope))
    except (NotImplementedError, ValueError, TypeError):
        pass
    function = sp.lambdify(radius, substituted, "mpmath")
    import mpmath

    with mpmath.workdps(40):
        small, smaller = mpmath.mpf("1e-20"), mpmath.mpf("1e-22")
        try:
            low, lower = function(small), function(smaller)
            return float((mpmath.log(low) - mpmath.log(lower)) / mpmath.log(small / smaller))
        except (ValueError, ZeroDivisionError, TypeError):
            return math.nan
