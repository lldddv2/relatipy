"""C source of an ``rp_external_acceleration_fn`` from a symbolic definition.

:func:`build_force_source` accepts a Newtonian potential or acceleration as
a Python callable (called with SymPy symbols) or as a string (parsed in a
restricted namespace), in Cartesian, cylindrical or spherical coordinates of
the black-hole frame, and returns one C translation unit that implements the
callback of ``native/src/force/external_force.h``::

    int name(const double position[3], const double velocity[3],
             const void *model, double acceleration[3]);

Declared parameters are read at run time from ``model``, interpreted as
``const double *``: ``model[i]`` is the ``i``-th declared parameter.  Changing
a parameter value therefore needs no recompilation.
"""

from __future__ import annotations

import ast
from collections.abc import Callable, Sequence
from dataclasses import dataclass

from ._sympy import require_sympy

sp = require_sympy()

from ._naming import C_RESERVED, IDENTIFIER  # noqa: E402
from ._printer import c_expression  # noqa: E402
from .optimizer import OptimizationReport, optimize, reduce_angle  # noqa: E402
from .special import C_HELPERS, Quad  # noqa: E402

__all__ = ["COORDINATES", "KINDS", "ForceSource", "build_force_source", "parse_expression"]

KINDS = ("potential", "force")
COORDINATES = ("cartesian", "cylindrical", "spherical")

COORDINATE_NAMES = {
    "cartesian": ("x", "y", "z"),
    "cylindrical": ("R", "phi", "z"),
    "spherical": ("r", "theta", "phi"),
}
VELOCITY_NAMES = {
    "cartesian": ("vx", "vy", "vz"),
    "cylindrical": ("v_R", "v_phi", "v_z"),
    "spherical": ("v_r", "v_theta", "v_phi"),
}
_SIGNATURE_NAMES = frozenset({"position", "velocity", "model", "acceleration", "parameters"})

_FUNCTIONS = {
    "sqrt": sp.sqrt, "exp": sp.exp, "log": sp.log,
    "sin": sp.sin, "cos": sp.cos, "tan": sp.tan,
    "asin": sp.asin, "acos": sp.acos, "atan": sp.atan, "atan2": sp.atan2,
    "sinh": sp.sinh, "cosh": sp.cosh, "tanh": sp.tanh,
    "abs": sp.Abs, "Abs": sp.Abs, "erf": sp.erf,
}
_CONSTANTS = {"pi": sp.pi}
_ALLOWED_NODES = (
    ast.Expression, ast.BinOp, ast.UnaryOp, ast.Call, ast.Name, ast.Load,
    ast.Constant, ast.Tuple, ast.List,
    ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow, ast.USub, ast.UAdd,
)


@dataclass(frozen=True)
class ForceSource:
    """Generated callback source.

    Attributes
    ----------
    source : str
        Complete C11 translation unit; it includes only ``<math.h>``.
    function_name : str
        Name of the exported ``rp_external_acceleration_fn``.
    parameters : tuple of str
        Parameter names; ``model`` must point to ``len(parameters)`` doubles
        in this order (it may be NULL when there are none).
    kind : str
        ``"potential"`` or ``"force"``.
    coords : str
        ``"cartesian"``, ``"cylindrical"`` or ``"spherical"``.
    acceleration : tuple
        The three Cartesian acceleration components as SymPy expressions of
        ``x, y, z, vx, vy, vz`` and the parameters.
    report : OptimizationReport
        Operation counts of the generated body.
    """

    source: str
    function_name: str
    parameters: tuple[str, ...]
    kind: str
    coords: str
    acceleration: tuple
    report: OptimizationReport


def parse_expression(text: str, namespace: dict):
    """Parse ``text`` with :func:`sympy.parse_expr` in a restricted namespace.

    Only numbers, the names in ``namespace``, the functions ``sqrt exp log
    sin cos tan asin acos atan atan2 sinh cosh tanh abs erf``, the constant
    ``pi``, the operators ``+ - * / **`` and tuples or lists are accepted.

    Raises
    ------
    ValueError
        If the text contains any other construct or name.
    """
    if not isinstance(text, str):
        raise TypeError("expression text must be a string")
    try:
        tree = ast.parse(text.strip(), mode="eval")
    except SyntaxError as exc:
        raise ValueError(f"invalid expression {text!r}: {exc.msg}") from exc
    allowed_names = set(namespace) | set(_FUNCTIONS) | set(_CONSTANTS)
    for node in ast.walk(tree):
        if isinstance(node, ast.BitXor):
            raise ValueError("use ** for powers; ^ is not accepted")
        if not isinstance(node, _ALLOWED_NODES):
            raise ValueError(f"unsupported syntax {type(node).__name__} in {text!r}")
        if isinstance(node, ast.Constant) and not (
            isinstance(node.value, (int, float)) and not isinstance(node.value, bool)
        ):
            raise ValueError(f"unsupported constant {node.value!r} in {text!r}")
        if isinstance(node, ast.Name) and node.id not in allowed_names:
            raise ValueError(f"unknown name {node.id!r} in {text!r}")
        if isinstance(node, ast.Call):
            if not isinstance(node.func, ast.Name) or node.func.id not in _FUNCTIONS:
                raise ValueError(f"unsupported function call in {text!r}")
            if node.keywords:
                raise ValueError(f"keyword arguments are not accepted in {text!r}")
    local = {**_FUNCTIONS, **_CONSTANTS, **namespace}
    global_namespace = {
        "__builtins__": {},
        "Integer": sp.Integer, "Float": sp.Float, "Rational": sp.Rational, "Symbol": sp.Symbol,
    }
    from sympy.parsing.sympy_parser import parse_expr, standard_transformations

    return parse_expr(
        text.strip(), local_dict=local, global_dict=global_namespace,
        transformations=standard_transformations, evaluate=True,
    )


def _user_symbols(coords: str):
    names = COORDINATE_NAMES[coords]
    positive = {"R", "r"}
    coordinates = tuple(
        sp.Symbol(name, positive=True) if name in positive else sp.Symbol(name, real=True)
        for name in names
    )
    velocities = tuple(sp.Symbol(name, real=True) for name in VELOCITY_NAMES[coords])
    return coordinates, velocities


def _validate_parameters(parameters: Sequence[str], coords: str) -> tuple[str, ...]:
    names = tuple(parameters)
    reserved = set(COORDINATE_NAMES[coords]) | set(VELOCITY_NAMES[coords]) | {"x", "y", "z", "vx", "vy", "vz"}
    for name in names:
        if not isinstance(name, str) or not IDENTIFIER.match(name):
            raise ValueError(f"parameter {name!r} is not a valid identifier")
        if name in C_RESERVED or name in _SIGNATURE_NAMES or name in reserved:
            raise ValueError(f"parameter name {name!r} is reserved")
        if name in _FUNCTIONS or name in _CONSTANTS:
            raise ValueError(f"parameter name {name!r} shadows a function or constant")
        if name.startswith("__"):
            raise ValueError(f"parameter name {name!r} must not start with '__'")
    if len(set(names)) != len(names):
        raise ValueError("parameter names must be unique")
    return names


def _definition_expressions(definition, kind, coordinates, velocities, parameter_symbols):
    if callable(definition):
        arguments = (*coordinates, *velocities, *parameter_symbols) if kind == "force" else (
            *coordinates, *parameter_symbols)
        value = definition(*arguments)
    else:
        namespace = {symbol.name: symbol for symbol in (*coordinates, *parameter_symbols)}
        if kind == "force":
            namespace |= {symbol.name: symbol for symbol in velocities}
        if isinstance(definition, str):
            value = parse_expression(definition, namespace)
        elif isinstance(definition, Sequence):
            value = [parse_expression(item, namespace) for item in definition]
        else:
            raise TypeError("definition must be a callable, a string or a sequence of strings")
    if kind == "potential":
        if isinstance(value, (sp.MatrixBase, sp.NDimArray, tuple, list)):
            raise ValueError("a potential must be one scalar expression")
        return sp.sympify(value)
    if isinstance(value, (sp.MatrixBase, sp.Tuple, tuple, list)):
        components = list(value)
    else:
        components = None
    if components is None or len(components) != 3:
        raise ValueError("a force must have three components")
    return [sp.sympify(component) for component in components]


def _to_cartesian(expr, coords, coordinates, cartesian):
    """Rewrite an expression of user coordinates in Cartesian coordinates."""
    x, y, z = cartesian
    rho = sp.sqrt(x**2 + y**2)
    sine_phi, cosine_phi = sp.Dummy("sin_phi"), sp.Dummy("cos_phi")
    if coords == "cartesian":
        return expr.xreplace(dict(zip(coordinates, cartesian)))
    if coords == "cylindrical":
        radius, phi, height = coordinates
        expr = reduce_angle(expr, phi, sine_phi, cosine_phi)
        expr = expr.xreplace({sine_phi: y / rho, cosine_phi: x / rho})
        return expr.xreplace({phi: sp.atan2(y, x), radius: rho, height: z})
    radius, theta, phi = coordinates
    distance = sp.sqrt(x**2 + y**2 + z**2)
    sine_theta, cosine_theta = sp.Dummy("sin_theta"), sp.Dummy("cos_theta")
    expr = reduce_angle(expr, phi, sine_phi, cosine_phi)
    expr = reduce_angle(expr, theta, sine_theta, cosine_theta)
    expr = expr.xreplace({
        sine_phi: y / rho, cosine_phi: x / rho,
        sine_theta: rho / distance, cosine_theta: z / distance,
    })
    return expr.xreplace({phi: sp.atan2(y, x), theta: sp.atan2(rho, z), radius: distance})


def _basis(coords, coordinates):
    """Orthonormal basis vectors (Cartesian components) and scale factors."""
    if coords == "cartesian":
        return (sp.Matrix([1, 0, 0]), sp.Matrix([0, 1, 0]), sp.Matrix([0, 0, 1])), (1, 1, 1)
    if coords == "cylindrical":
        radius, phi, _ = coordinates
        basis = (
            sp.Matrix([sp.cos(phi), sp.sin(phi), 0]),
            sp.Matrix([-sp.sin(phi), sp.cos(phi), 0]),
            sp.Matrix([0, 0, 1]),
        )
        return basis, (1, radius, 1)
    radius, theta, phi = coordinates
    basis = (
        sp.Matrix([sp.sin(theta) * sp.cos(phi), sp.sin(theta) * sp.sin(phi), sp.cos(theta)]),
        sp.Matrix([sp.cos(theta) * sp.cos(phi), sp.cos(theta) * sp.sin(phi), -sp.sin(theta)]),
        sp.Matrix([-sp.sin(phi), sp.cos(phi), 0]),
    )
    return basis, (1, radius, radius * sp.sin(theta))


def _extract_quadratures(acceleration, inputs, name, has_parameters):
    """Replace run-time integrals by array reads; return their C integrands and calls."""
    found = sorted(
        set().union(*(sp.sympify(item).atoms(Quad) for item in acceleration)),
        key=sp.default_sort_key,
    )
    if not found:
        return acceleration, [], []
    groups: dict = {}
    for quad in found:
        groups.setdefault(quad.args[1:], []).append(quad)
    replacements = {}
    sources, calls = [], []
    limit_names = {symbol: sp.Symbol(text) for symbol, text in inputs.items()}
    for index, (key, members) in enumerate(groups.items()):
        variable, lower, upper, singular, flags = key
        if len(members) > 16:
            raise ValueError("at most 16 integrals may share one variable and range")
        for k, quad in enumerate(members):
            if any(item.atoms(Quad) for item in quad.args[0].args):
                raise ValueError("nested integrals are not supported")
            replacements[quad] = sp.Symbol(f"rp_quad_{index}_{k}", real=True)
        integrand_name = f"{name}_integrand_{index}"
        integrand_inputs = dict(inputs) | {variable: "rp_s"}
        code = optimize(
            sp.Matrix([quad.args[0] for quad in members]), inputs=integrand_inputs,
            outputs="out", cse_modes=("basic",), transforms=("none",),
            reserved=("rp_s", "position", "velocity", "parameters", "out", integrand_name),
        )
        sources.append(code.c_function(
            integrand_name,
            ["double rp_s", "const double *position", "const double *velocity",
             "const double *parameters", "double *out"],
            static=True,
            prologue=("(void)rp_s;", "(void)position;", "(void)velocity;", "(void)parameters;"),
        ) + "\n")

        def limit(expr):
            return c_expression(sp.sympify(expr).xreplace(limit_names))

        calls.append(f"double rp_quad_{index}[{len(members)}];")
        calls.append(
            f"rp_jit_quad({integrand_name}, {len(members)}, {limit(lower)}, {limit(upper)}, "
            f"{limit(singular)}, {int(flags)}, position, velocity, "
            f"{'parameters' if has_parameters else '(const double *)0'}, rp_quad_{index});"
        )
        for k, quad in enumerate(members):
            inputs[replacements[quad]] = f"rp_quad_{index}[{k}]"
    acceleration = [sp.sympify(item).xreplace(replacements) for item in acceleration]
    return acceleration, sources, calls


def build_force_source(
    definition: Callable | str | Sequence[str],
    *,
    kind: str,
    coords: str = "cartesian",
    parameters: Sequence[str] = (),
    name: str = "rp_jit_external_acceleration",
) -> ForceSource:
    """Build the C source of an ``rp_external_acceleration_fn``.

    Parameters
    ----------
    definition
        For ``kind="potential"``, a scalar potential ``Phi``; for
        ``kind="force"``, three acceleration components in the orthonormal
        basis of ``coords``.  Either a callable or text: a callable receives
        the coordinate symbols, then (for a force) the three velocity
        components in the same orthonormal basis, then the parameter symbols,
        all positionally; text is parsed by :func:`parse_expression` with
        those names (one string, a tuple string ``"(a, b, c)"`` for a force,
        or a sequence of three strings).
    kind
        ``"potential"`` (the acceleration is ``-grad Phi``, obtained by the
        chain rule) or ``"force"``.
    coords
        ``"cartesian"`` with names ``x, y, z`` and velocities ``vx, vy, vz``;
        ``"cylindrical"`` with ``R, phi, z`` and ``v_R, v_phi, v_z``; or
        ``"spherical"`` with ``r, theta, phi`` and ``v_r, v_theta, v_phi``.
        Axes are those of the black-hole Cartesian frame (``z`` along the
        spin); ``R`` and ``r`` are the flat cylindrical and spherical radii
        of that frame.
    parameters
        Names of run-time parameters, read in this order from ``model``.
    name
        C name of the generated function.

    Returns
    -------
    ForceSource
        Source, layout and report.  The function returns ``0`` on success
        and ``1`` when an acceleration component is not finite.

    Raises
    ------
    ValueError
        For an unknown ``kind`` or ``coords``, an invalid name, a definition
        that uses undeclared names, or a component count other than three.
    """
    if kind not in KINDS:
        raise ValueError(f"kind must be one of {KINDS}, not {kind!r}")
    if coords not in COORDINATES:
        raise ValueError(f"coords must be one of {COORDINATES}, not {coords!r}")
    if not IDENTIFIER.match(name) or name in C_RESERVED:
        raise ValueError(f"{name!r} is not a usable C function name")
    parameter_names = _validate_parameters(parameters, coords)
    coordinates, velocities = _user_symbols(coords)
    parameter_symbols = tuple(sp.Symbol(item, real=True) for item in parameter_names)
    value = _definition_expressions(definition, kind, coordinates, velocities, parameter_symbols)

    basis, scale = _basis(coords, coordinates)
    if kind == "potential":
        components = [-sp.diff(value, q) / h for q, h in zip(coordinates, scale)]
    else:
        components = value
    allowed = set(coordinates) | set(parameter_symbols) | (set(velocities) if kind == "force" else set())
    unknown = set().union(*(sp.sympify(c).free_symbols for c in components)) - allowed
    if unknown:
        names = ", ".join(sorted(symbol.name for symbol in unknown))
        raise ValueError(f"definition uses undeclared names: {names}")

    cartesian = sp.symbols("x y z", real=True)
    cartesian_velocity = sp.symbols("vx vy vz", real=True)
    velocity_vector = sp.Matrix(cartesian_velocity)
    velocity_components = {
        symbol: _to_cartesian((basis[k].T * velocity_vector)[0], coords, coordinates, cartesian)
        for k, symbol in enumerate(velocities)
    }
    acceleration_vector = sp.zeros(3, 1)
    for component, vector in zip(components, basis):
        acceleration_vector += sp.sympify(component) * vector
    acceleration = []
    for entry in acceleration_vector:
        entry = entry.xreplace(velocity_components)
        acceleration.append(_to_cartesian(entry, coords, coordinates, cartesian))

    inputs = {symbol: f"position[{i}]" for i, symbol in enumerate(cartesian)}
    inputs |= {symbol: f"velocity[{i}]" for i, symbol in enumerate(cartesian_velocity)}
    inputs |= {symbol: f"parameters[{i}]" for i, symbol in enumerate(parameter_symbols)}
    acceleration, quadratures, quad_calls = _extract_quadratures(
        acceleration, inputs, name, bool(parameter_symbols))
    used_raw = set().union(*(sp.sympify(item).free_symbols for item in acceleration))
    code = optimize(
        sp.Matrix(acceleration), inputs=inputs, outputs="acceleration",
        reserved=("position", "velocity", "model", "parameters", "acceleration", name),
    )
    used = used_raw
    prologue = []
    if parameter_symbols:
        prologue.append("const double *parameters = (const double *)model;")
    else:
        prologue.append("(void)model;")
    if quad_calls:
        prologue.extend(quad_calls)
    else:
        if not used & set(cartesian):
            prologue.append("(void)position;")
        if not used & set(cartesian_velocity):
            prologue.append("(void)velocity;")
    layout = ", ".join(f"model[{i}] = {item}" for i, item in enumerate(parameter_names)) or "none"
    doc = (
        f"Run-time {kind} in {coords} coordinates of the black-hole frame.\n"
        f"Parameters read from model as const double: {layout}.\n"
        "Returns 0 on success and 1 when an acceleration component is not finite."
    )
    function = code.c_function(
        name,
        [
            "const double position[3]",
            "const double velocity[3]",
            "const void *model",
            "double acceleration[3]",
        ],
        return_type="int",
        doc=doc,
        prologue=prologue,
        epilogue=(
            "if (!isfinite(acceleration[0]) || !isfinite(acceleration[1])",
            "    || !isfinite(acceleration[2])) {",
            "    return 1;",
            "}",
            "return 0;",
        ),
    )
    header = (
        "/* Generated at run time by relatipy._codegen.forces; do not edit. */\n"
        "#include <math.h>\n\n"
    )
    if quadratures or "rp_jit_ellip" in function:
        header += C_HELPERS + "\n" + "".join(quadratures)
    return ForceSource(
        source=header + function,
        function_name=name,
        parameters=parameter_names,
        kind=kind,
        coords=coords,
        acceleration=tuple(acceleration),
        report=code.report,
    )
