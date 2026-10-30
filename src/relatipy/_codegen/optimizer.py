"""Symbolic optimization of SymPy expressions into semantic C99.

:func:`optimize` turns scalar expressions, matrices, tensors, or sequences
and mappings of them into straight-line C99 statements:

1. trigonometric functions of integer multiples of each declared angle are
   reduced to the sine and cosine of that angle, which are evaluated (or
   read) once; sums of the form ``A - A cos^2`` become ``A sin^2`` and
   ``A - A sin^2`` become ``A cos^2``;
2. a reproducible strategy search applies algebraic pre-transforms and
   common-subexpression elimination across all outputs, and keeps the
   candidate with the fewest operations, then the fewest temporaries, then
   the lexicographically smallest strategy name;
3. integer powers become multiplication chains that reuse named powers
   (``r2 = r * r``, ``r3 = r2 * r``, ``r4 = r2 * r2``), and half-integer
   powers reuse one ``sqrt`` per radicand;
4. every temporary gets a deterministic semantic name derived from its
   definition (never ``tmp`` or ``cse``), and statements are emitted as
   ``const double`` declarations grouped under comments.
"""

from __future__ import annotations

import heapq
import itertools
from collections import Counter
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field

from ._sympy import require_sympy

sp = require_sympy()

from ._naming import C_RESERVED, IDENTIFIER, Namer  # noqa: E402
from ._printer import c_expression, count_ops as _count_expr, wrap_statement  # noqa: E402

__all__ = [
    "Definition",
    "OperationCount",
    "OptimizationReport",
    "OptimizedCode",
    "count_operations",
    "optimize",
    "reduce_angle",
]

_TRIGONOMETRIC = (sp.sin, sp.cos, sp.tan, sp.cot, sp.sec, sp.csc)
_INTERNAL_PREFIX = "__cg_"


def _identity(expr):
    return expr


TRANSFORMS: dict[str, Callable] = {
    "none": _identity,
    "factor_terms": sp.factor_terms,
    "cancel": sp.cancel,
    "factor": sp.factor,
}
DEFAULT_TRANSFORMS = ("none", "factor_terms", "cancel")
CSE_MODES = {"raw": None, "basic": "basic"}


# ---------------------------------------------------------------------------
# Public result types
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class OperationCount:
    """C-level arithmetic operations of generated code.

    Attributes
    ----------
    add, mul, div, call : int
        Additions and subtractions, multiplications, divisions, and function
        calls (``sqrt``, ``sin``, ``pow`` and others).
    """

    add: int = 0
    mul: int = 0
    div: int = 0
    call: int = 0

    @classmethod
    def from_counter(cls, counter: Mapping[str, int]) -> OperationCount:
        return cls(**{key: int(counter.get(key, 0)) for key in ("add", "mul", "div", "call")})

    @property
    def total(self) -> int:
        """Sum of all operations."""
        return self.add + self.mul + self.div + self.call

    def __add__(self, other: OperationCount) -> OperationCount:
        return OperationCount(
            self.add + other.add, self.mul + other.mul,
            self.div + other.div, self.call + other.call,
        )

    def __str__(self) -> str:
        return (
            f"{self.total} (add {self.add}, mul {self.mul}, "
            f"div {self.div}, call {self.call})"
        )


@dataclass(frozen=True)
class CandidateSummary:
    """Cost of one strategy of the search."""

    name: str
    operations: int
    temporaries: int


@dataclass(frozen=True)
class OptimizationReport:
    """Operation counts before and after optimization.

    Attributes
    ----------
    strategy : str
        Name of the selected strategy, ``"<transform>+cse:<mode>"``.
    before : OperationCount
        Operations of the input expressions evaluated as written, with
        integer powers by binary exponentiation and every trigonometric
        function as a call.
    after : OperationCount
        Operations of the emitted statements.
    temporaries : int
        Number of ``const double`` temporaries other than input bindings.
    candidates : tuple of CandidateSummary
        Every evaluated strategy, in search order.
    """

    strategy: str
    before: OperationCount
    after: OperationCount
    temporaries: int
    candidates: tuple[CandidateSummary, ...]

    def __str__(self) -> str:
        lines = [
            f"operations before: {self.before}",
            f"operations after:  {self.after}",
            f"temporaries:       {self.temporaries}",
            f"strategy:          {self.strategy}",
            "candidates:",
        ]
        lines.extend(
            f"  {item.name:<32} ops {item.operations:>6}  temporaries {item.temporaries:>4}"
            for item in self.candidates
        )
        return "\n".join(lines)


@dataclass(frozen=True)
class Definition:
    """One emitted ``const double`` statement.

    Attributes
    ----------
    name : str
        C identifier of the temporary.
    kind : str
        ``"input"``, ``"angle"``, ``"power"`` or ``"auxiliary"``.
    expression : sympy.Expr
        Definition in terms of earlier definition names.
    c_text : str
        Right-hand side as printed C.
    group : str
        Comment heading of the group the statement is emitted under.
    """

    name: str
    kind: str
    expression: object
    c_text: str
    group: str


@dataclass(frozen=True)
class OptimizedCode:
    """Result of :func:`optimize`.

    Attributes
    ----------
    definitions : tuple of Definition
        Statements in emission order.
    outputs : tuple of tuple
        ``(lvalue, expression, c_text)`` for each output component.
    report : OptimizationReport
        Operation counts and the strategy search.
    originals : tuple
        The flattened input expressions, aligned with ``outputs``.
    """

    definitions: tuple[Definition, ...]
    outputs: tuple[tuple[str, object, str], ...]
    report: OptimizationReport
    originals: tuple
    _expansion: Mapping = field(repr=False, compare=False)
    _input_symbols: tuple = field(repr=False, compare=False)
    _angles: tuple = field(repr=False, compare=False, default=())

    @property
    def lvalues(self) -> tuple[str, ...]:
        """Output lvalues in emission order."""
        return tuple(item[0] for item in self.outputs)

    @property
    def angles(self) -> tuple:
        """Declared angle symbols."""
        return self._angles

    @property
    def input_symbols(self) -> tuple:
        """Input and angle symbols the expanded outputs may depend on."""
        return self._input_symbols

    def body(self, indent: str = "    ") -> str:
        """C99 statements: grouped declarations followed by output stores."""
        lines: list[str] = []
        group = None
        for definition in self.definitions:
            if definition.c_text == definition.name:
                continue
            if definition.group != group:
                if lines:
                    lines.append("")
                lines.append(f"{indent}/* {definition.group} */")
                group = definition.group
            lines.extend(
                wrap_statement(f"const double {definition.name} = ", definition.c_text, indent)
            )
        if lines:
            lines.append("")
        lines.append(f"{indent}/* Outputs. */")
        for lvalue, _, text in self.outputs:
            lines.extend(wrap_statement(f"{lvalue} = ", text, indent))
        return "\n".join(lines)

    def prototype(
        self, name: str, parameters: Sequence[str], *, return_type: str = "void",
        static: bool = False,
    ) -> str:
        """Function declarator without a trailing semicolon or body."""
        if not IDENTIFIER.match(name):
            raise ValueError(f"{name!r} is not a valid C identifier")
        storage = "static " if static else ""
        joined = ",\n".join(f"    {parameter}" for parameter in parameters) or "    void"
        return f"{storage}{return_type} {name}(\n{joined}\n)"

    def c_function(
        self, name: str, parameters: Sequence[str], *, return_type: str = "void",
        static: bool = False, doc: str | None = None,
        prologue: Sequence[str] = (), epilogue: Sequence[str] = (),
    ) -> str:
        """Complete C function wrapping :meth:`body`.

        Parameters
        ----------
        name
            C identifier of the function.
        parameters
            Parameter declarations, e.g. ``"const double u[4]"``.
        return_type
            C return type; a non-``void`` type needs an ``epilogue`` that
            returns a value.
        static
            Emit the function with internal linkage.
        doc
            Optional comment placed above the function.
        prologue, epilogue
            Extra statements emitted before and after the body.
        """
        lines = []
        if doc:
            lines.append("/**")
            lines.extend(f" * {line}".rstrip() for line in doc.strip().splitlines())
            lines.append(" */")
        lines.append(self.prototype(name, parameters, return_type=return_type, static=static))
        lines.append("{")
        lines.extend(f"    {statement}" for statement in prologue)
        if prologue:
            lines.append("")
        lines.append(self.body("    "))
        if epilogue:
            lines.append("")
        lines.extend(f"    {statement}" for statement in epilogue)
        lines.append("}")
        return "\n".join(lines) + "\n"

    def expanded(self) -> tuple:
        """Outputs rewritten in terms of the inputs and declared angles."""
        return tuple(sp.sympify(expr).xreplace(self._expansion) for _, expr, _ in self.outputs)

    def evaluate(self, values: Mapping) -> list:
        """Evaluate the statements in order with NumPy, mirroring the C code.

        Parameters
        ----------
        values
            Mapping from every input symbol (and every angle symbol whose sine
            or cosine is read from an input expression) to a float or array.

        Returns
        -------
        list of numpy.ndarray
            One array per output component.
        """
        import numpy as np

        environment: dict[str, object] = {}
        by_name = {symbol.name: value for symbol, value in values.items()}
        for symbol in self._input_symbols:
            if symbol.name in by_name:
                environment[symbol.name] = by_name[symbol.name]

        def run(expression, scope):
            expression = sp.sympify(expression)
            arguments = sorted(expression.free_symbols, key=lambda item: item.name)
            function = sp.lambdify(arguments, expression, modules="numpy")
            return function(*(scope[a.name] for a in arguments))

        for definition in self.definitions:
            if definition.kind == "input":
                continue
            if definition.kind == "angle":
                environment[definition.name] = run(definition.expression, by_name)
                continue
            environment[definition.name] = run(definition.expression, environment)
        return [
            np.asarray(run(expression, environment), dtype=float) * np.ones_like(
                np.asarray(next(iter(by_name.values()), 0.0), dtype=float)
            )
            for _, expression, _ in self.outputs
        ]


# ---------------------------------------------------------------------------
# Trigonometric reduction
# ---------------------------------------------------------------------------
def reduce_angle(expr, angle, sine, cosine):
    """Rewrite trigonometric functions of ``n * angle`` with ``sine, cosine``.

    Parameters
    ----------
    expr
        SymPy expression.
    angle
        Angle symbol.
    sine, cosine
        Symbols that stand for ``sin(angle)`` and ``cos(angle)``.

    Returns
    -------
    sympy.Expr
        Expression in which every ``sin, cos, tan, cot, sec, csc`` of an
        integer multiple of ``angle`` is a rational function of ``sine`` and
        ``cosine``, followed by the ``A - A cos^2 -> A sin^2`` contraction.
        Other functions of ``angle`` are left unchanged.
    """
    expr = sp.sympify(expr)
    base = {
        sp.sin(angle): sine,
        sp.cos(angle): cosine,
        sp.tan(angle): sine / cosine,
        sp.cot(angle): cosine / sine,
        sp.sec(angle): 1 / cosine,
        sp.csc(angle): 1 / sine,
    }
    mapping = {}
    for function in expr.atoms(*_TRIGONOMETRIC):
        argument = function.args[0]
        multiple = sp.simplify(argument / angle)
        if not multiple.is_Integer or multiple == 0:
            continue
        if multiple == 1:
            mapping[function] = base[type(function)(angle)]
            continue
        if isinstance(function, (sp.tan, sp.cot, sp.sec, sp.csc)):
            sine_n = sp.expand_trig(sp.sin(multiple * angle))
            cosine_n = sp.expand_trig(sp.cos(multiple * angle))
            rewritten = {
                sp.tan: sine_n / cosine_n, sp.cot: cosine_n / sine_n,
                sp.sec: 1 / cosine_n, sp.csc: 1 / sine_n,
            }[type(function)]
        else:
            rewritten = sp.expand_trig(type(function)(multiple * angle))
        mapping[function] = rewritten.xreplace(base)
    reduced = expr.xreplace(mapping) if mapping else expr
    return contract_pythagorean(reduced, sine, cosine)


def contract_pythagorean(expr, sine, cosine):
    """Replace ``A - A cosine^2`` by ``A sine^2`` (and vice versa) in sums.

    Only exact pairs of terms are combined, so the rewrite removes a
    cancellation without introducing one.
    """
    expr = sp.sympify(expr)
    if expr.is_Atom:
        return expr
    arguments = tuple(contract_pythagorean(argument, sine, cosine) for argument in expr.args)
    if arguments != expr.args:
        expr = expr.func(*arguments)
    if expr.is_Add:
        expr = _contract_sum(expr, sine, cosine)
    return expr


def _contract_sum(expr, sine, cosine):
    terms: dict[tuple, object] = {}
    order: list[tuple] = []
    for term in expr.args:
        coefficient, rest = term.as_coeff_Mul()
        powers = rest.as_powers_dict() if rest != 1 else {}
        powers = dict(powers)
        p = powers.pop(cosine, sp.S.Zero)
        q = powers.pop(sine, sp.S.Zero)
        remainder = sp.Mul(*[b**e for b, e in sorted(powers.items(), key=lambda item: sp.default_sort_key(item[0]))])
        key = (remainder, p, q)
        if key not in terms:
            order.append(key)
            terms[key] = sp.S.Zero
        terms[key] += coefficient
    changed = True
    modified = False
    while changed:
        changed = False
        for key in list(order):
            if key not in terms or terms[key] == 0:
                continue
            remainder, p, q = key
            for partner, target in (
                ((remainder, p + 2, q), (remainder, p, q + 2)),
                ((remainder, p, q + 2), (remainder, p + 2, q)),
            ):
                if partner in terms and terms[partner] != 0 and terms[partner] == -terms[key]:
                    value = terms[key]
                    terms[key] = sp.S.Zero
                    terms[partner] = sp.S.Zero
                    if target not in terms:
                        order.append(target)
                        terms[target] = sp.S.Zero
                    terms[target] += value
                    changed = modified = True
                    break
    if not modified:
        return expr
    return sp.Add(*[
        terms[key] * key[0] * cosine**key[1] * sine**key[2]
        for key in order if terms.get(key, 0) != 0
    ])


# ---------------------------------------------------------------------------
# Inputs and outputs
# ---------------------------------------------------------------------------
def _is_array(value) -> bool:
    return isinstance(value, (sp.MatrixBase, sp.NDimArray))


def _components(value, target) -> list[tuple[str | None, object]]:
    if _is_array(value):
        array = sp.Array(value) if isinstance(value, sp.MatrixBase) else value
        shape = array.shape
        if isinstance(value, sp.MatrixBase) and 1 in shape:
            shape = (max(shape),)
            array = sp.Array(list(value))
        indices = list(itertools.product(*(range(size) for size in shape)))
        values = [array[index] for index in indices]
        if isinstance(target, str):
            lvalues = [target + "".join(f"[{i}]" for i in index) for index in indices]
        else:
            lvalues = list(target)
            if len(lvalues) != len(values):
                raise ValueError(
                    f"{len(lvalues)} output lvalues given for {len(values)} components"
                )
        return list(zip(lvalues, (sp.sympify(item) for item in values)))
    if isinstance(target, str):
        return [(target, sp.sympify(value))]
    lvalues = list(target)
    if len(lvalues) != 1:
        raise ValueError("a scalar expression needs exactly one output lvalue")
    return [(lvalues[0], sp.sympify(value))]


def _flatten(expressions, outputs) -> list[tuple[str, object]]:
    if isinstance(expressions, Mapping):
        if outputs is None:
            targets = {key: str(key) for key in expressions}
        elif isinstance(outputs, Mapping):
            targets = dict(outputs)
        else:
            raise TypeError("outputs must be None or a mapping when expressions is a mapping")
        items = []
        for key, value in expressions.items():
            items.extend(_components(value, targets[key]))
    elif isinstance(expressions, sp.Basic) or _is_array(expressions) or isinstance(expressions, (int, float)):
        if outputs is None:
            raise TypeError("outputs is required unless expressions is a mapping")
        items = _components(expressions, outputs)
    elif isinstance(expressions, Sequence):
        if outputs is None:
            raise TypeError("outputs is required unless expressions is a mapping")
        if isinstance(outputs, str):
            if any(_is_array(item) for item in expressions):
                raise TypeError("a single output name needs a sequence of scalar expressions")
            items = _components(sp.Matrix(list(expressions)), outputs)
        else:
            targets = list(outputs)
            if len(targets) != len(expressions):
                raise ValueError("outputs must have one entry per expression")
            items = []
            for value, target in zip(expressions, targets):
                items.extend(_components(value, target))
    else:
        raise TypeError(f"unsupported expressions of type {type(expressions).__name__}")
    return [(lvalue, expr) for lvalue, expr in items if lvalue is not None]


@dataclass(frozen=True)
class _Angle:
    symbol: object
    sine: object
    cosine: object
    sine_c: str | None
    cosine_c: str | None


def _parse_inputs(inputs) -> dict:
    if isinstance(inputs, Mapping):
        parsed = {sp.sympify(key): value for key, value in inputs.items()}
    else:
        parsed = {sp.sympify(key): None for key in inputs}
    for symbol, value in parsed.items():
        if not isinstance(symbol, sp.Symbol):
            raise TypeError(f"input {symbol!r} is not a SymPy Symbol")
        if not IDENTIFIER.match(symbol.name) or symbol.name in C_RESERVED:
            raise ValueError(f"input name {symbol.name!r} is not a usable C identifier")
        if symbol.name.startswith(_INTERNAL_PREFIX):
            raise ValueError(f"input names must not start with {_INTERNAL_PREFIX!r}")
        if value is not None and not isinstance(value, str):
            raise TypeError(f"the C expression of input {symbol.name!r} must be a string")
    return parsed


def _parse_angles(angles, inputs: dict) -> list[_Angle]:
    parsed = []
    items = angles.items() if isinstance(angles, Mapping) else ((angle, None) for angle in angles)
    for angle, source in items:
        angle = sp.sympify(angle)
        if not isinstance(angle, sp.Symbol):
            raise TypeError(f"angle {angle!r} is not a SymPy Symbol")
        sine = sp.Symbol(f"{_INTERNAL_PREFIX}sin_{angle.name}")
        cosine = sp.Symbol(f"{_INTERNAL_PREFIX}cos_{angle.name}")
        if source is None:
            if angle not in inputs:
                raise ValueError(
                    f"angle {angle.name!r} must be an input unless its sine and cosine are given"
                )
            parsed.append(_Angle(angle, sine, cosine, None, None))
        else:
            sine_c, cosine_c = source
            parsed.append(_Angle(angle, sine, cosine, str(sine_c), str(cosine_c)))
    return parsed


# ---------------------------------------------------------------------------
# Candidate programs
# ---------------------------------------------------------------------------
class _Program:
    """Straight-line program over internal symbols (one strategy)."""

    def __init__(self, auxiliaries: dict, radicals: dict, outputs: list, protected: set):
        self.aux = dict(auxiliaries)  # symbol -> expression
        self.radicals = dict(radicals)  # sqrt symbol -> radicand (symbol or expression)
        self.outputs = list(outputs)
        self.protected = set(protected)  # inputs and angle symbols
        self.order = {symbol: index for index, symbol in enumerate(self.aux)}
        self.power_symbols: dict[tuple, object] = {}
        self.power_defs: dict = {}
        self.counter = itertools.count()
        self.base_aux: dict = {}

    # -- helpers ---------------------------------------------------------
    def _new(self, prefix: str):
        return sp.Symbol(f"{_INTERNAL_PREFIX}{prefix}{next(self.counter)}")

    def _substitute(self, mapping: dict) -> None:
        for symbol in list(self.aux):
            self.aux[symbol] = self.aux[symbol].xreplace(mapping)
        self.outputs = [expr.xreplace(mapping) for expr in self.outputs]
        for radical, radicand in list(self.radicals.items()):
            self.radicals[radical] = sp.sympify(radicand).xreplace(mapping)

    def _add_aux(self, expression):
        for symbol, value in self.aux.items():
            if value == expression and symbol not in self.radicals:
                return symbol
        symbol = self._new("b")
        self.aux[symbol] = expression
        self.order[symbol] = len(self.order) + 0.5
        return symbol

    def _power(self, base, exponent: int):
        key = (base, exponent)
        if key not in self.power_symbols:
            self.power_symbols[key] = self._new("p")
        return self.power_symbols[key]

    # -- passes ----------------------------------------------------------
    def inline_trivial(self) -> bool:
        mapping = {}
        for symbol, expr in list(self.aux.items()):
            if symbol in self.radicals:
                continue
            if expr.is_Symbol or expr.is_Number or (
                expr.is_Pow and expr.base.is_Symbol and expr.exp.is_Integer and expr.exp > 0
            ):
                mapping[symbol] = expr
        if not mapping:
            return False
        # Resolve chains such as a -> b -> c before substituting.
        for symbol in list(mapping):
            value = mapping[symbol]
            for _ in range(len(mapping)):
                replaced = value.xreplace(mapping)
                if replaced == value:
                    break
                value = replaced
            mapping[symbol] = value
        for symbol in mapping:
            del self.aux[symbol]
        self._substitute(mapping)
        return True

    def resolve_powers(self) -> bool:
        changed = False

        def resolve(expr):
            nonlocal changed
            if expr.is_Atom:
                return expr
            arguments = tuple(resolve(argument) for argument in expr.args)
            if arguments != expr.args:
                expr = expr.func(*arguments)
            if expr.is_Mul:
                factors = tuple(_resolve_pow(factor) for factor in expr.args)
                if factors != expr.args:
                    expr = sp.Mul(*factors)
                    if expr.is_Pow:
                        expr = _resolve_pow(expr)
            elif expr.is_Pow:
                expr = _resolve_pow(expr)
            return expr

        def _resolve_pow(expr):
            nonlocal changed
            if not expr.is_Pow or not expr.exp.is_Integer:
                return expr
            base, exponent = expr.args
            magnitude = abs(int(exponent))
            if magnitude < 2:
                return expr
            changed = True
            if base in self.radicals:
                pairs, odd = divmod(magnitude, 2)
                radicand = self._radicand_symbol(base)
                value = self._power(radicand, pairs) if pairs > 1 else radicand
                if odd:
                    value = value * base
            else:
                if not base.is_Symbol:
                    if base not in self.base_aux:
                        self.base_aux[base] = self._add_aux(base)
                    base = self.base_aux[base]
                value = self._power(base, magnitude)
            return value if exponent > 0 else 1 / value

        for symbol in list(self.aux):
            self.aux[symbol] = resolve(self.aux[symbol])
        self.outputs = [resolve(expr) for expr in self.outputs]
        return changed

    def _radicand_symbol(self, radical):
        radicand = self.radicals[radical]
        if radicand.is_Symbol:
            return radicand
        symbol = self._add_aux(radicand)
        self.radicals[radical] = symbol
        self.aux[radical] = sp.sqrt(symbol)
        self.order[symbol] = self.order.get(radical, 0) - 0.5
        return symbol

    def _uses(self) -> Counter:
        counts: Counter = Counter()
        expressions = list(self.aux.values()) + list(self.outputs)
        for expression in expressions:
            for node in sp.preorder_traversal(expression):
                if node.is_Symbol:
                    counts[node] += 1
        return counts

    def inline_single_use(self) -> bool:
        counts = self._uses()
        pinned = {base for base, _ in self.power_symbols} | set(self.radicals)
        pinned |= {value for value in self.radicals.values() if sp.sympify(value).is_Symbol}
        candidates = [
            symbol for symbol in self.aux
            if counts[symbol] <= 1 and symbol not in pinned
        ]
        if not candidates:
            return False
        symbol = candidates[0]
        expression = self.aux.pop(symbol)
        self._substitute({symbol: expression})
        return True

    def run(self) -> None:
        for _ in range(1000):
            changed = self.inline_trivial()
            changed |= self.resolve_powers()
            changed |= self.inline_single_use()
            if not changed:
                break
        self._build_power_chains()

    def _build_power_chains(self) -> None:
        bases: dict = {}
        for base, exponent in self.power_symbols:
            bases.setdefault(base, set()).add(exponent)
        for base in sorted(bases, key=sp.default_sort_key):
            available = {1: base}
            for exponent in sorted(bases[base]):
                self._chain(base, exponent, available)

    def _chain(self, base, exponent: int, available: dict):
        if exponent in available:
            return available[exponent]
        half = exponent // 2
        if exponent % 2 == 0 and half in available:
            left = right = available[half]
        else:
            pair = max(
                (k for k in available if k < exponent and exponent - k in available),
                default=None,
            )
            if pair is not None:
                left, right = available[pair], available[exponent - pair]
            elif exponent % 2 == 0:
                left = right = self._chain(base, half, available)
            else:
                left, right = self._chain(base, exponent - 1, available), base
        symbol = self._power(base, exponent)
        self.power_defs[symbol] = (left, right, base, exponent)
        available[exponent] = symbol
        return symbol

    # -- cost ------------------------------------------------------------
    def operation_count(self, angle_calls: int) -> OperationCount:
        counter: Counter = Counter()
        for expression in self.aux.values():
            counter.update(_count_expr(expression))
        for expression in self.outputs:
            counter.update(_count_expr(expression))
        counter["mul"] += len(self.power_defs)
        counter["call"] += angle_calls
        return OperationCount.from_counter(counter)

    @property
    def temporaries(self) -> int:
        return len(self.aux) + len(self.power_defs)


def _extract_radicals(expressions: list) -> tuple[list, dict]:
    """Replace half-integer powers by integer powers of sqrt symbols."""
    radicals: dict = {}  # radicand expression -> sqrt symbol
    counter = itertools.count()

    def visit(expr):
        if expr.is_Atom:
            return expr
        arguments = tuple(visit(argument) for argument in expr.args)
        if arguments != expr.args:
            expr = expr.func(*arguments)
        if expr.is_Pow and expr.exp.is_Rational and expr.exp.q == 2:
            base = expr.base
            if base not in radicals:
                radicals[base] = sp.Symbol(f"{_INTERNAL_PREFIX}q{next(counter)}")
            return radicals[base] ** int(expr.exp.p)
        return expr

    rewritten = [visit(expr) for expr in expressions]
    definitions = {symbol: sp.sqrt(base) for base, symbol in radicals.items()}
    return rewritten, definitions


def _polynomial_form(expr, generators, transform):
    try:
        combined = sp.together(expr)
        numerator, denominator = sp.fraction(combined)
        if denominator.has(*generators):
            return transform(expr)
        polynomial = sp.Poly(numerator, *generators)
    except sp.PolynomialError:
        return transform(expr)
    terms = []
    for monomial, coefficient in polynomial.terms():
        power = sp.Mul(*[g**k for g, k in zip(generators, monomial)])
        terms.append(transform(coefficient / denominator) * power)
    return sp.Add(*terms)


def _candidate(expressions, radical_defs, cse_mode, protected, angle_calls):
    pre_symbols = list(radical_defs)
    pre_expressions = [radical_defs[symbol] for symbol in pre_symbols]
    replacements, reduced = sp.cse(
        pre_expressions + list(expressions),
        symbols=sp.numbered_symbols(f"{_INTERNAL_PREFIX}x"),
        optimizations=CSE_MODES[cse_mode],
        order="canonical",
    )
    auxiliaries = dict(replacements)
    for symbol, expression in zip(pre_symbols, reduced[: len(pre_symbols)]):
        auxiliaries[symbol] = expression
    outputs = list(reduced[len(pre_symbols):])
    radicals = {
        symbol: auxiliaries[symbol].base for symbol in pre_symbols
        if auxiliaries[symbol].is_Pow and auxiliaries[symbol].exp == sp.S.Half
    }
    program = _Program(auxiliaries, radicals, outputs, protected)
    program.run()
    operations = program.operation_count(angle_calls)
    return program, operations


# ---------------------------------------------------------------------------
# Optimizer
# ---------------------------------------------------------------------------
def count_operations(expressions) -> OperationCount:
    """C-level operations of expressions evaluated as written."""
    counter: Counter = Counter()
    for expression in expressions:
        counter.update(_count_expr(expression))
    return OperationCount.from_counter(counter)


def optimize(
    expressions,
    *,
    inputs,
    outputs=None,
    angles=(),
    polynomial_in: Sequence = (),
    transforms: Sequence[str] = DEFAULT_TRANSFORMS,
    cse_modes: Sequence[str] = ("raw", "basic"),
    max_name_length: int = 40,
    reserved: Sequence[str] = (),
) -> OptimizedCode:
    """Optimize SymPy expressions into semantic C99 statements.

    Parameters
    ----------
    expressions
        A SymPy expression, ``Matrix`` or ``NDimArray``, a sequence of them,
        or a mapping from output names to them.
    inputs
        Input symbols.  A sequence makes each symbol a C variable of the same
        name (typically a function parameter); a mapping gives the C
        expression that reads it, emitted as ``const double <name> = <C>;``
        when the symbol is used.
    outputs
        Output lvalues.  A string ``"du"`` writes a scalar to ``du`` and an
        array to ``du[i]`` (``du[i][j]`` for rank two, and so on); a sequence
        gives one base or lvalue per expression, or one lvalue per component
        of a single array, ``None`` skipping a component.  With a mapping of
        expressions, ``None`` uses its keys and a mapping renames them.
    angles
        Angles whose trigonometric functions are reduced to one sine and one
        cosine.  A sequence of input symbols evaluates ``sin`` and ``cos``
        once; a mapping ``{angle: (sin_c, cos_c)}`` reads precomputed values
        from those C expressions.
    polynomial_in
        Symbols (for instance velocity components) in which the outputs are
        polynomials.  When given, the search adds strategies that collect
        each output by monomials of these symbols and transform each
        coefficient separately.
    transforms
        Pre-transforms tried by the search: names from ``"none"``,
        ``"factor_terms"``, ``"cancel"`` and ``"factor"``.
    cse_modes
        SymPy CSE optimization modes tried by the search: ``"raw"`` and
        ``"basic"``.
    max_name_length
        Longest structural temporary name before the category fallback.
    reserved
        Extra identifiers the generated temporaries must avoid.

    Returns
    -------
    OptimizedCode
        Statements, outputs and report of the selected strategy.

    Raises
    ------
    ValueError
        If an output depends on a symbol that is neither an input nor a
        declared angle, or if a name is not a usable C identifier.
    """
    targets = _flatten(expressions, outputs)
    input_map = _parse_inputs(inputs)
    angle_list = _parse_angles(angles, input_map)
    unknown_transforms = set(transforms) - set(TRANSFORMS)
    if unknown_transforms:
        raise ValueError(f"unknown transforms: {sorted(unknown_transforms)}")
    unknown_modes = set(cse_modes) - set(CSE_MODES)
    if unknown_modes:
        raise ValueError(f"unknown CSE modes: {sorted(unknown_modes)}")
    lvalues = [lvalue for lvalue, _ in targets]
    originals = tuple(expr for _, expr in targets)
    before = count_operations(originals)

    reduced = list(originals)
    for angle in angle_list:
        reduced = [reduce_angle(expr, angle.symbol, angle.sine, angle.cosine) for expr in reduced]
    allowed = set(input_map) | {a.sine for a in angle_list} | {a.cosine for a in angle_list}
    free = set().union(*(expr.free_symbols for expr in reduced)) if reduced else set()
    missing = sorted(symbol.name for symbol in free - allowed)
    if missing:
        raise ValueError(f"outputs depend on undeclared symbols: {', '.join(missing)}")
    protected = set(allowed)
    angle_calls_by_symbol = {
        a.sine: int(a.sine_c is None) for a in angle_list
    } | {a.cosine: int(a.cosine_c is None) for a in angle_list}

    generators = tuple(sp.sympify(symbol) for symbol in polynomial_in)
    contractions = [(a.sine, a.cosine) for a in angle_list]

    def contract(expr):
        for sine, cosine in contractions:
            expr = contract_pythagorean(expr, sine, cosine)
        return expr

    variants: list[tuple[str, list]] = []
    seen: list[list] = []
    for transform_name in transforms:
        transform = TRANSFORMS[transform_name]
        forms = [(transform_name, lambda e, t=transform: t(e))]
        if generators:
            forms.append((f"poly[{transform_name}]",
                          lambda e, t=transform: _polynomial_form(e, generators, t)))
        for form_name, function in forms:
            transformed = [contract(function(expr)) for expr in reduced]
            if any(transformed == previous for previous in seen):
                continue
            seen.append(transformed)
            variants.append((form_name, transformed))

    candidates = []
    best = None
    for form_name, transformed in variants:
        rewritten, radical_defs = _extract_radicals(transformed)
        for mode in cse_modes:
            name = f"{form_name}+cse:{mode}"
            used = set().union(*(expr.free_symbols for expr in rewritten)) if rewritten else set()
            used |= set().union(*(expr.free_symbols for expr in radical_defs.values())) if radical_defs else set()
            angle_calls = sum(count for symbol, count in angle_calls_by_symbol.items() if symbol in used)
            program, operations = _candidate(rewritten, radical_defs, mode, protected, angle_calls)
            candidates.append(CandidateSummary(name, operations.total, program.temporaries))
            key = (operations.total, program.temporaries, name)
            if best is None or key < best[0]:
                best = (key, name, program, operations)

    _, strategy, program, operations = best
    report = OptimizationReport(
        strategy=strategy,
        before=before,
        after=operations,
        temporaries=program.temporaries,
        candidates=tuple(candidates),
    )
    return _emit(program, report, lvalues, originals, input_map, angle_list,
                 max_name_length, reserved)


# ---------------------------------------------------------------------------
# Emission
# ---------------------------------------------------------------------------
def _emit(program, report, lvalues, originals, input_map, angle_list, max_name_length, reserved):
    used: set = set()
    for expression in list(program.aux.values()) + program.outputs:
        used |= expression.free_symbols
    for left, right, base, _ in program.power_defs.values():
        used |= {left, right, base}

    reserved_names = set(reserved)
    for lvalue in lvalues:
        reserved_names |= set(_identifiers(lvalue))
    for value in input_map.values():
        if value is not None:
            reserved_names |= set(_identifiers(value))
    for angle in angle_list:
        for value in (angle.sine_c, angle.cosine_c):
            if value is not None:
                reserved_names |= set(_identifiers(value))
    direct_angle_names = set()
    for angle in angle_list:
        for function, source in (("sin", angle.sine_c), ("cos", angle.cosine_c)):
            if source == f"{function}_{angle.symbol.name}":
                direct_angle_names.add(source)
    namer = Namer(
        max_length=max_name_length,
        reserved=reserved_names - {s.name for s in input_map} - direct_angle_names,
    )

    angle_inputs_needed = {a.symbol for a in angle_list if a.sine_c is None and (a.sine in used or a.cosine in used)}
    definitions: list[Definition] = []
    final: dict = {}
    expansion: dict = {}

    for symbol, source in input_map.items():
        namer.fix(symbol, symbol.name)
        final[symbol] = sp.Symbol(symbol.name)
        namer.names[final[symbol]] = symbol.name
        expansion[final[symbol]] = symbol
    for symbol, source in input_map.items():
        if source is None or source == symbol.name:
            continue
        if symbol in used or symbol in angle_inputs_needed:
            definitions.append(Definition(symbol.name, "input", symbol, source, "Inputs."))
    for angle in angle_list:
        namer.register_angle(angle.sine, angle.cosine, angle.symbol.name)
        for symbol in (angle.sine, angle.cosine):
            twin = sp.Symbol(namer.names[symbol])
            namer.names[twin] = namer.names[symbol]
            namer.angle_stems[twin] = namer.angle_stems[symbol]
        for symbol, source, function in (
            (angle.sine, angle.sine_c, "sin"), (angle.cosine, angle.cosine_c, "cos"),
        ):
            name = namer.names[symbol]
            final[symbol] = sp.Symbol(name)
            expansion[final[symbol]] = getattr(sp, function)(angle.symbol)
            if symbol not in used:
                continue
            text = source if source is not None else f"{function}({angle.symbol.name})"
            definitions.append(Definition(
                name, "angle", getattr(sp, function)(angle.symbol), text,
                f"Sine and cosine of {angle.symbol.name}.",
            ))

    # Dependency order: powers of inputs first, then auxiliaries in CSE order.
    nodes: dict = {}
    for symbol, expression in program.aux.items():
        nodes[symbol] = ("auxiliary", expression, set(expression.free_symbols))
    for symbol, (left, right, base, exponent) in program.power_defs.items():
        nodes[symbol] = ("power", (left, right, base, exponent), {left, right})
    protected = set(final)
    order_index = dict(program.order)
    for symbol, (_, _, base, exponent) in program.power_defs.items():
        if base in protected:
            order_index[symbol] = (-1.0, sp.default_sort_key(base), exponent)
        else:
            order_index[symbol] = (order_index.get(base, 0.0), "", exponent)

    def priority(symbol):
        value = order_index.get(symbol, float(len(order_index)))
        if isinstance(value, tuple):
            return value
        return (float(value), "", 0)

    dependencies = {
        symbol: {d for d in deps if d in nodes and d != symbol}
        for symbol, (_, _, deps) in nodes.items()
    }
    dependents: dict = {symbol: set() for symbol in nodes}
    for symbol, deps in dependencies.items():
        for dependency in deps:
            dependents[dependency].add(symbol)
    remaining = {symbol: len(deps) for symbol, deps in dependencies.items()}
    heap = [(priority(s), str(s), s) for s, n in remaining.items() if n == 0]
    heapq.heapify(heap)
    ordered = []
    while heap:
        _, _, symbol = heapq.heappop(heap)
        ordered.append(symbol)
        for dependent in dependents[symbol]:
            remaining[dependent] -= 1
            if remaining[dependent] == 0:
                heapq.heappush(heap, (priority(dependent), str(dependent), dependent))
    if len(ordered) != len(nodes):  # pragma: no cover - defensive
        raise RuntimeError("cyclic definitions in generated program")

    group_of_base: dict = {}
    for symbol in ordered:
        kind, payload, _ = nodes[symbol]
        if kind == "power":
            left, right, base, exponent = payload
            name = namer.assign(symbol, None, power_of=(base, exponent))
            final[symbol] = sp.Symbol(name)
            namer.names[final[symbol]] = name
            expression = final[left] * final[right]
            text = f"{final[left].name} * {final[right].name}"
            if base in protected:
                group = "Integer powers of the inputs."
            else:
                group = group_of_base.get(base, "Common subexpressions.")
            definitions.append(Definition(name, "power", expression, text, group))
            expansion[final[symbol]] = sp.sympify(expansion.get(final[base], final[base])) ** exponent
        else:
            expression = payload.xreplace(final)
            name = namer.assign(symbol, expression)
            final[symbol] = sp.Symbol(name)
            namer.names[final[symbol]] = name
            group = "Common subexpressions."
            group_of_base[symbol] = group
            definitions.append(Definition(name, "auxiliary", expression, c_expression(expression), group))
            expansion[final[symbol]] = expression.xreplace(expansion)

    outputs = []
    for lvalue, expression in zip(lvalues, program.outputs):
        expression = sp.sympify(expression).xreplace(final)
        outputs.append((lvalue, expression, c_expression(expression)))

    input_symbols = tuple(input_map) + tuple(a.symbol for a in angle_list if a.symbol not in input_map)
    return OptimizedCode(
        definitions=tuple(definitions),
        outputs=tuple(outputs),
        report=report,
        originals=tuple(originals),
        _expansion=expansion,
        _input_symbols=input_symbols,
        _angles=tuple(a.symbol for a in angle_list),
    )


def _identifiers(text: str) -> list[str]:
    import re

    return re.findall(r"[A-Za-z_][A-Za-z0-9_]*", text)
