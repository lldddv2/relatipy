"""Verification utilities for :func:`relatipy._codegen.optimize` results.

Three independent checks compare an :class:`OptimizedCode` with the
expressions it was built from:

- :func:`symbolic_residuals` proves equality symbolically when the
  difference is small enough to reduce;
- :func:`numeric_max_relative_error` evaluates the statements with NumPy at
  random points;
- :func:`compile_and_compare` compiles the statements with the system C
  compiler through :mod:`relatipy._jit` and calls them through ``ctypes``.

:func:`count_c_operations` counts the arithmetic operators of any C
fragment, so generated and hand-written kernels can be compared on the same
scale.
"""

from __future__ import annotations

import ctypes
import re
from collections import Counter
from collections.abc import Mapping

import numpy as np

from ._sympy import require_sympy

sp = require_sympy()

from .optimizer import OperationCount, OptimizedCode, reduce_angle  # noqa: E402

__all__ = [
    "compile_and_compare",
    "count_c_operations",
    "numeric_max_relative_error",
    "symbolic_residuals",
]


def symbolic_residuals(code: OptimizedCode, *, max_operations: int = 20000) -> list | None:
    """Symbolic differences between the original and optimized outputs.

    Each difference is written with ``S = sin`` and ``C = cos`` of every
    declared angle, put over a common denominator, and its numerator is
    reduced modulo ``S^2 + C^2 - 1``.

    Parameters
    ----------
    code
        Result of :func:`optimize`.
    max_operations
        Size bound (SymPy ``count_ops``) of a difference above which the
        check is declared intractable.

    Returns
    -------
    list or None
        One residual per output (``0`` when proven equal), or ``None`` when
        a difference exceeds ``max_operations``.
    """
    residuals = []
    trig = [
        (angle, sp.Symbol(f"__check_sin_{angle.name}"), sp.Symbol(f"__check_cos_{angle.name}"))
        for angle in code.angles
    ]
    for original, optimized in zip(code.originals, code.expanded()):
        difference = sp.sympify(original) - sp.sympify(optimized)
        if sp.count_ops(difference) > max_operations:
            return None
        for angle, sine, cosine in trig:
            difference = reduce_angle(difference, angle, sine, cosine)
        numerator, _ = sp.fraction(sp.together(difference))
        numerator = sp.expand(numerator)
        for _, sine, cosine in trig:
            if numerator.has(sine):
                numerator = sp.rem(
                    sp.Poly(numerator, sine), sp.Poly(sine**2 + cosine**2 - 1, sine)
                ).as_expr()
                numerator = sp.expand(numerator)
        residuals.append(sp.simplify(numerator) if numerator != 0 else sp.S.Zero)
    return residuals


def _sample_points(code: OptimizedCode, domain: Mapping, samples: int, seed: int) -> dict:
    generator = np.random.default_rng(seed)
    points = {}
    for symbol in code.input_symbols:
        low, high = domain.get(symbol, (0.5, 1.5))
        points[symbol] = generator.uniform(low, high, samples)
    return points


def _reference(code: OptimizedCode, points: dict) -> np.ndarray:
    symbols = list(points)
    values = [points[symbol] for symbol in symbols]
    rows = []
    for original in code.originals:
        function = sp.lambdify(symbols, original, modules="numpy")
        rows.append(np.broadcast_to(np.asarray(function(*values), dtype=float), values[0].shape))
    return np.array(rows)


def _relative_error(actual: np.ndarray, expected: np.ndarray) -> float:
    scale = np.max(np.abs(expected), axis=0)
    scale = np.where(scale > 0.0, scale, 1.0)
    return float(np.max(np.abs(actual - expected) / scale))


def numeric_max_relative_error(
    code: OptimizedCode, *, domain: Mapping | None = None, samples: int = 64, seed: int = 0,
) -> float:
    """Largest error of :meth:`OptimizedCode.evaluate` at random points.

    The error of a point is ``max_k |optimized_k - original_k|`` divided by
    ``max_k |original_k|`` (the infinity norm of the outputs at that point).

    Parameters
    ----------
    code
        Result of :func:`optimize`.
    domain
        ``{symbol: (low, high)}`` uniform sampling ranges; ``(0.5, 1.5)``
        for unlisted inputs.
    samples
        Number of random points.
    seed
        Seed of :func:`numpy.random.default_rng`.
    """
    points = _sample_points(code, domain or {}, samples, seed)
    expected = _reference(code, points)
    actual = np.array([
        np.broadcast_to(value, (samples,)) for value in code.evaluate(points)
    ])
    return _relative_error(actual, expected)


def check_source(code: OptimizedCode, name: str = "rp_codegen_check") -> str:
    """C translation unit ``void name(const double *in, double *out)``.

    Inputs are read from ``in`` in :attr:`OptimizedCode.input_symbols`
    order; outputs are written to ``out`` in output order.  Angles whose
    sine and cosine are read from C expressions in the optimized code are
    passed as angles and evaluated with ``sin`` and ``cos``.
    """
    lines = ["#include <math.h>", "", f"void {name}(const double *in, double *out)", "{"]
    for index, symbol in enumerate(code.input_symbols):
        lines.append(f"    const double {symbol.name} = in[{index}];")
        lines.append(f"    (void){symbol.name};")
    for definition in code.definitions:
        if definition.kind == "input":
            continue
        if definition.kind == "angle":
            function = "sin" if definition.expression.func == sp.sin else "cos"
            angle = definition.expression.args[0].name
            lines.append(f"    const double {definition.name} = {function}({angle});")
            continue
        lines.append(f"    const double {definition.name} = {definition.c_text};")
    for index, (_, _, text) in enumerate(code.outputs):
        lines.append(f"    out[{index}] = {text};")
    lines.append("}")
    return "\n".join(lines) + "\n"


def compile_and_compare(
    code: OptimizedCode, *, domain: Mapping | None = None, samples: int = 64, seed: int = 0,
    cache_dir=None,
) -> float:
    """Compile the statements, call them through ctypes and compare.

    The statements are wrapped by :func:`check_source`, compiled with
    :func:`relatipy._jit.compile_function` and evaluated at the random
    points of :func:`numeric_max_relative_error`; the return value is the
    same error measure against the original expressions evaluated with
    NumPy.

    Raises
    ------
    relatipy._jit.CompilerNotFoundError
        If no C compiler is available.
    """
    from .. import _jit

    function = _jit.compile_function(check_source(code), "rp_codegen_check", cache_dir=cache_dir)
    pointer = ctypes.POINTER(ctypes.c_double)
    call = function.as_ctypes(None, pointer, pointer)
    points = _sample_points(code, domain or {}, samples, seed)
    expected = _reference(code, points)
    actual = np.empty_like(expected)
    inputs = np.ascontiguousarray(np.array([points[s] for s in code.input_symbols]).T)
    out = np.empty(len(code.outputs))
    for index in range(samples):
        row = np.ascontiguousarray(inputs[index])
        call(row.ctypes.data_as(pointer), out.ctypes.data_as(pointer))
        actual[:, index] = out
    return _relative_error(actual, expected)


_COMMENT = re.compile(r"/\*.*?\*/|//[^\n]*", re.S)
_TOKEN = re.compile(
    r"\d+\.\d*(?:[eE][-+]?\d+)?|\.\d+(?:[eE][-+]?\d+)?|\d+(?:[eE][-+]?\d+)?"
    r"|->|[A-Za-z_][A-Za-z0-9_]*|[-+*/()\[\],;=]|\S"
)
_TYPE_WORDS = {"const", "double", "float", "int", "static", "void", "sizeof", "return"}


def count_c_operations(text: str) -> OperationCount:
    """Count arithmetic operators on the right-hand sides of C statements.

    Every binary ``+``, ``-``, ``*`` and ``/`` and every function call in the
    expression after ``=`` of each statement is counted; unary signs, array
    indices and member access are not.

    Parameters
    ----------
    text
        C statements or a complete function.
    """
    counter: Counter = Counter()
    text = _COMMENT.sub(" ", text)
    for statement in text.split(";"):
        if "=" not in statement:
            continue
        right = statement.split("=", 1)[1]
        tokens = _TOKEN.findall(right)
        previous = None
        depth = 0
        for index, token in enumerate(tokens):
            if token == "[":
                depth += 1
            elif token == "]":
                depth -= 1
            elif depth == 0 and token in {"+", "-", "*", "/"}:
                binary = previous is not None and (
                    previous in {")", "]"} or re.match(r"[A-Za-z0-9_.]", previous)
                )
                if binary:
                    counter[{"+": "add", "-": "add", "*": "mul", "/": "div"}[token]] += 1
            elif (
                depth == 0 and re.match(r"[A-Za-z_]", token) and token not in _TYPE_WORDS
                and index + 1 < len(tokens) and tokens[index + 1] == "("
            ):
                counter["call"] += 1
            if token != "->":
                previous = token
            else:
                previous = None
        # member access "a->b" leaves previous at "b", which is an operand
    return OperationCount.from_counter(counter)
