"""C99 printing, line wrapping and operation counting for generated code."""

from __future__ import annotations

import math
import re
from collections import Counter

from ._sympy import require_sympy

sp = require_sympy()
from sympy.printing.c import C99CodePrinter  # noqa: E402

LINE_WIDTH = 79

_BINARY_SPACING = re.compile(r"\s*([*/])\s*")
_PARENTHESIZED_LITERAL = re.compile(r"(?<![A-Za-z0-9_])\((\d+\.\d*(?:e[-+]?\d+)?)\)")
_DOUBLE_PARENTHESIZED_QUOTIENT = re.compile(r"\(\((\d+\.0 / \d+\.0)\)\)")


class GeneratedCPrinter(C99CodePrinter):
    """C99 printer with double literals and multiplication-only powers.

    Integers print as ``2.0``; dyadic rationals print as exact decimals
    (``0.5``) and other rationals as a constant quotient ``(1.0/3.0)`` that
    the compiler folds.  Integer powers that survive optimization print as
    explicit products, never as ``pow``.
    """

    def _print_Pi(self, expr):  # noqa: N802 - strict C11 has no M_PI
        return repr(math.pi)

    def _print_Exp1(self, expr):  # noqa: N802
        return repr(math.e)

    def _print_Integer(self, expr):  # noqa: N802 - SymPy printer protocol
        return f"{int(expr)}.0"

    def _print_Rational(self, expr):  # noqa: N802
        p, q = int(expr.p), int(expr.q)
        if q & (q - 1) == 0 and abs(p) < 2**53:
            return repr(p / q)
        return f"({p}.0/{q}.0)"

    def _print_Pow(self, expr):  # noqa: N802
        base, exponent = expr.args
        if exponent == -1:
            return f"1.0/{self.parenthesize(base, 50)}"
        if exponent == sp.S.Half:
            return f"sqrt({self._print(base)})"
        if exponent == -sp.S.Half:
            return f"1.0/sqrt({self._print(base)})"
        if exponent.is_Integer and int(exponent) > 0:
            factor = self.parenthesize(base, 50)
            return "(" + "*".join([factor] * int(exponent)) + ")"
        if exponent.is_Integer:
            factor = self.parenthesize(base, 50)
            return "1.0/(" + "*".join([factor] * -int(exponent)) + ")"
        return super()._print_Pow(expr)


_PRINTER = GeneratedCPrinter()


def c_expression(expr) -> str:
    """Print ``expr`` as a C99 expression with spaced binary operators."""
    text = _PARENTHESIZED_LITERAL.sub(r"\1", _PRINTER.doprint(expr))
    text = _BINARY_SPACING.sub(r" \1 ", text)
    return _DOUBLE_PARENTHESIZED_QUOTIENT.sub(r"(\1)", text)


def wrap_statement(head: str, text: str, indent: str, width: int = LINE_WIDTH) -> list[str]:
    """Wrap ``head + text + ';'`` at token boundaries.

    The printer separates tokens with spaces only between operators and
    operands, so a break at any space is valid C.  An overflowing line is
    broken before its last ``+`` or ``-`` when it has one; continuation lines
    get four extra spaces of indentation.
    """
    first = f"{indent}{head}{text};"
    if len(first) <= width:
        return [first]
    continuation = indent + "    "
    lines: list[str] = []
    prefix = f"{indent}{head}".rstrip()
    line: list[str] = []

    def render(tokens: list[str], lead: str) -> str:
        return lead + " ".join(tokens) if lead.strip() == "" else lead + " " + " ".join(tokens)

    lead = prefix
    for token in (text + ";").split(" "):
        if not line or len(render(line + [token], lead)) <= width:
            if not line and lead == prefix and len(render([token], lead)) > width:
                lines.append(prefix)
                lead = continuation
            line.append(token)
            continue
        split = max(
            (index for index, item in enumerate(line) if item in {"+", "-"} and index > 0),
            default=None,
        )
        if split is None or len(render(line[split:] + [token], continuation)) > width:
            split = len(line)
        lines.append(render(line[:split], lead))
        lead = continuation
        line = line[split:] + [token]
    lines.append(render(line, lead))
    return [item.rstrip() for item in lines]


def multiplication_chain_cost(exponent: int) -> int:
    """Multiplications of the binary method for a positive integer power."""
    exponent = abs(int(exponent))
    if exponent <= 1:
        return 0
    return exponent.bit_length() - 1 + bin(exponent).count("1") - 1


def count_ops(expr) -> Counter:
    """Count C-level operations of one expression.

    Additions and subtractions count as ``add``; products, numeric scaling
    and power chains as ``mul``; quotients as ``div``; ``sqrt``, ``pow`` and
    every other function application as ``call``.  Negations absorbed by an
    enclosing sum are free; a standalone negation counts as one ``add``.
    """
    counter: Counter = Counter()
    _count(sp.sympify(expr), counter, in_sum=False)
    return counter


def _count(expr, counter: Counter, in_sum: bool) -> None:
    if expr.is_Atom:
        return
    if expr.is_Add:
        counter["add"] += len(expr.args) - 1
        for term in expr.args:
            _count(term, counter, in_sum=True)
        return
    if expr.is_Mul:
        _count_mul(expr, counter, in_sum)
        return
    if expr.is_Pow:
        _count_pow(expr, counter)
        return
    if isinstance(expr, sp.Function):
        counter["call"] += 1
    for argument in expr.args:
        _count(argument, counter, in_sum=False)


def _count_mul(expr, counter: Counter, in_sum: bool) -> None:
    coefficient, factors = expr.as_coeff_mul()
    numerator = []
    denominator = []
    for factor in factors:
        if factor.is_Pow and factor.exp.is_Rational and factor.exp.is_negative:
            denominator.append(sp.Pow(factor.base, -factor.exp, evaluate=False))
        else:
            numerator.append(factor)
    scaled = abs(coefficient) != 1
    numerator_terms = len(numerator) + int(scaled)
    counter["mul"] += max(numerator_terms - 1, 0)
    if denominator:
        counter["div"] += 1
        counter["mul"] += len(denominator) - 1
    if coefficient == -1 and not in_sum:
        counter["add"] += 1
    for factor in numerator + denominator:
        if factor.is_Pow and factor.exp == 1:
            _count(factor.base, counter, in_sum=False)
        else:
            _count(factor, counter, in_sum=False)


def _count_pow(expr, counter: Counter) -> None:
    base, exponent = expr.args
    _count(base, counter, in_sum=False)
    if exponent.is_Integer:
        counter["mul"] += multiplication_chain_cost(int(exponent))
        if exponent.is_negative:
            counter["div"] += 1
        return
    if exponent.is_Rational and exponent.q == 2:
        counter["call"] += 1
        whole = (abs(int(exponent.p)) - 1) // 2
        if whole:
            counter["mul"] += multiplication_chain_cost(whole) + 1
        if exponent.is_negative:
            counter["div"] += 1
        return
    counter["call"] += 1
    _count(exponent, counter, in_sum=False)
