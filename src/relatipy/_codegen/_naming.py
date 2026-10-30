"""Deterministic semantic names for generated C temporaries.

Names are derived from the structure of a definition written in terms of
already named symbols: ``r * r`` gives ``r2``, ``spin**2 * r**2`` gives
``r2_spin2``, ``1 / sigma**2`` gives ``inv_sigma2``, ``r / sigma`` gives
``r_over_sigma`` and ``r**2 + spin**2`` gives ``r2_plus_spin2``.  A name
longer than the configured bound falls back to a readable category name such
as ``sum_r4`` or ``ratio_r``; collisions get a numeric suffix (``_2``,
``_3``, ...) in definition order, so the result depends only on the
expressions.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from ._sympy import require_sympy

sp = require_sympy()

IDENTIFIER = re.compile(r"[A-Za-z_][A-Za-z0-9_]*\Z")

C_RESERVED = frozenset(
    """
    auto break case char const continue default do double else enum extern
    float for goto if inline int long register restrict return short signed
    sizeof static struct switch typedef union unsigned void volatile while
    _Bool _Complex _Imaginary _Alignas _Alignof _Atomic _Generic _Noreturn
    _Static_assert _Thread_local
    acos acosh asin asinh atan atan2 atanh cbrt ceil copysign cos cosh erf
    erfc exp exp2 expm1 fabs fdim floor fma fmax fmin fmod frexp hypot ilogb
    ldexp lgamma log log10 log1p log2 logb lrint lround modf nan nearbyint
    nextafter pow remainder rint round scalbn sin sinh sqrt tan tanh tgamma
    trunc isfinite isinf isnan signbit fpclassify
    NULL M_PI INFINITY NAN HUGE_VAL main
    """.split()
)

_NUMBER_WORDS = {
    1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six",
    7: "seven", 8: "eight", 9: "nine", 10: "ten",
}
_FRACTION_WORDS = {2: "half", 3: "third", 4: "quarter"}


def number_name(value) -> str:
    """Readable identifier fragment of a rational or float constant."""
    value = sp.nsimplify(value) if isinstance(value, sp.Float) else sp.Rational(value)
    sign = "neg_" if value < 0 else ""
    value = abs(value)
    if value.is_Integer:
        integer = int(value)
        return sign + _NUMBER_WORDS.get(integer, f"c{integer}")
    if isinstance(value, sp.Rational):
        p, q = int(value.p), int(value.q)
        if p == 1 and q in _FRACTION_WORDS:
            return sign + _FRACTION_WORDS[q]
        return sign + f"c{p}_{q}"
    return sign + "c"


@dataclass
class Namer:
    """Assign unique semantic names in definition order.

    Parameters
    ----------
    max_length
        Upper bound of a structural name before the category fallback.
    reserved
        Identifiers that must not be produced (C keywords, math functions,
        caller parameters, output arrays).
    """

    max_length: int = 40
    reserved: set[str] = field(default_factory=set)
    names: dict = field(default_factory=dict)
    angle_stems: dict = field(default_factory=dict)
    used: set[str] = field(default_factory=set)

    def __post_init__(self) -> None:
        self.used |= set(C_RESERVED) | set(self.reserved)

    # -- registration ---------------------------------------------------
    def fix(self, symbol, name: str) -> str:
        """Register a symbol whose name is imposed (inputs and angles)."""
        if not IDENTIFIER.match(name):
            raise ValueError(f"{name!r} is not a valid C identifier")
        if name in C_RESERVED:
            raise ValueError(f"{name!r} is a reserved C identifier")
        self.names[symbol] = name
        self.used.add(name)
        return name

    def register_angle(self, sine, cosine, angle_name: str) -> None:
        """Register ``sin_<angle>`` and ``cos_<angle>`` with power stems."""
        self.fix(sine, self.unique(f"sin_{angle_name}"))
        self.fix(cosine, self.unique(f"cos_{angle_name}"))
        self.angle_stems[sine] = ("sin", angle_name)
        self.angle_stems[cosine] = ("cos", angle_name)

    def unique(self, candidate: str) -> str:
        if candidate[0].isdigit():
            candidate = "v" + candidate
        if candidate not in self.used:
            return candidate
        index = 2
        while f"{candidate}_{index}" in self.used:
            index += 1
        return f"{candidate}_{index}"

    def assign(self, symbol, expression, *, power_of=None) -> str:
        """Name ``symbol`` from its definition and register the name."""
        if power_of is not None:
            base, exponent = power_of
            candidate = self.power_name(base, exponent)
        else:
            candidate = self.structural(expression)
            if candidate is None or len(candidate) > self.max_length:
                candidate = self.fallback(expression)
        candidate = self.unique(candidate)
        self.names[symbol] = candidate
        self.used.add(candidate)
        return candidate

    # -- naming rules ---------------------------------------------------
    def power_name(self, base, exponent: int) -> str:
        if base in self.angle_stems:
            stem, angle = self.angle_stems[base]
            return f"{stem}{exponent}_{angle}"
        name = self.names[base]
        if "_" not in name and not name[-1].isdigit():
            return f"{name}{exponent}"
        if "_" not in name:
            return f"{name}_{exponent}"
        suffix = {2: "sq", 3: "cubed"}.get(exponent, f"pow{exponent}")
        return f"{name}_{suffix}"

    def structural(self, expr) -> str | None:
        expr = sp.sympify(expr)
        if expr.is_Symbol:
            return self.names.get(expr, None)
        if expr.is_Number:
            return number_name(expr)
        angle = self._angle_pattern(expr)
        if angle is not None:
            return angle
        if expr.is_Pow:
            return self._pow_name(expr)
        if expr.is_Mul:
            return self._mul_name(expr)
        if expr.is_Add:
            return self._add_name(expr)
        if isinstance(expr, sp.Function):
            parts = [self.structural(argument) for argument in expr.args]
            if any(part is None for part in parts):
                return None
            return "_".join([type(expr).__name__.lower(), *parts])
        return None

    def _angle_pattern(self, expr) -> str | None:
        sines = {stem[1]: symbol for symbol, stem in self.angle_stems.items() if stem[0] == "sin"}
        cosines = {stem[1]: symbol for symbol, stem in self.angle_stems.items() if stem[0] == "cos"}
        for angle, sine in sines.items():
            cosine = cosines[angle]
            if expr == 2 * sine * cosine:
                return f"sin_2{angle}"
            if expr == cosine**2 - sine**2:
                return f"cos_2{angle}"
            if expr == 8 * sine**2 * cosine**2:
                return f"one_minus_cos_4{angle}"
            if expr == cosine / sine:
                return f"cot_{angle}"
            if expr == sine / cosine:
                return f"tan_{angle}"
        return None

    def _pow_name(self, expr) -> str | None:
        base, exponent = expr.args
        base_name = self.structural(base)
        if base_name is None:
            return None
        if exponent == sp.S.Half:
            return f"sqrt_{base_name}"
        if exponent == -sp.S.Half:
            return f"inv_sqrt_{base_name}"
        if exponent == -1:
            return f"inv_{base_name}"
        if exponent.is_Integer and base.is_Symbol and base in self.names:
            power = self.power_name(base, abs(int(exponent)))
            return power if exponent > 0 else f"inv_{power}"
        return None

    def _factor_names(self, factors) -> list[str] | None:
        parts = []
        for factor in factors:
            part = self.structural(factor)
            if part is None:
                return None
            parts.append(part)
        return parts

    def _mul_name(self, expr) -> str | None:
        coefficient, factors = expr.as_coeff_mul()
        numerator = []
        denominator = []
        for factor in factors:
            if factor.is_Pow and factor.exp.is_Rational and factor.exp.is_negative:
                denominator.append(factor.base ** -factor.exp)
            else:
                numerator.append(factor)
        top = self._factor_names(numerator)
        bottom = self._factor_names(denominator)
        if top is None or bottom is None:
            return None
        prefix = []
        if coefficient != 1:
            prefix = [number_name(coefficient)] if coefficient != -1 else ["neg"]
        if not top and bottom:
            return "_".join(prefix + ["inv"] + bottom)
        name = "_".join(prefix + top)
        if bottom:
            name += "_over_" + "_".join(bottom)
        return name

    def _add_name(self, expr) -> str | None:
        name = ""
        for index, term in enumerate(expr.as_ordered_terms()):
            negative = term.could_extract_minus_sign()
            magnitude = -term if negative else term
            part = self.structural(magnitude)
            if part is None:
                return None
            if index == 0:
                name = ("neg_" if negative else "") + part
            else:
                name += ("_minus_" if negative else "_plus_") + part
        return name

    def fallback(self, expr) -> str:
        """Category name of an expression whose structural name is too long."""
        if expr.is_Add:
            kind = "sum"
            terms = expr.as_ordered_terms()
            lead_expr = terms[0]
        elif expr.is_Mul:
            numerator, denominator = sp.fraction(expr)
            kind = "ratio" if denominator != 1 else "prod"
            lead_expr = numerator if numerator != 1 else denominator
        elif expr.is_Pow:
            if expr.exp == sp.S.Half:
                kind = "sqrt"
            elif expr.exp == -sp.S.Half:
                kind = "inv_sqrt"
            else:
                kind = "inv" if expr.exp.is_negative else "pow"
            lead_expr = expr.base
        elif isinstance(expr, sp.Function):
            kind = type(expr).__name__.lower()
            lead_expr = expr.args[0]
        else:
            kind = "value"
            lead_expr = expr
        lead = self._lead(lead_expr)
        if not lead:
            return kind
        tokens = lead.split("_")
        name = kind
        for token in tokens:
            if len(name) + 1 + len(token) > max(self.max_length // 2, len(kind) + 2):
                break
            name += "_" + token
        return name if name != kind else f"{kind}_{tokens[0]}"

    def _lead(self, expr) -> str:
        """Short fragment naming the leading symbol content of ``expr``."""
        part = self.structural(expr)
        if part is not None and len(part) <= self.max_length // 2:
            return part
        while expr.is_Add or expr.is_Mul:
            terms = expr.as_ordered_terms() if expr.is_Add else expr.as_ordered_factors()
            terms = [item for item in terms if not item.is_Number] or terms
            expr = terms[0]
            part = self.structural(expr)
            if part is not None and len(part) <= self.max_length // 2:
                return part
        symbols = sorted(self.names[s] for s in expr.free_symbols if s in self.names)
        return "_".join(symbols[:2])
