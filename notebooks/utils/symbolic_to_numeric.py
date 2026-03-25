"""
Utilities to turn symbolic metric data into numeric-ready Python snippets.

This module loads text templates for metric and coordinate code generation,
renames SymPy symbols into valid Python identifiers, optimizes Christoffel
symbol tensors via substitution and symmetry, and fills ``$placeholder$``
templates with generated code fragments.

Notes
-----
The Christoffel optimization keeps only independent components with
:math:`\\nu \\leq \\sigma` in the last two indices and mirrors the rest from
symmetry. Repeated sub-expressions may be factored into auxiliary symbols.

Typical usage is from a notebook: pass a symbolic ``metric`` to
``export_christoffel_symbols_to_code`` and paste the returned string into a
numeric metric class.
"""
import os
import re
from itertools import product

import sympy as sp

this_dir = os.path.dirname(__file__)
template_path_metric = os.path.join(this_dir, "metric_template.txt")
template_path_coordinates = os.path.join(this_dir, "coordinate_template.txt")

with open(template_path_metric, encoding="utf-8") as _f:
    template_metric = _f.read()
with open(template_path_coordinates, encoding="utf-8") as _f:
    template_coordinates = _f.read()


def make_valid_varname(name):
    """
    Turn a LaTeX-like or symbolic string into a valid Python identifier.

    Replaces operators and parentheses with mnemonic tokens, maps patterns such
    as ``\\dot{x}^i`` to ``dxi_dt``, collapses adjacent powers, and replaces any
    remaining non-alphanumeric characters with underscores. If the result would
    start with a digit, a leading underscore is prepended.

    Parameters
    ----------
    name : str
        Raw symbol or expression string to sanitize.

    Returns
    -------
    str
        A string suitable for use as a Python variable name.

    Examples
    --------
    >>> make_valid_varname("x**2")
    'x_pow2'
    >>> make_valid_varname("1foo")
    '_1foo'
    """
    replaces = {
        r"\*\*": "_pow",
        r"\(": "_I",
        r"\)": "I_",
        r"/": "_over_",
        r"\+": "_plus_",
    }

    name = re.sub(r"\\dot\{x\}\^(-?\d+)", r"dx\1_dt", name)

    for pattern, replacement in replaces.items():
        name = re.sub(pattern, replacement, name)

    name = re.sub(
        r"_pow(-?\d+)_pow(-?\d+)",
        lambda m: f"_pow{int(m.group(1)) * int(m.group(2))}",
        name,
    )
    name = re.sub(r"[^0-9a-zA-Z_]", "_", name)
    if name and name[0].isdigit():
        name = "_" + name
    return name


def generate_new_symbols(original_symbols):
    """
    Map each SymPy symbol to a new symbol with a simplified string name.

    Coordinate-like names containing ``x^`` are shortened by removing the caret
    (e.g. ``x^i`` → ``xi``). All other symbols are passed through
    :func:`make_valid_varname`.

    Parameters
    ----------
    original_symbols : iterable
        Symbols to rename (e.g. ``expr.free_symbols``).

    Returns
    -------
    dict
        Mapping ``old_symbol -> new_symbol`` with string names adjusted as
        described above.

    Examples
    --------
    >>> x, y = sp.symbols("x y")
    >>> out = generate_new_symbols([x, y])
    >>> sorted(str(v) for v in out.values())
    ['x', 'y']
    """
    new_symbols = {}
    for symbol in original_symbols:
        str_symbol = str(symbol)
        if "x^" in str_symbol:
            str_symbol = str_symbol.replace("x^", "x")
        else:
            str_symbol = make_valid_varname(str_symbol)

        new_symbols[symbol] = sp.symbols(str_symbol)

    return new_symbols


def rename_symbols_of_expr(expr):
    """
    Substitute all free symbols in an expression with renamed counterparts.

    Parameters
    ----------
    expr : sympy.Expr
        Expression whose free symbols should be renamed.

    Returns
    -------
    new_expr : sympy.Expr
        Expression after substitution.
    coordinate_substitutions : dict
        Mapping from original symbols to new symbols (same as returned by
        :func:`generate_new_symbols`).

    Examples
    --------
    >>> t = sp.symbols("t")
    >>> e = t**2
    >>> new_e, subs = rename_symbols_of_expr(e)
    >>> new_e == t**2  # symbol object may differ but name matches
    True
    """
    original_symbols = expr.free_symbols
    coordinate_substitutions = generate_new_symbols(original_symbols)
    new_expr = expr.subs(coordinate_substitutions)

    return new_expr, coordinate_substitutions


def get_optimized_christoffel_symbols(metric):
    """
    Build a reduced Christoffel tensor and a history of symbolic substitutions.

    The Christoffel tensor :math:`\\Gamma^\\mu_{\\nu\\sigma}` is computed from
    ``metric``, symbols are renamed for valid Python identifiers, components with
    :math:`\\nu > \\sigma` are zeroed (symmetry in lower indices), and repeated
    sub-expressions of selected types are replaced by new symbols to shorten the
    representation.

    Parameters
    ----------
    metric : object
        Symbolic metric providing ``christoffel_symbols().tensor()`` (e.g. a
        Relatipy symbolic metric instance).

    Returns
    -------
    optimized_symbols : sympy.Array
        Christoffel tensor as a SymPy array of shape ``(4, 4, 4)`` after the
        optimizations above.
    substitution_history : list of dict
        List of substitution dictionaries: first entry maps non-coordinate
        renames; later entries map factored repeated expressions to new symbols.

    Notes
    -----
    Requires ``metric`` to expose ``christoffel_symbols().tensor()`` returning a
    SymPy array-like object that supports ``free_symbols`` and ``subs`` (e.g.
    ``sympy.tensor.array.ImmutableDenseNDimArray``).

    Examples
    --------
    >>> from unittest.mock import MagicMock
    >>> import sympy as sp
    >>> from sympy.tensor.array import ImmutableDenseNDimArray
    >>> arr = ImmutableDenseNDimArray([sp.Integer(0)] * 64, (4, 4, 4))
    >>> inner = MagicMock()
    >>> inner.tensor.return_value = arr
    >>> metric = MagicMock()
    >>> metric.christoffel_symbols.return_value = inner
    >>> ch, hist = get_optimized_christoffel_symbols(metric)
    >>> ch.shape
    (4, 4, 4)
    >>> isinstance(hist, list)
    True
    """
    christoffel_symbols = metric.christoffel_symbols().tensor()

    substitution_history = []

    optimized_symbols, coordinate_substitutions = rename_symbols_of_expr(
        christoffel_symbols
    )

    non_coordinate_subs = {
        symbol: replacement
        for symbol, replacement in coordinate_substitutions.items()
        if "x^" not in str(symbol)
    }
    substitution_history.append(non_coordinate_subs)

    symbol_array = sp.MutableDenseNDimArray(optimized_symbols.tolist(), (4, 4, 4))

    for mu, nu, sigma in product(range(4), repeat=3):
        if nu > sigma:
            symbol_array[mu, nu, sigma] = 0

    optimized_symbols = optimized_symbols.__class__(symbol_array.tolist(), (4, 4, 4))

    expression_types = [sp.Function, sp.Pow, sp.Pow, sp.Mul]

    for expr_type in expression_types:
        repeated_expressions = _find_repeated_expressions(optimized_symbols, expr_type)

        if repeated_expressions:
            expression_substitutions = generate_new_symbols(repeated_expressions)
            optimized_symbols = optimized_symbols.subs(expression_substitutions)
            substitution_history.append(expression_substitutions)

    return optimized_symbols, substitution_history


def _find_repeated_expressions(tensor, expression_type):
    """
    Collect sub-expressions of a given type that occur more than once in a tensor.

    Parameters
    ----------
    tensor : sympy.Array
        Array whose elements are searched for repeated atomic sub-expressions.
    expression_type : type
        SymPy class to match (e.g. ``sympy.Function``, ``sympy.Pow``, ``sympy.Mul``).

    Returns
    -------
    set
        Sub-expressions of ``expression_type`` that appear at least twice in
        ``tensor``.

    Examples
    --------
    >>> x = sp.symbols("x")
    >>> arr = sp.Array([x**2, x**2, 1])
    >>> sorted(_find_repeated_expressions(arr, sp.Pow), key=str)
    [x**2]
    """
    all_expressions = {expr for expr in tensor.atoms(expression_type)}

    repeated_expressions = {
        expr for expr in all_expressions if tensor.count(expr) > 1
    }

    return repeated_expressions


def fill_template(template: str, replacements: dict) -> str:
    """
    Replace ``$key$`` placeholders in a template string with provided values.

    Placeholders must appear alone on a line as ``$name$`` (after optional
    leading whitespace). For multi-line replacements, each non-empty line is
    prefixed with the leading whitespace of the placeholder line, with the
    first character of that whitespace run removed. Missing keys leave the
    placeholder unchanged.

    Parameters
    ----------
    template : str
        Template text containing ``$word$`` placeholders.
    replacements : dict
        Mapping from placeholder names (without ``$``) to replacement strings.

    Returns
    -------
    str
        Template with substitutions applied.

    Examples
    --------
    >>> fill_template("  $a$\\n", {"a": "line1\\nline2"})
    ' line1\\n line2\\n'
    """
    def replace(match):
        full_match = match.group(0)
        key = match.group(2)
        leading_spaces = match.group(1)[1:]

        if key not in replacements:
            return full_match

        content = replacements[key]

        if not content:
            return ""

        lines = content.split("\n")

        if len(lines) == 1:
            return content

        indented_lines = []

        for line in lines:
            if line.strip():
                indented_lines.append(leading_spaces + line)
            else:
                indented_lines.append("")

        return "\n".join(indented_lines)

    return re.sub(r"^(\s*)\$(\w+)\$$", replace, template, flags=re.MULTILINE)


def export_christoffel_symbols_to_code(metric):
    """
    Emit Python source text defining auxiliary variables and Christoffel entries.

    Uses :func:`get_optimized_christoffel_symbols` and :data:`template_metric`
    to build a string suitable for pasting into a numeric metric class: variable
    assignments for substituted symbols, then ``Gamma[mu, nu, rho] = ...`` lines
    with symmetry relations for :math:`\\nu > \\rho`.

    Parameters
    ----------
    metric : object
        Symbolic metric passed to :func:`get_optimized_christoffel_symbols`.

    Returns
    -------
    str
        Generated Python code as a single string.

    Examples
    --------
    >>> isinstance(fill_template("", {}), str)
    True
    """
    ch_subs, all_subs = get_optimized_christoffel_symbols(metric)

    auxiliary_variables = ""

    mega_dict = all_subs[0].copy()

    for key, value in all_subs[0].items():
        auxiliary_variables += f"{value} = self.{key}\n"
    auxiliary_variables += "\n"

    for i, sub_ in enumerate(all_subs[1:]):
        sub_ = dict(sorted(sub_.items(), key=lambda item: str(item[1])))
        for key, value in sub_.items():
            if value not in mega_dict.values():
                auxiliary_variables += f"{value} = {key}\n"
        i += 1
        mega_dict.update(sub_)
        auxiliary_variables += "\n"

    auxiliary_variables += "\n"

    gamma_txt = ""

    for mu in range(0, 4):
        gamma_txt += "\n"
        gamma_txt += f"# mu = {mu}\n"
        for nu, rho in product(range(4), repeat=2):
            if nu <= rho:
                if ch_subs[mu, nu, rho] != 0:
                    gamma_txt += f"Gamma[{mu}, {nu}, {rho}] = {ch_subs[mu, nu, rho]}\n"
        gamma_txt += "\n"
        for nu, rho in product(range(4), repeat=2):
            if nu > rho:
                if ch_subs[mu, rho, nu] != 0:
                    gamma_txt += (
                        f"Gamma[{mu}, {nu}, {rho}] = Gamma[{mu}, {rho}, {nu}]\n"
                    )

    txt = fill_template(
        template_metric,
        {
            "auxiliary_variables": auxiliary_variables,
            "gamma_txt": gamma_txt,
        },
    )

    return txt
