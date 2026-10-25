"""The private symbolic optimizer emits correct, semantic and compact C99."""

from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

sp = pytest.importorskip("sympy")

from relatipy import _jit  # noqa: E402
from relatipy._codegen import (  # noqa: E402
    compile_and_compare,
    count_c_operations,
    numeric_max_relative_error,
    optimize,
    symbolic_residuals,
)

ROOT = Path(__file__).resolve().parents[2]
r, th, spin, x, y = sp.symbols("r th spin x y", real=True)


def _has_compiler() -> bool:
    try:
        _jit.compiler_command()
    except _jit.CompilerNotFoundError:
        return False
    return True


needs_cc = pytest.mark.skipif(not _has_compiler(), reason="no C compiler")


def test_import_without_sympy_names_the_extra() -> None:
    script = (
        "import sys; sys.modules['sympy'] = None\n"
        "import relatipy._codegen as codegen\n"
        "try:\n"
        "    codegen.optimize\n"
        "except ImportError as error:\n"
        "    print(error)\n"
    )
    environment = dict(os.environ, PYTHONPATH=str(ROOT / "src"))
    completed = subprocess.run(
        [sys.executable, "-c", script], capture_output=True, text=True, env=environment, check=True,
    )
    assert "relatipy[symbolic]" in completed.stdout


def test_multiple_angles_reduce_to_one_sine_and_cosine() -> None:
    expressions = [
        sp.sin(2 * th) * r + sp.cos(2 * th),
        1 - sp.cos(4 * th),
        sp.cot(th) * r + sp.tan(th),
    ]
    code = optimize(expressions, inputs=(r,), angles={th: ("sin_value", "cos_value")}, outputs="out")
    body = code.body()
    assert "sin(" not in body and "cos(" not in body and "tan(" not in body
    assert code.report.after.call == 0
    assert "const double sin_th = sin_value;" in body
    assert symbolic_residuals(code) == [0, 0, 0]
    # 1 - cos(4 th) = 8 sin^2 cos^2 needs no addition.
    assert code.outputs[1][2] == "8.0 * cos2_th * sin2_th"


def test_computed_angle_calls_sine_and_cosine_once() -> None:
    code = optimize(
        [sp.sin(th) ** 3 + sp.sin(3 * th), sp.cos(th) * sp.sin(2 * th)],
        inputs=(r, th), angles=(th,), outputs="out",
    )
    body = code.body()
    assert body.count("sin(th)") == 1 and body.count("cos(th)") == 1
    assert code.report.after.call == 2
    assert numeric_max_relative_error(code, domain={th: (0.1, 3.0)}) < 1e-12


def test_pythagorean_pairs_contract_without_cancellation() -> None:
    expression = r * spin**2 - r * spin**2 * sp.cos(th) ** 2
    code = optimize(expression, inputs=(r, spin), angles={th: ("s", "c")}, outputs="out")
    assert "cos2_th" not in code.body()
    assert symbolic_residuals(code) == [0]


def test_integer_powers_are_multiplication_chains() -> None:
    code = optimize([r**k for k in range(2, 8)], inputs=(r,), outputs="out")
    powers = {d.name: d.c_text for d in code.definitions if d.kind == "power"}
    assert powers == {
        "r2": "r * r", "r3": "r2 * r", "r4": "r2 * r2",
        "r5": "r4 * r", "r6": "r3 * r3", "r7": "r6 * r",
    }
    assert "pow(" not in code.body()
    assert code.report.after.mul == 6


def test_shared_powers_reuse_named_products() -> None:
    expressions = [r**2 * spin**2, r**4 + spin**4, r**6]
    code = optimize(expressions, inputs=(r, spin), outputs="out")
    powers = {d.name: d.c_text for d in code.definitions if d.kind == "power"}
    assert powers["r6"] == "r4 * r2"
    assert powers["spin4"] == "spin2 * spin2"


def test_reciprocal_powers_and_square_roots_are_shared() -> None:
    sigma = r**2 + spin**2 * sp.cos(th) ** 2
    expressions = [r / sigma**2, spin / sigma**2, (x**2 + y**2) ** sp.Rational(3, 2),
                   1 / sp.sqrt(x**2 + y**2)]
    code = optimize(expressions, inputs=(r, spin, x, y), angles={th: ("s", "c")}, outputs="out")
    assert code.report.after.call == 1
    assert code.body().count("sqrt(") == 1
    assert code.report.after.div <= 2
    assert any(d.name.startswith("inv_") for d in code.definitions)
    assert numeric_max_relative_error(code, domain={th: (0.2, 2.9)}) < 1e-14


def test_common_subexpressions_are_shared_across_outputs() -> None:
    shared = sp.exp(r * spin) / (1 + r**2)
    code = optimize([shared * x, shared * y], inputs=(r, spin, x, y), outputs="out")
    assert code.body().count("exp(") == 1
    assert code.report.after.total < code.report.before.total


def test_names_are_semantic_bounded_and_deterministic() -> None:
    sigma = r**2 + spin**2 * sp.cos(th) ** 2
    expressions = [
        r * spin / sigma + sp.sin(2 * th) / sigma**2,
        (r**2 + spin**2) / sigma + r / sigma,
        (r**3 + 3 * spin * x + y * x**2 + r * spin * x * y + 7) ** 3 / sigma,
    ]
    arguments = dict(inputs=(r, spin, x, y), angles={th: ("s", "c")}, outputs="out",
                     max_name_length=24)
    first = optimize(expressions, **arguments)
    second = optimize(expressions, **arguments)
    assert first.body() == second.body()
    names = [d.name for d in first.definitions]
    assert len(names) == len(set(names))
    for name in names:
        assert re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", name)
        assert "tmp" not in name and "cse" not in name and "__" not in name
        assert len(name) <= 24 + 3
    assert "r2_plus_spin2" in names or "inv_sigma" not in names


def test_strategy_search_selects_the_cheapest_candidate() -> None:
    expressions = [r**2 / (r + spin) + r * spin / (r + spin), (r + spin) ** 2]
    code = optimize(expressions, inputs=(r, spin), outputs="out",
                    transforms=("none", "factor_terms", "cancel", "factor"))
    report = code.report
    assert report.candidates
    best = min(report.candidates, key=lambda c: (c.operations, c.temporaries, c.name))
    assert report.strategy == best.name
    assert report.after.total == best.operations
    assert report.temporaries == best.temporaries
    assert "operations before" in str(report)


def test_matrix_tensor_mapping_and_sequence_outputs() -> None:
    matrix = sp.Matrix([[r, spin], [r * spin, 0]])
    code = optimize(matrix, inputs=(r, spin), outputs="g")
    assert code.lvalues == ("g[0][0]", "g[0][1]", "g[1][0]", "g[1][1]")
    assert code.outputs[3][2] == "0.0"

    tensor = sp.MutableDenseNDimArray.zeros(2, 2, 2)
    tensor[1, 0, 1] = r**2
    code = optimize(tensor, inputs=(r,), outputs="gamma")
    assert code.lvalues[5] == "gamma[1][0][1]"

    code = optimize({"dx": sp.Matrix([r, spin]), "value": r * spin}, inputs=(r, spin))
    assert code.lvalues == ("dx[0]", "dx[1]", "value")

    code = optimize([r, sp.Matrix([spin, r])], inputs=(r, spin), outputs=["*first", "rest"])
    assert code.lvalues == ("*first", "rest[0]", "rest[1]")

    code = optimize(sp.Matrix([r, spin]), inputs=(r, spin), outputs=["a", None])
    assert code.lvalues == ("a",)


def test_input_bindings_and_reserved_names() -> None:
    code = optimize([r**2, spin * r**2], inputs={r: "state->x[1]", spin: None}, outputs="out",
                    reserved=("r2",))
    body = code.body()
    assert "const double r = state->x[1];" in body
    assert "const double r2 " not in body and "const double r2_2 = r * r;" in body


def test_invalid_requests_raise() -> None:
    with pytest.raises(ValueError, match="undeclared symbols: spin"):
        optimize(r * spin, inputs=(r,), outputs="out")
    with pytest.raises(ValueError, match="C identifier"):
        optimize(r, inputs=(sp.Symbol("double"),), outputs="out")
    with pytest.raises(ValueError, match="unknown transforms"):
        optimize(r, inputs=(r,), outputs="out", transforms=("simplify",))
    with pytest.raises(ValueError, match="must be an input"):
        optimize(sp.sin(th), inputs=(r,), angles=(th,), outputs="out")
    with pytest.raises(TypeError, match="outputs is required"):
        optimize(r, inputs=(r,))


def test_symbolic_residual_reports_intractable_differences() -> None:
    code = optimize(sp.sin(2 * th) * r**3 + sp.cos(3 * th), inputs=(r, th), angles=(th,),
                    outputs="out")
    assert symbolic_residuals(code) == [0]
    assert symbolic_residuals(code, max_operations=1) is None


@needs_cc
def test_emitted_function_compiles_and_matches(tmp_path) -> None:
    sigma = r**2 + spin**2 * sp.cos(th) ** 2
    expressions = [r * sp.sin(2 * th) / sigma**2, sp.sqrt(r**2 + spin**2) * sp.cot(th),
                   (r**2 + spin**2) ** 2 / sigma]
    code = optimize(expressions, inputs=(r, spin, th), angles=(th,), outputs="out")
    domain = {r: (2.0, 30.0), spin: (-1.0, 1.0), th: (0.1, 3.0)}
    assert compile_and_compare(code, domain=domain, cache_dir=tmp_path) < 1e-14

    source = "#include <math.h>\n" + code.c_function(
        "evaluate", ["double r", "double spin", "double th", "double out[3]"], doc="Test kernel."
    )
    path = tmp_path / "kernel.c"
    path.write_text(source)
    compiler = _jit.compiler_command()
    for standard in ("-std=c99", "-std=c11"):
        subprocess.run(
            [*compiler, standard, "-Wall", "-Wextra", "-Wpedantic", "-Werror", "-c", str(path),
             "-o", str(tmp_path / "kernel.o")],
            check=True, capture_output=True,
        )


def test_count_c_operations_matches_the_internal_count() -> None:
    code = optimize([r**3 / (1 + spin**2) - r, sp.sqrt(r) * spin], inputs=(r, spin), outputs="o")
    assert count_c_operations(code.body()) == code.report.after
    counted = count_c_operations("const double a = -b * c[1] + f(d) / e - state->x[2];")
    assert (counted.add, counted.mul, counted.div, counted.call) == (2, 1, 1, 1)


def test_evaluate_matches_the_original_expressions() -> None:
    code = optimize([r**3 * sp.cos(2 * th), spin / (r + 1)], inputs=(r, spin, th), angles=(th,),
                    outputs="out")
    values = {r: np.array([1.5, 2.0]), spin: np.array([0.3, -0.4]), th: np.array([0.4, 2.0])}
    expected = [1.5**3 * np.cos(0.8), 0.3 / 2.5]
    actual = [row[0] for row in code.evaluate(values)]
    np.testing.assert_allclose(actual, expected, rtol=1e-14)
