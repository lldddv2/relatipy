"""Acceptance: the optimizer reproduces the native Kerr geodesic kernels.

The Kerr timelike right-hand side ``-Gamma^mu_(alpha beta) u^alpha u^beta``
is derived symbolically from the metric of ``native/README.md`` (signature
``(-,+,+,+)``, ``Delta = r^2 - 2 M r + a^2``), optimized, compiled and
compared with the library's C kernels:

- the Levi--Civita contraction ``rp_kerr_geodesic_physic_evaluate_corrected_rhs``
  at every spin;
- the hand-tuned legacy-compatible ``rp_kerr_geodesic_physic_evaluate_rhs``
  at ``a = 0`` only, because for ``a != 0`` that kernel contracts the frozen
  legacy connection, which is not the connection of this metric;
- the same hand-tuned kernel at every spin when the optimizer is given that
  legacy connection, transcribed from ``rp_kerr_physic_evaluate_christoffel``.
"""

from __future__ import annotations

import ctypes
import re
from pathlib import Path

import numpy as np
import pytest

sp = pytest.importorskip("sympy")

from relatipy import _jit  # noqa: E402
from relatipy._codegen import count_c_operations, optimize, symbolic_residuals  # noqa: E402
from relatipy._codegen.kerr import (  # noqa: E402
    block_determinant_residual,
    geodesic_acceleration,
    kerr_symbols,
    simplified_geodesic_acceleration,
)

ROOT = Path(__file__).resolve().parents[2]
NATIVE = ROOT / "native"
SPINS = (-0.9, -0.5, 0.0, 0.3, 0.7, 0.99)
SIGNATURE = [
    "double mass", "double spin", "const struct rp_kerr_geodesic_state *state",
    "double dx[RP_KERR_DIM]", "double du[RP_KERR_DIM]",
]
WRAPPER = """
int rp_test_rhs(int which, double mass, double spin, const double x[4],
                const double u[4], double du[4])
{
    struct rp_kerr_geodesic_state state;
    double dx[4];
    int i;
    for (i = 0; i < 4; ++i) {
        state.x[i] = x[i];
        state.u[i] = u[i];
    }
    if (rp_kerr_geodesic_physic_prepare_state(mass, spin, &state) != RP_KERR_STATUS_OK) {
        return 1;
    }
    switch (which) {
    case 0: rp_kerr_geodesic_physic_evaluate_rhs(mass, spin, &state, dx, du); break;
    case 1: rp_kerr_geodesic_physic_evaluate_corrected_rhs(mass, spin, &state, dx, du); break;
    default: candidate_rhs(mass, spin, &state, dx, du); break;
    }
    for (i = 0; i < 4; ++i) {
        if (dx[i] != u[i]) {
            return 2;
        }
    }
    return 0;
}
"""


def _has_compiler() -> bool:
    try:
        _jit.compiler_command()
    except _jit.CompilerNotFoundError:
        return False
    return True


pytestmark = pytest.mark.skipif(not _has_compiler(), reason="no C compiler")


def _prepared_inputs(s):
    return {
        s.mass: None, s.spin: None, s.r: "state->x[1]",
        s.sigma: "state->sigma", s.delta: "state->delta",
        **{u: f"state->u[{i}]" for i, u in enumerate(s.u)},
    }


def _optimize_rhs(s, acceleration, **options):
    return optimize(
        {"dx": sp.Matrix(s.u), "du": sp.Matrix(acceleration)},
        inputs=_prepared_inputs(s),
        angles={s.theta: ("state->sin_theta", "state->cos_theta")},
        polynomial_in=s.u,
        **options,
    )


class NativeRhs:
    """The library kernels and a candidate compiled into one library."""

    def __init__(self, candidate_code, cache_dir) -> None:
        sources = ("geodesic/physic/kerr.c", "utils/numeric.c", "utils/tensor.c")
        text = "".join(f'#include "{(NATIVE / "src" / item).as_posix()}"\n' for item in sources)
        text += candidate_code.c_function("candidate_rhs", SIGNATURE, static=True) + WRAPPER
        self.function = _jit.compile_function(
            text, "rp_test_rhs", cache_dir=cache_dir,
            include_dirs=(NATIVE / "include", NATIVE / "src"),
        )
        double = ctypes.c_double
        pointer = ctypes.POINTER(double)
        self.call = self.function.as_ctypes(
            ctypes.c_int, ctypes.c_int, double, double, pointer, pointer, pointer,
        )

    def __call__(self, which: int, spin: float, x, u):
        x = np.ascontiguousarray(x, dtype=float)
        u = np.ascontiguousarray(u, dtype=float)
        out = np.zeros(4)
        pointer = ctypes.POINTER(ctypes.c_double)
        status = self.call(which, 1.0, spin, x.ctypes.data_as(pointer), u.ctypes.data_as(pointer),
                           out.ctypes.data_as(pointer))
        assert status == 0
        return out


def _regular_points(spin: float, count: int, seed: int):
    generator = np.random.default_rng(seed)
    for _ in range(count):
        r = generator.uniform(3.0, 60.0)
        theta = generator.uniform(0.15, np.pi - 0.15)
        x = np.array([0.0, r, theta, generator.uniform(0.0, 2.0 * np.pi)])
        u = np.array([1.0 + generator.random(), generator.normal(scale=0.3),
                      generator.normal(scale=0.3) / r, generator.normal(scale=0.3) / r])
        yield x, u


def _max_relative(actual, expected) -> float:
    return float(np.max(np.abs(actual - expected) / np.abs(expected)))


def test_block_determinant_identity() -> None:
    assert block_determinant_residual(kerr_symbols()) == 0


@pytest.mark.parametrize("form", ["prepared", "simplified"])
def test_metric_derived_rhs_matches_the_corrected_kernel(tmp_path, form) -> None:
    s = kerr_symbols()
    if form == "prepared":
        acceleration = geodesic_acceleration(s, prepared=True)
    else:
        acceleration = simplified_geodesic_acceleration(s)
    code = _optimize_rhs(s, acceleration)
    assert code.report.after.total < code.report.before.total
    native = NativeRhs(code, tmp_path)
    for index, spin in enumerate(SPINS):
        for x, u in _regular_points(spin, 40, seed=index):
            expected = native(1, spin, x, u)
            actual = native(2, spin, x, u)
            assert _max_relative(actual, expected) < 1e-12
            if spin == 0.0:
                legacy = native(0, spin, x, u)
                assert _max_relative(actual, legacy) < 1e-12


def test_explicit_metric_rhs_matches_symbolically_and_numerically(tmp_path) -> None:
    s = kerr_symbols()
    acceleration = geodesic_acceleration(s, prepared=False)
    code = optimize(
        sp.Matrix(acceleration), inputs=(s.mass, s.spin, s.r, *s.u, s.theta), angles=(s.theta,),
        outputs="du", polynomial_in=s.u, transforms=("none", "factor_terms"),
    )
    from relatipy._codegen import compile_and_compare

    domain = {s.mass: (1.0, 1.0), s.spin: (-0.95, 0.95), s.r: (3.0, 60.0),
              s.theta: (0.15, 3.0), s.u[0]: (1.0, 2.0), s.u[1]: (-0.3, 0.3),
              s.u[2]: (-0.01, 0.01), s.u[3]: (-0.01, 0.01)}
    assert compile_and_compare(code, domain=domain, samples=200, cache_dir=tmp_path) < 1e-12


def _legacy_connection(s):
    """Symbolic transcription of ``rp_kerr_physic_evaluate_christoffel``."""
    text = (NATIVE / "src" / "metric" / "physic" / "kerr.c").read_text(encoding="utf-8")
    start = text.index("void rp_kerr_physic_evaluate_christoffel(")
    body = text[text.index("{", start) + 1: text.index("\n}\n", start)]
    body = re.sub(r"/\*.*?\*/", " ", body, flags=re.S)
    replacements = {
        "point->x[1]": "r", "point->cos_theta": "cos(th)", "point->sin_theta": "sin(th)",
    }
    environment = {
        "r": s.r, "th": s.theta, "mass": s.mass, "spin": s.spin, "sin": sp.sin, "cos": sp.cos,
    }
    gamma = sp.MutableDenseNDimArray.zeros(4, 4, 4)
    for statement in body.split(";"):
        statement = " ".join(statement.split())
        if "=" not in statement:
            continue
        left, right = (part.strip() for part in statement.split("=", 1))
        for old, new in replacements.items():
            right = right.replace(old, new)
        right = re.sub(r"\b(\d+)\.0\b", r"\1", right)
        right = re.sub(r"christoffel\[(\d)\]\[(\d)\]\[(\d)\]", r"gamma_\1\2\3", right)
        value = sp.parse_expr(right, local_dict=environment)
        if left.startswith("const double "):
            environment[left.removeprefix("const double ").strip()] = value
            continue
        indices = tuple(int(item) for item in re.findall(r"\[(\d)\]", left))
        assert left.startswith("christoffel") and len(indices) == 3
        gamma[indices] = value
        environment["gamma_" + "".join(map(str, indices))] = value
    return gamma


def test_legacy_connection_reproduces_the_hand_tuned_kernel(tmp_path) -> None:
    s = kerr_symbols()
    gamma = _legacy_connection(s)
    nonzero = [index for index in np.ndindex(4, 4, 4) if gamma[index] != 0]
    assert len(nonzero) == 32
    acceleration = [
        -sum(gamma[mu, a, b] * s.u[a] * s.u[b] for a in range(4) for b in range(4))
        for mu in range(4)
    ]
    code = _optimize_rhs(s, acceleration, transforms=("none", "factor_terms"))
    assert symbolic_residuals(code) == [0] * 8
    native = NativeRhs(code, tmp_path)
    for index, spin in enumerate(SPINS):
        for x, u in _regular_points(spin, 40, seed=100 + index):
            assert _max_relative(native(2, spin, x, u), native(0, spin, x, u)) < 1e-12

    hand = (NATIVE / "src" / "geodesic" / "physic" / "kerr.c").read_text(encoding="utf-8")
    start = hand.index("void rp_kerr_geodesic_physic_evaluate_rhs(")
    hand_body = hand[hand.index("{", start): hand.index("\n}\n", start)]
    assert code.report.after.total <= count_c_operations(hand_body).total
