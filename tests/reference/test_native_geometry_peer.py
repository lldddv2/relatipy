"""Compare the experimental C Kerr geometry with independent references.

The test compiles the standalone C reference into a temporary shared library;
it does not select a production build system or expose a Python binding.
PyGRO is an optional, independent oracle for the metric. The connection
oracle is the Levi-Civita connection of the Boyer-Lindquist Kerr metric,
written out and differentiated here with SymPy, independently of
``relatipy._codegen``, and evaluated with 40-digit mpmath.
"""

from __future__ import annotations

import ctypes
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from .peer_oracles import make_pygro_engine


DIMENSION = 4
ROOT = Path(__file__).resolve().parents[2]

CASES = (
    (1.0, 0.5, (0.0, 8.0, 1.1, 0.3)),
    (2.0, -1.4, (0.2, 15.0, 0.8, -0.4)),
    (0.75, 0.675, (-1.0, 6.0, 1.35, 2.0)),
)
TANGENTS = (
    (1.1, -0.01, 0.002, 0.02),
    (1.0, 0.3, -0.01, 0.015),
    (0.9, 0.0, 0.0, 0.05),
)


@pytest.fixture(scope="module")
def native_geometry(tmp_path_factory: pytest.TempPathFactory):
    """Compile and load the experimental C geometry from a temporary path."""

    compiler = shutil.which("cc")
    if compiler is None:
        pytest.skip("a C compiler is required for the optional native peer check")

    output = tmp_path_factory.mktemp("native-kerr") / "librelatipy-kerr-reference.so"
    subprocess.run(
        [
            compiler,
            "-std=c11",
            "-shared",
            "-fPIC",
            "-Wall",
            "-Wextra",
            "-Wpedantic",
            "-Werror",
            f"-I{ROOT / 'native' / 'include'}",
            str(ROOT / "native" / "src" / "utils" / "numeric.c"),
            str(ROOT / "native" / "src" / "utils" / "tensor.c"),
            str(ROOT / "native" / "src" / "metric" / "physic" / "kerr.c"),
            str(ROOT / "native" / "src" / "metric" / "kerr.c"),
            str(ROOT / "native" / "src" / "geodesic" / "physic" / "kerr.c"),
            str(ROOT / "native" / "src" / "geodesic" / "physic" / "kerr_rhs_generated.c"),
            str(ROOT / "native" / "src" / "geodesic" / "kerr.c"),
            "-lm",
            "-o",
            str(output),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    library = ctypes.CDLL(str(output))
    pointer = ctypes.POINTER(ctypes.c_double)
    library.rp_kerr_metric.argtypes = [ctypes.c_double, ctypes.c_double, pointer, pointer]
    library.rp_kerr_metric.restype = ctypes.c_int
    library.rp_kerr_null_geodesic_rhs.argtypes = [
        ctypes.c_double, ctypes.c_double, pointer, pointer, pointer, pointer
    ]
    library.rp_kerr_null_geodesic_rhs.restype = ctypes.c_int
    return library


def _call_native(library, name: str, mass: float, spin: float, point, shape):
    """Evaluate one flat C output buffer and reshape it as a NumPy array."""

    coordinates = (ctypes.c_double * DIMENSION)(*point)
    size = int(np.prod(shape))
    output = (ctypes.c_double * size)()
    status = getattr(library, name)(mass, spin, coordinates, output)
    assert status == 0
    return np.ctypeslib.as_array(output).copy().reshape(shape)


@pytest.mark.parametrize(("mass", "spin", "point"), CASES)
def test_native_metric_matches_pygro(
    native_geometry, mass: float, spin: float, point: tuple[float, ...]
):
    """Match the complete C metric at regular Kerr points."""

    pytest.importorskip("pygro", reason="optional peer library PyGRO is not installed")
    pygro_metric, _engine = make_pygro_engine()
    pygro_metric.set_constant(m=mass, a=spin)

    c_metric = _call_native(
        native_geometry, "rp_kerr_metric", mass, spin, point, (4, 4)
    )
    peer_metric = np.asarray(pygro_metric._g_f(point), dtype=float)

    np.testing.assert_allclose(c_metric, peer_metric, rtol=2e-14, atol=2e-14)


@pytest.fixture(scope="module")
def levi_civita_acceleration():
    """``-Gamma^mu_(alpha beta) u^alpha u^beta`` of the Kerr metric at 40 digits."""

    sp = pytest.importorskip("sympy")
    mpmath = pytest.importorskip("mpmath")
    m, a, r, theta = sp.symbols("m a r theta", real=True)
    u = sp.symbols("u0:4", real=True)
    coordinates = (sp.Symbol("t"), r, theta, sp.Symbol("phi"))
    sigma = r**2 + a**2 * sp.cos(theta) ** 2
    delta = r**2 - 2 * m * r + a**2
    sin2 = sp.sin(theta) ** 2
    g = sp.zeros(4, 4)
    g[0, 0] = -(1 - 2 * m * r / sigma)
    g[0, 3] = g[3, 0] = -2 * m * a * r * sin2 / sigma
    g[1, 1] = sigma / delta
    g[2, 2] = sigma
    g[3, 3] = (r**2 + a**2 + 2 * m * a**2 * r * sin2 / sigma) * sin2
    inverse = g.inv()
    derivative = [[[sp.diff(g[i, j], coordinates[k]) for k in range(4)]
                   for j in range(4)] for i in range(4)]
    acceleration = []
    for mu in range(4):
        total = 0
        for alpha in range(4):
            for beta in range(4):
                gamma = sum(
                    inverse[mu, nu] * (derivative[nu][alpha][beta] + derivative[nu][beta][alpha]
                                       - derivative[alpha][beta][nu])
                    for nu in range(4)
                ) / 2
                total -= gamma * u[alpha] * u[beta]
        acceleration.append(total)
    function = sp.lambdify((m, a, r, theta, *u), acceleration, modules="mpmath")

    def evaluate(mass, spin, point, tangent):
        with mpmath.workdps(40):
            values = function(*(mpmath.mpf(value) for value in (mass, spin, point[1], point[2])),
                              *(mpmath.mpf(value) for value in tangent))
            return np.array([float(value) for value in values])

    return evaluate


@pytest.mark.parametrize(("mass", "spin", "point"), CASES)
@pytest.mark.parametrize("tangent", TANGENTS)
def test_native_rhs_matches_the_levi_civita_connection(
    native_geometry, levi_civita_acceleration, mass, spin, point, tangent
):
    """The C geodesic RHS contracts the Levi-Civita connection of the metric."""

    pointer = ctypes.POINTER(ctypes.c_double)
    coordinates = np.array(point, dtype=float)
    velocity = np.array(tangent, dtype=float)
    dx, du = np.zeros(4), np.zeros(4)
    status = native_geometry.rp_kerr_null_geodesic_rhs(
        mass, spin, *(array.ctypes.data_as(pointer) for array in (coordinates, velocity, dx, du))
    )
    assert status == 0
    np.testing.assert_array_equal(dx, velocity)
    expected = levi_civita_acceleration(mass, spin, point, tangent)
    np.testing.assert_allclose(du, expected, rtol=1e-13, atol=1e-14 * np.max(np.abs(expected)))
