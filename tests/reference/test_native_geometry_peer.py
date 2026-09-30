"""Compare the experimental C Kerr geometry with frozen and peer references.

The test compiles the standalone C reference into a temporary shared library;
it does not select a production build system or expose a Python binding.  The
Christoffel oracle is the literal output of the frozen RelatiPy legacy
``Kerr._get_christoffel_symbols`` implementation requested for compatibility.
PyGRO remains an optional, independent oracle for the metric only.
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

# Frozen with the 2026-09-01 legacy snapshot.  For each case the legacy
# constructor received dimensionless spin ``spin / mass``, which stores the C
# API's geometric spin length as ``self.a``.  This fixture intentionally tests
# literal compatibility and does not assert that these values are the
# Levi-Civita connection of the separately tested corrected metric.
LEGACY_CHRISTOFFEL_CASES = (
    (
        1.0,
        0.5,
        (0.0, 8.0, 1.1, 0.3),
        [
            [[0.0, 0.020756247543677078, -0.00039413983281485126, 0.0], [0.020756247543677078, 0.0, 0.0, -0.02469076273038043], [-0.00039413983281485126, 0.0, 0.0, 0.00015652289119530666], [0.0, -0.02469076273038043, 0.00015652289119530666, 0.0]],
            [[0.011620304889560888, 0.0, 0.0, -0.004614716824978808], [0.0, -0.02169724155453523, -0.001577826425797284, 0.0], [0.0, -0.001577826425797284, -5.9639567157709426, 0.0], [-0.004614716824978808, 0.0, 0.0, -4.735043332424259]],
            [[-6.153489274526096e-06, 0.0, 0.0, 0.0007907233717766035], [0.0, 3.304348535701119e-05, 0.12489961708420821, 0.0], [0.0, 0.12489961708420821, -0.001577826425797284, 0.0], [0.0007907233717766035, 0.0, 0.0, -0.40612845345216936]],
            [[0.0, 0.0001615272182387321, -0.000992482355935277, 0.0], [0.0001615272182387321, 0.0, 0.0, 0.12432147266357679], [-0.000992482355935277, 0.0, 0.0, 0.5093622450718792], [0.0, 0.12432147266357679, 0.5093622450718792, 0.0]],
        ],
    ),
    (
        2.0,
        -1.4,
        (0.2, 15.0, 0.8, -0.4),
        [
            [[0.0, 0.0119310628525823, -0.0011512299811903667, 0.0], [0.0119310628525823, 0.0, 0.0, 0.025783109802546084], [-0.0011512299811903667, 0.0, 0.0, -0.000829391742690033], [0.0, 0.025783109802546084, -0.000829391742690033, 0.0]],
            [[0.00633317380482608, 0.0, 0.0, 0.0045626696182046265], [0.0, -0.013349072449853517, -0.004335366801519993, 0.0], [0.0, -0.004335366801519993, -10.823567227775882, 0.0], [0.0045626696182046265, 0.0, 0.0, -5.5665179820373805]],
            [[-5.09503397777441e-06, 0.0, 0.0, -0.0008259777939969145], [0.0, 2.6590816986751676e-05, 0.06638596189754589, 0.0], [0.0, 0.06638596189754589, -0.004335366801519993, 0.0], [-0.0008259777939969145, 0.0, 0.0, -0.5032052700803743]],
            [[0.0, -7.359661611568215e-05, 0.001597954743669163, 0.0], [-7.359661611568215e-05, 0.0, 0.0, 0.06593189833572628], [0.001597954743669163, 0.0, 0.0, 0.9723658306316648], [0.0, 0.06593189833572628, 0.9723658306316648, 0.0]],
        ],
    ),
    (
        0.75,
        0.675,
        (-1.0, 6.0, 1.35, 2.0),
        [
            [[0.0, 0.027612208840615205, -0.0006753081944310627, 0.0], [0.027612208840615205, 0.0, 0.0, -0.052831868149926374], [-0.0006753081944310627, 0.0, 0.0, 0.0004339694880985594], [0.0, -0.052831868149926374, 0.0004339694880985594, 0.0]],
            [[0.01532407904006055, 0.0, 0.0, -0.009847626300758312], [0.0, -0.031216461830180325, -0.0027028725434599314, 0.0], [0.0, -0.0027028725434599314, -4.421378531006508, 0.0], [-0.009847626300758312, 0.0, 0.0, -4.202983520670714]],
            [[-1.874718060273453e-05, 0.0, 0.0, 0.0010125039790526878], [0.0, 0.00010182468200739069, 0.16656555413365384, 0.0], [0.0, 0.16656555413365384, -0.0027028725434599314, 0.0], [0.0010125039790526878, 0.0, 0.0, -0.21755674974308947]],
            [[0.0, 0.0005112583028658887, -0.0010508599566847194, 0.0], [0.0005112583028658887, 0.0, 0.0, 0.16360543781649062], [-0.0010508599566847194, 0.0, 0.0, 0.22513102644315708], [0.0, 0.16360543781649062, 0.22513102644315708, 0.0]],
        ],
    ),
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
    for name in ("rp_kerr_metric", "rp_kerr_christoffel"):
        function = getattr(library, name)
        function.argtypes = [ctypes.c_double, ctypes.c_double, pointer, pointer]
        function.restype = ctypes.c_int
    return library


def _call_native(library, name: str, mass: float, spin: float, point, shape):
    """Evaluate one flat C output buffer and reshape it as a NumPy array."""

    coordinates = (ctypes.c_double * DIMENSION)(*point)
    size = int(np.prod(shape))
    output = (ctypes.c_double * size)()
    status = getattr(library, name)(mass, spin, coordinates, output)
    assert status == 0
    return np.ctypeslib.as_array(output).copy().reshape(shape)


@pytest.mark.parametrize(
    ("mass", "spin", "point"),
    [
        (1.0, 0.5, (0.0, 8.0, 1.1, 0.3)),
        (2.0, -1.4, (0.2, 15.0, 0.8, -0.4)),
        (0.75, 0.675, (-1.0, 6.0, 1.35, 2.0)),
    ],
)
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


@pytest.mark.parametrize(
    ("mass", "spin", "point", "expected"), LEGACY_CHRISTOFFEL_CASES
)
def test_native_connection_matches_literal_legacy_snapshot(
    native_geometry,
    mass: float,
    spin: float,
    point: tuple[float, ...],
    expected: list[list[list[float]]],
):
    """Match all 64 C connection components to the frozen legacy output."""

    c_christoffel = _call_native(
        native_geometry, "rp_kerr_christoffel", mass, spin, point, (4, 4, 4)
    )
    np.testing.assert_allclose(
        c_christoffel, np.asarray(expected), rtol=5e-13, atol=5e-14
    )
    np.testing.assert_array_equal(c_christoffel, c_christoffel.swapaxes(1, 2))
