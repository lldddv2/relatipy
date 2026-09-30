"""Compare private native integrators with SciPy using the same Kerr RHS.

Both sides evaluate the corrected C Kerr RHS, so this checks numerical
integration behavior and is not an independent oracle for Kerr physics.
Separate native and peer-validation tests check the physical equations.
"""

from __future__ import annotations

import ctypes
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[2]
DIMENSION = 8


class IntegratorConfig(ctypes.Structure):
    _fields_ = [
        ("method", ctypes.c_int),
        ("dimension", ctypes.c_size_t),
        ("relative_tolerance", ctypes.c_double),
        ("absolute_tolerance", ctypes.c_double),
        ("initial_step", ctypes.c_double),
        ("maximum_step", ctypes.c_double),
        ("maximum_steps", ctypes.c_size_t),
        ("absolute_tolerances", ctypes.POINTER(ctypes.c_double)),
        ("step_observer", ctypes.c_void_p),
        ("step_observer_context", ctypes.c_void_p),
        ("jacobian", ctypes.c_void_p),
    ]


class IntegratorStats(ctypes.Structure):
    _fields_ = [
        ("accepted_steps", ctypes.c_size_t),
        ("rejected_steps", ctypes.c_size_t),
        ("rhs_evaluations", ctypes.c_size_t),
        ("jacobian_evaluations", ctypes.c_size_t),
        ("linear_solves", ctypes.c_size_t),
        ("final_independent_variable", ctypes.c_double),
    ]


class KerrContext(ctypes.Structure):
    _fields_ = [("mass", ctypes.c_double), ("spin", ctypes.c_double)]


@pytest.fixture(scope="module")
def native_integrators(tmp_path_factory: pytest.TempPathFactory):
    """Compile the standalone native prototype into a temporary library."""

    compiler = shutil.which("cc")
    if compiler is None:
        pytest.skip("a C compiler is required for the native integrator check")
    output = tmp_path_factory.mktemp("native-integrators") / "libintegrators.so"
    sources = [
        "native/src/utils/numeric.c",
        "native/src/utils/tensor.c",
        "native/src/metric/physic/kerr.c",
        "native/src/metric/kerr.c",
        "native/src/geodesic/physic/kerr.c",
        "native/src/geodesic/kerr.c",
        "native/src/geodesic/integrators/integrator.c",
        "native/src/geodesic/integrators/dop853.c",
        "native/src/geodesic/integrators/dp45.c",
        "native/src/geodesic/integrators/radau.c",
        "native/src/geodesic/integrators/kerr.c",
    ]
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
            f"-I{ROOT / 'native' / 'src'}",
            *(str(ROOT / source) for source in sources),
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
    rhs_type = ctypes.CFUNCTYPE(
        ctypes.c_int,
        ctypes.c_double,
        pointer,
        pointer,
        ctypes.c_void_p,
    )
    library.rp_kerr_four_velocity.argtypes = [
        ctypes.c_double,
        ctypes.c_double,
        pointer,
        pointer,
        pointer,
    ]
    library.rp_kerr_four_velocity.restype = ctypes.c_int
    library.rp_kerr_null_geodesic_rhs.argtypes = [
        ctypes.c_double,
        ctypes.c_double,
        pointer,
        pointer,
        pointer,
        pointer,
    ]
    library.rp_kerr_null_geodesic_rhs.restype = ctypes.c_int
    library.rp_integrator_integrate.argtypes = [
        ctypes.POINTER(IntegratorConfig),
        rhs_type,
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_double,
        pointer,
        ctypes.POINTER(IntegratorStats),
    ]
    library.rp_integrator_integrate.restype = ctypes.c_int
    callback = rhs_type(("rp_kerr_integrator_rhs", library))
    return library, callback


def _initial_state(library) -> np.ndarray:
    coordinates = (ctypes.c_double * 4)(0.0, 8.0, 1.1, 0.3)
    coordinate_velocity = (ctypes.c_double * 3)(-0.01, 0.002, 0.02)
    four_velocity = (ctypes.c_double * 4)()
    assert library.rp_kerr_four_velocity(
        1.0, 0.5, coordinates, coordinate_velocity, four_velocity
    ) == 0
    return np.array([*coordinates, *four_velocity], dtype=float)


def _integrate_native(library, callback, method: int, initial: np.ndarray):
    config = IntegratorConfig(method, DIMENSION, 1e-10, 1e-12, 0.0, 0.02, 100000)
    context = KerrContext(1.0, 0.5)
    state = (ctypes.c_double * DIMENSION)(*initial)
    stats = IntegratorStats()
    status = library.rp_integrator_integrate(
        ctypes.byref(config),
        callback,
        ctypes.cast(ctypes.byref(context), ctypes.c_void_p),
        0.0,
        0.2,
        state,
        ctypes.byref(stats),
    )
    assert status == 0
    assert stats.accepted_steps > 0
    return np.ctypeslib.as_array(state).copy()


def test_native_dop853_and_radau_match_scipy(native_integrators):
    """Match SciPy integration of the identical corrected C RHS."""

    scipy_integrate = pytest.importorskip("scipy.integrate")
    library, callback = native_integrators
    initial = _initial_state(library)

    def rhs(_affine_parameter: float, state: np.ndarray) -> np.ndarray:
        coordinates = (ctypes.c_double * 4)(*state[:4])
        four_velocity = (ctypes.c_double * 4)(*state[4:])
        coordinate_derivative = (ctypes.c_double * 4)()
        velocity_derivative = (ctypes.c_double * 4)()
        status = library.rp_kerr_null_geodesic_rhs(
            1.0,
            0.5,
            coordinates,
            four_velocity,
            coordinate_derivative,
            velocity_derivative,
        )
        if status != 0:
            raise RuntimeError(f"native Kerr RHS failed with status {status}")
        return np.array([*coordinate_derivative, *velocity_derivative])

    reference = scipy_integrate.solve_ivp(
        rhs,
        (0.0, 0.2),
        initial,
        method="DOP853",
        rtol=1e-13,
        atol=1e-15,
    )
    assert reference.success
    expected = reference.y[:, -1]
    dop853 = _integrate_native(library, callback, 0, initial)
    radau = _integrate_native(library, callback, 1, initial)
    np.testing.assert_allclose(dop853, expected, rtol=2e-9, atol=2e-11)
    np.testing.assert_allclose(radau, expected, rtol=2e-9, atol=2e-11)
