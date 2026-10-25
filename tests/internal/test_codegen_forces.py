"""Run-time forces compile into rp_external_acceleration_fn callbacks."""

from __future__ import annotations

import ctypes

import numpy as np
import pytest

sp = pytest.importorskip("sympy")

from relatipy import _jit  # noqa: E402
from relatipy._codegen import build_force_source  # noqa: E402
from relatipy._codegen.forces import parse_expression  # noqa: E402


def _has_compiler() -> bool:
    try:
        _jit.compiler_command()
    except _jit.CompilerNotFoundError:
        return False
    return True


needs_cc = pytest.mark.skipif(not _has_compiler(), reason="no C compiler")
_DOUBLES = ctypes.POINTER(ctypes.c_double)
_CALLBACK = ctypes.CFUNCTYPE(ctypes.c_int, _DOUBLES, _DOUBLES, ctypes.c_void_p, _DOUBLES)


class Callback:
    """ctypes view of a compiled rp_external_acceleration_fn."""

    def __init__(self, force, cache_dir) -> None:
        self.force = force
        self.function = _jit.compile_function(force.source, force.function_name, cache_dir=cache_dir)
        self.call = _CALLBACK(self.function.address)

    def __call__(self, position, velocity, parameters=()):
        position = np.ascontiguousarray(position, dtype=float)
        velocity = np.ascontiguousarray(velocity, dtype=float)
        model = np.ascontiguousarray(parameters, dtype=float)
        out = np.zeros(3)
        status = self.call(
            position.ctypes.data_as(_DOUBLES), velocity.ctypes.data_as(_DOUBLES),
            model.ctypes.data if model.size else None, out.ctypes.data_as(_DOUBLES),
        )
        return status, out


def _points(count: int = 32, seed: int = 4):
    generator = np.random.default_rng(seed)
    positions = generator.normal(scale=5.0, size=(count, 3))
    positions[0] = (0.0, 0.0, 2.0)  # on the spin axis
    positions[1] = (3.0, -1.0, 0.0)  # in the equatorial plane
    velocities = generator.normal(scale=0.3, size=(count, 3))
    return positions, velocities


def _plummer(position, gm, b):
    return -gm * position / (position @ position + b * b) ** 1.5


@needs_cc
def test_uniform_field_reads_parameters_at_run_time(tmp_path) -> None:
    callback = Callback(
        build_force_source(["0", "0", "-g"], kind="force", parameters=("g",)), tmp_path,
    )
    positions, velocities = _points()
    for g in (0.5, 2.0e-3):
        for position, velocity in zip(positions, velocities):
            status, acceleration = callback(position, velocity, [g])
            assert status == 0
            np.testing.assert_array_equal(acceleration, [0.0, 0.0, -g])
    assert "(void)position;" in callback.force.source and "(void)velocity;" in callback.force.source


@needs_cc
def test_uniform_field_as_potential(tmp_path) -> None:
    callback = Callback(
        build_force_source("g*z", kind="potential", parameters=("g",)), tmp_path,
    )
    status, acceleration = callback([1.0, 2.0, 3.0], [0.0, 0.0, 0.0], [0.25])
    assert status == 0
    np.testing.assert_array_equal(acceleration, [0.0, 0.0, -0.25])


@needs_cc
@pytest.mark.parametrize(
    ("definition", "coords"),
    [
        ("-GM/sqrt(R**2 + z**2 + b**2)", "cylindrical"),
        (lambda r, theta, phi, GM, b: -GM / sp.sqrt(r**2 + b**2), "spherical"),
        ("-GM/sqrt(x**2 + y**2 + z**2 + b**2)", "cartesian"),
    ],
)
def test_plummer_sphere_matches_numpy(tmp_path, definition, coords) -> None:
    force = build_force_source(definition, kind="potential", coords=coords, parameters=("GM", "b"))
    assert force.parameters == ("GM", "b")
    callback = Callback(force, tmp_path)
    positions, velocities = _points()
    for gm, b in ((1.0e-2, 0.5), (3.0, 2.0)):
        for position, velocity in zip(positions, velocities):
            status, acceleration = callback(position, velocity, [gm, b])
            assert status == 0
            np.testing.assert_allclose(acceleration, _plummer(position, gm, b), rtol=1e-14, atol=0)


@needs_cc
def test_plummer_potential_with_azimuthal_and_polar_terms(tmp_path) -> None:
    # A potential written with explicit angles: Phi = -GM / sqrt(r^2 + b^2)
    # + eps * r^2 * sin(theta)^2 * cos(2 phi), a quadrupolar perturbation.
    force = build_force_source(
        "-GM/sqrt(r**2 + b**2) + eps*r**2*sin(theta)**2*cos(2*phi)",
        kind="potential", coords="spherical", parameters=("GM", "b", "eps"),
    )
    callback = Callback(force, tmp_path)
    positions, velocities = _points()
    gm, b, eps = 1.0, 0.7, 1.0e-3
    for position, velocity in zip(positions, velocities):
        x, y, _ = position
        expected = _plummer(position, gm, b) - eps * np.array([2.0 * x, -2.0 * y, 0.0])
        status, acceleration = callback(position, velocity, [gm, b, eps])
        assert status == 0
        np.testing.assert_allclose(acceleration, expected, rtol=1e-13, atol=1e-16)


@needs_cc
@pytest.mark.parametrize(
    ("definition", "coords"),
    [
        ("(-gamma*vx, -gamma*vy, -gamma*vz)", "cartesian"),
        (["-gamma*v_R", "-gamma*v_phi", "-gamma*v_z"], "cylindrical"),
        (lambda r, theta, phi, v_r, v_theta, v_phi, gamma: (-gamma * v_r, -gamma * v_theta,
                                                            -gamma * v_phi), "spherical"),
    ],
)
def test_velocity_dependent_drag_matches_numpy(tmp_path, definition, coords) -> None:
    callback = Callback(
        build_force_source(definition, kind="force", coords=coords, parameters=("gamma",)), tmp_path,
    )
    positions, velocities = _points()
    for position, velocity in zip(positions[1:], velocities[1:]):
        status, acceleration = callback(position, velocity, [0.3])
        assert status == 0
        np.testing.assert_allclose(acceleration, -0.3 * velocity, rtol=1e-13, atol=1e-16)


@needs_cc
def test_radial_drag_in_spherical_coordinates(tmp_path) -> None:
    callback = Callback(
        build_force_source("(-k*v_r, 0, 0)", kind="force", coords="spherical", parameters=("k",)),
        tmp_path,
    )
    positions, velocities = _points()
    for position, velocity in zip(positions, velocities):
        unit = position / np.linalg.norm(position)
        status, acceleration = callback(position, velocity, [2.0])
        assert status == 0
        np.testing.assert_allclose(acceleration, -2.0 * (velocity @ unit) * unit, rtol=1e-13,
                                   atol=1e-16)


@needs_cc
def test_non_finite_acceleration_returns_nonzero(tmp_path) -> None:
    callback = Callback(build_force_source("1/x", kind="potential"), tmp_path)
    status, _ = callback([0.0, 1.0, 1.0], [0.0, 0.0, 0.0])
    assert status == 1
    status, acceleration = callback([2.0, 1.0, 1.0], [0.0, 0.0, 0.0])
    assert status == 0
    np.testing.assert_allclose(acceleration, [0.25, 0.0, 0.0])


@pytest.mark.parametrize(
    "text",
    [
        "__import__('os').system('true')",
        "x.real",
        "(lambda: 1)()",
        "'text'",
        "x ^ 2",
        "y + w",
        "sqrt(x, evaluate=False)",
        "[i for i in (1, 2)]",
        "x if y else z",
    ],
)
def test_restricted_parser_rejects_other_constructs(text) -> None:
    x, y, z = sp.symbols("x y z")
    with pytest.raises(ValueError):
        parse_expression(text, {"x": x, "y": y, "z": z})


def test_restricted_parser_accepts_math() -> None:
    x = sp.Symbol("x")
    value = parse_expression("2*atan2(x, 1) + exp(-x**2)/sqrt(pi) - 0.5", {"x": x})
    assert value == 2 * sp.atan2(x, 1) + sp.exp(-x**2) / sp.sqrt(sp.pi) - sp.Float("0.5")


def test_invalid_definitions_raise() -> None:
    with pytest.raises(ValueError, match="kind"):
        build_force_source("x", kind="field")
    with pytest.raises(ValueError, match="coords"):
        build_force_source("x", kind="potential", coords="polar")
    with pytest.raises(ValueError, match="reserved"):
        build_force_source("x", kind="potential", parameters=("model",))
    with pytest.raises(ValueError, match="reserved"):
        build_force_source("x", kind="potential", coords="cylindrical", parameters=("R",))
    with pytest.raises(ValueError, match="three components"):
        build_force_source(["x", "y"], kind="force")
    with pytest.raises(ValueError, match="unknown name 'vx'"):
        build_force_source("vx", kind="potential")
    with pytest.raises(ValueError, match="undeclared names: w"):
        build_force_source(lambda x, y, z: sp.Symbol("w") * x, kind="potential")
    with pytest.raises(ValueError, match="scalar"):
        build_force_source("(x, y, z)", kind="potential")


def test_source_documents_the_parameter_layout() -> None:
    force = build_force_source("-GM/R", kind="potential", coords="cylindrical",
                               parameters=("GM",), name="rp_test_force")
    assert force.function_name == "rp_test_force"
    assert "model[0] = GM" in force.source
    assert "int rp_test_force(" in force.source
    assert "const double *parameters = (const double *)model;" in force.source
    assert force.report.after.total <= force.report.before.total
