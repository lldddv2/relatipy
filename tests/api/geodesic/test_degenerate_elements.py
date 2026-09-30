"""Degenerate osculating elements preserve valid physical trajectories."""

from dataclasses import FrozenInstanceError

import numpy as np
import pytest
from astropy import units as u
from astropy.constants import c

from relatipy.coordinates import OrbitalElements
from relatipy.geodesic import Solution
from relatipy.geodesic.orbit import _state_from_native
from relatipy.metrics.kerr import Kerr


_RADIAL_DEFINED = np.array([True, True, False, False, False, False])
_PHYSICAL_FIELDS = (
    "tau", "t", "xyz", "vxyz", "uxyz", "r", "theta", "phi", "R",
    "Theta", "Phi", "vr", "vtheta", "vphi", "vR", "vTheta", "vPhi",
    "ut", "ur", "utheta", "uphi", "uR", "uTheta", "uPhi",
)


def _radial_orbit():
    metric = Kerr(mass=1 * u.Msun, spin=0)
    return metric.orbit(x=20 * metric.r_g, vx=-0.05 * c)


def _assert_radial(state) -> None:
    for name in _PHYSICAL_FIELDS:
        assert np.all(np.isfinite(getattr(state, name).value)), name
    elements = state.orbital_elements()
    expected = np.broadcast_to(_RADIAL_DEFINED, elements.a.shape + (6,))
    np.testing.assert_array_equal(elements.defined, expected)
    for name in ("inc", "Omega", "omega", "f"):
        assert np.all(np.isnan(getattr(elements, name).value)), name
    assert np.all(np.isfinite(elements.a.value))
    np.testing.assert_allclose(elements.e, 1.0, rtol=0, atol=1e-14)
    with pytest.raises(ValueError):
        elements.defined.setflags(write=True)


def test_defined_mask_is_strongly_immutable_and_keeps_six_argument_constructor():
    elements = OrbitalElements(
        np.inf * u.km, 1.0, np.nan * u.rad, np.nan * u.rad,
        np.nan * u.rad, np.nan * u.rad,
    )
    assert elements.defined.dtype == np.bool_
    assert elements.defined.shape == (6,)
    np.testing.assert_array_equal(elements.defined, _RADIAL_DEFINED)
    with pytest.raises(ValueError):
        elements.defined[0] = False
    with pytest.raises(ValueError):
        elements.defined.setflags(write=True)
    with pytest.raises((FrozenInstanceError, AttributeError, TypeError)):
        elements.defined = np.ones(6, dtype=bool)


def test_defined_mask_reports_each_series_field_without_reclassifying_data():
    elements = OrbitalElements(
        [10, np.inf, -np.inf, np.nan] * u.km,
        [0.2, 1, np.inf, np.nan],
        [0, np.nan, np.inf, 0.1] * u.rad,
        [0, np.nan, 0, 0.2] * u.rad,
        [0, np.nan, 0, 0.3] * u.rad,
        [0, np.nan, 0, 0.4] * u.rad,
    )
    np.testing.assert_array_equal(
        elements.defined,
        [
            [True] * 6,
            _RADIAL_DEFINED,
            [False, False, False, True, True, True],
            [False, False, True, True, True, True],
        ],
    )
    with pytest.raises(ValueError):
        elements.defined[:, :2].setflags(write=True)


@pytest.mark.parametrize("method", ["radau", "dop853", "dp45"])
def test_radial_schwarzschild_initial_integrate_solve_copy_and_reset(method):
    orbit = _radial_orbit()
    initial_xyz = orbit.xyz.copy()
    _assert_radial(orbit.initial)
    orbit.integrate(0.04 * orbit._metric._time_scale, method=method)
    _assert_radial(orbit)
    assert orbit.x < initial_xyz[0]
    clone = orbit.copy()
    _assert_radial(clone)
    clone.reset()
    _assert_radial(clone.initial)
    np.testing.assert_array_equal(clone.xyz, initial_xyz)
    assert orbit.tau != clone.tau
    orbit.reset()
    np.testing.assert_array_equal(orbit.xyz, initial_xyz)
    times = np.linspace(0, 0.04, 5) * orbit._metric._time_scale
    solution = orbit.solve(tau_eval=times, method=method)
    assert solution.success
    assert solution.status == 0
    _assert_radial(solution)
    _assert_radial(solution[2])
    _assert_radial(solution[[4, 0, 2]])
    _assert_radial(solution.at(tau=times[2]))
    _assert_radial(solution.at(tau=[times[2].value] * times.unit))
    assert orbit.tau == orbit.initial.tau


def test_radial_solution_interpolation_preserves_missing_angles_and_physical_views():
    orbit = _radial_orbit()
    scale = orbit._metric._time_scale
    solution = orbit.solve(tau_eval=np.linspace(0, 0.04, 5) * scale, method="dp45")
    selected = solution.at(tau=np.array([0.02, 0.015, 0.015, 0]) * scale)
    _assert_radial(selected)
    np.testing.assert_array_equal(selected.xyz[[0, 3]], solution.xyz[[2, 0]])
    np.testing.assert_array_equal(selected.xyz[1], selected.xyz[2])
    _assert_radial(solution.at(tau=0.015 * scale))
    _assert_radial(solution.at(t=0.5 * (solution.t[1] + solution.t[2])))
    np.testing.assert_allclose(
        selected.uxyz.to_value(u.m / u.s),
        selected.ut.value[:, None] * selected.vxyz.to_value(u.m / u.s),
        rtol=1e-13,
    )


def test_near_radial_initial_state_keeps_finite_physics_and_undefined_angles():
    metric = Kerr(mass=1 * u.Msun, spin=0)
    orbit = metric.orbit(x=20 * metric.r_g, vx=-0.05 * c, vy=1e-17 * c)
    _assert_radial(orbit.initial)
    orbit.integrate(0.01 * metric._time_scale, method="dp45")
    _assert_radial(orbit)


def test_native_mixed_rows_keep_scalar_series_and_exact_selection_masks():
    from relatipy import _core

    metric = Kerr(mass=1 * u.Msun, spin=0)
    inputs = (
        np.array([0, 20, 0, 0, 0, 0.1, 0]),
        np.array([0, 20, 0, 0, -0.05, 0, 0]),
        np.array([0, 8, 0, 0, 0, 0.5, 0]),
    )
    canonical = np.stack([_core.initial_kerr(0, "cartesian", row) for row in inputs])
    state = _state_from_native(
        metric, canonical, np.array([0., 1., 2.]),
        length_unit=u.km, time_unit=u.s,
    )
    expected = np.array([[True] * 6, _RADIAL_DEFINED, [True] * 6])
    np.testing.assert_array_equal(state.orbital_elements().defined, expected)
    assert np.isposinf(state.orbital_elements().a.value[2])
    source = _radial_orbit().solve(
        tau_eval=np.array([0, 0.01]) * metric._time_scale, method="dp45"
    )
    solution = Solution(
        state=state, integration=source.integration, status=0, message="mixed elements",
        _metric=metric,
    )
    np.testing.assert_array_equal(solution[1].orbital_elements().defined, expected[1])
    np.testing.assert_array_equal(
        solution[[2, 1, 0]].orbital_elements().defined, expected[[2, 1, 0]]
    )
    np.testing.assert_array_equal(
        solution.at(tau=[2, 1, 2] * metric._time_scale).orbital_elements().defined,
        expected[[2, 1, 2]],
    )
    selected = solution.at(tau=np.array([2, 0.5, 1, 1.5, 0]) * metric._time_scale)
    np.testing.assert_array_equal(
        selected.orbital_elements().defined[[0, 2, 4]], expected[[2, 1, 0]]
    )
    assert np.all(selected.orbital_elements().defined[[1, 3]])
    assert np.isposinf(selected.orbital_elements().a.value[0])
    for name in _PHYSICAL_FIELDS:
        assert np.all(np.isfinite(getattr(selected, name).value)), name
    with pytest.raises(ValueError):
        selected.orbital_elements().defined.setflags(write=True)


@pytest.mark.parametrize("method", ["radau", "dp45"])
def test_exactly_parabolic_initial_state_does_not_block_physical_evolution(method):
    metric = Kerr(mass=1 * u.Msun, spin=0)
    orbit = metric.orbit(x=8 * metric.r_g, vy=0.5 * c)
    assert np.isposinf(orbit.orbital_elements().a.value)
    times = np.linspace(0, 0.04, 5) * metric._time_scale
    solution = orbit.solve(tau_eval=times, method=method)
    assert solution.success
    assert np.isposinf(solution.orbital_elements().a.value[0])
    assert np.all(solution.orbital_elements().defined)
    selected = solution.at(tau=0.015 * metric._time_scale)
    assert np.all(selected.orbital_elements().defined)
    for name in _PHYSICAL_FIELDS:
        assert np.all(np.isfinite(getattr(solution, name).value)), name
        assert np.all(np.isfinite(getattr(selected, name).value)), name
    orbit.integrate(times[-1], method=method)
    assert np.all(orbit.orbital_elements().defined)
    clone = orbit.copy()
    clone.reset()
    assert np.isposinf(clone.orbital_elements().a.value)
    assert np.all(clone.orbital_elements().defined)
    orbit.reset()
    assert np.isposinf(orbit.orbital_elements().a.value)


@pytest.mark.parametrize("family, column, value", [
    ("cartesian", 0, np.nan), ("spherical", 6, np.inf),
    ("elements", 0, np.nan), ("elements", 0, -np.inf),
    ("elements", 1, np.nan), ("elements", 1, np.inf),
    ("elements", 2, np.inf), ("elements", 5, -np.inf),
])
def test_invalid_native_values_cannot_enter_orbit_or_interpolated_state(
    monkeypatch, family, column, value,
):
    from relatipy import _core

    metric = Kerr(mass=1 * u.Msun, spin=0)
    orbit = metric.orbit(
        R=20 * metric.r_g, Theta=90 * u.deg, Phi=0 * u.rad,
        vR=-0.05 * c,
    )
    scale = orbit._metric._time_scale
    solution = orbit.solve(tau_eval=np.linspace(0, 0.04, 5) * scale, method="dp45")
    reconstruct = _core.reconstruct_canonical_family_batch

    def corrupt(*args):
        rows, status = reconstruct(*args)
        rows = rows.copy()
        if args[2] == family:
            rows[0, column] = value
        return rows, status

    monkeypatch.setattr(_core, "reconstruct_canonical_family_batch", corrupt)
    accessor = {
        "cartesian": lambda state: state.txyz,
        "spherical": lambda state: state.trqp,
        "elements": lambda state: state.orbital_elements(),
    }[family]
    with pytest.raises(ValueError, match="non-finite"):
        accessor(solution)


def test_invalid_cartesian_native_values_remain_rejected(monkeypatch):
    from relatipy import _core

    reconstruct = _core.reconstruct_canonical_family_batch

    def corrupt(*args):
        rows, status = reconstruct(*args)
        rows = rows.copy()
        rows[0, 1] = np.nan
        return rows, status

    monkeypatch.setattr(_core, "reconstruct_canonical_family_batch", corrupt)
    metric = Kerr(mass=1 * u.Msun, spin=0)
    orbit = metric.orbit(R=20 * metric.r_g, Theta=90 * u.deg, Phi=0 * u.rad)
    with pytest.raises(ValueError, match="non-finite"):
        _ = orbit.xyz


def test_preview_rejects_radial_and_exactly_parabolic_states_explicitly():
    radial = _radial_orbit()
    with pytest.raises(ValueError, match="near-radial states have undefined angles"):
        radial.preview()
    metric = Kerr(mass=1 * u.Msun, spin=0)
    parabolic = metric.orbit(x=8 * metric.r_g, vy=0.5 * c)
    assert np.isposinf(parabolic.orbital_elements().a.value)
    np.testing.assert_array_equal(parabolic.orbital_elements().defined, [True] * 6)
    with pytest.raises(ValueError, match="exactly parabolic states have a = \\+inf"):
        parabolic.preview()
