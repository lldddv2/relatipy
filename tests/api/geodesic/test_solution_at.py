"""Contract tests for exact selection and post-integration interpolation."""

from __future__ import annotations

import sys
from concurrent.futures import ThreadPoolExecutor
from types import SimpleNamespace

import numpy as np
import pytest
from astropy import units as u
from astropy.constants import G, c

from relatipy.geodesic import Solution
from relatipy.metrics import Kerr
from relatipy.geodesic.interpolation import interpolate_cartesian

from support.builders import make_series_state, make_solution


@pytest.mark.parametrize(
    ("query", "expected"),
    [
        ({"tau": 2 * u.s}, [2]),
        ({"t": 2 * u.s}, [2]),
        ({"tau": [3, 1, 1] * u.s}, [3, 1, 1]),
        ({"t": [3, 1, 1] * u.s}, [3, 1, 1]),
    ],
)
def test_exact_queries_select_original_samples(query, expected) -> None:
    solution = make_solution()
    selected = solution.at(**query)
    original = solution[expected[0] if len(expected) == 1 else expected]
    assert selected.tau.shape == (() if len(expected) == 1 else (len(expected),))
    for name in ("xyz", "vxyz", "R", "uR", "ut"):
        np.testing.assert_array_equal(getattr(selected, name), getattr(original, name))
    np.testing.assert_array_equal(
        selected.orbital_elements().a, original.orbital_elements().a
    )
    with pytest.raises(ValueError):
        selected.xyz.value.flags.writeable = True


def test_scalar_and_single_element_array_keep_distinct_shapes() -> None:
    solution = make_solution()
    assert solution.at(tau=1 * u.s).xyz.shape == (3,)
    assert solution.at(tau=[1] * u.s).xyz.shape == (1, 3)


@pytest.mark.parametrize(
    "query, error",
    [
        ({}, ValueError),
        ({"tau": 1 * u.s, "t": 1 * u.s}, ValueError),
        ({"tau": 1.0}, TypeError),
        ({"tau": [[1]] * u.s}, ValueError),
        ({"tau": np.nan * u.s}, ValueError),
        ({"tau": -1 * u.s}, ValueError),
        ({"tau": [1, 5] * u.s}, ValueError),
        ({"t": 5 * u.s}, ValueError),
    ],
)
def test_invalid_queries_fail_as_whole(query, error) -> None:
    with pytest.raises(error):
        make_solution().at(**query)


def test_between_samples_needs_geometry_context() -> None:
    with pytest.raises(ValueError, match="requires the Kerr metric"):
        make_solution().at(tau=0.5 * u.s)


def test_single_sample_allows_only_exact_queries() -> None:
    source = make_solution()
    solution = Solution(
        state=make_series_state(1),
        integration=source.integration,
        status=0,
        message="single sample",
    )
    assert solution.at(tau=0 * u.s).tau.shape == ()
    assert solution.at(t=[0] * u.s).tau.shape == (1,)
    with pytest.raises(ValueError, match="outside"):
        solution.at(tau=0.1 * u.s)


@pytest.mark.parametrize("metric", [object(), 1 * u.M_sun, 0.5])
def test_invalid_private_metric_is_rejected(metric) -> None:
    source = make_solution()
    with pytest.raises(TypeError, match="_metric must be a Kerr metric"):
        Solution(
            state=source._state,
            integration=source.integration,
            status=0,
            message="invalid context",
            _metric=metric,
        )


def test_private_metric_is_stored_by_reference() -> None:
    source = make_solution()
    metric = Kerr(mass=1 * u.M_sun, spin=0.2)
    solution = Solution(
        state=source._state,
        integration=source.integration,
        status=0,
        message="metric context",
        _metric=metric,
    )
    assert solution._metric is metric
    with pytest.raises(AttributeError):
        solution._metric = Kerr(mass=2 * u.M_sun, spin=0.2)


def test_cartesian_position_and_velocity_are_interpolated_separately() -> None:
    solution = make_solution()
    data = interpolate_cartesian(solution._state, tau=[2.5, 0.5, 2.5] * u.s)
    np.testing.assert_allclose(data.xyz[:, 0], [3.5, 1.5, 3.5])
    np.testing.assert_allclose(data.vxyz[:, 0], [6.5, 4.5, 6.5])
    np.testing.assert_allclose(data.t, [2.5, 0.5, 2.5])
    np.testing.assert_array_equal(data.exact_indices, [-1, -1, -1])


def test_coordinate_time_inverse_uses_monotone_pchip() -> None:
    solution = make_solution()
    t_samples = np.array([0.0, 0.5, 2.0, 5.0])
    source = solution._state
    from relatipy.coordinates import (
        BoyerLindquistCoordinates,
        CartesianCoordinates,
        SphericalCoordinates,
    )
    from relatipy.geodesic import State

    coordinate_time = t_samples * u.s
    state = State(
        tau=source.tau,
        txyz=CartesianCoordinates(coordinate_time, source.x, source.y, source.z),
        trqp=SphericalCoordinates(coordinate_time, source.r, source.theta, source.phi),
        tRQP=BoyerLindquistCoordinates(
            coordinate_time, source.R, source.Theta, source.Phi
        ),
        vxyz=source.vxyz,
        vrqp=source.vrqp,
        vRQP=source.vRQP,
        ut=source.ut,
        uxyz=source.uxyz,
        urqp=source.urqp,
        uRQP=source.uRQP,
        orbital_elements=source.orbital_elements(),
    )
    values = interpolate_cartesian(state, t=[1.0, 0.5, 1.0] * u.s)
    np.testing.assert_allclose(values.t, [1.0, 0.5, 1.0], rtol=0, atol=1e-12)
    np.testing.assert_array_equal(values.exact_indices, [-1, 1, -1])
    assert values.tau[0] == values.tau[2]
    assert 1.0 < values.tau[0] < 2.0


def test_one_native_batch_preserves_exact_rows(monkeypatch) -> None:
    source = make_solution()
    solution = Solution(
        state=source._state,
        integration=source.integration,
        status=source.status,
        message=source.message,
        _metric=Kerr(mass=1 * u.M_sun, spin=0.2),
    )
    calls = []

    def reconstruct_batch(spin, cartesian):
        calls.append((spin, cartesian.copy()))
        return (
            np.ones((cartesian.shape[0], 29)),
            np.zeros(cartesian.shape[0], dtype=np.int32),
        )

    monkeypatch.setitem(
        sys.modules,
        "relatipy._core",
        SimpleNamespace(reconstruct_batch=reconstruct_batch),
    )
    monkeypatch.setattr(sys.modules["relatipy"], "_core", sys.modules["relatipy._core"], raising=False)
    selected = solution.at(tau=[2, 0.5, 2, 1.5] * u.s)
    assert len(calls) == 1
    assert calls[0][1].shape == (2, 7)
    np.testing.assert_array_equal(selected.xyz[[0, 2]], source.xyz[[2, 2]])
    np.testing.assert_array_equal(selected.R[[0, 2]], source.R[[2, 2]])
    np.testing.assert_array_equal(
        selected.orbital_elements().a[[0, 2]], source.orbital_elements().a[[2, 2]]
    )
    scale = (G * (1 * u.M_sun) / c**2).to(u.km).value
    np.testing.assert_allclose(selected.R[[1, 3]].to_value(u.km), scale)
    with pytest.raises(ValueError):
        selected.xyz.value.flags.writeable = True


def test_native_failure_reports_original_query_row(monkeypatch) -> None:
    source = make_solution()
    solution = Solution(
        state=source._state,
        integration=source.integration,
        status=0,
        message="test",
        _metric=Kerr(mass=1 * u.M_sun, spin=0),
    )

    def fail(spin, cartesian):
        return np.zeros((1, 29)), np.array([7], dtype=np.int32)

    monkeypatch.setitem(
        sys.modules, "relatipy._core", SimpleNamespace(reconstruct_batch=fail)
    )
    monkeypatch.setattr(sys.modules["relatipy"], "_core", sys.modules["relatipy._core"], raising=False)
    with pytest.raises(ValueError, match="query row 1.*native status 7"):
        solution.at(tau=[2, 0.5] * u.s)


def test_compiled_reconstruction_preserves_cartesian_relations() -> None:
    pytest.importorskip("relatipy._core")
    source = make_solution()
    solution = Solution(
        state=source._state,
        integration=source.integration,
        status=0,
        message="native reconstruction",
        _metric=Kerr(mass=1e29 * u.kg, spin=0.2),
    )
    selected = solution.at(tau=[0.5, 1.5] * u.s)
    assert selected.xyz.shape == (2, 3)
    assert selected.uxyz.shape == (2, 3)
    np.testing.assert_allclose(
        selected.r.to_value(u.km),
        np.linalg.norm(selected.xyz.to_value(u.km), axis=1),
        rtol=1e-13,
    )
    np.testing.assert_allclose(
        selected.uxyz.to_value(u.km / u.s),
        selected.ut.to_value(u.one)[:, None] * selected.vxyz.to_value(u.km / u.s),
        rtol=1e-13,
    )
    assert np.all(selected.ut > 1 * u.one)


def test_compiled_reconstruction_is_reentrant_across_threads() -> None:
    native = pytest.importorskip("relatipy._core")
    batches = []
    for radius, velocity in ((20.0, 0.15), (25.0, 0.14)):
        rows = np.empty((512, 7), dtype=np.float64)
        rows[:, 0] = np.linspace(0, 1, 512)
        rows[:, 1] = radius + np.linspace(0, 0.1, 512)
        rows[:, 2] = 3.0
        rows[:, 3] = 2.0
        rows[:, 4] = 0.01
        rows[:, 5] = velocity
        rows[:, 6] = 0.02
        batches.append(rows)

    sequential = [native.reconstruct_batch(0.2, rows) for rows in batches]
    with ThreadPoolExecutor(max_workers=2) as pool:
        parallel = list(
            pool.map(lambda rows: native.reconstruct_batch(0.2, rows), batches)
        )
    for (expected, expected_status), (actual, actual_status) in zip(
        sequential, parallel, strict=True
    ):
        np.testing.assert_array_equal(expected_status, np.zeros(512, dtype=np.int32))
        np.testing.assert_array_equal(actual_status, expected_status)
        np.testing.assert_array_equal(actual, expected)


def test_interpolated_azimuth_keeps_unwrapped_sample_branch() -> None:
    pytest.importorskip("relatipy._core")
    from relatipy.coordinates import (
        BoyerLindquistCoordinates,
        CartesianCoordinates,
        SphericalCoordinates,
    )
    from relatipy.geodesic import State

    source = make_solution()
    state = source._state
    phase = np.array([3.0, 3.1, 3.2, 3.3])
    coordinate_time = state.t
    custom = State(
        tau=state.tau,
        txyz=CartesianCoordinates(
            coordinate_time,
            20 * np.cos(phase) * u.km,
            20 * np.sin(phase) * u.km,
            np.full(4, 2.0) * u.km,
        ),
        trqp=SphericalCoordinates(
            coordinate_time, state.r, state.theta, phase * u.rad
        ),
        tRQP=BoyerLindquistCoordinates(
            coordinate_time, state.R, state.Theta, phase * u.rad
        ),
        vxyz=np.tile([0.0, 10.0, 1.0], (4, 1)) * u.km / u.s,
        vrqp=state.vrqp,
        vRQP=state.vRQP,
        ut=state.ut,
        uxyz=state.uxyz,
        urqp=state.urqp,
        uRQP=state.uRQP,
        orbital_elements=state.orbital_elements(),
    )
    solution = Solution(
        state=custom,
        integration=source.integration,
        status=0,
        message="unwrapped phase",
        _metric=Kerr(mass=1e29 * u.kg, spin=0.2),
    )
    selected = solution.at(tau=[1, 1.5, 2] * u.s)
    np.testing.assert_allclose(selected.Phi.to_value(u.rad)[[0, 2]], [3.1, 3.2])
    np.testing.assert_allclose(selected.phi.to_value(u.rad)[[0, 2]], [3.1, 3.2])
    assert 3.1 < selected.Phi.to_value(u.rad)[1] < 3.2
    assert 3.1 < selected.phi.to_value(u.rad)[1] < 3.2
