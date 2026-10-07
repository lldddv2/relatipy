"""Null states and solutions: immutable views, selection and queries in t."""

from __future__ import annotations

import numpy as np
import pytest
from astropy import units as u
from astropy.constants import c
from scipy.optimize import brentq

from relatipy import Kerr


_COMPONENTS = (
    "t", "x", "y", "z", "r", "theta", "phi", "R", "Theta", "Phi",
    "vx", "vy", "vz", "vr", "vtheta", "vphi", "vR", "vTheta", "vPhi",
)
_VECTORS = ("xyz", "vxyz")
_NULL_ABSENT = (
    "tau", "ut", "uxyz", "urqp", "uRQP", "ux", "uy", "uz", "ur",
    "utheta", "uphi", "uR", "uTheta", "uPhi", "orbital_elements", "preview",
)


@pytest.fixture(scope="module")
def radial_null():
    metric = Kerr(mass=1 * u.Msun, spin=0)
    return metric.null(
        R=10 * metric.r_g, Theta=np.pi / 2 * u.rad, Phi=0 * u.rad, vR=c,
    )


@pytest.fixture(scope="module")
def radial_solution(radial_null):
    time_scale = (radial_null.R / (10 * c)).to(u.s)
    samples = np.linspace(0, 2, 41) * time_scale
    solution = radial_null.solve(t_eval=samples, method="dp45")
    assert solution.status == 0
    return solution


def _quantity_views(state):
    """Include component arrays inside each public coordinate container."""
    result = {
        name: getattr(state, name) for name in (*_COMPONENTS, *_VECTORS)
    }
    for container_name, components in (
        ("txyz", ("t", "x", "y", "z", "xyz")),
        ("trqp", ("t", "r", "theta", "phi")),
        ("tRQP", ("t", "R", "Theta", "Phi")),
        ("vrqp", ("vr", "vtheta", "vphi")),
        ("vRQP", ("vR", "vTheta", "vPhi")),
        ("state_vector", ("xyz", "vxyz")),
    ):
        container = getattr(state, container_name)
        for name in components:
            result[f"{container_name}.{name}"] = getattr(container, name)
    return result


def test_null_exposes_scalar_coordinate_views(radial_null) -> None:
    for name in _COMPONENTS:
        quantity = getattr(radial_null, name)
        assert isinstance(quantity, u.Quantity), name
        assert quantity.shape == (), name
    for name in _VECTORS:
        assert getattr(radial_null, name).shape == (3,), name
    np.testing.assert_array_equal(radial_null.txyz.xyz, radial_null.xyz)
    np.testing.assert_array_equal(radial_null.state_vector.xyz, radial_null.xyz)
    np.testing.assert_array_equal(radial_null.state_vector.vxyz, radial_null.vxyz)
    assert radial_null.txyz.t == radial_null.trqp.t == radial_null.tRQP.t


@pytest.mark.parametrize("name", _NULL_ABSENT)
def test_null_does_not_expose_orbit_only_api(radial_null, name) -> None:
    assert not hasattr(radial_null, name)
    with pytest.raises(AttributeError):
        getattr(radial_null, name)


@pytest.mark.parametrize("name", ["b", "eta"])
def test_null_invariants_are_read_only(radial_null, name) -> None:
    invariant = getattr(radial_null, name)
    assert invariant.shape == ()
    assert invariant.value.flags.writeable is False
    with pytest.raises(AttributeError):
        setattr(radial_null, name, invariant)


@pytest.mark.parametrize(
    "kind", ["null", "initial", "solution", "indexed", "sliced", "interpolated"]
)
def test_null_public_views_are_read_only(radial_null, radial_solution, kind) -> None:
    states = {
        "null": radial_null,
        "initial": radial_null.initial,
        "solution": radial_solution,
        "indexed": radial_solution[1],
        "sliced": radial_solution[1:4],
        "interpolated": radial_solution.at(
            (radial_solution.t[1] + radial_solution.t[2]) / 2
        ),
    }
    for name, quantity in _quantity_views(states[kind]).items():
        assert isinstance(quantity, u.Quantity), name
        assert quantity.value.flags.writeable is False, name
    with pytest.raises(ValueError):
        states[kind].xyz[...] = 0 * states[kind].xyz.unit


def test_null_solution_series_and_scalar_index_shapes(radial_solution) -> None:
    assert len(radial_solution) == 41
    assert np.all(np.diff(radial_solution.t.to_value(u.s)) > 0)
    for name in _COMPONENTS:
        assert getattr(radial_solution, name).shape == (41,), name
        assert getattr(radial_solution[-1], name).shape == (), name
    for name in _VECTORS:
        assert getattr(radial_solution, name).shape == (41, 3), name
        assert getattr(radial_solution[-1], name).shape == (3,), name
    assert radial_solution[-1].t == radial_solution.t[-1]


@pytest.mark.parametrize(
    ("index", "expected_indices"),
    [
        (slice(1, 4), [1, 2, 3]),
        (np.array([3, 1, 1]), [3, 1, 1]),
        (np.arange(41) % 10 == 0, [0, 10, 20, 30, 40]),
    ],
)
def test_null_solution_selection_keeps_numpy_order(
    radial_solution, index, expected_indices,
) -> None:
    selected = radial_solution[index]
    assert selected.t.shape == (len(expected_indices),)
    assert selected.xyz.shape == (len(expected_indices), 3)
    for name in (*_COMPONENTS, *_VECTORS):
        np.testing.assert_array_equal(
            getattr(selected, name), getattr(radial_solution, name)[expected_indices]
        )


@pytest.mark.parametrize("name", ["tau", "orbital_elements", "plot_evol"])
def test_null_solution_has_no_timelike_only_api(radial_solution, name) -> None:
    assert not hasattr(radial_solution, name)
    with pytest.raises(AttributeError):
        getattr(radial_solution, name)


@pytest.mark.parametrize("vector_query", [False, True])
def test_null_solution_at_exact_samples_is_identical(
    radial_solution, vector_query,
) -> None:
    index = np.array([3, 1, 3]) if vector_query else 3
    selected = radial_solution.at(radial_solution.t[index])
    assert selected.t.shape == ((3,) if vector_query else ())
    assert selected.xyz.shape == ((3, 3) if vector_query else (3,))
    for name in (*_COMPONENTS, *_VECTORS):
        np.testing.assert_array_equal(
            getattr(selected, name), getattr(radial_solution, name)[index]
        )


def test_null_solution_at_preserves_single_element_query_shape(radial_solution) -> None:
    query = (radial_solution.t[1] + radial_solution.t[2]) / 2
    scalar = radial_solution.at(query)
    vector = radial_solution.at(query.reshape(1))
    assert scalar.t.shape == ()
    assert scalar.xyz.shape == (3,)
    assert vector.t.shape == (1,)
    assert vector.xyz.shape == (1, 3)
    for name in (*_COMPONENTS, *_VECTORS):
        np.testing.assert_array_equal(getattr(scalar, name), getattr(vector, name)[0])


def test_null_solution_at_radial_schwarzschild_matches_analytic_time(
    radial_null, radial_solution,
) -> None:
    # Outgoing radial null ray, r in r_g and t in T0:
    # t - t0 = r - r0 + 2 log((r - 2)/(r0 - 2)).
    # Dense stored samples and 1e-6 relative tolerance allow C Hermite error.
    query = (radial_solution.t[:-1] + radial_solution.t[1:]) / 2
    interpolated = radial_solution.at(query)
    radius_scale = radial_null.R / 10
    time_scale = (radius_scale / c).to(u.s)
    elapsed = ((query - radial_null.t) / time_scale).to_value(u.one)
    expected = np.array([
        brentq(
            lambda radius: radius - 10 + 2 * np.log((radius - 2) / 8) - time,
            10, 10 + time,
        )
        for time in elapsed
    ])
    radius = (interpolated.R / radius_scale).to_value(u.one)
    np.testing.assert_allclose(radius, expected, rtol=1e-6, atol=0)
    np.testing.assert_array_equal(interpolated.t, query)
    for name, quantity in _quantity_views(interpolated).items():
        assert np.all(np.isfinite(quantity.value)), name


@pytest.mark.parametrize("which", ["before", "after", "mixed"])
def test_null_solution_at_rejects_outside_domain(radial_solution, which) -> None:
    spacing = radial_solution.t[1] - radial_solution.t[0]
    queries = {
        "before": radial_solution.t[0] - spacing,
        "after": radial_solution.t[-1] + spacing,
        "mixed": u.Quantity([radial_solution.t[1], radial_solution.t[-1] + spacing]),
    }
    with pytest.raises(ValueError):
        radial_solution.at(queries[which])


def test_null_solution_at_requires_time_units(radial_solution) -> None:
    with pytest.raises(TypeError):
        radial_solution.at(radial_solution.t[1].to_value(u.s))


def test_t_eval_in_other_unit_than_initial_time_is_accepted() -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0.6)
    time_scale = (metric.r_g / c).to(u.s)
    photon = metric.null(
        R=8 * metric.r_g, Theta=1.0 * u.rad, Phi=0 * u.rad,
        vR=-0.1 * c, vPhi=0.02 * u.rad / time_scale, t=1 * u.s,
    )
    seconds = np.array([1.0, 1.0 + 5 * time_scale.value, 1.0 + 10 * time_scale.value])
    for unit in (u.ms, u.us):
        solution = photon.solve(t_eval=(seconds * u.s).to(unit), method="dp45")
        assert solution.status == 0
        assert len(solution) == 3
    end = (1.0 + 10 * time_scale.value) * u.s
    spanned = photon.solve(
        t_span=(1000 * u.ms, end), t_eval=[1000.0, end.to_value(u.ms)] * u.ms,
        method="dp45",
    )
    assert spanned.status == 0


def test_null_states_survive_copy_and_pickle(radial_solution) -> None:
    import copy
    import pickle

    state = radial_solution[3]
    for clone in (copy.copy(state), copy.deepcopy(state), pickle.loads(pickle.dumps(state))):
        assert clone.t == state.t
        np.testing.assert_array_equal(clone.xyz.value, state.xyz.value)
        np.testing.assert_array_equal(clone.vRQP.vR.value, state.vRQP.vR.value)
    series = pickle.loads(pickle.dumps(radial_solution[1:4]))
    np.testing.assert_array_equal(series.t.value, radial_solution.t[1:4].value)
