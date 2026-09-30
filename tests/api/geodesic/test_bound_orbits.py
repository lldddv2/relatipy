"""Public stable bound-orbit inputs against frozen independent KerrGeoPy states."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import FrozenInstanceError
import json

import numpy as np
import pytest
from astropy import units as u
from astropy.constants import c

from relatipy.coordinates import KerrOrbitalElements
from relatipy.metrics.kerr import Kerr

from support import FIXTURES_DIR


_REFERENCE = json.loads(
    (FIXTURES_DIR / "bound_orbit_reference.json").read_text(encoding="utf-8")
)


def _canonical(state, metric: Kerr) -> np.ndarray:
    """Convert public Boyer--Lindquist fields to G=c=M=1 state order."""
    return np.stack((
        (state.t / metric._time_scale).to_value(u.one),
        (state.R / metric.r_g).to_value(u.one),
        state.Theta.to_value(u.rad),
        state.Phi.to_value(u.rad),
        state.ut.to_value(u.one),
        (state.uR / c).to_value(u.one),
        (state.uTheta * metric._time_scale).to_value(u.rad),
        (state.uPhi * metric._time_scale).to_value(u.rad),
    ), axis=-1)


def _invariants(spin: float, state: np.ndarray) -> dict[str, np.ndarray]:
    """Use Schmidt (2002), Eqs. (3)--(8), independent of production C."""
    radius, theta = state[..., 1], state[..., 2]
    sine2 = np.sin(theta)**2
    sigma = radius**2 + spin**2 * np.cos(theta)**2
    delta = radius**2 - 2 * radius + spin**2
    inverse = np.zeros((*state.shape[:-1], 4, 4))
    inverse[..., 0, 0] = -((radius**2 + spin**2)**2 - delta * spin**2 * sine2) / (delta * sigma)
    inverse[..., 0, 3] = inverse[..., 3, 0] = -2 * spin * radius / (delta * sigma)
    inverse[..., 1, 1] = delta / sigma
    inverse[..., 2, 2] = 1 / sigma
    inverse[..., 3, 3] = (delta - spin**2 * sine2) / (delta * sigma * sine2)
    momentum = np.einsum("...ij,...j->...i", np.linalg.inv(inverse), state[..., 4:])
    energy, angular_momentum = -momentum[..., 0], momentum[..., 3]
    return {
        "norm": np.einsum("...i,...i->...", momentum, state[..., 4:]),
        "energy": energy,
        "Lz": angular_momentum,
        "CarterQ": momentum[..., 2]**2 + np.cos(theta)**2 * (
            spin**2 * (1 - energy**2) + angular_momentum**2 / sine2
        ),
    }


def _raw_orbit(case: dict, *, mass: u.Quantity = 1 * u.Msun):
    parameters = case["parameters"]
    metric = Kerr(mass=mass, spin=parameters["spin"])
    orbit = metric.orbit(
        p=parameters["p_over_M"] * metric.r_g,
        e=parameters["e"], x=parameters["x"],
        q_r0=parameters["q_r0_rad"] * u.rad,
        q_theta0=parameters["q_theta0_rad"] * u.rad,
        q_phi0=parameters["q_phi0_rad"] * u.rad,
    )
    return metric, orbit


@pytest.mark.parametrize("case", _REFERENCE["cases"], ids=lambda case: case["id"])
def test_bound_initial_state_matches_independent_kerrgeopy(case: dict) -> None:
    metric, orbit = _raw_orbit(case)
    expected = np.asarray(case["initial_x_u"])
    actual = _canonical(orbit.initial, metric)
    if case["id"] == "circular_inclined":
        # KerrGeoPy leaves a ~1e-8 radial derivative at an exact circle.
        expected = expected.copy()
        expected[5] = 0.0
        assert actual[5] == 0.0
    np.testing.assert_allclose(actual, expected, rtol=2e-11, atol=2e-12)
    for name, value in case["constants"].items():
        np.testing.assert_allclose(_invariants(metric.spin, actual)[name], value,
                                   rtol=3e-10, atol=3e-10)
    np.testing.assert_allclose(_invariants(metric.spin, actual)["norm"], -1,
                               rtol=0, atol=3e-11)


def test_bound_object_units_and_default_phases() -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0.7)
    parameters = _REFERENCE["cases"][0]["parameters"]
    elements = KerrOrbitalElements(
        (parameters["p_over_M"] * metric.r_g).to(u.km),
        parameters["e"] * u.one, parameters["x"] * u.one,
        parameters["q_r0_rad"] * u.rad.to(u.deg) * u.deg,
        parameters["q_theta0_rad"] * u.rad,
        parameters["q_phi0_rad"] * u.rad,
    )
    from_object = metric.orbit(elements=elements, t=3 * u.s, tau=5 * u.s)
    from_raw = metric.orbit(
        p=elements.p, e=elements.e, x=elements.x,
        q_r0=elements.q_r0, q_theta0=elements.q_theta0,
        q_phi0=elements.q_phi0, t=3 * u.s, tau=5 * u.s,
    )
    np.testing.assert_allclose(_canonical(from_object.initial, metric),
                               _canonical(from_raw.initial, metric), rtol=0, atol=3e-13)
    assert from_object.t == 3 * u.s
    assert from_object.tau == 5 * u.s
    assert from_object.initial.tau == 5 * u.s
    implicit = metric.orbit(p=12 * metric.r_g, e=0.3, x=0.7)
    explicit = metric.orbit(p=12 * metric.r_g, e=0.3, x=0.7,
                            q_r0=0 * u.rad, q_theta0=0 * u.rad,
                            q_phi0=0 * u.rad)
    np.testing.assert_allclose(_canonical(implicit.initial, metric),
                               _canonical(explicit.initial, metric), rtol=0, atol=3e-13)
    with pytest.raises((FrozenInstanceError, AttributeError)):
        elements.e = 0.8
    with pytest.raises((TypeError, ValueError)):
        elements.p.value[...] = 1.0


def test_phase_periods_preserve_state_and_unwrapped_azimuth() -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0.7)
    kwargs = dict(p=12 * metric.r_g, e=0.3, x=0.7,
                  q_r0=0.4 * u.rad, q_theta0=1.2 * u.rad,
                  q_phi0=0.2 * u.rad)
    reference = _canonical(metric.orbit(**kwargs).initial, metric)
    wrapped = _canonical(metric.orbit(
        **{**kwargs, "q_r0": (0.4 + 2*np.pi) * u.rad,
           "q_theta0": (1.2 + 2*np.pi) * u.rad,
           "q_phi0": (0.2 + 2*np.pi) * u.rad}
    ).initial, metric)
    np.testing.assert_allclose(wrapped[[0, 1, 2, 4, 5, 6, 7]],
                               reference[[0, 1, 2, 4, 5, 6, 7]],
                               rtol=2e-11, atol=2e-11)
    np.testing.assert_allclose(wrapped[3] - reference[3], 2*np.pi,
                               rtol=0, atol=2e-11)


def test_cartesian_x_remains_cartesian_and_families_do_not_mix() -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0.5)
    old = metric.orbit(x=12 * metric.r_g, vy=0.1 * c)
    assert old.xyz[0] == 12 * metric.r_g
    with pytest.raises(ValueError, match="families"):
        metric.orbit(p=12 * metric.r_g, e=0.3, x=0.7, y=1 * metric.r_g)
    with pytest.raises(ValueError, match="families"):
        metric.orbit(p=12 * metric.r_g, e=0.3, x=0.7, a=12 * metric.r_g)
    elements = KerrOrbitalElements(12 * metric.r_g, 0.3, 0.7)
    with pytest.raises(ValueError, match="families"):
        metric.orbit(elements=elements, p=12 * metric.r_g)


def test_bound_input_rejects_missing_and_invalid_values() -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0.5)
    for kwargs in ({"p": 12 * metric.r_g, "e": 0.3},
                   {"p": 12 * metric.r_g, "x": 0.7},
                   {"q_r0": 0.1 * u.rad},
                   {"p": 12 * metric.r_g, "e": 0.3, "x": 0.7,
                    "q_theta0": 1 * u.s}):
        with pytest.raises((ValueError, TypeError, u.UnitConversionError)):
            metric.orbit(**kwargs)
    for p, e, x in ((0 * metric.r_g, 0.3, 0.7),
                    (12 * metric.r_g, -0.1, 0.7),
                    (12 * metric.r_g, 1.0, 0.7),
                    (12 * metric.r_g, 0.3, 1.1),
                    (3 * metric.r_g, 0.2, 0.7)):
        with pytest.raises(ValueError):
            metric.orbit(p=p, e=e, x=x)
    with pytest.raises(u.UnitConversionError):
        metric.orbit(p=12 * u.s, e=0.3, x=0.7)
    with pytest.raises((TypeError, u.UnitConversionError)):
        metric.orbit(p=12 * metric.r_g, e=0.3 * u.s, x=0.7)


def test_bound_container_rejects_invalid_shapes_values_and_types() -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0.5)
    for p in (np.array([12.0, 13.0]) * metric.r_g,
              np.nan * metric.r_g):
        with pytest.raises(ValueError):
            KerrOrbitalElements(p, 0.3, 0.7)
    with pytest.raises(u.UnitConversionError):
        KerrOrbitalElements(12 * metric.r_g, 0.3, 0.7, q_r0=1 * u.s)
    with pytest.raises(ValueError):
        KerrOrbitalElements(12 * metric.r_g, 0.3, 0.7,
                            q_phi0=np.inf * u.rad)
    with pytest.raises(TypeError):
        metric.orbit(elements=object())


@pytest.mark.parametrize(
    ("kwargs", "error", "message"),
    [
        ({"p": (12 + 1j) * u.km}, TypeError, "p must be real"),
        ({"q_theta0": 1j * u.rad}, TypeError, "q_theta0 must be real"),
        ({"p": 0 * u.km}, ValueError, "p must be finite and positive"),
        ({"p": np.inf * u.km}, ValueError, "p must be finite and positive"),
        ({"e": -0.1}, ValueError, r"e must be finite and in \[0, 1\]"),
        ({"e": np.nan}, ValueError, r"e must be finite and in \[0, 1\]"),
        ({"e": 1.0}, ValueError, "bound eccentricity e must satisfy"),
        ({"x": 1.5}, ValueError, r"x must be finite and in \[-1, 1\]"),
        ({"q_r0": np.nan * u.rad}, ValueError, "q_r0 must be finite"),
    ],
)
def test_bound_container_reports_each_invalid_field(kwargs, error, message) -> None:
    fields = {"p": 12 * u.km, "e": 0.3, "x": 0.7} | kwargs
    with pytest.raises(error, match=message):
        KerrOrbitalElements(**fields)


def test_bound_container_stores_dimensionless_fields_as_float() -> None:
    elements = KerrOrbitalElements(12 * u.km, 0.3 * u.one, -1 * u.one)
    assert type(elements.e) is float and elements.e == 0.3
    assert type(elements.x) is float and elements.x == -1.0


def test_polar_bound_orbit_needs_off_axis_initial_phase() -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0.5)
    orbit = metric.orbit(p=12 * metric.r_g, e=0.2, x=0,
                         q_theta0=0.5 * u.rad)
    assert 0 * u.rad < orbit.Theta < np.pi * u.rad
    np.testing.assert_allclose(_invariants(metric.spin, _canonical(orbit, metric))["norm"],
                               -1, rtol=0, atol=3e-11)
    with pytest.raises(ValueError, match="Boyer-Lindquist"):
        metric.orbit(p=12 * metric.r_g, e=0.2, x=0,
                     q_theta0=0 * u.rad)


def test_bound_solve_sampling_and_mutable_orbit_lifecycle() -> None:
    metric, orbit = _raw_orbit(_REFERENCE["cases"][0])
    initial = _canonical(orbit.initial, metric)
    copy = orbit.copy()
    evaluation = np.array([0.0, 0.1, 0.2]) * metric._time_scale
    solution = orbit.solve(tau_eval=evaluation, method="dp45")
    assert solution.success and len(solution) == len(evaluation)
    np.testing.assert_array_equal(solution.tau.to_value(u.s), evaluation.to_value(u.s))
    np.testing.assert_array_equal(_canonical(orbit, metric), initial)
    orbit.integrate(evaluation[-1], method="dp45")
    assert orbit.tau == evaluation[-1]
    np.testing.assert_array_equal(_canonical(copy, metric), initial)
    orbit.reset()
    np.testing.assert_array_equal(_canonical(orbit, metric), initial)


@pytest.mark.parametrize("method", ["radau", "dop853", "dp45"])
def test_bound_initial_state_integrates_and_preserves_invariants(method: str) -> None:
    metric, orbit = _raw_orbit(_REFERENCE["cases"][0])
    end = 2 * metric._time_scale
    # tau_span retains accepted native steps; tau_eval interpolates state
    # components and is therefore unsuitable for a conservation check.
    solution = orbit.solve(tau_span=(0 * u.s, end), method=method,
                           rtol=1e-10, atol=1e-12)
    assert solution.success
    assert solution.tau[-1] == end
    invariants = _invariants(metric.spin, _canonical(solution, metric))
    np.testing.assert_allclose(invariants["norm"], -1, rtol=0, atol=2e-9)
    for name in ("energy", "Lz", "CarterQ"):
        np.testing.assert_allclose(invariants[name], invariants[name][0],
                                   rtol=2e-9, atol=2e-9)


def test_concurrent_bound_initialization_is_deterministic() -> None:
    cases = _REFERENCE["cases"][:2]

    def initialize(case: dict) -> np.ndarray:
        metric, orbit = _raw_orbit(case)
        return _canonical(orbit.initial, metric)

    serial = [initialize(case) for case in cases]
    with ThreadPoolExecutor(max_workers=2) as pool:
        concurrent = list(pool.map(initialize, cases))
    for actual, expected in zip(concurrent, serial, strict=True):
        np.testing.assert_array_equal(actual, expected)


def test_extremal_spin_bound_initial_state_is_timelike() -> None:
    metric = Kerr(mass=1 * u.Msun, spin=1.0)
    orbit = metric.orbit(p=12 * metric.r_g, e=0.2, x=0.8,
                         q_r0=0.4 * u.rad, q_theta0=0.6 * u.rad)
    actual = _canonical(orbit.initial, metric)
    assert np.all(np.isfinite(actual))
    np.testing.assert_allclose(_invariants(metric.spin, actual)["norm"], -1,
                               rtol=0, atol=3e-11)
