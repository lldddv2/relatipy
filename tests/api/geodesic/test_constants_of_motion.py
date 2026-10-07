"""Public E, Lz and Carter Q accessors of Orbit and Solution."""

from __future__ import annotations

import json

import numpy as np
import pytest
from astropy import units as u
from astropy.constants import c

from relatipy.geodesic import IntegrationInfo, Solution, State
from relatipy.metrics.kerr import Kerr

from support import FIXTURES_DIR


_REFERENCE = json.loads(
    (FIXTURES_DIR / "bound_orbit_reference.json").read_text(encoding="utf-8")
)
_NAMES = (("get_E", "energy"), ("get_Lz", "Lz"), ("get_Q", "CarterQ"))


def _bound_orbit(case: dict):
    parameters = case["parameters"]
    metric = Kerr(mass=4e6 * u.Msun, spin=parameters["spin"])
    return metric.orbit(
        p=parameters["p_over_M"] * metric.r_g,
        e=parameters["e"], x=parameters["x"],
        q_r0=parameters["q_r0_rad"] * u.rad,
        q_theta0=parameters["q_theta0_rad"] * u.rad,
        q_phi0=parameters["q_phi0_rad"] * u.rad,
    )


@pytest.mark.parametrize("case", _REFERENCE["cases"], ids=lambda case: case["id"])
def test_orbit_constants_match_independent_kerrgeopy(case: dict) -> None:
    orbit = _bound_orbit(case)
    for method, name in _NAMES:
        value = getattr(orbit, method)()
        assert type(value) is float
        np.testing.assert_allclose(value, case["constants"][name],
                                   rtol=3e-10, atol=3e-10)


def test_solution_constants_are_conserved_read_only_series() -> None:
    orbit = _bound_orbit(_REFERENCE["cases"][0])
    end = (5000 * orbit._metric._time_scale).to(u.s)
    solution = orbit.solve(tau_span=(0 * u.s, end))
    for method, _ in _NAMES:
        series = getattr(solution, method)()
        assert series.shape == (len(solution),)
        assert series.dtype == np.float64
        assert not series.flags.writeable
        np.testing.assert_allclose(series[0], getattr(orbit, method)(), rtol=0, atol=0)
        # Automatic tolerances keep the drift near 1e-8 over a few orbits.
        np.testing.assert_allclose(series, series[0], rtol=1e-7, atol=1e-8)
    sampled = orbit.solve(tau_eval=np.linspace(0, end.value, 500) * end.unit)
    assert sampled.get_E().shape == (500,)


def test_orbit_constants_follow_current_state() -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0.5)
    orbit = metric.orbit(x=20 * metric.r_g, vy=0.2 * c, vz=0.05 * c)
    initial = (orbit.get_E(), orbit.get_Lz(), orbit.get_Q())
    assert initial[0] < 1 and initial[2] > 0
    orbit.integrate(500 * metric._time_scale)
    current = (orbit.get_E(), orbit.get_Lz(), orbit.get_Q())
    assert current != initial
    np.testing.assert_allclose(current, initial, rtol=1e-8)
    orbit.reset()
    assert (orbit.get_E(), orbit.get_Lz(), orbit.get_Q()) == initial


def test_equatorial_orbit_has_zero_carter_constant() -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0.9)
    orbit = metric.orbit(x=15 * metric.r_g, vy=0.25 * c)
    assert abs(orbit.get_Q()) <= 1e-24


def test_solution_without_native_samples_rejects_constants() -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0.5)
    native = metric.orbit(x=20 * metric.r_g, vy=0.2 * c).solve(
        tau_span=(0 * u.s, 1e-4 * u.s)
    )
    plain = State(
        tau=native.tau, txyz=native.txyz, trqp=native.trqp, tRQP=native.tRQP,
        vxyz=native.vxyz, vrqp=native.vrqp, vRQP=native.vRQP, ut=native.ut,
        uxyz=native.uxyz, urqp=native.urqp, uRQP=native.uRQP,
        orbital_elements=native.orbital_elements(),
    )
    solution = Solution(state=plain, integration=native.integration,
                        status=native.status, message=native.message)
    assert isinstance(solution.integration, IntegrationInfo)
    with pytest.raises(ValueError, match="native Kerr samples"):
        solution.get_Q()
