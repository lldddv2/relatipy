"""Public Radau trajectory regressions after enabling the native Jacobian."""

from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest
from astropy import units as u
from astropy.constants import c

from relatipy.metrics.kerr import Kerr


def _solve_regular(spin: float, method: str = "radau"):
    metric = Kerr(mass=1 * u.Msun, spin=spin)
    orbit = metric.orbit(
        R=8 * metric.r_g, Theta=1.1 * u.rad, Phi=0.3 * u.rad,
        vR=-0.01 * c, vTheta=0.002 * u.rad / metric._time_scale,
        vPhi=0.02 * u.rad / metric._time_scale,
    )
    return orbit.solve(
        tau_eval=np.linspace(0, 5, 11) * metric._time_scale,
        method=method, rtol=1e-10, atol=1e-13,
        first_step=0.01 * metric._time_scale,
        max_step=0.1 * metric._time_scale,
    )


@pytest.mark.parametrize("spin", [0.0, 0.5, 1.0])
def test_radau_matches_tight_explicit_trajectory(spin):
    radau = _solve_regular(spin)
    reference = _solve_regular(spin, "dop853")
    assert radau.success and reference.success
    for name in ("t", "xyz", "vxyz", "ut"):
        actual = getattr(radau, name)
        expected = getattr(reference, name).to_value(actual.unit)
        np.testing.assert_allclose(actual.value, expected, rtol=2e-9, atol=1e-9)
    assert radau.integration.nfev > 0


@pytest.mark.parametrize("pole", ["north", "south"])
@pytest.mark.parametrize("distance", [1e-4, 1e-8, 1e-10, 128 * np.finfo(float).eps])
def test_radau_preserves_valid_near_polar_chart(pole, distance):
    metric = Kerr(mass=1 * u.Msun, spin=0)
    theta = distance if pole == "north" else np.pi - distance
    orbit = metric.orbit(R=8 * metric.r_g, Theta=theta * u.rad, Phi=0 * u.rad)
    solution = orbit.solve(
        tau_eval=np.array([0, 1e-4, 1e-3]) * metric._time_scale,
        method="radau", rtol=1e-9, atol=1e-12,
        first_step=1e-4 * metric._time_scale,
        max_step=1e-4 * metric._time_scale,
    )
    assert solution.success
    np.testing.assert_array_equal(solution.Theta.to_value(u.rad), [theta] * 3)
    assert np.all(np.isfinite(solution.xyz.value))
    assert np.all(np.isfinite(solution.vxyz.value))
    assert orbit.tau == orbit.initial.tau


def test_concurrent_radau_calls_keep_context_and_counters_independent():
    spins = [0.0, 1.0]
    serial = [_solve_regular(spin) for spin in spins]
    with ThreadPoolExecutor(max_workers=2) as pool:
        concurrent = list(pool.map(_solve_regular, spins))
    for actual, expected in zip(concurrent, serial, strict=True):
        np.testing.assert_array_equal(actual.xyz.value, expected.xyz.value)
        assert actual.integration.nfev == expected.integration.nfev
        assert actual.integration.n_steps == expected.integration.n_steps
