"""Conversion regressions: saved canonical velocity must not be normalized again."""

import numpy as np
import pytest
from astropy import units as u
from astropy.constants import c

from relatipy import _core
from relatipy.metrics.kerr import Kerr
from relatipy.geodesic.orbit import _state_from_native


def test_stored_canonical_output_preserves_drift_and_azimuth_branch():
    metric = Kerr(mass=1 * u.Msun, spin=0.5)
    source = np.array([[0.0, 8.0, 1.1, 20.0, 1.3, -0.01, 0.002, 0.02]])
    rows, cartesian, statuses = _core.reconstruct_canonical_batch(metric.spin, source)
    np.testing.assert_array_equal(statuses, [0])
    np.testing.assert_array_equal(rows[:, :4], source[:, :4])
    np.testing.assert_array_equal(rows[:, 7:11], source[:, 4:])
    np.testing.assert_array_equal(rows[:, 4:7], source[:, 5:] / source[:, 4:5])
    state = _state_from_native(
        metric, source, np.array([0.0]), length_unit=u.m, time_unit=u.s,
    )
    np.testing.assert_array_equal(state.ut.to_value(u.one), source[:, 4])
    np.testing.assert_allclose((state.uR / c).to_value(u.one), source[:, 5], rtol=1e-15)
    np.testing.assert_allclose((state.uxyz / c).to_value(u.one),
                               source[:, 4:5] * cartesian[:, 4:], rtol=1e-15)
    np.testing.assert_array_equal(state.Phi.to_value(u.rad), source[:, 3])


@pytest.mark.parametrize("method", ["radau", "dop853", "dp45"])
def test_public_solve_keeps_integrated_velocity_on_native_steps(monkeypatch, method):
    metric = Kerr(mass=1 * u.Msun, spin=0.5)
    orbit = metric.orbit(
        R=8 * metric.r_g, Theta=1.1 * u.rad, Phi=0.2 * u.rad,
        vR=-0.01 * c, vTheta=0.002 * u.rad / metric._time_scale,
        vPhi=0.02 * u.rad / metric._time_scale,
    )
    captured = []
    original = _core.integrate_kerr

    def capture(*args, **kwargs):
        result = original(*args, **kwargs)
        captured.append(result[1].copy())
        return result

    monkeypatch.setattr(_core, "integrate_kerr", capture)
    result = orbit.solve(tau_span=(0 * u.s, 10 * metric._time_scale),
                         method=method, rtol=1e-6, atol=1e-8)
    assert result.success and len(captured) == 1
    np.testing.assert_array_equal(result.ut.to_value(u.one), captured[0][:, 4])
    np.testing.assert_allclose((result.uR / c).to_value(u.one), captured[0][:, 5], rtol=1e-15)
    np.testing.assert_array_equal(result.Phi.to_value(u.rad), captured[0][:, 3])
