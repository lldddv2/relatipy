"""Integration options share validation and failed solves never invent samples."""

import numpy as np
import pytest
from astropy import units as u

from relatipy import IntegrationError, IntegrationInfo, Kerr
from relatipy.geodesic._native import describe_integrator_status


def _plunging_orbit():
    metric = Kerr(mass=1e6 * u.M_sun, spin=0.0)
    orbit = metric.orbit(
        R=6 * metric.r_g, Theta=np.pi / 2 * u.rad, Phi=0 * u.rad,
        vR=-3e4 * u.km / u.s,
    )
    return metric, orbit


def test_projection_failure_has_distinct_native_diagnostic() -> None:
    assert describe_integrator_status(10) == "invariant projection failed; native status 10"


@pytest.mark.parametrize(
    ("options", "error", "match"),
    [
        ({"rtol": "1e-3"}, TypeError, "rtol must be a real number"),
        ({"rtol": True}, TypeError, "rtol must be a real number"),
        ({"rtol": -1.0}, ValueError, "rtol must be finite and non-negative"),
        ({"first_step": -1 * u.s}, ValueError, "first_step must be finite and positive"),
        ({"max_step": np.inf * u.s}, ValueError, "max_step must be finite and positive"),
    ],
)
def test_solve_and_integration_info_share_option_rules(options, error, match) -> None:
    _, orbit = _plunging_orbit()
    settings = {"rtol": 1e-3, "first_step": None, "max_step": None} | options
    with pytest.raises(error, match=match):
        orbit.solve(tau_span=(0 * u.s, 1 * u.s), **options)
    with pytest.raises(error, match=match):
        IntegrationInfo(
            "radau", settings["rtol"], 1e-6, settings["max_step"],
            settings["first_step"], 0, 0,
        )


def test_failure_before_every_tau_eval_sample_raises() -> None:
    metric, orbit = _plunging_orbit()
    full = orbit.solve(tau_span=(0 * u.s, 1e5 * u.s))
    assert full.status == -1
    assert "right-hand side evaluation failed; native status 4" in full.message
    late = full.tau[-1] + np.array([1.0, 2.0]) * metric._time_scale
    with pytest.raises(IntegrationError, match="no tau_eval sample was reached"):
        orbit.solve(tau_eval=late)


def test_failure_keeps_only_reached_tau_eval_samples() -> None:
    _, orbit = _plunging_orbit()
    sampled = orbit.solve(tau_eval=[0, 1e5] * u.s)
    assert sampled.status == -1
    assert sampled.tau.to_value(u.s).tolist() == [0.0]


def test_termination_before_every_tau_eval_sample_raises(monkeypatch) -> None:
    from relatipy import IntegrationTerminated
    from relatipy.geodesic.orbit import Orbit

    _, orbit = _plunging_orbit()

    def stop_at_start(self, y0, tau0, tau1, options, *, store_steps):
        stats = {"n_steps": 0, "nfev": 0, "crossed_outer_horizon": True}
        return np.array([tau0]), np.array([y0]), stats, 1

    monkeypatch.setattr(Orbit, "_call_native", stop_at_start)
    with pytest.raises(IntegrationTerminated) as caught:
        orbit.solve(tau_eval=[1, 2] * u.s)
    assert caught.value.reason == "outer_horizon"
    assert caught.value.tau == orbit.initial.tau
    assert "before the first tau_eval sample" in str(caught.value)
