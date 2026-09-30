"""Contract checks for the reusable Kerr observable evaluator."""

import numpy as np
import pytest
from astropy import constants as const
from astropy import units as u

from relatipy import Kerr, KerrMcmcModel


MASS = 4e6 * u.M_sun
DISTANCE = 8 * u.kpc
ORBIT = {
    "a": 1000 * u.au,
    "e": 0.5,
    "inc": 0.6 * u.rad,
    "Omega": 0.4 * u.rad,
    "omega": 0.3 * u.rad,
}


def _sample_kerr(spin=0.4, vec=(0.0, 2.0, 0.0)):
    return {"mass": MASS, "spin": spin, "vec": vec}


def test_general_evaluator_sorts_times_normalizes_spin_and_restores_order(monkeypatch):
    from relatipy import _mcmc_core

    captured = {}

    def fake_evaluate(*args):
        captured["args"] = args
        return np.array([[1e-6, 2e-6, 0.01], [3e-6, 4e-6, 0.02]])

    monkeypatch.setattr(_mcmc_core, "evaluate_general", fake_evaluate, raising=False)
    model = KerrMcmcModel()
    t0 = 2002 * u.yr
    dt = (const.G * MASS / const.c**3).to(u.s)
    epochs = u.Quantity([t0 + 2 * dt, t0 - 3 * dt])
    alpha, delta, velocity = model.get_ra_dec_vr(
        kerr=_sample_kerr(), orbit=ORBIT | {"t": t0},
        sol={"t_obs": epochs}, distance=DISTANCE,
    )

    spin, rotation, family, components, times, tau0, kind, scale, *solver = captured["args"]
    assert spin == pytest.approx(0.4)
    np.testing.assert_allclose(rotation[:, 2], [0, 1, 0], atol=1e-15)
    assert family == "elements"
    assert components[0] == 0.0
    np.testing.assert_allclose(times, [-3, 2], atol=1e-6)
    assert tau0 == 0.0
    assert kind == 2
    assert scale == pytest.approx((const.G * MASS / const.c**2 / DISTANCE).to_value(u.one))
    assert solver == [0, 1e-9, 1e-12, 100000]
    np.testing.assert_allclose(alpha, np.array([3e-6, 1e-6]) * u.rad.to(u.arcsec))
    np.testing.assert_allclose(delta, np.array([4e-6, 2e-6]) * u.rad.to(u.arcsec))
    np.testing.assert_allclose(
        velocity, np.array([0.02, 0.01]) * const.c.to_value(u.km / u.s)
    )
    assert not alpha.flags.writeable
    assert not delta.flags.writeable
    assert not velocity.flags.writeable


@pytest.mark.parametrize(
    ("orbit", "family"),
    [
        (ORBIT, "elements"),
        ({"x": 100 * u.au, "y": 1 * u.au, "z": 0 * u.au}, "cartesian"),
        ({"r": 100 * u.au, "theta": 1 * u.rad, "phi": 0.2 * u.rad}, "spherical"),
        ({"R": 100 * u.au, "Theta": 1 * u.rad, "Phi": 0.2 * u.rad}, "bl"),
        ({"p": 100 * u.au, "e": 0.2, "x": 0.4}, "bound"),
    ],
)
def test_all_kerr_orbit_input_families_reach_native(monkeypatch, orbit, family):
    from relatipy import _mcmc_core

    seen = []

    def fake_evaluate(*args):
        seen.append((args[2], args[3].copy()))
        return np.zeros((1, 3))

    monkeypatch.setattr(_mcmc_core, "evaluate_general", fake_evaluate, raising=False)
    KerrMcmcModel().get_ra_dec_vr(
        kerr=_sample_kerr(spin=0, vec=(0, 0, 0)), orbit=orbit,
        sol={"t_eval": [0] * u.s}, distance=DISTANCE,
    )
    assert seen[0][0] == family
    assert seen[0][1].shape == (7,)


def test_proper_time_uses_proper_epoch_and_solver_override(monkeypatch):
    from relatipy import _mcmc_core

    captured = {}

    def fake_evaluate(*args):
        captured["args"] = args
        return np.zeros((2, 3))

    monkeypatch.setattr(_mcmc_core, "evaluate_general", fake_evaluate, raising=False)
    model = KerrMcmcModel()
    model.set_solver(method="radau", rtol=1e-8, atol=1e-10)
    dt = (const.G * MASS / const.c**3).to(u.s)
    model.get_ra_dec_vr(
        kerr=_sample_kerr(), orbit=ORBIT | {"t": 2000 * u.yr, "tau": 7 * u.s},
        sol={"tau_eval": u.Quantity([7 * u.s + dt, 7 * u.s - 2 * dt]), "atol": 1e-11},
        distance=DISTANCE,
    )
    args = captured["args"]
    np.testing.assert_allclose(args[4], [-2, 1], atol=1e-13)
    assert args[6] == 1
    assert args[8:] == (1, 1e-8, 1e-11, 100000)


@pytest.mark.parametrize(
    ("kerr", "orbit", "sol", "distance", "error"),
    [
        (_sample_kerr(vec=(0, 0, 0)), ORBIT, {"t_eval": [0] * u.s}, DISTANCE, ValueError),
        (_sample_kerr(), ORBIT, {"t_eval": [0] * u.s, "t_obs": [0] * u.s}, DISTANCE, ValueError),
        (_sample_kerr(), ORBIT, {"t_eval": [] * u.s}, DISTANCE, ValueError),
        (_sample_kerr(), ORBIT, {"t_eval": [0.0]}, DISTANCE, TypeError),
        (_sample_kerr(), ORBIT, {"t_eval": [0] * u.s}, 8.0, TypeError),
        (_sample_kerr(), ORBIT | {"bad": 1}, {"t_eval": [0] * u.s}, DISTANCE, ValueError),
    ],
)
def test_invalid_general_inputs_fail_before_native(kerr, orbit, sol, distance, error):
    with pytest.raises(error):
        KerrMcmcModel().get_ra_dec_vr(
            kerr=kerr, orbit=orbit, sol=sol, distance=distance,
        )


def test_legacy_call_needs_legacy_epochs():
    model = KerrMcmcModel()
    with pytest.raises(ValueError, match="legacy evaluation requires"):
        model(np.zeros(13))


def test_general_initial_projection_matches_kerr_orbit():
    initial = Kerr(mass=MASS, spin=0).orbit(**ORBIT).initial
    alpha, delta, velocity = KerrMcmcModel().get_ra_dec_vr(
        kerr=_sample_kerr(spin=0, vec=(0, 0, 1)),
        orbit=ORBIT,
        sol={"t_eval": [0] * u.s},
        distance=DISTANCE,
    )
    expected_alpha = (initial.y / DISTANCE).to_value(
        u.arcsec, equivalencies=u.dimensionless_angles()
    )
    expected_delta = (initial.x / DISTANCE).to_value(
        u.arcsec, equivalencies=u.dimensionless_angles()
    )
    expected_velocity = (
        initial.ut.to_value(u.one) + (initial.uz / const.c).to_value(u.one) - 1
    ) * const.c.to_value(u.km / u.s)
    assert alpha[0] == pytest.approx(expected_alpha, abs=1e-10)
    assert delta[0] == pytest.approx(expected_delta, abs=1e-10)
    assert velocity[0] == pytest.approx(expected_velocity, abs=1e-6)


def test_three_time_coordinates_meet_at_the_same_initial_state():
    t0, tau0 = 12345 * u.s, 678 * u.s
    initial = Kerr(mass=MASS, spin=0).orbit(
        **ORBIT, t=t0, tau=tau0
    ).initial
    arrival = t0 + initial.z / const.c
    model = KerrMcmcModel()
    outputs = [
        model.get_ra_dec_vr(
            kerr=_sample_kerr(spin=0, vec=(0, 0, 1)),
            orbit=ORBIT | {"t": t0, "tau": tau0},
            sol={key: u.Quantity([epoch])},
            distance=DISTANCE,
        )
        for key, epoch in (
            ("t_eval", t0), ("tau_eval", tau0), ("t_obs", arrival)
        )
    ]
    for alternative in outputs[1:]:
        for expected, actual in zip(outputs[0], alternative):
            np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-7)


def test_rotated_spin_elements_and_cartesian_have_same_geodesic():
    # A Schwarzschild initial state supplies the Cartesian Kepler state in
    # observer axes; its position and velocity are independent of the spin.
    initial = Kerr(mass=MASS, spin=0).orbit(**ORBIT).initial
    cartesian = {
        key: getattr(initial, key)
        for key in ("x", "y", "z", "vx", "vy", "vz")
    }
    model = KerrMcmcModel()
    sol = {"t_eval": u.Quantity([0, 1e5, 2e5], u.s)}
    kerr = _sample_kerr(spin=0.6, vec=(1, 2, 3))
    by_elements = model.get_ra_dec_vr(
        kerr=kerr, orbit=ORBIT, sol=sol, distance=DISTANCE
    )
    by_cartesian = model.get_ra_dec_vr(
        kerr=kerr, orbit=cartesian, sol=sol, distance=DISTANCE
    )
    for expected, actual in zip(by_elements, by_cartesian):
        np.testing.assert_allclose(actual, expected, rtol=0, atol=2e-7)


def test_s2_like_periapsis_observed_time_matches_initial_projection():
    distance = 8.33 * u.kpc
    mass = 4.35e6 * u.M_sun
    t0 = 2002.33 * u.yr
    a = (0.1255 * u.arcsec).to_value(u.rad) * distance
    orbit = {
        "a": a,
        "e": 0.8839,
        "inc": np.deg2rad(134.18) * u.rad,
        "Omega": np.deg2rad(226.94) * u.rad,
        "omega": np.deg2rad(65.51) * u.rad,
        "f": 0 * u.rad,
        "t": t0,
    }
    initial = Kerr(mass=mass, spin=0).orbit(**orbit).initial
    arrival = t0 + initial.z / const.c
    epochs = u.Quantity([arrival - 0.05 * u.yr, arrival, arrival + 0.05 * u.yr])
    alpha, delta, velocity = KerrMcmcModel().get_ra_dec_vr(
        kerr={"mass": mass, "spin": 0, "vec": (0, 0, 0)},
        orbit=orbit,
        sol={"t_obs": epochs},
        distance=distance,
    )
    assert np.all(np.isfinite(alpha))
    assert np.all(np.isfinite(delta))
    assert np.all(np.isfinite(velocity))
    expected_alpha = (initial.y / distance).to_value(
        u.arcsec, equivalencies=u.dimensionless_angles()
    )
    expected_delta = (initial.x / distance).to_value(
        u.arcsec, equivalencies=u.dimensionless_angles()
    )
    expected_velocity = (
        initial.ut.to_value(u.one) + (initial.uz / const.c).to_value(u.one) - 1
    ) * const.c.to_value(u.km / u.s)
    assert alpha[1] == pytest.approx(expected_alpha, abs=1e-8)
    assert delta[1] == pytest.approx(expected_delta, abs=1e-8)
    assert velocity[1] == pytest.approx(expected_velocity, abs=1e-4)
