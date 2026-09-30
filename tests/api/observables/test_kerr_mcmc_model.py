"""End-to-end checks for direct Kerr MCMC observables."""

import numpy as np
import pytest
from astropy import constants as const
from astropy import units as u

from relatipy import Kerr, KerrMcmcModel


def sample_parameters():
    return np.array(
        [8.33, 4.35, 2002.33, 0.1255, 0.8839,
         np.deg2rad(134.18), np.deg2rad(226.94), np.deg2rad(65.51),
         0.001, -0.002, 0.0001, -0.0002, 12.0],
        dtype=float,
    )


def test_initial_observables_match_kerr_orbit_elements():
    params = sample_parameters()
    params[8:] = 0.0
    distance_kpc, mass_million, periapsis_year, semi_major_arcsec = params[:4]
    distance = distance_kpc * u.kpc
    mass = mass_million * 1e6 * u.M_sun
    semi_major = (semi_major_arcsec * u.arcsec).to_value(u.rad) * distance
    initial = Kerr(mass=mass, spin=0.0).orbit(
        a=semi_major, e=params[4], inc=params[5] * u.rad,
        Omega=params[6] * u.rad, omega=params[7] * u.rad,
        f=0.0 * u.rad,
    ).initial
    observed_periapsis = periapsis_year + (
        initial.z / const.c
    ).to_value(u.yr)
    model = KerrMcmcModel(
        [observed_periapsis], [observed_periapsis],
        reference_epoch=periapsis_year,
    )
    model.set_solver(method="dop853", rtol=1e-9, atol=1e-12)

    alpha, delta, velocity = model(params)
    expected_alpha = (initial.y / distance).to_value(u.arcsec, equivalencies=u.dimensionless_angles())
    expected_delta = (initial.x / distance).to_value(u.arcsec, equivalencies=u.dimensionless_angles())
    expected_velocity = (
        initial.ut.to_value(u.one) + (initial.uz / const.c).to_value(u.one) - 1.0
    ) * const.c.to_value(u.km / u.s)
    assert alpha[0] == pytest.approx(expected_alpha, abs=1e-10)
    assert delta[0] == pytest.approx(expected_delta, abs=1e-10)
    assert velocity[0] == pytest.approx(expected_velocity, abs=1e-6)


def test_periapsis_projection_and_observation_order():
    p = sample_parameters()
    distance, mass, t_peri, semi_major, eccentricity = p[:5]
    inclination, ascending_node, periapsis_argument = p[5:8]
    gravitational_radius = (
        const.G * mass * 1e6 * const.M_sun / const.c**2
    ).to_value(u.au)
    distance_au = (distance * u.kpc).to_value(u.au)
    radius = semi_major * (1 * u.arcsec).to_value(u.rad)
    radius *= distance_au / gravitational_radius * (1 - eccentricity)
    x = radius * (
        np.cos(ascending_node) * np.cos(periapsis_argument)
        - np.sin(ascending_node) * np.cos(inclination)
        * np.sin(periapsis_argument)
    )
    y = radius * (
        np.sin(ascending_node) * np.cos(periapsis_argument)
        + np.cos(ascending_node) * np.cos(inclination)
        * np.sin(periapsis_argument)
    )
    z = radius * np.sin(inclination) * np.sin(periapsis_argument)
    observed_periapsis = t_peri + z * gravitational_radius / const.c.to_value(u.au / u.yr)
    epochs = np.array([observed_periapsis + 0.1, observed_periapsis,
                       observed_periapsis - 0.1])
    model = KerrMcmcModel(epochs, epochs[[2, 1]], reference_epoch=2009.2)
    model.set_solver(method="dop853", rtol=1e-9, atol=1e-12)
    alpha, delta, velocity = model(p)
    scale = gravitational_radius / distance_au * (1 * u.rad).to_value(u.arcsec)
    expected_alpha = y * scale + p[8] + p[10] * (observed_periapsis - 2009.2)
    expected_delta = x * scale + p[9] + p[11] * (observed_periapsis - 2009.2)

    assert alpha[1] == pytest.approx(expected_alpha, abs=2e-9)
    assert delta[1] == pytest.approx(expected_delta, abs=2e-9)
    assert alpha.shape == delta.shape == (3,)
    assert velocity.shape == (2,)
    assert not alpha.flags.writeable
    assert not delta.flags.writeable
    assert not velocity.flags.writeable

    sorted_model = KerrMcmcModel(
        np.sort(epochs), np.sort(epochs[[2, 1]]), reference_epoch=2009.2
    )
    sorted_model.set_solver(method="dop853", rtol=1e-9, atol=1e-12)
    sorted_alpha, sorted_delta, sorted_velocity = sorted_model(p)
    np.testing.assert_allclose(alpha[np.argsort(epochs)], sorted_alpha, atol=2e-9)
    np.testing.assert_allclose(delta[np.argsort(epochs)], sorted_delta, atol=2e-9)
    np.testing.assert_allclose(velocity[np.argsort(epochs[[2, 1]])], sorted_velocity, atol=1e-5)


def test_spin_changes_prediction_and_radau_runs():
    p = sample_parameters()
    epochs = np.linspace(1998.0, 2014.0, 15)
    zero_spin = KerrMcmcModel(epochs, epochs, reference_epoch=2000.0)
    zero_spin.set_solver(method="dop853", rtol=1e-8, atol=1e-10)
    rotating = KerrMcmcModel(
        epochs, epochs, reference_epoch=2000.0,
        spin_vector=(0.2, 0.3, 0.4),
    )
    rotating.set_solver(method="radau", rtol=1e-8, atol=1e-10)
    zero = zero_spin(p)
    nonzero = rotating(p)
    assert all(np.all(np.isfinite(series)) for series in nonzero)
    assert np.max(np.abs(zero[0] - nonzero[0])) > 1e-9


def test_spin_can_change_per_proposal_without_changing_model_default():
    params = sample_parameters()
    epochs = np.linspace(1998.0, 2014.0, 15)
    default_spin = (0.0, 0.0, 0.5)
    model = KerrMcmcModel(
        epochs, epochs[::3], reference_epoch=2000.0,
        spin_vector=default_spin,
    )
    model.set_solver(method="dop853", rtol=1e-9, atol=1e-12)

    baseline = model(params)
    same = model(params, spin_vector=default_spin)
    changed = model(params, spin_vector=(0.0, 0.5, 0.0))
    restored = model(params)

    for expected, actual in zip(baseline, same):
        np.testing.assert_array_equal(actual, expected)
    for expected, actual in zip(baseline, restored):
        np.testing.assert_array_equal(actual, expected)
    assert np.max(np.abs(changed[0] - baseline[0])) > 1e-10

    zero = KerrMcmcModel(epochs, epochs[::3], reference_epoch=2000.0)
    zero.set_solver(method="dop853", rtol=1e-9, atol=1e-12)
    for expected, actual in zip(zero(params), model(params, spin_vector=(0, 0, 0))):
        np.testing.assert_array_equal(actual, expected)

    with pytest.raises(ValueError, match="magnitude"):
        model(params, spin_vector=(0.0, 0.0, 1.1))
    with pytest.raises(ValueError, match="three finite"):
        model(params, spin_vector=(0.0, float("nan"), 0.0))


def test_dop853_observables_converge_for_rotating_case():
    p = sample_parameters()
    epochs = np.linspace(1993.0, 2018.0, 35)
    medium = KerrMcmcModel(
        epochs, epochs[::3], reference_epoch=2000.0,
        spin_vector=(0.2, 0.3, 0.4),
    )
    fine = KerrMcmcModel(
        epochs, epochs[::3], reference_epoch=2000.0,
        spin_vector=(0.2, 0.3, 0.4),
    )
    medium.set_solver(method="dop853", rtol=1e-9, atol=1e-12)
    fine.set_solver(method="dop853", rtol=1e-11, atol=1e-14)
    medium_output = medium(p)
    fine_output = fine(p)
    np.testing.assert_allclose(medium_output[0], fine_output[0], atol=1e-7)
    np.testing.assert_allclose(medium_output[1], fine_output[1], atol=1e-7)
    np.testing.assert_allclose(medium_output[2], fine_output[2], atol=1e-3)


def test_invalid_input_is_rejected():
    p = sample_parameters()
    model = KerrMcmcModel([2002.0], [2002.0], reference_epoch=2000.0)
    with pytest.raises(ValueError, match="set_solver"):
        model(p)
    with pytest.raises(ValueError, match="method"):
        model.set_solver(method="RK45", rtol=1e-6, atol=1e-9)
    for method in ("DOP853", "Radau", "projection_radau", "dp45"):
        with pytest.raises(ValueError, match="method"):
            model.set_solver(method=method, rtol=1e-6, atol=1e-9)
    model.set_solver(method="dop853", rtol=1e-6, atol=1e-9)
    with pytest.raises(ValueError, match="thirteen"):
        model(p[:-1])
    p[4] = 1.0
    with pytest.raises(ValueError, match="e must"):
        model(p)
    with pytest.raises(ValueError, match="magnitude"):
        KerrMcmcModel(
            [2002.0], [], reference_epoch=2000.0,
            spin_vector=(0.0, 0.0, 1.1),
        )


def test_reference_epoch_is_explicit_and_reanchors_frame_drift():
    with pytest.raises(TypeError, match="reference_epoch"):
        KerrMcmcModel([2002.0], [2002.0])

    params = sample_parameters()
    epochs = np.array([2001.0, 2002.0, 2003.0])
    earlier = KerrMcmcModel(epochs, epochs, reference_epoch=2000.0)
    later = KerrMcmcModel(epochs, epochs, reference_epoch=2010.0)
    earlier.set_solver(method="dop853", rtol=1e-8, atol=1e-10)
    later.set_solver(method="dop853", rtol=1e-8, atol=1e-10)
    early_alpha, early_delta, early_velocity = earlier(params)
    late_alpha, late_delta, late_velocity = later(params)
    np.testing.assert_allclose(early_alpha - late_alpha, 10.0 * params[10], atol=1e-10)
    np.testing.assert_allclose(early_delta - late_delta, 10.0 * params[11], atol=1e-10)
    np.testing.assert_allclose(early_velocity, late_velocity, atol=1e-9)


@pytest.mark.parametrize(
    ("options", "error", "match"),
    [
        ({"rtol": "x"}, TypeError, "rtol must be a real number"),
        ({"rtol": True}, TypeError, "rtol must be a real number"),
        ({"atol": 0.0}, ValueError, "atol must be positive"),
        ({"max_steps": 5.0}, TypeError, "max_steps must be an integer"),
        ({"max_steps": 0}, ValueError, "max_steps must be a positive integer"),
    ],
)
def test_solver_options_are_validated_with_shared_rules(options, error, match):
    model = KerrMcmcModel([2002.0], [2002.0], reference_epoch=2000.0)
    settings = {"method": "dop853", "rtol": 1e-8, "atol": 1e-10} | options
    with pytest.raises(error, match=match):
        model.set_solver(**settings)


def test_reference_epoch_accepts_years_and_rejects_non_numbers():
    params = sample_parameters()
    by_number = KerrMcmcModel([2002.0], [2002.0], reference_epoch=2000.0)
    by_quantity = KerrMcmcModel([2002.0], [2002.0], reference_epoch=2000.0 * u.yr)
    for model in (by_number, by_quantity):
        model.set_solver(method="dop853", rtol=1e-8, atol=1e-10)
    np.testing.assert_array_equal(by_number(params)[0], by_quantity(params)[0])
    for invalid, error in (("2000", TypeError), (True, TypeError), (np.nan, ValueError)):
        with pytest.raises(error, match="reference_epoch"):
            KerrMcmcModel([2002.0], [2002.0], reference_epoch=invalid)


def test_missing_extension_is_not_an_integration_error(monkeypatch):
    import builtins

    model = KerrMcmcModel([2002.0], [2002.0], reference_epoch=2000.0)
    model.set_solver(method="dop853", rtol=1e-8, atol=1e-10)
    real_import = builtins.__import__

    def blocked(name, globals=None, locals=None, fromlist=(), level=0):
        if fromlist and "_mcmc_core" in fromlist:
            raise ImportError("blocked for test")
        return real_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", blocked)
    with pytest.raises(ImportError, match="build relatipy first"):
        model(sample_parameters())
