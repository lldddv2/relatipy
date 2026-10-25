"""Contract checks for the perturbed multi-star model and its likelihood."""

import math
import multiprocessing
import pickle

import numpy as np
import pytest
from astropy import units as u

from relatipy import (
    KerrMcmcGaussianLikelihood, KerrMcmcMultiModel, PerturbedKerrMcmcModel,
    PerturbedKerrMcmcMultiModel,
)
from relatipy.forces import Plummer

SHARED = {
    "D": 8.2759, "M": 4.2996, "xS0": -7.2e-4, "yS0": -1.2e-4,
    "vxS0": 7.8e-5, "vyS0": 4.1e-5, "v_LSR": -2.0,
}
SPIN = {"spin_dec": 0.2, "spin_ra": -0.3, "spin_away": 0.5}
ELEMENTS = {
    "S2": (2018.3789, 0.12502, 0.88444, 134.67, 228.21, 66.279),
    "S38": (2022.6843, 0.14249, 0.81807, 168.69, 122.43, 40.065),
}
SOLVER = {"method": "dop853", "rtol": 1e-10, "atol": 1e-12}
FORCE = Plummer(mass=0 * u.M_sun, radius=1 * u.mpc)
EPOCHS = {
    "S2": (np.linspace(2024.0, 2004.0, 12), np.linspace(2010.0, 2022.0, 5)),
    "S38": (np.linspace(2005.0, 2023.0, 10), np.array([2012.0, 2021.0])),
}


def _model(cls=PerturbedKerrMcmcMultiModel, stars=("S2", "S38"), **kwargs):
    model = cls(
        {star: EPOCHS[star][0] for star in stars},
        {star: EPOCHS[star][1] for star in stars},
        reference_epoch=2000.0, **kwargs,
    )
    model.set_solver(**SOLVER)
    return model


def _named(stars=("S2", "S38"), force_mass=0.0, force_radius=1.0):
    values = SHARED | SPIN
    for star in stars:
        t_p, a, e, inc, node, peri = ELEMENTS[star]
        values |= {
            f"{star}.t_p": t_p, f"{star}.a": a, f"{star}.e": e,
            f"{star}.inc": np.deg2rad(inc), f"{star}.Omega": np.deg2rad(node),
            f"{star}.omega": np.deg2rad(peri),
        }
    return values, {"force.mass": force_mass, "force.radius": force_radius}


def _vector(model, **kwargs):
    values, force = _named(model.stars, **kwargs)
    return model.pack(values | force)


def _assert_equal(first, second):
    assert first.keys() == second.keys()
    for star in first:
        for one, other in zip(first[star], second[star]):
            np.testing.assert_array_equal(one, other)


def test_layout_appends_the_force_block_last() -> None:
    model = _model(force=FORCE)
    homogeneous = _model(KerrMcmcMultiModel)
    assert model.n_parameters == homogeneous.n_parameters + 2
    assert model.parameter_names[:-2] == homogeneous.parameter_names
    assert model.parameter_names[-2:] == ("force.mass", "force.radius")
    assert model.parameter_units["force.mass"] == u.M_sun
    assert model.parameter_units["force.radius"] == u.mpc
    assert model.parameter_slices["force"] == slice(
        model.n_parameters - 2, model.n_parameters
    )
    assert model.force is FORCE
    assert model.pack(
        _named()[0] | {"force.mass": 1 * u.M_sun, "force.radius": 2 * u.au}
    )[-1] == pytest.approx((2 * u.au).to_value(u.mpc))


def test_zero_force_reproduces_the_homogeneous_model() -> None:
    model = _model(force=FORCE)
    homogeneous = _model(KerrMcmcMultiModel)
    values, force = _named()
    perturbed = model(model.pack(values | force))
    reference = homogeneous(homogeneous.pack(values))
    for star in reference:
        for one, other in zip(perturbed[star], reference[star]):
            np.testing.assert_allclose(
                one, other, rtol=0, atol=1e-10 * np.max(np.abs(other))
            )


def test_one_star_equals_the_single_star_legacy_call() -> None:
    model = _model(stars=("S2",), force=FORCE)
    values, force = _named(("S2",), force_mass=4e3)
    vector = model.pack(values | force)
    single = PerturbedKerrMcmcModel(
        *EPOCHS["S2"], reference_epoch=2000.0,
        spin_vector=tuple(SPIN.values()), force=FORCE,
    )
    single.set_solver(**SOLVER)
    legacy = np.array([
        values["D"], values["M"],
        *(values[f"S2.{key}"] for key in ("t_p", "a", "e", "inc", "Omega", "omega")),
        values["xS0"], values["yS0"], values["vxS0"], values["vyS0"], values["v_LSR"],
        4e3, 1.0,
    ])
    _assert_equal(model(vector), {"S2": single(legacy)})


def test_in_place_updates_equal_fresh_models() -> None:
    model = _model(force=FORCE)
    first = _vector(model, force_mass=4e3)
    second = _vector(model, force_mass=1e3, force_radius=0.5)
    second[1] = 4.1
    result_first = model(first)
    result_second = model(second)
    _assert_equal(model(first), result_first)
    _assert_equal(_model(force=FORCE)(second), result_second)
    heavier = FORCE.with_parameters(mass=1e3 * u.M_sun, radius=0.5 * u.mpc)
    _assert_equal(_model(force=heavier)(second), result_second)
    assert np.max(np.abs(result_first["S2"][0] - result_second["S2"][0])) > 1e-4


def test_force_block_is_validated_before_native_evaluation() -> None:
    model = _model(force=FORCE)
    likelihood_data = {
        star: {
            "alpha": np.zeros(EPOCHS[star][0].size),
            "delta": np.zeros(EPOCHS[star][0].size),
            "v_los": np.zeros(EPOCHS[star][1].size),
            "alpha_err": np.ones(EPOCHS[star][0].size),
            "delta_err": np.ones(EPOCHS[star][0].size),
            "v_los_err": np.ones(EPOCHS[star][1].size),
        }
        for star in model.stars
    }
    likelihood = KerrMcmcGaussianLikelihood(model, likelihood_data)
    with pytest.raises(ValueError, match="mass must not be negative"):
        model(_vector(model, force_mass=-1.0))
    with pytest.raises(ValueError, match="mass must not be negative"):
        likelihood(_vector(model, force_mass=-1.0))
    with pytest.raises(ValueError, match="set_solver"):
        PerturbedKerrMcmcMultiModel(
            {"S2": EPOCHS["S2"][0]}, {"S2": EPOCHS["S2"][1]},
            reference_epoch=2000.0, force=FORCE,
        )(np.zeros(21))
    with pytest.raises(TypeError, match="ExternalForce"):
        PerturbedKerrMcmcMultiModel(
            {"S2": EPOCHS["S2"][0]}, {"S2": EPOCHS["S2"][1]},
            reference_epoch=2000.0, force=None,
        )
    with pytest.raises(ValueError, match="reserved"):
        PerturbedKerrMcmcMultiModel(
            {"force": EPOCHS["S2"][0]}, {"force": EPOCHS["S2"][1]},
            reference_epoch=2000.0, force=FORCE,
        )


def _likelihood(force_mass=4e3):
    model = _model(force=FORCE)
    truth = _vector(model, force_mass=force_mass)
    data = {}
    for star, (alpha, delta, v_los) in model(truth).items():
        data[star] = {
            "alpha": alpha, "delta": delta, "v_los": v_los,
            "alpha_err": np.full(alpha.size, 5e-5),
            "delta_err": np.full(delta.size, 5e-5),
            "v_los_err": np.full(v_los.size, 10.0),
        }
    return KerrMcmcGaussianLikelihood(model, data), truth


def test_likelihood_peaks_at_the_true_force_and_fails_soft() -> None:
    likelihood, truth = _likelihood()
    best = likelihood(truth)
    assert math.isfinite(best)
    lighter = truth.copy()
    lighter[-2] = 0.0
    assert likelihood(lighter) < best - 10
    unbound = truth.copy()
    unbound[likelihood.model.parameter_slices["S2"].start + 1] = 1e-9
    assert likelihood(unbound) == -math.inf


def test_model_and_likelihood_pickle_and_run_in_a_pool() -> None:
    likelihood, truth = _likelihood()
    expected = likelihood(truth)
    restored = pickle.loads(pickle.dumps(likelihood))
    assert restored(truth) == expected
    model = likelihood.model
    _assert_equal(pickle.loads(pickle.dumps(model))(truth), model(truth))
    lighter = truth.copy()
    lighter[-2] = 1e3
    context = multiprocessing.get_context("spawn")
    with context.Pool(2) as pool:
        results = pool.map(likelihood, [truth, lighter])
    assert results[0] == expected
    assert results[1] == likelihood(lighter)
