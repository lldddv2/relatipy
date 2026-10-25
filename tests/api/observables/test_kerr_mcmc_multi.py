"""Contract checks for the multi-star Kerr MCMC model and its likelihood."""

import math
import pickle

import numpy as np
import pytest
from astropy import units as u

from relatipy import KerrMcmcGaussianLikelihood, KerrMcmcModel, KerrMcmcMultiModel


SHARED = {
    "D": 8.2759, "M": 4.2996, "xS0": -7.2e-4, "yS0": -1.2e-4,
    "vxS0": 7.8e-5, "vyS0": 4.1e-5, "v_LSR": -2.0,
}
SPIN = {"spin_dec": 0.2, "spin_ra": -0.3, "spin_away": 0.5}
ELEMENTS = {
    "S2": (2018.3789, 0.12502, 0.88444, 134.67, 228.21, 66.279),
    "S38": (2022.6843, 0.14249, 0.81807, 168.69, 122.43, 40.065),
}
STAR_KEYS = ("t_p", "a", "e", "inc", "Omega", "omega")
SOLVER = {"method": "dop853", "rtol": 1e-9, "atol": 1e-12}
REFERENCE = 2000.0


def _epochs():
    return {
        "S2": (np.linspace(2024.0, 2004.0, 12), np.linspace(2010.0, 2022.0, 5)),
        "S38": (np.linspace(2005.0, 2023.0, 10), np.array([2012.0, 2021.0])),
    }


def _model(stars=("S2", "S38")):
    epochs = _epochs()
    model = KerrMcmcMultiModel(
        {star: epochs[star][0] for star in stars},
        {star: epochs[star][1] for star in stars},
        reference_epoch=REFERENCE,
    )
    model.set_solver(**SOLVER)
    return model


def _named(stars=("S2", "S38")):
    values = SHARED | SPIN
    for star in stars:
        t_p, a, e, inc, node, peri = ELEMENTS[star]
        values |= {
            f"{star}.t_p": t_p, f"{star}.a": a, f"{star}.e": e,
            f"{star}.inc": np.deg2rad(inc), f"{star}.Omega": np.deg2rad(node),
            f"{star}.omega": np.deg2rad(peri),
        }
    return values


def _legacy(values, star):
    return np.array([
        values["D"], values["M"],
        *(values[f"{star}.{key}"] for key in STAR_KEYS),
        values["xS0"], values["yS0"], values["vxS0"], values["vyS0"],
        values["v_LSR"],
    ])


def _spin(values):
    return (values["spin_dec"], values["spin_ra"], values["spin_away"])


def test_one_star_reproduces_kerr_mcmc_model_exactly():
    model = _model(("S2",))
    values = _named(("S2",))
    astrometry, spectroscopy = _epochs()["S2"]
    single = KerrMcmcModel(
        astrometry, spectroscopy, reference_epoch=REFERENCE,
        spin_vector=_spin(values),
    )
    single.set_solver(**SOLVER)

    result = model(model.pack(values))
    assert list(result) == ["S2"]
    for expected, actual in zip(single(_legacy(values, "S2")), result["S2"]):
        np.testing.assert_array_equal(actual, expected)
        assert not actual.flags.writeable
    assert [series.shape for series in result["S2"]] == [(12,), (12,), (5,)]


def test_each_star_matches_its_own_kerr_model_with_shared_parameters():
    model = _model()
    values = _named()
    result = model(model.pack(values))
    assert list(result) == ["S2", "S38"]
    for star in model.stars:
        astrometry, spectroscopy = _epochs()[star]
        single = KerrMcmcModel(astrometry, spectroscopy, reference_epoch=REFERENCE)
        single.set_solver(**SOLVER)
        expected = single(_legacy(values, star), spin_vector=_spin(values))
        for reference, actual in zip(expected, result[star]):
            np.testing.assert_array_equal(actual, reference)


@pytest.mark.parametrize(
    ("name", "step"),
    [("D", 0.05), ("M", 0.02), ("xS0", 1e-3), ("vyS0", 1e-4), ("v_LSR", 5.0),
     ("spin_ra", 0.4)],
)
def test_shared_parameters_change_every_star(name, step):
    model = _model()
    values = _named()
    baseline = model(model.pack(values))
    changed = model(model.pack(values | {name: values[name] + step}))
    for star in model.stars:
        difference = max(
            np.max(np.abs(new - old), initial=0.0)
            for new, old in zip(changed[star], baseline[star])
        )
        assert difference > 1e-12, (name, star)


def test_star_parameters_change_only_their_star():
    model = _model()
    values = _named()
    baseline = model(model.pack(values))
    changed = model(model.pack(values | {"S38.a": values["S38.a"] * 1.01}))
    for new, old in zip(changed["S2"], baseline["S2"]):
        np.testing.assert_array_equal(new, old)
    assert np.max(np.abs(changed["S38"][0] - baseline["S38"][0])) > 1e-6


def test_parameter_names_layout_and_packing():
    model = _model()
    shared = ("D", "M", "xS0", "yS0", "vxS0", "vyS0", "v_LSR")
    spin = ("spin_dec", "spin_ra", "spin_away")
    stars = tuple(f"{star}.{key}" for star in ("S2", "S38") for key in STAR_KEYS)
    assert model.stars == ("S2", "S38")
    assert model.parameter_names == shared + spin + stars
    assert model.n_parameters == 22
    assert dict(model.parameter_slices) == {
        "shared": slice(0, 7), "spin": slice(7, 10),
        "S2": slice(10, 16), "S38": slice(16, 22),
    }
    assert model.parameter_units["M"] == u.Unit(1e6 * u.M_sun)
    assert model.parameter_units["S2.inc"] == u.rad
    assert model.parameter_units["vxS0"] == u.arcsec / u.yr

    values = _named()
    vector = model.pack(values)
    assert vector.dtype == np.float64
    for block, names in (("shared", shared), ("spin", spin)):
        assert tuple(np.asarray(model.parameter_names)[model.parameter_slices[block]]) == names
    np.testing.assert_array_equal(
        vector[model.parameter_slices["S38"]],
        [values[f"S38.{key}"] for key in STAR_KEYS],
    )
    assert model.unpack(vector) == pytest.approx(values)
    assert list(model.unpack(vector)) == list(model.parameter_names)

    with_units = values | {
        "D": 8275.9 * u.pc, "M": 4.2996e6 * u.M_sun,
        "S2.inc": 134.67 * u.deg, "v_LSR": -2000.0 * u.m / u.s,
    }
    np.testing.assert_allclose(model.pack(with_units), vector, rtol=1e-15)


def test_epochs_accept_quantities_and_are_read_only():
    model = KerrMcmcMultiModel(
        {"S2": [2010.0, 2011.0] * u.yr}, {"S2": [365.25 * 2012] * u.day},
        reference_epoch=2000.0 * u.yr,
    )
    np.testing.assert_allclose(model.astrometry_times["S2"], [2010.0, 2011.0])
    np.testing.assert_allclose(model.spectroscopy_times["S2"], [2012.0])
    assert model.reference_epoch == 2000.0
    assert not model.astrometry_times["S2"].flags.writeable
    with pytest.raises(TypeError):
        model.astrometry_times["S2"] = np.zeros(2)


def _synthetic(model, vector, scale=1.0):
    errors = {"alpha": 5e-5, "delta": 6e-5, "v_los": 10.0}
    data = {}
    for star, (alpha, delta, v_los) in model(vector).items():
        data[star] = {
            "alpha": alpha, "delta": delta, "v_los": v_los,
            **{f"{key}_err": np.full(series.size, value * scale)
               for (key, value), series in zip(errors.items(), (alpha, delta, v_los))},
        }
    return data


def test_likelihood_is_maximal_at_truth_for_noise_free_data():
    model = _model()
    truth = model.pack(_named())
    data = _synthetic(model, truth)
    likelihood = KerrMcmcGaussianLikelihood(model, data)
    count = sum(series.size for star in data.values() for key, series in star.items()
                if not key.endswith("_err"))
    sigma_logs = sum(np.sum(np.log(star[f"{key}_err"])) for star in data.values()
                     for key in ("alpha", "delta", "v_los"))
    at_truth = likelihood(truth)
    assert likelihood.model is model
    assert at_truth == pytest.approx(-sigma_logs - 0.5 * count * math.log(2 * math.pi))

    steps = {"D": 1e-3, "M": 1e-4, "xS0": 1e-5, "v_LSR": 1.0, "spin_away": 0.3,
             "S2.t_p": 1e-4, "S2.e": 1e-5, "S38.Omega": 1e-4, "S38.a": 1e-5}
    index = {name: position for position, name in enumerate(model.parameter_names)}
    for name, step in steps.items():
        for sign in (-1.0, 1.0):
            moved = truth.copy()
            moved[index[name]] += sign * step
            assert likelihood(moved) < at_truth, name

    doubled = KerrMcmcGaussianLikelihood(model, _synthetic(model, truth, scale=2.0))
    assert doubled(truth) == pytest.approx(at_truth - count * math.log(2.0))


def test_likelihood_accepts_quantities():
    model = _model(("S2",))
    truth = model.pack(_named(("S2",)))
    data = _synthetic(model, truth)
    star = data["S2"]
    converted = {
        "alpha": star["alpha"] * u.arcsec, "delta": (star["delta"] * u.arcsec).to(u.mas),
        "v_los": (star["v_los"] * u.km / u.s).to(u.m / u.s),
        "alpha_err": star["alpha_err"] * u.arcsec,
        "delta_err": (star["delta_err"] * u.arcsec).to(u.mas),
        "v_los_err": star["v_los_err"] * u.km / u.s,
    }
    plain = KerrMcmcGaussianLikelihood(model, data)
    with_units = KerrMcmcGaussianLikelihood(model, {"S2": converted})
    assert with_units(truth) == pytest.approx(plain(truth), rel=1e-12)


def test_likelihood_returns_minus_infinity_on_native_failures():
    model = _model()
    truth = model.pack(_named())
    likelihood = KerrMcmcGaussianLikelihood(model, _synthetic(model, truth))
    rejected = truth.copy()
    rejected[model.parameter_names.index("S38.a")] = 1e-8
    assert likelihood(rejected) == -math.inf

    model.set_solver(method="dop853", rtol=1e-9, atol=1e-12, max_steps=3)
    assert likelihood(truth) == -math.inf
    with pytest.raises(ValueError, match="e must"):
        likelihood(np.where(np.arange(truth.size) == 12, 1.0, truth))


def test_model_and_likelihood_are_picklable():
    model = _model()
    truth = model.pack(_named())
    likelihood = KerrMcmcGaussianLikelihood(model, _synthetic(model, truth, 3.0))
    restored = pickle.loads(pickle.dumps(likelihood))
    assert restored(truth) == likelihood(truth)
    restored_model = pickle.loads(pickle.dumps(model))
    assert restored_model.parameter_names == model.parameter_names
    for star in model.stars:
        for expected, actual in zip(model(truth)[star], restored_model(truth)[star]):
            np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize(
    ("astrometry", "spectroscopy", "error", "match"),
    [
        ([2010.0], {"S2": [2010.0]}, TypeError, "mappings"),
        ({}, {}, ValueError, "at least one star"),
        ({"S2": [2010.0]}, {"S29": [2010.0]}, ValueError, "same stars"),
        ({1: [2010.0]}, {1: [2010.0]}, TypeError, "strings"),
        ({"S.2": [2010.0]}, {"S.2": [2010.0]}, ValueError, "without dots"),
        ({"": [2010.0]}, {"": [2010.0]}, ValueError, "nonempty"),
        ({"spin": [2010.0]}, {"spin": [2010.0]}, ValueError, "reserved"),
        ({"S2": [[2010.0]]}, {"S2": [2010.0]}, ValueError, "one-dimensional"),
        ({"S2": [np.nan]}, {"S2": [2010.0]}, ValueError, "finite"),
        ({"S2": []}, {"S2": []}, ValueError, "at least one observation"),
    ],
)
def test_construction_is_validated(astrometry, spectroscopy, error, match):
    with pytest.raises(error, match=match):
        KerrMcmcMultiModel(astrometry, spectroscopy, reference_epoch=REFERENCE)


def test_reference_epoch_and_solver_are_validated():
    with pytest.raises(TypeError, match="reference_epoch"):
        KerrMcmcMultiModel({"S2": [2010.0]}, {"S2": []}, reference_epoch="2000")
    model = KerrMcmcMultiModel({"S2": [2010.0]}, {"S2": []}, reference_epoch=REFERENCE)
    with pytest.raises(ValueError, match="set_solver"):
        model(np.zeros(model.n_parameters))
    with pytest.raises(ValueError, match="method"):
        model.set_solver(method="RK45", rtol=1e-6, atol=1e-9)
    with pytest.raises(ValueError, match="atol must be positive"):
        model.set_solver(method="dop853", rtol=1e-6, atol=0.0)
    with pytest.raises(ValueError, match="set_solver"):
        model(np.zeros(model.n_parameters))


@pytest.mark.parametrize(
    ("name", "value", "match"),
    [
        ("D", 0.0, "D and M"),
        ("M", -1.0, "D and M"),
        ("S38.a", 0.0, "a must"),
        ("S2.e", 1.0, "e must"),
        ("S38.e", -0.1, "e must"),
        ("spin_away", 1.01, "magnitude"),
        ("xS0", np.nan, "finite"),
    ],
)
def test_proposal_is_validated(name, value, match):
    model = _model()
    vector = model.pack(_named())
    vector[model.parameter_names.index(name)] = value
    with pytest.raises(ValueError, match=match):
        model(vector)


def test_wrong_length_and_bad_packing_are_rejected():
    model = _model()
    values = _named()
    with pytest.raises(ValueError, match="22 finite"):
        model(model.pack(values)[:-1])
    with pytest.raises(ValueError, match="22 numbers"):
        model.unpack(np.zeros(21))
    with pytest.raises(TypeError, match="mapping"):
        model.pack(list(values.values()))
    with pytest.raises(ValueError, match="missing parameters: S38.omega"):
        model.pack({key: value for key, value in values.items() if key != "S38.omega"})
    with pytest.raises(ValueError, match="unknown parameters: S29.a"):
        model.pack(values | {"S29.a": 0.39})
    with pytest.raises(TypeError, match="real number"):
        model.pack(values | {"D": True})
    with pytest.raises(ValueError, match="scalar"):
        model.pack(values | {"D": [8.0, 8.1]})
    with pytest.raises(ValueError, match="finite"):
        model.pack(values | {"D": np.inf})
    with pytest.raises(u.UnitConversionError):
        model.pack(values | {"D": 1.0 * u.s})


def test_likelihood_data_are_validated():
    model = _model()
    data = _synthetic(model, model.pack(_named()))
    with pytest.raises(TypeError, match="multi-star"):
        KerrMcmcGaussianLikelihood(KerrMcmcModel(), data)
    with pytest.raises(TypeError, match="mapping"):
        KerrMcmcGaussianLikelihood(model, [data])
    with pytest.raises(ValueError, match="missing stars: S38"):
        KerrMcmcGaussianLikelihood(model, {"S2": data["S2"]})
    with pytest.raises(ValueError, match="unknown stars: S29"):
        KerrMcmcGaussianLikelihood(model, data | {"S29": data["S2"]})
    without = {key: value for key, value in data["S2"].items() if key != "v_los_err"}
    with pytest.raises(ValueError, match="missing v_los_err"):
        KerrMcmcGaussianLikelihood(model, data | {"S2": without})
    with pytest.raises(ValueError, match="unknown keys: rv"):
        KerrMcmcGaussianLikelihood(model, data | {"S2": data["S2"] | {"rv": [1.0]}})
    with pytest.raises(ValueError, match="must have 12 values"):
        bad = data["S2"] | {"alpha": data["S2"]["alpha"][:-1]}
        KerrMcmcGaussianLikelihood(model, data | {"S2": bad})
    with pytest.raises(ValueError, match="positive"):
        bad = data["S2"] | {"delta_err": np.zeros(12)}
        KerrMcmcGaussianLikelihood(model, data | {"S2": bad})
    with pytest.raises(ValueError, match="finite"):
        bad = data["S2"] | {"v_los": np.full(5, np.nan)}
        KerrMcmcGaussianLikelihood(model, data | {"S2": bad})
    with pytest.raises(u.UnitConversionError):
        bad = data["S2"] | {"v_los": data["S2"]["v_los"] * u.arcsec}
        KerrMcmcGaussianLikelihood(model, data | {"S2": bad})
