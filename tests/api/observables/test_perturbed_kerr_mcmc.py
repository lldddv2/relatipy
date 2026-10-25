"""Contract checks for the perturbed single-star Kerr observable model."""

import pickle

import numpy as np
import pytest
from astropy import units as u

from relatipy import IntegrationError, KerrMcmcModel, PerturbedKerrMcmcModel
from relatipy.forces import ThinDisk

MASS = 4.3e6 * u.M_sun
DISTANCE = 8.2 * u.kpc
KERR = {"mass": MASS, "spin": 0.5, "vec": (0.2, -0.3, 0.9)}
ORBIT = {
    "a": 1000 * u.au, "e": 0.88, "inc": 2.3 * u.rad,
    "Omega": 4.0 * u.rad, "omega": 1.2 * u.rad, "t": 2018.4 * u.yr,
}
TIMES = {"t_obs": np.linspace(2024.0, 2004.0, 15) * u.yr}
ZERO = ThinDisk(mass=0 * u.M_sun, r_out=1362.5 * u.au)
DISK = ThinDisk(mass=4e3 * u.M_sun, r_out=1362.5 * u.au)
# (disk_mass, normal_theta, normal_phi)
ZERO_BLOCK = [0.0, 0.0, 0.0]
DISK_BLOCK = [4e3, 0.0, 0.0]
# (D, M, t_p, a, e, inc, Omega, omega, xS0, yS0, vxS0, vyS0, v_LSR)
LEGACY = np.array([
    8.2759, 4.2996, 2018.3789, 0.12502, 0.88444, 2.3505, 3.9830, 1.1568,
    -7.2e-4, -1.2e-4, 7.8e-5, 4.1e-5, -2.0,
])
SPIN = (0.2, -0.3, 0.5)
EPOCHS = (np.linspace(2024.0, 2004.0, 12), np.linspace(2010.0, 2022.0, 5))


def _legacy(model_type, **kwargs):
    model = model_type(*EPOCHS, reference_epoch=2000.0, spin_vector=SPIN, **kwargs)
    model.set_solver(method="dop853", rtol=1e-10, atol=1e-12)
    return model


def _assert_close(first, second, rtol=1e-10):
    for one, other in zip(first, second):
        scale = np.max(np.abs(other))
        np.testing.assert_allclose(one, other, rtol=0, atol=rtol * scale)


@pytest.mark.parametrize("method", ["dop853", "radau"])
@pytest.mark.parametrize(
    ("kind", "rtol"),
    # Proper-time sampling amplifies the roundoff of the two right-hand
    # sides; that mode's own integration error at this tolerance is larger.
    [("t_obs", 1e-10), ("t_eval", 1e-10), ("tau_eval", 1e-6)],
)
def test_zero_force_reproduces_the_general_kerr_evaluator(method, kind, rtol) -> None:
    sol = {kind: TIMES["t_obs"], "method": method, "rtol": 1e-12, "atol": 1e-14}
    reference = KerrMcmcModel().get_ra_dec_vr(
        kerr=KERR, orbit=ORBIT, sol=sol, distance=DISTANCE
    )
    perturbed = PerturbedKerrMcmcModel(force=ZERO).get_ra_dec_vr(
        kerr=KERR, orbit=ORBIT, sol=sol, distance=DISTANCE
    )
    _assert_close(perturbed, reference, rtol)
    assert all(not value.flags.writeable for value in perturbed)


def test_force_argument_overrides_the_model_force() -> None:
    model = PerturbedKerrMcmcModel(force=ZERO)
    plain = model.get_ra_dec_vr(kerr=KERR, orbit=ORBIT, sol=TIMES, distance=DISTANCE)
    forced = model.get_ra_dec_vr(
        kerr=KERR, orbit=ORBIT, sol=TIMES, distance=DISTANCE, force=DISK
    )
    same = PerturbedKerrMcmcModel(force=DISK).get_ra_dec_vr(
        kerr=KERR, orbit=ORBIT, sol=TIMES, distance=DISTANCE
    )
    for one, other in zip(forced, same):
        np.testing.assert_array_equal(one, other)
    assert np.max(np.abs(forced[0] - plain[0])) > 1e-4
    again = model.get_ra_dec_vr(kerr=KERR, orbit=ORBIT, sol=TIMES, distance=DISTANCE)
    for one, other in zip(again, plain):
        np.testing.assert_array_equal(one, other)


def test_legacy_call_appends_the_force_block() -> None:
    reference = _legacy(KerrMcmcModel)(LEGACY)
    model = _legacy(PerturbedKerrMcmcModel, force=ZERO)
    assert model.n_force_parameters == 3
    _assert_close(model(np.append(LEGACY, ZERO_BLOCK)), reference)
    forced = model(np.append(LEGACY, DISK_BLOCK))
    same = _legacy(PerturbedKerrMcmcModel, force=DISK)(np.append(LEGACY, DISK_BLOCK))
    for one, other in zip(forced, same):
        np.testing.assert_array_equal(one, other)
    assert np.max(np.abs(forced[0] - reference[0])) > 1e-4


def test_legacy_call_validates_vector_and_force_block() -> None:
    model = _legacy(PerturbedKerrMcmcModel, force=ZERO)
    with pytest.raises(ValueError, match="16 finite numbers"):
        model(LEGACY)
    with pytest.raises(ValueError, match="mass must not be negative"):
        model(np.append(LEGACY, [-1.0, 0.0, 0.0]))
    unset = PerturbedKerrMcmcModel(*EPOCHS, reference_epoch=2000.0, force=ZERO)
    with pytest.raises(ValueError, match="set_solver"):
        unset(np.append(LEGACY, ZERO_BLOCK))
    with pytest.raises(ValueError, match="observation epochs"):
        PerturbedKerrMcmcModel(force=ZERO)(np.append(LEGACY, ZERO_BLOCK))


def test_inputs_are_validated_like_the_kerr_model() -> None:
    with pytest.raises(TypeError, match="ExternalForce"):
        PerturbedKerrMcmcModel(force="thin_disk")
    model = PerturbedKerrMcmcModel(force=ZERO)
    with pytest.raises(TypeError, match="ExternalForce"):
        model.get_ra_dec_vr(kerr=KERR, orbit=ORBIT, sol=TIMES, distance=DISTANCE,
                            force=1.0)
    with pytest.raises(ValueError, match="kerr requires"):
        model.get_ra_dec_vr(kerr={"mass": MASS}, orbit=ORBIT, sol=TIMES,
                            distance=DISTANCE)
    with pytest.raises(ValueError, match="method"):
        model.get_ra_dec_vr(kerr=KERR, orbit=ORBIT, sol=TIMES | {"method": "dp45"},
                            distance=DISTANCE)


def test_integration_failure_raises_integration_error() -> None:
    model = PerturbedKerrMcmcModel(force=DISK)
    with pytest.raises(IntegrationError):
        model.get_ra_dec_vr(kerr=KERR, orbit=ORBIT, sol=TIMES | {"max_steps": 3},
                            distance=DISTANCE)


def test_model_pickles_without_its_native_handle() -> None:
    model = _legacy(PerturbedKerrMcmcModel, force=DISK)
    vector = np.append(LEGACY, DISK_BLOCK)
    expected = model(vector)
    restored = pickle.loads(pickle.dumps(model))
    assert restored.force == DISK
    for one, other in zip(restored(vector), expected):
        np.testing.assert_array_equal(one, other)
