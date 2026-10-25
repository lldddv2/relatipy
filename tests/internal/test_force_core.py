"""Direct checks of relatipy._force_core and the private bound-force handles."""

import numpy as np
import pytest
from astropy import constants as const
from astropy import units as u

from relatipy import _force_core, _mcmc_core
from relatipy.forces import ThinDisk
from relatipy.forces._bound import BoundForce, geometric_conversion

SPIN = 0.6
POSITIONS = np.array([[8.0, -3.0, 1.5], [-120.0, 40.0, -300.0]])
VELOCITIES = np.zeros((2, 3))
COMPONENTS = np.array([0.0, 4000.0, 0.8, 0.7, 0.3, 0.2, 0.0])
DISK = ThinDisk(mass=4e3 * u.M_sun, r_out=1 * u.mpc)


def _general(method=1, time_kind=2, times=np.array([1.0])):
    return (np.eye(3), "elements", COMPONENTS, times, 0.0, time_kind, 1e-8,
            method, 1e-10, 1e-12, 1000)


def test_null_force_is_bitwise_the_homogeneous_observables() -> None:
    rotation = np.eye(3)
    times = np.linspace(-3e6, 3e6, 9)
    args = (rotation, 4000.0, 0.8, 0.7, 0.3, 0.2, times, 1e-8,
            1e-9, -2e-9, 3e-12, 4e-12, 1e-5, 0.5, 0, 1e-10, 1e-12, 100000)
    reference = _mcmc_core.evaluate(SPIN, *args)
    np.testing.assert_array_equal(_force_core.evaluate(SPIN, 0, 0, *args), reference)
    for method in (0, 1):
        for time_kind in (0, 1, 2):
            general = (rotation, "elements", COMPONENTS, np.sort(times), 0.0,
                       time_kind, 1e-8, method, 1e-10, 1e-12, 100000)
            np.testing.assert_array_equal(
                _force_core.evaluate_general(SPIN, 0, 0, *general),
                _mcmc_core.evaluate_general(SPIN, *general),
            )


def test_disk_callback_reads_its_placed_model() -> None:
    bound = DISK._bind(4e6 * u.M_sun)
    first = _force_core.acceleration(bound.callback, bound.model, POSITIONS, VELOCITIES)
    assert first.shape == (2, 3) and np.all(np.isfinite(first))
    heavier = DISK.with_parameters(disk_mass=8e3 * u.M_sun)._bind(4e6 * u.M_sun)
    np.testing.assert_allclose(
        _force_core.acceleration(heavier.callback, heavier.model, POSITIONS, VELOCITIES),
        2 * first, rtol=1e-14,
    )


def test_failing_callback_maps_to_the_existing_errors() -> None:
    # The thin-disk callback fails on a NULL model.
    address = DISK._bind(4e6 * u.M_sun).callback
    invalid = 0
    with pytest.raises(RuntimeError, match="callback failed at row 0"):
        _force_core.acceleration(address, invalid, POSITIONS, VELOCITIES)
    rotation = np.eye(3)
    with pytest.raises(RuntimeError, match="callback failure"):
        _force_core.evaluate(
            SPIN, address, invalid, rotation, 4000.0, 0.8, 0.7, 0.3,
            0.2, np.array([1.0]), 1e-8, 0, 0, 0, 0, 0, 0, 0, 1e-10, 1e-12, 1000,
        )
    with pytest.raises(RuntimeError, match="callback failure"):
        _force_core.evaluate_general(SPIN, address, invalid, *_general(method=0))


def test_methods_and_arguments_are_validated() -> None:
    with pytest.raises(ValueError, match="projection_radau"):
        _force_core.evaluate_general(SPIN, 0, 0, *_general(method=3))
    with pytest.raises(ValueError, match="method"):
        _force_core.evaluate_general(SPIN, 0, 0, *_general(method=7))
    with pytest.raises(ValueError, match="projection_radau"):
        _force_core.evaluate(
            SPIN, 0, 0, np.eye(3), 4000.0, 0.8, 0.7, 0.3, 0.2, np.array([1.0]),
            1e-8, 0, 0, 0, 0, 0, 0, 3, 1e-10, 1e-12, 1000,
        )
    with pytest.raises(ValueError, match="time_kind"):
        _force_core.evaluate_general(SPIN, 0, 0, *_general(time_kind=5))
    with pytest.raises(ValueError, match="at least one"):
        _force_core.evaluate_general(SPIN, 0, 0, *_general(times=np.empty(0)))
    with pytest.raises(ValueError, match="nonzero"):
        _force_core.acceleration(0, 0, POSITIONS, VELOCITIES)
    with pytest.raises(ValueError, match="shape"):
        _force_core.acceleration(1, 0, POSITIONS[:, :2].copy(), VELOCITIES)
    with pytest.raises(ValueError, match="normal"):
        DISK._bind(4e6 * u.M_sun)._buffer.set(1.0, 1.0, [0.0, 1.0])
    with pytest.raises(ValueError, match="invalid placed force model"):
        DISK._bind(4e6 * u.M_sun)._buffer.set(1.0, -1.0, [0.0, 0.0, 1.0])


def test_geometric_conversion_uses_length_time_and_mass_dimensions() -> None:
    mass = 4e6 * u.M_sun
    r_g = const.G * mass / const.c**2
    t_g = r_g / const.c
    cases = [
        (u.pc, (1 * u.pc / r_g).decompose()),
        (u.M_sun, (1 * u.M_sun / mass).decompose()),
        (u.km / u.s, (1 * u.km / u.s / const.c).decompose()),
        (u.yr, (1 * u.yr / t_g).decompose()),
        (u.pc**3 / u.yr**2, (1 * u.pc**3 / u.yr**2 / (r_g**3 / t_g**2)).decompose()),
        (u.deg, (1 * u.deg).to(u.rad)),
        (u.one, 1.0 * u.one),
    ]
    mass_kg = mass.to_value(u.kg)
    for unit, expected in cases:
        factor, power = geometric_conversion(unit, "p")
        assert factor * mass_kg ** (-power) == pytest.approx(expected.value, rel=1e-12)
    with pytest.raises(u.UnitConversionError):
        geometric_conversion(u.K, "p")


def test_in_place_update_equals_a_fresh_binding() -> None:
    bound = DISK._bind(4.3e6 * u.M_sun)
    assert isinstance(bound, BoundForce)
    address = bound.model
    other = DISK.with_parameters(disk_mass=1e3 * u.M_sun, normal_theta=0.5, normal_phi=1.0)
    bound.update((4.1e6 * u.M_sun).to_value(u.kg), np.array([1e3, 0.5, 1.0]))
    assert bound.model == address
    np.testing.assert_array_equal(
        bound.model_data, other._bind(4.1e6 * u.M_sun).model_data
    )
    r_g = const.G * 4.1e6 * u.M_sun / const.c**2
    normal = [np.sin(0.5) * np.cos(1.0), np.sin(0.5) * np.sin(1.0), np.cos(0.5)]
    np.testing.assert_allclose(
        bound.model_data, [1e3 / 4.1e6, (1 * u.mpc / r_g).decompose(), *normal],
        rtol=1e-12,
    )
    assert bound.callback == DISK._native_callback()[1].callback
    assert not bound.model_data.flags.writeable
