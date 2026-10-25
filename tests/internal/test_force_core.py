"""Direct checks of relatipy._force_core and the private bound-force handles."""

import numpy as np
import pytest
from astropy import constants as const
from astropy import units as u

from relatipy import _core, _force_core, _jit, _mcmc_core
from relatipy.forces import Potential, ThinDisk
from relatipy.forces._bound import BoundForce, geometric_conversion

SPIN = 0.6
STATE = _core.initial_kerr(
    SPIN, "elements", np.array([0.0, 200.0, 0.5, 0.7, 0.3, 0.2, 0.1])
)
SPAN = (0.0, 2e4)
POSITIONS = np.array([[8.0, -3.0, 1.5], [-120.0, 40.0, -300.0]])
VELOCITIES = np.zeros((2, 3))


DISK = ThinDisk(mass=4e3 * u.M_sun, r_out=1 * u.mpc)


def _compiles() -> bool:
    try:
        import sympy  # noqa: F401

        _jit.compiler_command()
    except (ImportError, _jit.CompilerNotFoundError):
        return False
    return True


needs_compiled = pytest.mark.skipif(not _compiles(), reason="needs SymPy and a C compiler")


def _plummer_potential() -> Potential:
    return Potential("-m / sqrt(x**2 + y**2 + z**2 + b**2)",
                     parameters={"m": 4e3 * u.M_sun, "b": 1 * u.mpc})


@pytest.mark.parametrize("method", ["dop853", "dp45", "radau"])
def test_null_force_is_bitwise_the_homogeneous_solve(method) -> None:
    reference = _core.integrate_kerr(SPIN, STATE, *SPAN, method, 1e-10, 1e-12)
    forced = _force_core.integrate_kerr(SPIN, 0, 0, STATE, *SPAN, method, 1e-10, 1e-12)
    np.testing.assert_array_equal(forced[0], reference[0])
    np.testing.assert_array_equal(forced[1], reference[1])
    assert forced[3] == reference[3] == 0
    assert forced[2]["force_active"] is False
    assert {key: forced[2][key] for key in reference[2]} == reference[2]


def test_null_force_is_bitwise_the_homogeneous_observables() -> None:
    rotation = np.eye(3)
    times = np.linspace(-3e6, 3e6, 9)
    args = (rotation, 4000.0, 0.8, 0.7, 0.3, 0.2, times, 1e-8,
            1e-9, -2e-9, 3e-12, 4e-12, 1e-5, 0.5, 0, 1e-10, 1e-12, 100000)
    reference = _mcmc_core.evaluate(SPIN, *args)
    np.testing.assert_array_equal(_force_core.evaluate(SPIN, 0, 0, *args), reference)
    components = np.array([0.0, 4000.0, 0.8, 0.7, 0.3, 0.2, 0.0])
    general = (rotation, "elements", components, np.sort(times), 0.0, 2, 1e-8,
               1, 1e-10, 1e-12, 100000)
    np.testing.assert_array_equal(
        _force_core.evaluate_general(SPIN, 0, 0, *general),
        _mcmc_core.evaluate_general(SPIN, *general),
    )


@needs_compiled
def test_compiled_callback_matches_the_formula() -> None:
    sphere = Potential("-m / sqrt(x**2 + y**2 + z**2 + b**2)",
                       parameters={"m": 0.03, "b": 4.0})
    bound = sphere._bind(4e6 * u.M_sun)
    d2 = np.sum(POSITIONS**2, axis=1) + 16.0
    np.testing.assert_allclose(
        _force_core.acceleration(bound.callback, bound.model, POSITIONS, VELOCITIES),
        -0.03 * POSITIONS / d2[:, None] ** 1.5, rtol=1e-15,
    )


def test_bl_to_cartesian_uses_the_oblate_black_hole_frame() -> None:
    rows = np.array([[0.0, 10.0, 1.1, 0.3], [5.0, 3.0, 0.2, -2.0]])
    position = _force_core.bl_to_cartesian(SPIN, rows)
    rho = np.sqrt(rows[:, 1] ** 2 + SPIN**2) * np.sin(rows[:, 2])
    np.testing.assert_allclose(position[:, 0], rho * np.cos(rows[:, 3]), rtol=1e-15)
    np.testing.assert_allclose(position[:, 1], rho * np.sin(rows[:, 3]), rtol=1e-15)
    np.testing.assert_allclose(position[:, 2], rows[:, 1] * np.cos(rows[:, 2]), rtol=1e-15)
    with pytest.raises(ValueError, match="spin"):
        _force_core.bl_to_cartesian(1.5, rows)


def test_failing_callback_maps_to_the_existing_errors() -> None:
    # The thin-disk callback fails on a NULL model.
    address = DISK._bind(4e6 * u.M_sun).callback
    invalid = 0
    with pytest.raises(RuntimeError, match="callback failed at row 0"):
        _force_core.acceleration(address, invalid, POSITIONS, VELOCITIES)
    with pytest.raises(RuntimeError, match="external force callback failed"):
        _force_core.four_acceleration(SPIN, address, invalid,
                                      STATE[:4], STATE[4:])
    assert _force_core.CALLBACK_FAILURE == 32
    taus, states, stats, status = _force_core.integrate_kerr(
        SPIN, address, invalid, STATE, *SPAN, "dop853", 1e-10, 1e-12
    )
    assert status == -1 and stats["native_status"] == 4 and len(taus) == 1
    rotation = np.eye(3)
    with pytest.raises(RuntimeError, match="callback failure"):
        _force_core.evaluate(
            SPIN, address, invalid, rotation, 4000.0, 0.8, 0.7, 0.3,
            0.2, np.array([1.0]), 1e-8, 0, 0, 0, 0, 0, 0, 0, 1e-10, 1e-12, 1000,
        )
    with pytest.raises(RuntimeError, match="callback failure"):
        _force_core.evaluate_general(
            SPIN, address, invalid, rotation, "elements",
            np.array([0.0, 4000.0, 0.8, 0.7, 0.3, 0.2, 0.0]), np.array([1.0]),
            0.0, 0, 1e-8, 0, 1e-10, 1e-12, 1000,
        )


def test_methods_and_arguments_are_validated() -> None:
    with pytest.raises(ValueError, match="projection_radau"):
        _force_core.integrate_kerr(SPIN, 0, 0, STATE, *SPAN, "projection_radau")
    with pytest.raises(ValueError, match="method"):
        _force_core.integrate_kerr(SPIN, 0, 0, STATE, *SPAN, "rk4")
    with pytest.raises(ValueError, match="y0"):
        _force_core.integrate_kerr(SPIN, 0, 0, STATE[:7], *SPAN)
    with pytest.raises(ValueError, match="projection_radau"):
        _force_core.evaluate_general(
            SPIN, 0, 0, np.eye(3), "elements",
            np.array([0.0, 4000.0, 0.8, 0.7, 0.3, 0.2, 0.0]), np.array([1.0]),
            0.0, 0, 1e-8, 3, 1e-10, 1e-12, 1000,
        )
    with pytest.raises(ValueError, match="nonzero"):
        _force_core.acceleration(0, 0, POSITIONS, VELOCITIES)
    with pytest.raises(OverflowError):
        _force_core.integrate_kerr(SPIN, -1, 0, STATE, *SPAN)


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


@needs_compiled
def test_in_place_update_of_parameter_storage_equals_a_fresh_binding() -> None:
    sphere = _plummer_potential()
    bound = sphere._bind(4.3e6 * u.M_sun)
    address = bound.model
    other = sphere.with_parameters(m=1e3 * u.M_sun, b=2 * u.mpc)
    bound.update((4.1e6 * u.M_sun).to_value(u.kg), np.array([1e3, 2.0]))
    assert bound.model == address
    np.testing.assert_array_equal(
        bound.model_data, other._bind(4.1e6 * u.M_sun).model_data
    )
    np.testing.assert_allclose(
        bound.model_data,
        [1e3 / 4.1e6, (2 * u.mpc / (const.G * 4.1e6 * u.M_sun / const.c**2)).decompose()],
        rtol=1e-12,
    )
    assert bound.callback == sphere._native_callback()[0]
    assert not bound.model_data.flags.writeable
