"""Run-time Force and Potential compile to native callbacks and pickle by source."""

import pickle
import sys

import numpy as np
import pytest
from astropy import units as u

from relatipy import _force_core, _jit
from relatipy.forces import Force, Potential


def _has_compiler() -> bool:
    try:
        _jit.compiler_command()
    except _jit.CompilerNotFoundError:
        return False
    return True


pytest.importorskip("sympy")
pytestmark = pytest.mark.skipif(not _has_compiler(), reason="no C compiler")

POSITIONS = np.array([[10.0, 5.0, 3.0], [-40.0, 2.0, -7.0]])
VELOCITIES = np.array([[0.01, -0.02, 0.003], [0.0, 0.05, -0.01]])


@pytest.fixture(autouse=True)
def _cache(tmp_path_factory, monkeypatch) -> None:
    monkeypatch.setenv("RELATIPY_CACHE_DIR", str(tmp_path_factory.getbasetemp() / "jit"))


def _acceleration(force, mass=4e6 * u.M_sun) -> np.ndarray:
    bound = force._bind(mass)
    return _force_core.acceleration(bound.callback, bound.model, POSITIONS, VELOCITIES)


def test_cylindrical_potential_string_gives_minus_its_gradient() -> None:
    potential = Potential(
        "-k / sqrt(R**2 + z**2 + b**2)", coords="cylindrical",
        parameters={"k": 1e-3, "b": 50.0},
    )
    d2 = np.sum(POSITIONS**2, axis=1) + 2500.0
    expected = -1e-3 * POSITIONS / d2[:, None] ** 1.5
    np.testing.assert_allclose(_acceleration(potential), expected, rtol=1e-13)


def test_callable_force_in_cylindrical_components() -> None:
    force = Force(
        lambda rho, phi, z, vrho, vphi, vz, k: (-k * rho, -k * vphi, -k * z),
        coords="cylindrical", parameters={"k": 2.0},
    )
    rho = np.hypot(POSITIONS[:, 0], POSITIONS[:, 1])
    radial = POSITIONS[:, :2] / rho[:, None]
    azimuthal = np.column_stack((-radial[:, 1], radial[:, 0]))
    vphi = np.sum(VELOCITIES[:, :2] * azimuthal, axis=1)
    expected = np.column_stack((
        -2.0 * rho[:, None] * radial - 2.0 * vphi[:, None] * azimuthal,
        -2.0 * POSITIONS[:, 2],
    ))
    np.testing.assert_allclose(_acceleration(force), expected, rtol=1e-13, atol=1e-15)


def test_quantity_parameters_are_converted_at_bind_time() -> None:
    mass = 4e6 * u.M_sun
    r_g = (mass * 6.67430e-11 * u.m**3 / u.kg / u.s**2 / (299792458 * u.m / u.s) ** 2).to(u.pc)
    potential = Potential(
        "-m / sqrt(x**2 + y**2 + z**2 + b**2)",
        parameters={"m": 4e4 * u.M_sun, "b": 1000 * r_g},
    )
    np.testing.assert_allclose(
        potential._bind(mass).model_data, [0.01, 1000.0], rtol=1e-12
    )


def test_parameter_change_reuses_the_compiled_callback() -> None:
    potential = Potential("-k / sqrt(x**2 + y**2 + z**2 + 1)", parameters={"k": 1.0})
    first = potential._bind(1 * u.M_sun)
    stronger = potential.with_parameters(k=3.0)
    second = stronger._bind(1 * u.M_sun)
    assert second.callback == first.callback
    np.testing.assert_allclose(
        _force_core.acceleration(second.callback, second.model, POSITIONS, VELOCITIES),
        3 * _force_core.acceleration(first.callback, first.model, POSITIONS, VELOCITIES),
        rtol=1e-15,
    )


def test_pickle_uses_source_and_recompiles_from_cache_without_sympy(monkeypatch) -> None:
    potential = Potential(lambda x, y, z, k: -k / (x**2 + y**2 + z**2 + 4) ** 0.5,
                          parameters={"k": 0.5})
    expected = _acceleration(potential)
    payload = pickle.dumps(potential)
    monkeypatch.setitem(sys.modules, "sympy", None)
    restored = pickle.loads(payload)
    assert restored == potential and hash(restored) == hash(potential)
    assert "compiled source" in repr(restored)
    np.testing.assert_array_equal(_acceleration(restored), expected)


def test_missing_sympy_raises_an_import_error_naming_the_extra(monkeypatch) -> None:
    monkeypatch.setitem(sys.modules, "sympy", None)
    potential = Potential("-1 / sqrt(x**2 + y**2 + z**2 + 1)")
    with pytest.raises(ImportError, match=r"relatipy\[symbolic\]"):
        potential._bind(1 * u.M_sun)
    with pytest.raises(ImportError, match=r"relatipy\[symbolic\]"):
        potential.source


def test_missing_compiler_raises_a_runtime_error(monkeypatch, tmp_path) -> None:
    potential = Potential("-2 / sqrt(x**2 + y**2 + z**2 + 9)")
    potential.source
    monkeypatch.setenv("RELATIPY_CACHE_DIR", str(tmp_path))
    monkeypatch.setenv("CC", "relatipy-no-such-compiler")
    with pytest.raises(RuntimeError, match="no C compiler found"):
        potential._bind(1 * u.M_sun)


def test_invalid_definition_raises_value_error() -> None:
    with pytest.raises(ValueError):
        Potential("-k / r", coords="spherical").source


def test_zero_potential_reproduces_kerr_orbits_and_observables() -> None:
    from astropy.constants import c

    from relatipy import Kerr, KerrMcmcModel, PerturbedKerr, PerturbedKerrMcmcModel

    kerr = Kerr(mass=4e6 * u.M_sun, spin=0.6)
    orbit = {"a": 200 * kerr.r_g, "e": 0.5, "inc": 0.7 * u.rad}
    span = (0 * u.s, (2e4 * kerr.r_g / c).to(u.s))
    reference = kerr.orbit(**orbit).solve(tau_span=span, method="dop853")
    perturbed = PerturbedKerr(mass=kerr.mass, spin=0.6, force=Potential("0")).orbit(
        **orbit
    ).solve(tau_span=span, method="dop853")
    assert perturbed.tau[-1] == reference.tau[-1]
    np.testing.assert_allclose(
        perturbed.x[-1].to_value(u.m), reference.x[-1].to_value(u.m),
        rtol=1e-12, atol=1e-12 * kerr.r_g.to_value(u.m),
    )
    inputs = {
        "kerr": {"mass": kerr.mass, "spin": 0.6, "vec": (0.0, 0.3, 1.0)},
        "orbit": {"a": 1000 * u.au, "e": 0.8, "inc": 1.0 * u.rad},
        "sol": {"t_obs": np.linspace(0.0, 30.0, 9) * u.yr, "rtol": 1e-12, "atol": 1e-14},
        "distance": 8 * u.kpc,
    }
    expected = KerrMcmcModel().get_ra_dec_vr(**inputs)
    result = PerturbedKerrMcmcModel(force=Potential("0")).get_ra_dec_vr(**inputs)
    for one, other in zip(result, expected):
        np.testing.assert_allclose(one, other, rtol=0, atol=1e-10 * np.max(np.abs(other)))


def test_compiled_axisymmetric_potential_conserves_e_prime_and_l_prime() -> None:
    from astropy.constants import c

    from relatipy import Kerr, PerturbedKerr
    from relatipy.geodesic.solution import constants_of_motion

    kerr = Kerr(mass=1e6 * u.M_sun, spin=0.7)
    # Miyamoto-Nagai disk aligned with the spin, written in cylindrical
    # coordinates of the black-hole frame.
    disk = Potential(
        lambda R, phi, z, m, a, b: -m / (R**2 + (a + (z**2 + b**2) ** 0.5) ** 2) ** 0.5,
        coords="cylindrical",
        parameters={"m": 0.2 * kerr.mass, "a": 40 * kerr.r_g, "b": 10 * kerr.r_g},
    )
    metric = PerturbedKerr(mass=kerr.mass, spin=0.7, force=disk)
    solution = metric.orbit(
        a=120 * kerr.r_g, e=0.4, inc=0.9 * u.rad, Omega=0.4 * u.rad
    ).solve(tau_span=(0 * u.s, (1e5 * kerr.r_g / c).to(u.s)), method="dop853",
            rtol=1e-12)
    state = solution[:]
    constants = constants_of_motion(state)
    m, a, b = metric._binding().model_data
    rows = np.column_stack((
        np.zeros(len(solution)), (state.R / kerr.r_g).to_value(u.one),
        state.Theta.to_value(u.rad), state.Phi.to_value(u.rad),
    ))
    x, y, z = _force_core.bl_to_cartesian(0.7, np.ascontiguousarray(rows)).T
    weight = np.exp(-m / np.sqrt(x**2 + y**2 + (a + np.sqrt(z**2 + b**2)) ** 2))
    assert np.ptp(constants[:, 0]) / constants[0, 0] > 1e-4
    for column in (0, 1):
        invariant = constants[:, column] * weight
        assert np.max(np.abs(invariant / invariant[0] - 1)) < 1e-9
