"""Public null initial conditions, units, validation and type boundaries."""

import numpy as np
import pytest
from astropy import units as u
from astropy.constants import c

import relatipy
from relatipy.metrics import Kerr


def _position(metric, family, *, unit=u.m):
    radius = (10 * metric.r_g).to(unit)
    if family == "cartesian":
        return {"x": radius}
    if family == "spherical":
        return {"r": radius, "theta": 90 * u.deg, "phi": 0 * u.rad}
    return {"R": radius, "Theta": 90 * u.deg, "Phi": 0 * u.rad}


def _constants(metric):
    # Stay clear of radial turning points, whose potential can round negative.
    return {
        "b": 3 * metric.r_g,
        "eta": 5 * metric.r_g**2,
        "radial_sign": 1,
        "polar_sign": -1,
    }


def _direction(metric, family):
    if family == "cartesian":
        return {"vy": 0.1 * c}
    angular = (0.1 * c / (10 * metric.r_g)) * u.rad
    return {"vphi" if family == "spherical" else "vPhi": angular}


@pytest.mark.parametrize("family", ["cartesian", "spherical", "bl"])
@pytest.mark.parametrize("mode", ["direction", "constants"])
@pytest.mark.parametrize("time", [None, 3 * u.ms])
def test_null_families_modes_and_initial_time(family, mode, time) -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0)
    options = _position(metric, family, unit=u.km)
    options.update(
        _direction(metric, family) if mode == "direction" else _constants(metric)
    )
    if time is not None:
        options["t"] = time
    photon = metric.null(**options)
    expected_time = 0 * u.s if time is None else time
    assert isinstance(photon, relatipy.Null)
    assert photon.t == expected_time
    assert photon.initial.t == expected_time
    assert photon.t.unit == expected_time.unit
    assert photon.initial.t.unit == expected_time.unit
    for name in ("x", "y", "z", "r", "R", "b"):
        assert getattr(photon, name).unit == u.km
    assert photon.eta.unit == u.km**2
    np.testing.assert_allclose((photon.R / metric.r_g).to_value(u.one), 10, rtol=1e-12)
    if family == "cartesian":
        # Canonical BL reconstruction can leave roundoff from sin/cos(pi/2).
        for name in ("y", "z"):
            np.testing.assert_allclose(
                (getattr(photon, name) / metric.r_g).to_value(u.one),
                0, rtol=0, atol=1e-14,
            )
    if mode == "constants":
        np.testing.assert_allclose((photon.b / metric.r_g).to_value(u.one), 3, rtol=1e-12)
        np.testing.assert_allclose((photon.eta / metric.r_g**2).to_value(u.one), 5, rtol=1e-12)
        assert photon.vR > 0 * u.km / u.s
        assert photon.vTheta < 0 * u.rad / u.s


def test_null_uses_first_supplied_cartesian_position_unit() -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0)
    photon = metric.null(y=(10 * metric.r_g).to(u.km), z=1 * u.m, vx=c, t=2 * u.ms)
    np.testing.assert_allclose(
        (photon.x / metric.r_g).to_value(u.one), 0, rtol=0, atol=1e-14,
    )
    assert photon.x.unit == u.km
    assert photon.z.unit == u.km
    assert photon.vx.unit == u.km / u.ms
    assert photon.t.unit == u.ms


def test_null_schwarzschild_photon_sphere_impact_parameter() -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0)
    photon = metric.null(
        R=3 * metric.r_g, Theta=np.pi / 2 * u.rad, Phi=0 * u.rad,
        vPhi=1 * u.rad / u.s,
    )
    # Analytic equatorial Schwarzschild null circle: b/r_g = 3 sqrt(3), Q = 0.
    np.testing.assert_allclose(
        (photon.b / metric.r_g).to_value(u.one), 3 * np.sqrt(3), rtol=1e-10,
    )
    np.testing.assert_allclose(
        (photon.eta / metric.r_g**2).to_value(u.one), 0, rtol=0, atol=1e-12,
    )


def test_null_schwarzschild_radial_coordinate_speed() -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0)
    photon = metric.null(**_position(metric, "bl"), vR=0.01 * c)
    # Outgoing radial null condition at R/r_g=10: dR/dt = (1-2/10)c.
    np.testing.assert_allclose((photon.vR / c).to_value(u.one), 0.8, rtol=1e-12)
    assert photon.vTheta == 0 * u.rad / u.s
    assert photon.vPhi == 0 * u.rad / u.s


def test_null_schwarzschild_tangential_coordinate_speed() -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0)
    photon = metric.null(**_position(metric, "bl"), vPhi=1 * u.rad / u.s)
    # Equatorial tangential null condition: (R dPhi/dt)^2 = (1-2/R)c^2.
    tangential = photon.R * photon.vPhi / u.rad
    np.testing.assert_allclose((tangential / c).to_value(u.one), np.sqrt(0.8), rtol=1e-12)
    assert photon.vR == 0 * u.m / u.s
    assert photon.vTheta == 0 * u.rad / u.s


def test_null_direction_constants_round_trip_outside_equatorial_plane() -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0.7)
    radius = 10 * metric.r_g
    position = {"R": radius, "Theta": 0.9 * u.rad, "Phi": 0.4 * u.rad}
    direction = metric.null(
        **position, vR=0.15 * c,
        vTheta=(-0.08 * c / radius) * u.rad,
        vPhi=(0.2 * c / radius) * u.rad,
    )
    rebuilt = metric.null(
        **position, b=direction.b, eta=direction.eta,
        radial_sign=int(np.sign(direction.vR.value)),
        polar_sign=int(np.sign(direction.vTheta.value)),
    )
    for name in ("vR", "vTheta", "vPhi", "vx", "vy", "vz"):
        expected = getattr(direction, name)
        np.testing.assert_allclose(
            getattr(rebuilt, name).to_value(expected.unit), expected.value,
            rtol=1e-9, atol=0,
        )


@pytest.mark.parametrize("options", [{}, {"vx": 1 * u.m / u.s}, {"b": 1 * u.m}])
def test_null_requires_position_family(options) -> None:
    with pytest.raises(ValueError, match="position"):
        Kerr(mass=1 * u.Msun, spin=0).null(**options)


@pytest.mark.parametrize("options", [
    {"x": 1 * u.m, "R": 1 * u.m},
    {"x": 1 * u.m, "r": 1 * u.m},
    {"r": 1 * u.m, "R": 1 * u.m},
])
def test_null_rejects_mixed_position_families(options) -> None:
    with pytest.raises(ValueError, match="incompatible"):
        Kerr(mass=1 * u.Msun, spin=0).null(**options)


@pytest.mark.parametrize("family", ["spherical", "bl"])
@pytest.mark.parametrize("missing_index", [0, 1, 2])
def test_null_requires_all_spherical_and_bl_position_components(family, missing_index) -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0)
    options = _position(metric, family)
    options.pop(tuple(options)[missing_index])
    options.update(_direction(metric, family))
    with pytest.raises(ValueError, match="requires"):
        metric.null(**options)


@pytest.mark.parametrize("family, velocity", [
    ("cartesian", "vR"), ("cartesian", "vr"),
    ("spherical", "vx"), ("spherical", "vR"),
    ("bl", "vx"), ("bl", "vr"),
])
def test_null_rejects_velocity_from_another_family(family, velocity) -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0)
    with pytest.raises(ValueError, match="incompatible"):
        metric.null(**_position(metric, family), **{velocity: c})


@pytest.mark.parametrize("key", ["b", "eta", "radial_sign", "polar_sign"])
def test_null_rejects_mixing_direction_and_constants(key) -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0)
    with pytest.raises(ValueError, match="either"):
        metric.null(x=10 * metric.r_g, vx=c, **{key: _constants(metric)[key]})


@pytest.mark.parametrize("missing", ["b", "eta", "radial_sign", "polar_sign"])
def test_null_requires_complete_constants(missing) -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0)
    options = _constants(metric)
    options.pop(missing)
    with pytest.raises(ValueError, match="requires"):
        metric.null(**_position(metric, "bl"), **options)


@pytest.mark.parametrize("key", ["radial_sign", "polar_sign"])
@pytest.mark.parametrize("value", [0, 2, 1.5])
def test_null_rejects_sign_outside_integer_plus_or_minus_one(key, value) -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0)
    options = _constants(metric) | {key: value}
    with pytest.raises(ValueError, match="sign"):
        metric.null(**_position(metric, "bl"), **options)


@pytest.mark.parametrize("key", ["radial_sign", "polar_sign"])
@pytest.mark.parametrize("value", [True, "1"])
def test_null_rejects_boolean_and_string_signs(key, value) -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0)
    with pytest.raises(TypeError):
        metric.null(**_position(metric, "bl"), **(_constants(metric) | {key: value}))


@pytest.mark.parametrize("options", [
    {"x": 1, "vx": c},
    {"x": 10 * u.km, "vx": 1},
    {"x": 10 * u.km, "vx": c, "t": 0},
    {"r": 10 * u.km, "theta": 1, "phi": 0 * u.rad, "vr": c},
])
def test_null_direction_inputs_require_units(options) -> None:
    with pytest.raises(TypeError):
        Kerr(mass=1 * u.Msun, spin=0).null(**options)


@pytest.mark.parametrize("key", ["b", "eta"])
def test_null_constants_require_units(key) -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0)
    with pytest.raises(TypeError):
        metric.null(**_position(metric, "bl"), **(_constants(metric) | {key: 1}))


@pytest.mark.parametrize("options", [
    {"x": 1 * u.kg, "vx": c},
    {"x": 10 * u.km, "vx": 1 * u.s},
    {"x": 10 * u.km, "vx": c, "t": 1 * u.kg},
])
def test_null_direction_inputs_reject_incompatible_units(options) -> None:
    with pytest.raises(u.UnitConversionError):
        Kerr(mass=1 * u.Msun, spin=0).null(**options)


@pytest.mark.parametrize("key, value", [("b", 1 * u.s), ("eta", 1 * u.m)])
def test_null_constants_reject_incompatible_units(key, value) -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0)
    with pytest.raises(u.UnitConversionError):
        metric.null(**_position(metric, "bl"), **(_constants(metric) | {key: value}))


@pytest.mark.parametrize("key, unit", [
    ("R", u.m), ("Theta", u.rad), ("vR", u.m / u.s), ("t", u.s),
    ("b", u.m), ("eta", u.m**2),
])
@pytest.mark.parametrize("value", [np.nan, np.inf, np.array([1.0])])
def test_null_rejects_nonfinite_or_nonscalar_input(key, unit, value) -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0)
    options = _position(metric, "bl")
    options.update(_constants(metric) if key in ("b", "eta") else {"vR": c})
    options[key] = value * unit
    with pytest.raises(ValueError):
        metric.null(**options)


@pytest.mark.parametrize("explicit_zero", [False, True])
def test_null_rejects_zero_or_omitted_direction(explicit_zero) -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0)
    options = {"x": 10 * metric.r_g}
    if explicit_zero:
        options.update(vx=0 * c, vy=0 * c, vz=0 * c)
    with pytest.raises(ValueError, match="direction"):
        metric.null(**options)


@pytest.mark.parametrize("factor", [1.0, 1 + 1e-7])
@pytest.mark.parametrize("mode", ["direction", "constants"])
def test_null_rejects_horizon_and_margin_within_valid_bl_chart(factor, mode) -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0.7)
    options = _position(metric, "bl") | {"R": factor * metric.horizons.event}
    options.update({"vR": c} if mode == "direction" else _constants(metric))
    with pytest.raises(ValueError, match="horizon"):
        metric.null(**options)


def test_null_rejects_ambiguous_future_roots_in_ergoregion() -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0.9)
    with pytest.raises(ValueError, match="ergoregion"):
        metric.null(
            R=1.8 * metric.r_g, Theta=np.pi / 2 * u.rad, Phi=0 * u.rad,
            vPhi=1 * u.rad / u.s,
        )


def test_null_rejects_direction_without_real_future_tangent() -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0.9)
    with pytest.raises(ValueError, match="null"):
        metric.null(
            R=1.8 * metric.r_g, Theta=np.pi / 2 * u.rad, Phi=0 * u.rad, vR=c,
        )


def test_null_rejects_constants_without_real_tangent() -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0)
    with pytest.raises(ValueError, match="null"):
        metric.null(
            **_position(metric, "bl"),
            **(_constants(metric) | {"b": 100 * metric.r_g}),
        )


@pytest.mark.parametrize("theta", [0, np.pi])
def test_null_rejects_polar_coordinate_singularity(theta) -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0)
    with pytest.raises(ValueError):
        metric.null(**(_position(metric, "bl") | {"Theta": theta * u.rad}), vR=c)


@pytest.mark.parametrize("key, value", [
    ("tau", 0 * u.s), ("a", 10 * u.km), ("e", 0.1),
    ("p", 10 * u.km), ("elements", object()),
])
def test_null_rejects_orbit_only_keywords(key, value) -> None:
    metric = Kerr(mass=1 * u.Msun, spin=0)
    with pytest.raises(TypeError):
        metric.null(x=10 * metric.r_g, vx=c, **{key: value})


def test_null_public_reexports_are_distinct_from_timelike_types() -> None:
    assert relatipy.Null is relatipy.geodesic.Null
    assert relatipy.NullSolution is relatipy.geodesic.NullSolution
    assert not issubclass(relatipy.Null, relatipy.Orbit)
    assert not issubclass(relatipy.NullSolution, relatipy.Solution)
