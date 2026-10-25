"""ThinDisk and Spheroid: validation, shared tables, orientation and physics."""

import pickle

import numpy as np
import pytest
from astropy import units as u
from astropy.constants import c

from relatipy import (
    Kerr, KerrMcmcMultiModel, PerturbedKerr, PerturbedKerrMcmcMultiModel,
    Spheroid, ThinDisk, _force_core,
)
from relatipy.forces import ExternalForce
from relatipy.geodesic.solution import constants_of_motion

MASS = 4.3e6 * u.M_sun
FIDUCIAL = ThinDisk(mass=1e3 * u.M_sun, r_out=1362.5 * u.au)


def _acceleration(bound, positions):
    return _force_core.acceleration(
        bound.callback, bound.model, np.ascontiguousarray(positions),
        np.zeros_like(positions),
    )


def test_parameters_units_and_value_semantics() -> None:
    assert isinstance(FIDUCIAL, ExternalForce)
    assert FIDUCIAL.parameter_names == ("disk_mass", "normal_theta", "normal_phi")
    assert dict(FIDUCIAL.parameter_units) == {
        "disk_mass": u.M_sun, "normal_theta": u.rad, "normal_phi": u.rad,
    }
    tilted = ThinDisk(mass=1e3 * u.M_sun, r_out=1362.5 * u.au, normal=(0, 2, 2))
    np.testing.assert_allclose(tilted.normal, [0, 2**-0.5, 2**-0.5], atol=1e-15)
    assert tilted.parameters["normal_theta"].to_value(u.rad) == pytest.approx(np.pi / 4)
    assert tilted.parameters["normal_phi"].to_value(u.rad) == pytest.approx(np.pi / 2)
    assert tilted != FIDUCIAL
    assert FIDUCIAL == ThinDisk(mass=1e3 * u.M_sun, r_out=1362.5 * u.au)
    assert FIDUCIAL != ThinDisk(mass=1e3 * u.M_sun, r_out=1362.5 * u.au, slope=-0.5)
    assert repr(FIDUCIAL).startswith("ThinDisk(mass=<Quantity 1000. solMass>, r_out=")
    heavier = FIDUCIAL.with_parameters(disk_mass=2e3 * u.M_sun, normal_theta=0.3)
    assert heavier.parameters["normal_theta"].value == 0.3
    assert heavier == pickle.loads(pickle.dumps(heavier))
    sphere = Spheroid(mass=1e3 * u.M_sun, s_out=1000 * u.au)
    assert sphere.parameter_names == ("mass", "normal_theta", "normal_phi")
    assert pickle.loads(pickle.dumps(sphere)) == sphere


@pytest.mark.parametrize(
    ("factory", "error"),
    [
        (lambda: ThinDisk(mass=1.0, r_out=1 * u.au), TypeError),
        (lambda: ThinDisk(mass=1 * u.M_sun, r_out=1.0), TypeError),
        (lambda: ThinDisk(mass=1 * u.M_sun, r_out=1 * u.s), u.UnitConversionError),
        (lambda: ThinDisk(mass=-1 * u.M_sun, r_out=1 * u.au), ValueError),
        (lambda: ThinDisk(mass=1 * u.M_sun, r_out=0 * u.au), ValueError),
        (lambda: ThinDisk(mass=1 * u.M_sun, r_out=1 * u.au, r_in=1 * u.au), ValueError),
        (lambda: ThinDisk(mass=1 * u.M_sun, r_out=1 * u.au, slope=-2), ValueError),
        (lambda: ThinDisk(mass=1 * u.M_sun, r_out=1 * u.au, softening=0), ValueError),
        (lambda: ThinDisk(mass=1 * u.M_sun, r_out=1 * u.au, normal=(0, 0, 0)), ValueError),
        (lambda: ThinDisk(mass=1 * u.M_sun, r_out=1 * u.au, normal=(0, 1)), ValueError),
        (lambda: Spheroid(mass=1 * u.M_sun, s_out=1 * u.au, axis_ratio=0), ValueError),
        (lambda: Spheroid(mass=1 * u.M_sun, s_out=1 * u.au, slope=-3), ValueError),
    ],
)
def test_invalid_shapes_are_rejected(factory, error) -> None:
    with pytest.raises(error):
        factory()


def test_tables_are_shared_per_shape_and_rebuilt_lazily_after_pickle() -> None:
    first = FIDUCIAL._bind(MASS)
    other_radius = ThinDisk(mass=5e2 * u.M_sun, r_out=2000 * u.au)._bind(MASS)
    restored = pickle.loads(pickle.dumps(FIDUCIAL))
    assert restored._native is None
    again = restored._bind(MASS)
    assert first._keepalive is other_radius._keepalive is again._keepalive
    assert first.model != again.model
    np.testing.assert_array_equal(first.model_data, again.model_data)


def test_spheroid_sphere_matches_the_enclosed_mass_field() -> None:
    sphere = Spheroid(mass=1e3 * u.M_sun, s_out=1000 * u.au, slope=-2.0)
    bound = sphere._bind(MASS)
    mass, scale = bound.model_data[:2]
    directions = np.array([[1.0, 0, 0], [0.3, -0.5, 0.8], [0, 0, -1], [0.6, 0.8, 0]])
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    for radius in scale * np.array([1e-3, 0.05, 0.4, 0.9, 1.5, 8.0]):
        positions = radius * directions
        enclosed = mass * min(1.0, radius / scale)
        expected = -enclosed * positions / radius**3
        np.testing.assert_allclose(
            _acceleration(bound, positions), expected, rtol=1e-10,
            atol=1e-10 * enclosed / radius**2,
        )


def test_normal_rotates_the_field_and_observer_axes_map_through_the_spin() -> None:
    aligned = Spheroid(mass=1e3 * u.M_sun, s_out=1000 * u.au, axis_ratio=0.5)
    normal = np.array([0.3, -0.4, 0.866])
    normal /= np.linalg.norm(normal)
    tilted = aligned.with_parameters(
        normal_theta=np.arccos(normal[2]), normal_phi=np.arctan2(normal[1], normal[0])
    )
    flat, rotated = aligned._bind(MASS), tilted._bind(MASS)
    rotation = np.array(_rotation_from_normal(normal))
    body = np.array([[2000.0, 500.0, 300.0], [-800.0, 100.0, -50.0]])
    np.testing.assert_allclose(
        _acceleration(rotated, body),
        _acceleration(flat, body @ rotation.T) @ rotation, rtol=1e-12, atol=1e-25,
    )
    spin_frame = np.linalg.qr(np.array([[1.0, 0.2, 0.1], [0.1, 1.0, 0.3], [0.2, 0.1, 1.0]]))[0]
    observer = tilted._bind(MASS, rotation=spin_frame)
    np.testing.assert_allclose(observer.model_data[2:], spin_frame.T @ normal, rtol=1e-12)


def _rotation_from_normal(normal):
    """Rows: two axes orthogonal to ``normal``, then ``normal``."""
    helper = np.array([1.0, 0.0, 0.0]) if abs(normal[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    first = np.cross(helper, normal)
    first /= np.linalg.norm(first)
    return first, np.cross(normal, first), normal


def test_thin_disk_far_field_tends_to_the_monopole() -> None:
    bound = FIDUCIAL._bind(MASS)
    mass, scale = bound.model_data[:2]
    positions = 30 * scale * np.array([[1.0, 0, 0], [0, 0.6, 0.8]])
    expected = -mass * positions / (30 * scale) ** 3
    np.testing.assert_allclose(_acceleration(bound, positions), expected, rtol=2e-3)
    with pytest.raises(RuntimeError, match="callback failed"):
        _acceleration(bound, np.zeros((1, 3)))


def test_thin_disk_orbit_conserves_e_prime_to_the_interpolation_accuracy() -> None:
    kerr = Kerr(mass=1e6 * u.M_sun, spin=0.7)
    disk = ThinDisk(mass=1e4 * u.M_sun, r_out=1000 * kerr.r_g)
    metric = PerturbedKerr(mass=kerr.mass, spin=0.7, force=disk)
    solution = metric.orbit(
        a=300 * kerr.r_g, e=0.5, inc=0.6 * u.rad, Omega=0.4 * u.rad
    ).solve(tau_span=(0 * u.s, (2e5 * kerr.r_g / c).to(u.s)), method="dop853",
            rtol=1e-12)
    assert solution.status == 0
    state = solution[:]
    constants = constants_of_motion(state)
    rows = np.column_stack((
        np.zeros(len(solution)), (state.R / kerr.r_g).to_value(u.one),
        state.Theta.to_value(u.rad), state.Phi.to_value(u.rad),
    ))
    positions = _force_core.bl_to_cartesian(0.7, np.ascontiguousarray(rows))
    potential = metric._binding()._buffer.potential(positions)
    weight = np.exp(potential)
    drift_kerr = np.ptp(constants[:, 0]) / constants[0, 0]
    for column in (0, 1):
        invariant = constants[:, column] * weight
        drift = np.max(np.abs(invariant / invariant[0] - 1))
        # Phi and the acceleration are separate interpolants (about 1e-5
        # apart in the softened plane), so the bound is loose.
        assert drift < 1e-6
        assert drift < 1e-2 * drift_kerr


def _multi(force):
    astrometry = {"S2": np.linspace(2004.0, 2024.0, 10)}
    spectroscopy = {"S2": np.array([2010.0, 2018.4, 2022.0])}
    if force is None:
        model = KerrMcmcMultiModel(astrometry, spectroscopy, reference_epoch=2000.0)
    else:
        model = PerturbedKerrMcmcMultiModel(
            astrometry, spectroscopy, reference_epoch=2000.0, force=force
        )
    model.set_solver(method="dop853", rtol=1e-10, atol=1e-12)
    values = {
        "D": 8.2759, "M": 4.2996, "xS0": 0.0, "yS0": 0.0, "vxS0": 0.0,
        "vyS0": 0.0, "v_LSR": 0.0, "spin_dec": 0.2, "spin_ra": -0.3,
        "spin_away": 0.5, "S2.t_p": 2018.3789, "S2.a": 0.12502, "S2.e": 0.88444,
        "S2.inc": 2.3505, "S2.Omega": 3.9830, "S2.omega": 1.1568,
    }
    return model, values


def test_disk_force_block_in_the_multi_star_model() -> None:
    model, values = _multi(FIDUCIAL)
    assert model.parameter_names[-3:] == (
        "force.disk_mass", "force.normal_theta", "force.normal_phi"
    )
    homogeneous, _ = _multi(None)
    reference = homogeneous(homogeneous.pack(values))["S2"]
    zero = model(model.pack(values | {
        "force.disk_mass": 0.0, "force.normal_theta": 0.4, "force.normal_phi": 1.0,
    }))["S2"]
    for one, other in zip(zero, reference):
        np.testing.assert_allclose(one, other, rtol=0, atol=1e-10 * np.max(np.abs(other)))
    face_on = model.pack(values | {
        "force.disk_mass": 1e4, "force.normal_theta": 0.0, "force.normal_phi": 0.0,
    })
    edge_on = face_on.copy()
    edge_on[-2] = np.pi / 2
    first, second = model(face_on)["S2"], model(edge_on)["S2"]
    assert np.max(np.abs(first[0] - second[0])) > 1e-6
    fresh, _ = _multi(FIDUCIAL)
    for one, other in zip(fresh(face_on)["S2"], first):
        np.testing.assert_array_equal(one, other)
    with pytest.raises(ValueError, match="mass must not be negative"):
        bad = face_on.copy()
        bad[-3] = -1.0
        model(bad)
