"""ThinDisk: validation, value semantics, shared tables, orientation and field."""

import copy
import pickle

import numpy as np
import pytest
from astropy import units as u

import relatipy as rp
from relatipy import ThinDisk, _force_core
from relatipy.forces import ExternalForce

MASS = 4.3e6 * u.M_sun
FIDUCIAL = ThinDisk(mass=1e3 * u.M_sun, r_out=1362.5 * u.au)


def _disk(mass=1 * u.M_sun, r_out=1 * u.au, **kwargs):
    return ThinDisk(mass=mass, r_out=r_out, **kwargs)


def _acceleration(bound, positions):
    return _force_core.acceleration(
        bound.callback, bound.model, np.ascontiguousarray(positions),
        np.zeros_like(positions),
    )


def test_package_exports_the_thin_disk() -> None:
    assert rp.ThinDisk is ThinDisk and rp.forces.ThinDisk is ThinDisk
    assert {"ThinDisk", "forces"} <= set(rp.__all__)
    for name in ("Force", "Potential", "Plummer", "Spheroid", "PerturbedKerr"):
        assert not hasattr(rp, name) and not hasattr(rp.forces, name)


def test_thin_disk_stores_parameters_with_their_units() -> None:
    disk = ThinDisk(mass=1e4 * u.M_sun, r_out=0.01 * u.pc)
    assert isinstance(disk, ExternalForce)
    assert disk.parameter_names == ("disk_mass", "normal_theta", "normal_phi")
    assert dict(disk.parameter_units) == {
        "disk_mass": u.M_sun, "normal_theta": u.rad, "normal_phi": u.rad,
    }
    assert disk.parameters["disk_mass"] == 1e4 * u.M_sun
    assert not disk.parameters["disk_mass"].flags.writeable
    assert repr(disk) == (
        "ThinDisk(mass=<Quantity 10000. solMass>, r_out=<Quantity 0.01 pc>, "
        "r_in=<Quantity 0. pc>, slope=-1.0, softening=0.01, normal=(0.0, 0.0, 1.0))"
    )


def test_parameters_units_and_value_semantics() -> None:
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


@pytest.mark.parametrize(
    ("factory", "error"),
    [
        (lambda: ThinDisk(mass=1.0, r_out=1 * u.au), TypeError),
        (lambda: ThinDisk(mass=1 * u.M_sun, r_out=1.0), TypeError),
        (lambda: ThinDisk(mass=1 * u.M_sun, r_out=1 * u.s), u.UnitConversionError),
        (lambda: ThinDisk(mass=1 * u.s, r_out=1 * u.au), u.UnitConversionError),
        (lambda: ThinDisk(mass=-1 * u.M_sun, r_out=1 * u.au), ValueError),
        (lambda: ThinDisk(mass=np.inf * u.M_sun, r_out=1 * u.au), ValueError),
        (lambda: ThinDisk(mass=[1, 2] * u.M_sun, r_out=1 * u.au), ValueError),
        (lambda: ThinDisk(mass=1 * u.M_sun, r_out=0 * u.au), ValueError),
        (lambda: ThinDisk(mass=1 * u.M_sun, r_out=1 * u.au, r_in=1 * u.au), ValueError),
        (lambda: ThinDisk(mass=1 * u.M_sun, r_out=1 * u.au, slope=-2), ValueError),
        (lambda: ThinDisk(mass=1 * u.M_sun, r_out=1 * u.au, slope=True), TypeError),
        (lambda: ThinDisk(mass=1 * u.M_sun, r_out=1 * u.au, softening=0), ValueError),
        (lambda: ThinDisk(mass=1 * u.M_sun, r_out=1 * u.au, normal=(0, 0, 0)), ValueError),
        (lambda: ThinDisk(mass=1 * u.M_sun, r_out=1 * u.au, normal=(0, 1)), ValueError),
    ],
)
def test_invalid_shapes_are_rejected(factory, error) -> None:
    with pytest.raises(error):
        factory()


def test_forces_are_immutable() -> None:
    disk = _disk()
    with pytest.raises(AttributeError, match="immutable"):
        disk._parameters = {}
    with pytest.raises(AttributeError, match="immutable"):
        del disk._native


def test_equality_and_hash_follow_type_and_parameters() -> None:
    disk = _disk()
    assert disk == _disk()
    assert hash(disk) == hash(_disk())
    assert disk != _disk(mass=2 * u.M_sun)
    assert disk != _disk(r_out=2 * u.au)
    assert disk != "ThinDisk"


def test_with_parameters_returns_a_validated_copy() -> None:
    disk = _disk()
    heavier = disk.with_parameters(disk_mass=3 * u.M_sun)
    assert heavier == _disk(mass=3 * u.M_sun)
    assert disk.parameters["disk_mass"] == 1 * u.M_sun
    with pytest.raises(ValueError, match="unknown parameters: b"):
        disk.with_parameters(b=1 * u.au)
    with pytest.raises(ValueError, match="mass must not be negative"):
        disk.with_parameters(disk_mass=-1 * u.M_sun)


@pytest.mark.parametrize(
    "roundtrip",
    [copy.copy, copy.deepcopy, lambda value: pickle.loads(pickle.dumps(value))],
)
def test_thin_disk_survives_copy_and_pickle(roundtrip) -> None:
    disk = _disk(normal=(0.0, 1.0, 1.0))
    assert roundtrip(disk) == disk


def test_tables_are_shared_per_shape_and_rebuilt_lazily_after_pickle() -> None:
    first = FIDUCIAL._bind(MASS)
    other_radius = ThinDisk(mass=5e2 * u.M_sun, r_out=2000 * u.au)._bind(MASS)
    restored = pickle.loads(pickle.dumps(FIDUCIAL))
    assert restored._native is None
    again = restored._bind(MASS)
    assert first._keepalive is other_radius._keepalive is again._keepalive
    assert first.model != again.model
    np.testing.assert_array_equal(first.model_data, again.model_data)
    assert first._keepalive.kind == "thin_disk" and first._keepalive.nbytes > 0


def test_bound_handles_pickle_by_rebinding() -> None:
    bound = FIDUCIAL._bind(MASS)
    restored = pickle.loads(pickle.dumps(bound))
    assert restored.force == FIDUCIAL and restored.mass == bound.mass
    assert restored.model != bound.model
    np.testing.assert_array_equal(restored.model_data, bound.model_data)


def test_normal_rotates_the_field_and_observer_axes_map_through_the_spin() -> None:
    aligned = ThinDisk(mass=1e3 * u.M_sun, r_out=1000 * u.au, slope=-0.5)
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


def test_thin_disk_field_is_plane_symmetric_and_scales_with_mass() -> None:
    bound = FIDUCIAL._bind(MASS)
    heavier = FIDUCIAL.with_parameters(disk_mass=3e3 * u.M_sun)._bind(MASS)
    scale = bound.model_data[1]
    above = scale * np.array([[0.3, 0.1, 0.02], [1.5, -0.4, 0.3]])
    below = above * np.array([1.0, 1.0, -1.0])
    first = _acceleration(bound, above)
    np.testing.assert_allclose(
        _acceleration(bound, below), first * np.array([1.0, 1.0, -1.0]), rtol=1e-13
    )
    np.testing.assert_allclose(_acceleration(heavier, above), 3 * first, rtol=1e-13)
    assert np.all(first[:, 2] < 0)
