"""Forces are immutable, validated, hashable and picklable values."""

import copy
import pickle

import numpy as np
import pytest
from astropy import units as u

import relatipy as rp
from relatipy.forces import ExternalForce, Force, Potential, ThinDisk


def _disk(mass=1 * u.M_sun, r_out=1 * u.au, **kwargs):
    return ThinDisk(mass=mass, r_out=r_out, **kwargs)


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


def test_package_exports_the_force_layer() -> None:
    assert rp.Force is Force and rp.Potential is Potential and rp.ThinDisk is ThinDisk
    assert rp.forces.ThinDisk is ThinDisk
    assert {"Force", "Potential", "ThinDisk", "PerturbedKerr", "forces"} <= set(rp.__all__)


@pytest.mark.parametrize(
    ("kwargs", "error"),
    [
        ({"mass": 1.0, "r_out": 1 * u.pc}, TypeError),
        ({"mass": 1 * u.M_sun, "r_out": 1 * u.s}, u.UnitConversionError),
        ({"mass": -1 * u.M_sun, "r_out": 1 * u.pc}, ValueError),
        ({"mass": 1 * u.M_sun, "r_out": 0 * u.pc}, ValueError),
        ({"mass": np.inf * u.M_sun, "r_out": 1 * u.pc}, ValueError),
        ({"mass": [1, 2] * u.M_sun, "r_out": 1 * u.pc}, ValueError),
    ],
)
def test_thin_disk_rejects_invalid_parameters(kwargs, error) -> None:
    with pytest.raises(error):
        ThinDisk(**kwargs)


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


def test_compiled_force_constructors_validate_without_sympy() -> None:
    potential = Potential("-k / r", coords="spherical", parameters={"k": 1.0})
    assert potential.coords == "spherical"
    assert potential.parameter_names == ("k",)
    assert potential.parameter_units["k"] == u.one
    assert repr(potential) == (
        "Potential('-k / r', coords='spherical', parameters={'k': <Quantity 1.>})"
    )
    with pytest.raises(ValueError, match="coords"):
        Potential("0", coords="polar")
    with pytest.raises(TypeError, match="definition"):
        Potential(("x", "y", "z"))
    with pytest.raises(TypeError, match="definition"):
        Force(3.0)
    with pytest.raises(TypeError, match="parameters must be a mapping"):
        Potential("0", parameters=[("k", 1.0)])
    with pytest.raises(ValueError, match="identifier"):
        Potential("0", parameters={"not a name": 1.0})
    with pytest.raises(TypeError, match="real number"):
        Potential("0", parameters={"k": "1"})
    with pytest.raises(u.UnitConversionError, match="length, time, mass and angle"):
        Potential("0", parameters={"T": 1 * u.K})
    assert Force(("0", "0", "0")).coords == "cartesian"
