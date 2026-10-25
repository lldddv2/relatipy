"""Forces are immutable, validated, hashable and picklable values."""

import copy
import pickle

import numpy as np
import pytest
from astropy import units as u

import relatipy as rp
from relatipy.forces import ExternalForce, Force, Plummer, Potential


def test_plummer_stores_parameters_with_their_units() -> None:
    sphere = Plummer(mass=1e4 * u.M_sun, radius=0.01 * u.pc)
    assert isinstance(sphere, ExternalForce)
    assert sphere.parameter_names == ("mass", "radius")
    assert dict(sphere.parameter_units) == {"mass": u.M_sun, "radius": u.pc}
    assert sphere.parameters["mass"] == 1e4 * u.M_sun
    assert not sphere.parameters["radius"].flags.writeable
    assert repr(sphere) == (
        "Plummer(mass=<Quantity 10000. solMass>, radius=<Quantity 0.01 pc>)"
    )


def test_package_exports_the_force_layer() -> None:
    assert rp.Force is Force and rp.Potential is Potential and rp.Plummer is Plummer
    assert rp.forces.Plummer is Plummer
    assert {"Force", "Potential", "Plummer", "PerturbedKerr", "forces"} <= set(rp.__all__)


@pytest.mark.parametrize(
    ("kwargs", "error"),
    [
        ({"mass": 1.0, "radius": 1 * u.pc}, TypeError),
        ({"mass": 1 * u.M_sun, "radius": 1 * u.s}, u.UnitConversionError),
        ({"mass": -1 * u.M_sun, "radius": 1 * u.pc}, ValueError),
        ({"mass": 1 * u.M_sun, "radius": 0 * u.pc}, ValueError),
        ({"mass": np.inf * u.M_sun, "radius": 1 * u.pc}, ValueError),
        ({"mass": [1, 2] * u.M_sun, "radius": 1 * u.pc}, ValueError),
    ],
)
def test_plummer_rejects_invalid_parameters(kwargs, error) -> None:
    with pytest.raises(error):
        Plummer(**kwargs)


def test_forces_are_immutable() -> None:
    sphere = Plummer(mass=1 * u.M_sun, radius=1 * u.au)
    with pytest.raises(AttributeError, match="immutable"):
        sphere._parameters = {}
    with pytest.raises(AttributeError, match="immutable"):
        del sphere._native


def test_equality_and_hash_follow_type_and_parameters() -> None:
    sphere = Plummer(mass=1 * u.M_sun, radius=1 * u.au)
    assert sphere == Plummer(mass=1 * u.M_sun, radius=1 * u.au)
    assert hash(sphere) == hash(Plummer(mass=1 * u.M_sun, radius=1 * u.au))
    assert sphere != Plummer(mass=2 * u.M_sun, radius=1 * u.au)
    assert sphere != "Plummer"


def test_with_parameters_returns_a_validated_copy() -> None:
    sphere = Plummer(mass=1 * u.M_sun, radius=1 * u.au)
    heavier = sphere.with_parameters(mass=3 * u.M_sun)
    assert heavier == Plummer(mass=3 * u.M_sun, radius=1 * u.au)
    assert sphere.parameters["mass"] == 1 * u.M_sun
    with pytest.raises(ValueError, match="unknown parameters: b"):
        sphere.with_parameters(b=1 * u.au)
    with pytest.raises(ValueError, match="radius must be positive"):
        sphere.with_parameters(radius=-1 * u.au)


@pytest.mark.parametrize(
    "roundtrip",
    [copy.copy, copy.deepcopy, lambda value: pickle.loads(pickle.dumps(value))],
)
def test_plummer_survives_copy_and_pickle(roundtrip) -> None:
    sphere = Plummer(mass=1 * u.M_sun, radius=1 * u.au)
    assert roundtrip(sphere) == sphere


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
