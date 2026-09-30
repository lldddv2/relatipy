"""Coordinate classes keep one identity across package and module imports."""

from importlib import import_module

import pytest


@pytest.mark.parametrize(
    ("name", "module_name"),
    [
        ("CartesianCoordinates", "cartesian"),
        ("SphericalCoordinates", "spherical"),
        ("BoyerLindquistCoordinates", "boyer_lindquist"),
        ("SphericalVelocity", "spherical"),
        ("BoyerLindquistVelocity", "boyer_lindquist"),
        ("SphericalFourVelocity", "spherical"),
        ("BoyerLindquistFourVelocity", "boyer_lindquist"),
        ("CartesianStateVector", "cartesian"),
        ("OrbitalElements", "orbital_elements"),
        ("KerrOrbitalElements", "kerr_orbital_elements"),
    ],
)
def test_coordinate_imports_share_class_identity(name: str, module_name: str) -> None:
    package = import_module("relatipy.coordinates")
    module = import_module(f"relatipy.coordinates.{module_name}")

    assert getattr(package, name) is getattr(module, name)
