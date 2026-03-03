"""
Coordinate systems module for relatipy.

This module contains implementations of different coordinate systems
used in general relativity, including Cartesian, spherical,
and Boyer-Lindquist coordinates, and cylindrical coordinates.
"""

from .base import CoordinateBase
from .cartesian import Cartesian
from .spherical import Spherical
from .boyer_lindquist import BoyerLindquist
from .cilindrical import Cylindrical
from .orbital_elements import OrbitalElements

# Dictionary mapping coordinate system names to their classes
coordinate_systems = {
    "CoordinateBase": CoordinateBase,
    "Cartesian": Cartesian,
    "Spherical": Spherical,
    "BoyerLindquist": BoyerLindquist,
    "Cylindrical": Cylindrical,
    "OrbitalElements": OrbitalElements,
}

__all__ = [
    "CoordinateBase",
    "Cartesian",
    "Spherical",
    "BoyerLindquist",
    "Cylindrical",
    "OrbitalElements",
    "coordinate_systems",
]
