"""
Coordinate representations for the relatipy numeric stack.

This subpackage provides coordinate classes for flat and curved spacetime
settings used in general relativity and orbital mechanics: Cartesian,
spherical, cylindrical, Boyer-Lindquist (Kerr), and Keplerian orbital elements.

In standard spherical coordinates, a spatial position is written as

.. math::

    (x, y, z) = (r \\sin\\theta \\cos\\phi,\\;
                 r \\sin\\theta \\sin\\phi,\\;
                 r \\cos\\theta),

with polar angle :math:`\\theta` and azimuth :math:`\\phi`.

Attributes
----------
coordinate_systems : dict of str to type
    Registry mapping human-readable class names to coordinate classes
    (``"Cartesian"`` → :class:`~relatipy.numeric.coordinates.cartesian.Cartesian`,
    etc.) for discovery and dynamic construction.

Notes
-----
The implementation module for cylindrical coordinates is named
``cilindrical`` (historical spelling); the public class is
:class:`~relatipy.numeric.coordinates.cilindrical.Cylindrical`.

See Also
--------
relatipy.numeric.coordinates.base.CoordinateBase :
    Abstract base for coordinate implementations.

Examples
--------
>>> from relatipy.numeric.coordinates import Cartesian, coordinate_systems
>>> "Cartesian" in coordinate_systems
True
>>> coordinate_systems["Cartesian"] is Cartesian
True
"""

from .base import CoordinateBase
from .cartesian import Cartesian
from .spherical import Spherical
from .boyer_lindquist import BoyerLindquist
from .cilindrical import Cylindrical
from .orbital_elements import OrbitalElements
from .proper_orbital_elements import ProperOrbitalElements

coordinate_systems = {
    "CoordinateBase": CoordinateBase,
    "Cartesian": Cartesian,
    "Spherical": Spherical,
    "BoyerLindquist": BoyerLindquist,
    "Cylindrical": Cylindrical,
    "OrbitalElements": OrbitalElements,
    "ProperOrbitalElements": ProperOrbitalElements,
}

__all__ = [
    "CoordinateBase",
    "Cartesian",
    "Spherical",
    "BoyerLindquist",
    "Cylindrical",
    "OrbitalElements",
    "ProperOrbitalElements"
    "coordinate_systems",
]
