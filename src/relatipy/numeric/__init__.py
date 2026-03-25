"""
Numerical routines for general relativity in relatipy.

This subpackage exposes implementations used for coordinate charts, spacetime
metrics, geodesic integration, and physical constants in geometric or SI-related
forms. Submodules are re-exported at this level for a single import path.

Notes
-----
The speed of light :math:`c` and gravitational constant :math:`G` in geometric
units are available as ``_c`` and ``_G`` (both unity in the current setup).
For SI-scale reference values, see ``relatipy.numeric.constants``.

Examples
--------
>>> from relatipy.numeric import metrics, geodesic
>>> from relatipy.numeric import _c, _G
>>> _c
1.0
"""

from . import constants
from . import coordinates
from . import metrics
from . import geodesic

from .constants import _c, _G

__all__ = [
    "constants",
    "coordinates",
    "metrics",
    "geodesic",
    "_c",
    "_G",
]
