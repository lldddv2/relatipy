"""
Kerrpy numeric module.

This module provides numerical implementations for general relativity calculations,
including coordinate systems, metrics, geodesics, and physical constants.
"""

# Import submódulos
from . import constants
from . import coordinates
from . import metrics
from . import geodesic

# Import constants at module level for convenience
from .constants import _c, _G

__all__ = [
    # Submódulos
    "constants",
    "coordinates",
    "metrics",
    "geodesic",
    # Constants for convenience
    "_c",
    "_G",
]
