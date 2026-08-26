"""
Symbolic spacetime metrics and coordinate utilities for relatipy.

This subpackage provides SymPy-based symbolic representations of metrics and
coordinate systems. It is intended to complement the numerical routines in
:mod:`relatipy.numeric`, where tensors are evaluated on arrays.

Submodules
----------
metrics
    Symbolic metric tensors (for example, Kerr in Boyer–Lindquist form).
coordinates
    Symbolic coordinate definitions and related helpers.

See Also
--------
relatipy.numeric : Numerical metrics, coordinates, and geodesic integration.

Examples
--------
>>> import relatipy.symbolic as sym
>>> sorted(sym.__all__)
['coordinates', 'metrics']
>>> sym.metrics.__name__
'relatipy.symbolic.metrics'
"""

from . import metrics
from . import coordinates

__all__ = [
    "metrics",
    "coordinates",
]
