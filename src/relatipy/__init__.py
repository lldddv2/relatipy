"""
Relatipy: numeric, symbolic, and visualization tools for relativistic geometry.

This package provides submodules for numerical computations, symbolic
manipulation, and plotting related to black-hole spacetimes and related
physics. The public API exposes the package version and the top-level
subpackages ``numeric``, ``symbolic``, and ``visualization``.

Notes
-----
The version string is available as :data:`relatipy.__version__` and is
re-exported in :data:`__all__` together with the subpackage names.

Examples
--------
>>> import relatipy
>>> hasattr(relatipy, "numeric")
True
>>> from relatipy import numeric, symbolic, visualization
"""

__version__ = "0.3.0"

# Expose subpackages for convenience
from . import numeric
from . import symbolic
from . import visualization

__all__ = ["__version__", "numeric", "symbolic", "visualization"]
