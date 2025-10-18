"""Top-level package for relatipy (KerrPy rename).

Exports package version and convenience imports if needed.
"""

__version__ = "0.0.0"

# Expose subpackages for convenience
from . import numeric
from . import symbolic
from . import visualization

__all__ = ["__version__", "numeric", "symbolic", "visualization"]
