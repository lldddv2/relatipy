"""
Symbolic spacetime metrics for general relativity.

This subpackage exposes SymPy-based metric tensors and related geometric
quantities (for example Christoffel symbols) built on ``einsteinpy.symbolic``.
The public API is limited to the symbols listed in ``__all__``.

Notes
-----
Concrete metric implementations live in sibling modules (for example
``kerr_metric``) and are re-exported here for a stable import path.

Examples
--------
>>> from relatipy.symbolic.metrics import Kerr
>>> k = Kerr()
>>> g = k.metric()
>>> g.tensor().shape == (4, 4)
True
"""

from .kerr_metric import Kerr

__all__ = [
    "Kerr",
]
