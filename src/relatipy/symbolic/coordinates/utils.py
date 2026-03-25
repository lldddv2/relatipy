"""
Shared imports for symbolic coordinate helpers.

This module centralizes common dependencies used by :mod:`relatipy.symbolic.coordinates`
submodules: SymPy for symbolic algebra, :func:`itertools.product` for Cartesian
products of index ranges, and NumPy for array-oriented interoperability where
needed.

At present the module exposes only these third-party objects; add functions or
constants here when shared coordinate utilities are factored out of other
modules.

Notes
-----
Sibling modules such as :mod:`relatipy.symbolic.coordinates.base` define chart
classes and symbolic variables; this file is the natural place for small
cross-chart helpers that do not belong on a single class.

Examples
--------
>>> from relatipy.symbolic.coordinates import utils
>>> x = utils.sp.Symbol("x")
>>> x
x
>>> list(utils.product([0, 1], repeat=2))
[(0, 0), (0, 1), (1, 0), (1, 1)]
"""

import sympy as sp
from itertools import product
import numpy
