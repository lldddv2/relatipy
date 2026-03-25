"""
Geodesic integration for test particles in a given spacetime metric.

This subpackage exposes :class:`~relatipy.numeric.geodesic.geodesic.Geodesic`,
which integrates the geodesic equation using the metric's Christoffel symbols.
It is the public entry point for constructing a geodesic solver bound to a
:class:`~relatipy.numeric.metrics.base.BaseMetric` instance (metrics typically
attach a ``Geodesic`` helper internally as well).

The second-order geodesic equation in affine parameter :math:`\\lambda` is

.. math::

    \\frac{\\mathrm{d}^2 x^\\sigma}{\\mathrm{d}\\lambda^2}
    + \\Gamma^\\sigma_{\\mu\\nu}
    \\frac{\\mathrm{d} x^\\mu}{\\mathrm{d}\\lambda}
    \\frac{\\mathrm{d} x^\\nu}{\\mathrm{d}\\lambda} = 0,

where :math:`\\Gamma^\\sigma_{\\mu\\nu}` are the Christoffel symbols of the
metric.

Notes
-----
Only symbols listed in ``__all__`` are part of the stable public API for this
package. Implementation details live in :mod:`relatipy.numeric.geodesic.geodesic`.

Examples
--------
>>> from relatipy.numeric.geodesic import Geodesic
>>> from relatipy.numeric.metrics import Schwarzschild
>>> geo = Geodesic(Schwarzschild(5.972e24))
>>> geo.metric is not None
True
"""

from .geodesic import Geodesic

__all__ = ["Geodesic"]
