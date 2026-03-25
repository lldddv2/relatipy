"""
Black-hole and vacuum spacetime metrics for numerical relativity.

This subpackage provides the public metric classes used by ``relatipy``'s
numeric stack: a shared :class:`~relatipy.numeric.metrics.base.BaseMetric`
interface plus concrete stationary solutions (non-rotating Schwarzschild
and rotating Kerr). Each implementation supplies the covariant metric
:math:`g_{\\mu\\nu}` in its natural chart and shares helpers for Christoffel
symbols and geodesic integration.

The Schwarzschild line element in spherical coordinates
:math:`(t,r,\\theta,\\phi)` (geometric units with Schwarzschild radius
:math:`R_s = 2GM/c^2`) is

.. math::

    ds^2 = \\left(1 - \\frac{R_s}{r}\\right)\\,c^2\\,dt^2
    - \\left(1 - \\frac{R_s}{r}\\right)^{-1} dr^2
    - r^2\\left(d\\theta^2 + \\sin^2\\theta\\, d\\phi^2\\right).

The Kerr implementation uses Boyer--Lindquist coordinates; see
:class:`~relatipy.numeric.metrics.kerr_metric.Kerr` for spin parameter
conventions.

Notes
-----
:class:`~relatipy.numeric.metrics.base.BaseMetric` attaches a
:class:`~relatipy.numeric.geodesic.geodesic.Geodesic` instance for tracing
test-particle trajectories in the chosen spacetime.

Exports are limited to ``__all__`` below; additional metric modules may
exist alongside this package but are not part of the stable public API
unless listed here.

References
----------
Schwarzschild, K. (1916). Über das Gravitationsfeld eines Massenpunktes
nach der Einsteinschen Theorie. *Sitzungsber. Preuss. Akad. Wiss.*
1916, 189--196.

Kerr, R. P. (1963). Gravitational field of a spinning mass as an example
of algebraically special metrics. *Phys. Rev. Lett.* **11**, 237--238.

Examples
--------
>>> import numpy as np
>>> from relatipy.numeric.metrics import Schwarzschild, Kerr, BaseMetric
>>> xs = np.array([0.0, 1e9, np.pi / 2, 0.0])
>>> g = Schwarzschild(5.972e24).metric(xs)
>>> g.shape
(4, 4)
>>> g_kerr = Kerr(5.972e24, 0.5).metric(xs)
>>> g_kerr.shape
(4, 4)
"""

from .base import BaseMetric
from .schwarzschild_metric import Schwarzschild
from .kerr_metric import Kerr

__all__ = [
    "BaseMetric",
    "Schwarzschild",
    "Kerr",
]
