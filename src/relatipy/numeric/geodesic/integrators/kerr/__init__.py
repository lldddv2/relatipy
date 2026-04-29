"""
Kerr-metric geodesic integrators (Boyer-Lindquist coordinates).

This subpackage collects all integrators specialized for the Kerr metric.
All integrators accept the standard 8-state convention
:math:`[q^0, q^1, q^2, q^3, u^0, u^1, u^2, u^3]` (position + contravariant
four-velocity in Boyer-Lindquist coordinates) and return the same format.

Available integrators
---------------------
Yoshida6 (symplectic, 6th-order)
    Long-baseline integrations; uses a C backend when available.
Radau2 (implicit Runge-Kutta, order 5)
    Fixed-step Radau IIA with Kerr constraint projection.
Mino (explicit RK4 in Mino time)
    Decouples radial and polar motion via the Mino time parameter
    :math:`\\mathrm{d}\\tau/\\mathrm{d}\\lambda = \\Sigma(r,\\theta)`.

Examples
--------
>>> from relatipy.numeric.geodesic.integrators.kerr import Yoshida6Integrator
>>> from relatipy.numeric.geodesic.integrators.kerr import _integrate_kerr_radau2
>>> from relatipy.numeric.geodesic.integrators.kerr import _integrate_kerr_mino
"""

from .yoshida6 import Yoshida6Integrator, yoshida6_integrate_geodesic
from .radau import _integrate_kerr_radau2, project_kerr_trajectory
from .mino import _integrate_kerr_mino

__all__ = [
    "Yoshida6Integrator",
    "yoshida6_integrate_geodesic",
    "_integrate_kerr_radau2",
    "project_kerr_trajectory",
    "_integrate_kerr_mino",
]
