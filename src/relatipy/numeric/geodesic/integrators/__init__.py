"""
Public API for geodesic time-stepping integrators.

This subpackage re-exports high-order symplectic integrators for the geodesic
equation in Hamiltonian form. The state convention matches the rest of the
geodesic machinery: an 8-vector :math:`[q^0, q^1, q^2, q^3, u^0, u^1, u^2, u^3]`
(position and contravariant four-velocity).

The Yoshida sixth-order scheme composes Störmer--Verlet (drift--kick--drift)
steps with fixed weights to achieve sixth-order accuracy in the affine
parameter while preserving the symplectic structure of Hamiltonian dynamics.

Notes
-----
Only names listed in ``__all__`` are guaranteed as the stable public surface
of this package. Other integrator implementations may exist as submodules
(e.g. specialized Radau solvers) without being imported here.

References
----------
Yoshida, H. (1990). Construction of higher order symplectic integrators.
Physics Letters A, 150(5-7), 262-268.

Examples
--------
>>> from relatipy.numeric.geodesic.integrators import Yoshida6Integrator
>>> from relatipy.numeric.metrics import Schwarzschild
>>> metric = Schwarzschild(5.972e24)
>>> integrator = Yoshida6Integrator(metric)
>>> integrator  # doctest: +ELLIPSIS
<relatipy.numeric.geodesic.integrators.kerr.yoshida6.Yoshida6Integrator object at ...>
"""

from .kerr.yoshida6 import Yoshida6Integrator, yoshida6_integrate_geodesic

__all__ = [
    "Yoshida6Integrator",
    "yoshida6_integrate_geodesic",
]
