"""External forces and potentials acting on a perturbed Kerr geodesic.

Every force is an immutable value with named run-time parameters. Bound to a
black hole, it becomes one native callback evaluated in the black-hole
Cartesian frame in geometric units; no Python code runs inside an
integrator. Use a force with :class:`~relatipy.metrics.PerturbedKerr` or the
perturbed observable models of :mod:`relatipy.observables`.

- :class:`Force` and :class:`Potential` compile a user definition at run
  time (``relatipy[symbolic]`` and a C compiler are required).
- :class:`ThinDisk` is the tabulated power-law thin disk of Piran et al.
  (2026).
"""

from .axisymmetric import ThinDisk
from .compiled import Force, Potential
from .core import ExternalForce

__all__ = ["ExternalForce", "Force", "Potential", "ThinDisk"]
