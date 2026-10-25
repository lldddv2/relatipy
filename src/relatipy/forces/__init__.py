"""External forces and potentials acting on a perturbed Kerr geodesic.

Every force is an immutable value with named run-time parameters. Bound to a
black hole, it becomes one native callback evaluated in the black-hole
Cartesian frame in geometric units; no Python code runs inside an
integrator. Use a force with :class:`~relatipy.metrics.PerturbedKerr` or the
perturbed observable models of :mod:`relatipy.observables`.

- :class:`Force` and :class:`Potential` compile a user definition at run
  time (``relatipy[symbolic]`` and a C compiler are required).
- :class:`Plummer` is a built-in Plummer sphere.
- :class:`ThinDisk` and :class:`Spheroid` are the tabulated power-law disk
  and spheroid of Piran et al. (2026).
"""

from .axisymmetric import Spheroid, ThinDisk
from .compiled import Force, Potential
from .core import ExternalForce
from .plummer import Plummer

__all__ = ["ExternalForce", "Force", "Plummer", "Potential", "Spheroid", "ThinDisk"]
