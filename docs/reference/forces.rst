.. _forces-reference:

External forces
===============

.. py:module:: relatipy.forces

:mod:`relatipy.forces` defines the external Newtonian forces that perturb a
Kerr geodesic in the perturbed observable models of
:mod:`relatipy.observables`. It provides the tabulated thin disk of Piran et
al. (2026), a tabulated thin disk of any axisymmetric surface density, and
forces and potentials written as expressions and compiled to C at run time.
:func:`get_force` builds the last three from text. See
:doc:`/user-guide/forces` for a guide and
:doc:`/development/external-forces` for the equation of motion and the
native contract.

Every force is an immutable value with named run-time parameters
(:attr:`~ExternalForce.parameter_names`, :attr:`~ExternalForce.parameter_units`).
It acts in the black-hole Cartesian frame, ``z`` along the spin, in geometric
units of the black hole. Each parameter keeps the unit it was given in and is
converted when the force is bound to a black-hole mass: a quantity of length
dimension ``L``, time ``T`` and mass ``M`` is divided by
``r_g**L * t_g**T * M**M``; plain numbers are dimensionless.

.. list-table::
   :header-rows: 1
   :widths: 22 38 40

   * - Force
     - Run-time parameters
     - Native resources
   * - :class:`ThinDisk`
     - ``disk_mass``, ``normal_theta``, ``normal_phi``
     - table per (slope, ``r_in / r_out``, softening), about 1 s and 11 MB
   * - :class:`SurfaceDensityDisk`
     - ``disk_mass``, ``normal_theta``, ``normal_phi``
     - compiled ``Sigma(R)``; table per (density and its parameters,
       ``r_in / r_out``, softening)
   * - :class:`Force`, :class:`Potential`
     - the keys of ``parameters``, in order
     - shared library compiled from the generated C source, cached on disk

In the perturbed MCMC models the disk normal is read in observer axes
``(Dec, RA, away)``. Tabulated disks are picklable by their constructor
arguments; tables are rebuilt on first use in each process. :class:`Force`
and :class:`Potential` pickle their generated C source and parameters.

:class:`SurfaceDensityDisk`, :class:`Force`, :class:`Potential` and
:func:`get_force` need SymPy (``pip install sympy``) to generate code and a
C11 compiler (``$CC`` or ``cc``) to compile it. Compiled libraries are cached
in ``$RELATIPY_CACHE_DIR`` or, when it is unset, ``~/.cache/relatipy/jit``.

API
---

.. autoclass:: relatipy.forces.ExternalForce
   :members: parameters, parameter_names, parameter_units, with_parameters

.. autoclass:: relatipy.forces.ThinDisk
   :members: normal
   :inherited-members: ExternalForce
   :no-show-inheritance:

.. autoclass:: relatipy.forces.SurfaceDensityDisk
   :members: normal
   :inherited-members: ExternalForce
   :no-show-inheritance:

.. autofunction:: relatipy.forces.get_force

.. autoclass:: relatipy.forces.Force
   :members: coords, source
   :inherited-members: ExternalForce
   :no-show-inheritance:

.. autoclass:: relatipy.forces.Potential
   :members: coords, source
   :inherited-members: ExternalForce
   :no-show-inheritance:

Building blocks of callable definitions
---------------------------------------

These functions return SymPy expressions for the callable definitions of
:class:`Force` and :class:`Potential`; they are not available in string
definitions.

.. autofunction:: relatipy.forces.elliptic_k

.. autofunction:: relatipy.forces.elliptic_e

.. autofunction:: relatipy.forces.integral
