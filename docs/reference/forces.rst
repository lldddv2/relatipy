.. _forces-reference:

External forces
===============

.. py:module:: relatipy.forces

:mod:`relatipy.forces` defines the external Newtonian forces that perturb a
Kerr geodesic in the perturbed observable models of
:mod:`relatipy.observables`. The library ships one force, the tabulated thin
disk of Piran et al. (2026). See :doc:`/user-guide/forces` for a guide and
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

In the perturbed MCMC models the disk normal is read in observer axes
``(Dec, RA, away)``. Forces are picklable by their constructor arguments;
tables are rebuilt on first use in each process.

API
---

.. autoclass:: relatipy.forces.ExternalForce
   :members: parameters, parameter_names, parameter_units, with_parameters

.. autoclass:: relatipy.forces.ThinDisk
   :members: normal
   :no-show-inheritance:
