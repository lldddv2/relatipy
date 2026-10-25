API reference
=============

The public Python interface is grouped by responsibility. The
:doc:`/tutorials/index` show these objects in use; the pages below give
signatures, units, array shapes, and raised exceptions.

Main classes
------------

.. autosummary::

   relatipy.metrics.Kerr
   relatipy.geodesic.Orbit
   relatipy.geodesic.Solution
   relatipy.geodesic.State
   relatipy.coordinates.KerrOrbitalElements
   relatipy.observables.KerrMcmcModel
   relatipy.observables.PerturbedKerrMcmcMultiModel
   relatipy.forces.ThinDisk

Modules
-------

.. list-table::
   :widths: 30 70

   * - :doc:`metrics`
     - Kerr metric: mass, spin, horizons, and characteristic radii.
   * - :doc:`forces`
     - External forces on a Kerr geodesic: the tabulated thin disk.
   * - :doc:`geodesic`
     - Orbits, integration controls, states, solutions, and integration outcomes.
   * - :doc:`coordinates`
     - Coordinate, velocity, and orbital-element value objects for initial conditions.
   * - :doc:`observables`
     - Astrometry and line-of-sight velocity predictions for sampler workflows,
       for one or several stars, with or without an external force.
   * - :doc:`plotting`
     - Static and interactive figures from stored solutions and observables.

.. toctree::
   :hidden:
   :maxdepth: 2

   metrics
   forces
   geodesic
   coordinates
   observables
   plotting

The native C interfaces under ``native/`` and the Cython modules under
``bindings/`` are implementation details. See :doc:`/user-guide/scope-and-limits`
for the implemented numerical scope and current limitations.
