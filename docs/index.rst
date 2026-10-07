RelatiPy
========

RelatiPy provides numerical tools for relativistic geometry: timelike
geodesics in the Kerr spacetime, with a Python interface built on Astropy
units and a native C core. Python validates inputs, handles physical units,
and returns immutable result objects; the physics and the numerical
integration run in C without Python callbacks.

Install it from PyPI:

.. code-block:: console

   $ python -m pip install relatipy

Integrate one orbit around a spinning black hole:

.. code-block:: python

   from astropy import units as u
   from astropy.constants import c
   from relatipy import Kerr

   bh = Kerr(mass=1 * u.Msun, spin=0.5)
   orb = bh.orbit(x=12 * bh.r_g, vy=0.1 * c)
   sol = orb.solve(tau_span=(0 * u.s, 2e-8 * u.s), method="dp45")

   mid = sol.at(tau=sol.tau[-1] / 2)   # interpolated state, no reintegration
   elements = mid.orbital_elements()   # instantaneous Kepler conic

Read :doc:`/user-guide/scope-and-limits` before relying on results.

Where to go next
----------------

:doc:`/tutorials/index`
   Executable notebooks, each with an *Open in Colab* button: first orbit,
   initial conditions, integration controls, solutions, plots, and MCMC
   observables.

:doc:`/user-guide/index`
   Installation, a narrative guide to orbits and solutions, observables for
   MCMC, and the current numerical scope.

:doc:`/reference/index`
   Signatures, units, shapes, and exceptions of the public modules.

:doc:`/development/index`
   Architecture, contribution workflow, and scientific validation.

.. toctree::
   :hidden:
   :maxdepth: 2

   user-guide/index
   tutorials/index
   reference/index
   development/index

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
