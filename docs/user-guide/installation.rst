Installation
============

RelatiPy requires Python 3.11 or later. Install it from PyPI:

.. code-block:: console

   $ python -m pip install relatipy

This also installs its runtime dependencies:

- NumPy;
- Astropy;
- SciPy;
- Matplotlib, for static figures;
- Plotly with nbformat, for interactive figures, which are the default for
  ``projection="3d"`` in :meth:`~relatipy.geodesic.Orbit.preview` and
  :meth:`~relatipy.geodesic.Solution.plot`.

Prebuilt wheels are published for Linux x86_64 and CPython 3.11 to 3.14. On
other platforms pip builds the package from the source distribution, which
requires a C11 compiler and the standard C math library.

Samplers
--------

:class:`~relatipy.observables.KerrMcmcModel` returns model predictions
(astrometry and line-of-sight velocity); the likelihood and the sampler are
left to the user, so RelatiPy does not depend on a specific sampler. The
example notebooks use `emcee <https://emcee.readthedocs.io/>`_:

.. code-block:: console

   $ python -m pip install emcee

Development installation
------------------------

Contributors working from a checkout of the repository use
`uv <https://docs.astral.sh/uv/>`_; see :doc:`/development/contributing`.
