Installation
============

Requirements
------------

RelatiPy requires Python 3.11 or later, a C11 compiler, and the standard C
math library. Source installation builds the private Cython extensions used by
the public Python interface.

Install from source
-------------------

Clone the repository and install the package from its root directory:

.. code-block:: console

   $ git clone https://github.com/lldddv2/relatipy.git
   $ cd relatipy
   $ python -m pip install .

Install optional dependencies with extras:

.. code-block:: console

   $ python -m pip install ".[plot]"
   $ python -m pip install ".[interactive]"
   $ python -m pip install ".[notebook]"

``plot`` installs Matplotlib for static figures. ``interactive`` installs
Plotly and nbformat for interactive figures, which are the default for
``projection="3d"`` in :meth:`~relatipy.geodesic.Orbit.preview` and
:meth:`~relatipy.geodesic.Solution.plot`. ``notebook`` installs Matplotlib,
Plotly, pandas, the Jupyter kernel and client packages, and emcee.

The ``mcmc`` extra installs ``emcee`` for sampler workflows. The
:class:`~relatipy.observables.KerrMcmcModel` is included in the base package.

For an editable development installation, use
`uv <https://docs.astral.sh/uv/>`_:

.. code-block:: console

   $ uv sync --group dev

Build the documentation
-----------------------

Install the documentation dependency group and run a strict Sphinx build:

.. code-block:: console

   $ uv sync --group docs
   $ uv run --group docs sphinx-build -W --keep-going -n -b html docs docs/_build/html

Open ``docs/_build/html/index.html`` in a browser to inspect the result.
