Contributing
============

RelatiPy separates the public Python API from private Cython bindings and the
native C implementation. Keep physics and numerical integration in the native
layer. Do not expose private extension modules as public API.

Development environment
-----------------------

Clone the repository and create the environment with
`uv <https://docs.astral.sh/uv/>`_. Building from source requires a C11
compiler and the standard C math library.

.. code-block:: console

   $ git clone https://github.com/lldddv2/relatipy.git
   $ cd relatipy
   $ uv sync --group dev

Add ``--group notebook`` for the example notebooks (pandas, ipykernel,
nbclient, emcee) and ``--group docs`` for the documentation build.

Documentation
-------------

Documentation source is reStructuredText in ``docs/``. Public Python APIs use
NumPy-style docstrings, rendered by Sphinx through Napoleon. Describe only
implemented, tested behavior. Keep units, shapes, exceptions, and numerical
limits explicit.

Before submitting documentation changes, run the strict HTML build from the repository root:

.. code-block:: console

   $ uv run --group docs sphinx-build -W --keep-going -n -b html docs docs/_build/html

Warnings fail the build locally and on Read the Docs. Documentation changes
must meet these requirements:

- The configuration sets ``nitpicky = True``, so every Python cross-reference
  must resolve, either to a RelatiPy object or through intersphinx to Python,
  NumPy, Astropy, or Matplotlib.
- New pages must be included in a ``toctree``.
- Public objects are referred to by their package paths, such as
  ``relatipy.geodesic.Orbit``.
- API names stay identical to the implementation under ``src/relatipy/``.

Open
``docs/_build/html/index.html`` in a browser to inspect the built pages.

Use Sphinx roles for cross-references, for example
``:py:mod:`relatipy.coordinates``` and ``:doc:`/user-guide/installation```.
Do not change runtime behavior as part of a documentation contribution.
Report an ambiguous or incorrect API contract for a code owner to resolve.

Tutorial notebooks
------------------

Tutorials live in ``docs/tutorials/`` and are rendered by MyST-NB. The
documentation build does not execute them, so commit each notebook with its
outputs. A tutorial notebook must:

- start with a Markdown cell holding the title and an *Open in Colab* badge
  that points to its own path on the ``main`` branch, for example
  ``https://colab.research.google.com/github/lldddv2/relatipy/blob/main/docs/tutorials/01-getting-started.ipynb``;
- install the latest RelatiPy in its second cell with
  ``%pip install -Uq relatipy``;
- use static Matplotlib figures, run in about a minute, and stay small;
- be listed in ``docs/tutorials/index.rst``.

Re-execute notebooks after changing them with the helper script, which skips
the install cell so the development installation is not replaced by the PyPI
release:

.. code-block:: console

   $ uv run --group notebook python docs/tools/execute_notebooks.py docs/tutorials/*.ipynb

Tests
-----

Examples in public docstrings are doctests. Run them, together with the test
suite, before submitting a change:

.. code-block:: console

   $ MPLBACKEND=Agg uv run --with pytest pytest --doctest-modules src/relatipy
   $ uv run --with pytest pytest

See :doc:`/development/architecture` for the layer boundaries and the test
layout.
