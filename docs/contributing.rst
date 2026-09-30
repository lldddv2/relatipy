Contributing
============

RelatiPy separates the public Python API from private Cython bindings and the
native C implementation. Keep physics and numerical integration in the native
layer. Do not expose private extension modules as public API.

Documentation source is reStructuredText in ``docs/``. Public Python APIs use
NumPy-style docstrings, rendered by Sphinx through Napoleon. Describe only
implemented, tested behavior. Keep units, shapes, exceptions, and numerical
limits explicit.

Before submitting documentation changes, run the strict HTML build from the repository root:

.. code-block:: console

   $ uv run --group docs sphinx-build -W --keep-going -n -b html docs docs/_build/html

Warnings fail the build locally and on Read the Docs. The configuration sets
``nitpicky = True``, so every Python cross-reference must resolve, either to a
RelatiPy object or through intersphinx to Python, NumPy, Astropy, or
Matplotlib. New pages must be included in a ``toctree``. Refer to public
objects by their package paths, such as ``relatipy.geodesic.Orbit``, and keep
API names identical to the implementation under ``src/relatipy/``.

Examples in public docstrings are doctests. Run them, together with the test
suite, before submitting a change:

.. code-block:: console

   $ MPLBACKEND=Agg uv run --with pytest pytest --doctest-modules src/relatipy
   $ uv run --with pytest pytest

See :doc:`developer-architecture` for the layer boundaries and the test
layout.

Use Sphinx roles for cross-references, for example
``:py:mod:`relatipy.coordinates``` and ``:doc:`installation```. Do not change
runtime behavior as part of a documentation contribution. Report an ambiguous
or incorrect API contract for a code owner to resolve.
