.. _developer-architecture:

Architecture and development workflow
=====================================

RelatiPy keeps its public API in Python, uses Cython only as a private
binding, and places Kerr physics and numerical integration in C. Public values
are validated and expressed with Astropy quantities in Python. The native
layer uses geometric units and returns status codes; it does not depend on
Python or Astropy. A binding maps validated contiguous arrays and native status
codes to Python results and exceptions. No Python callback runs inside an
integrator.

Private extensions
------------------

The package builds two private extension modules:

* ``relatipy._core`` supports native trajectory solving, state
  reconstruction, and the Kerr radii and surface profiles used by
  :class:`~relatipy.metrics.Kerr`, :class:`~relatipy.geodesic.Orbit`, and
  :class:`~relatipy.geodesic.Solution`.
* ``relatipy._mcmc_core`` evaluates the dedicated observable model of
  :class:`~relatipy.observables.KerrMcmcModel`.

Neither extension is public API or a stable C ABI. Arrays copied from native
output have their own NumPy storage and are made read-only before they reach
public objects. A native buffer is never exposed after it is freed.

Source layout
-------------

``native/include/relatipy/`` contains internal C headers for local Kerr
geometry and geodesic right-hand sides. ``native/src/`` separates:

- shared numeric utilities (``utils/``);
- metric routines (``metric/``);
- initial-state conversion (``geodesic/initial/``);
- solution reconstruction (``geodesic/solution/``);
- adaptive integrators (``geodesic/integrators/``).

The native integrator names are ``dop853``, ``dp45``, ``radau``, and
``projection_radau``. None of these headers is installed with the package.

``bindings/cython/_core.pyx`` and ``bindings/cython/_mcmc_core.pyx``:

- declare the private C interfaces;
- validate boundary arrays;
- release the GIL only around native calls;
- construct NumPy results.

They contain no Kerr formulas and no per-step Python callbacks. Generated
``.c`` files beside the Cython sources are build artifacts.

``setup.py`` is the active extension definition. The build backend is
``setuptools.build_meta`` with Cython and NumPy as build requirements. Both
extensions:

- are compiled from the ``.pyx`` and C sources listed in ``setup.py``;
- include the NumPy and native headers;
- link ``libm``;
- pass ``-std=c11`` and ``-O3``.

The release workflow in
``.github/workflows/publish.yml`` runs ``uv build --no-sources``. A migration
to CMake and scikit-build-core has been proposed but is not the current build
path.

Native conventions
------------------

The native geometry uses:

- Boyer--Lindquist component order ``(t, r, theta, phi)``;
- metric signature ``(-,+,+,+)``;
- geometric units ``G = c = 1``.

The solver normalizes the mass to one, so lengths are in
``G M / c^2`` and times in ``G M / c^3``; the Python frontend converts these
to and from physical quantities. Local geometry operations write into
caller-owned fixed-size buffers and return typed status codes. The header
comments state the ownership, buffer sizes, and status codes of each native
entry point.

Testing workflow
----------------

``pyproject.toml`` configures pytest (``testpaths = ["tests"]`` and
``pythonpath = ["tests"]``), but no dependency group installs pytest. Run the
suite from the repository root with a pytest available in the environment,
for example:

.. code-block:: console

   $ uv run --with pytest pytest

Python tests are organized by purpose, as described in ``tests/README.md``:

* ``tests/api/`` checks public contracts, grouped like ``src/relatipy/``:
  coordinates, the Kerr metric, states, orbit solving, solution lookup,
  observables, and plotting.
* ``tests/internal/`` checks private helpers and direct calls to the private
  extensions.
* ``tests/reference/`` checks native computations against independent
  references, invariants, and exact solutions; see :doc:`/development/peer-validation`.
* ``tests/fixtures/`` holds frozen JSON references and validation reports, and
  ``tests/support/`` holds shared test builders.

Tests that need optional peer libraries skip when those libraries are not
installed. Native C tests live under ``native/tests/`` as standalone C
programs grouped under ``unit/``, ``initial/``, and ``integrators/``. Build
each one with the source list and compiler flags given in
``native/README.md``; the source dependencies differ by test. That file also
lists optional AddressSanitizer and UndefinedBehaviorSanitizer commands.

A native, binding, or public-contract change needs tests at the changed layer
and at the exposed Python boundary. Report a failing test or build as a
failure.

Documentation uses reStructuredText and NumPy-style public docstrings. Build
it with warnings treated as errors and nitpicky cross-reference checking:

.. code-block:: console

   $ uv run --group docs sphinx-build -W --keep-going -n -b html docs docs/_build/html

Change boundaries
-----------------

Each layer keeps its own responsibilities:

- C keeps physical equations, native state transitions, integration, events,
  and native memory ownership.
- Python keeps public signatures, units, object creation, and
  post-integration interpolation.
- Cython is limited to adapting the two layers.

A runtime change must preserve this boundary and add tests proportional to
its public behavior.

The following are open decisions:

- supported platforms;
- a required minimum C standard;
- a build-system migration.

Do not infer them from the current compiler flags, the native test commands,
or the private extension names.
