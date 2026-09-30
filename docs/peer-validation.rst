.. _peer-validation:

Scientific validation
=====================

The ``tests/reference/`` suite compares RelatiPy with independent references:
exact Schwarzschild solutions, the analytic Kerr trajectories of KerrGeoPy,
the numerical geodesics of PyGRO, and conserved Kerr quantities computed
outside the native core. Frozen reference data and the recorded reports live
in ``tests/fixtures/``.

.. warning::

   These checks cover the fixed cases and tolerances listed on this page. They
   do not establish accuracy for other spins, orbits, durations, or chart
   regions. Independent scientific review of the references and of the
   recorded reports is still pending.

Conventions
-----------

All reference comparisons use geometric units ``G = c = M = 1``, metric
signature ``(-,+,+,+)``, Boyer--Lindquist component order
``(t, r, theta, phi)``, timelike proper time ``tau``, and the contravariant
four-velocity ``u^mu = dx^mu/dtau``. Compared states have the order
``(t, r, theta, phi, u^t, u^r, u^theta, u^phi)``. Equivalent azimuth branches
separated by an integer multiple of ``2*pi`` are aligned before subtraction.
Error limits are maximum absolute componentwise limits in these units.

KerrGeoPy produces its trajectories in Mino time ``lambda``. The reference
exporter maps them to proper time with

.. math::

   \frac{d\tau}{d\lambda} = r(\lambda)^2 + a^2 \cos^2\!\left(\theta(\lambda)\right),

integrated with SciPy DOP853 at ``rtol=2.3e-14`` and ``atol=1e-13``. PyGRO
integrates the same initial position and four-velocity with its ``dp853``
integrator and an initial step of ``0.05 M``. These settings are recorded in
``tests/fixtures/peer_validation_cases.json``. KerrGeoPy and PyGRO are test
oracles only; neither is a RelatiPy dependency.

Test suites
-----------

.. list-table:: Reference test modules in ``tests/reference/``
   :header-rows: 1
   :widths: 34 66

   * - Module
     - What it checks
   * - ``test_relatipy_peer_contract.py``
     - Public :class:`~relatipy.geodesic.Orbit` endpoints and
       :meth:`~relatipy.geodesic.Solution.at` trajectories against frozen
       exact, KerrGeoPy, and PyGRO states for ``radau``, ``dop853``, and
       ``dp45``. Runs without the peer libraries.
   * - ``test_orbit_invariant_drift.py``
     - Drift of the timelike norm, energy, axial angular momentum, and Carter
       constant over several orbital periods, computed with an independent
       metric implementation on accepted native steps.
   * - ``test_projection_radau.py``
     - Invariant preservation of ``projection_radau`` and its domain and
       horizon-failure contract.
   * - ``test_orbit_domain_extremes.py``
     - Near-axis and near-horizon states, extreme spin, exact Schwarzschild
       radial infall, and rejected initial conditions.
   * - ``test_peer_library_validation.py``
     - The KerrGeoPy and PyGRO oracles themselves against exact circular
       states and each other. Skips when either library is missing.
   * - ``test_native_geometry_peer.py``
     - The native metric against PyGRO and the native connection against a
       frozen compatibility snapshot. Needs a C compiler (``cc``).
   * - ``test_native_integrators_scipy.py``
     - The native DOP853 and Radau integrators against SciPy on the same C
       right-hand side. This checks time stepping only; it is not an
       independent physics oracle. Needs a C compiler (``cc``).

Public orbit against peer references
------------------------------------

The frozen references in ``tests/fixtures/orbit_peer_reference.json`` were
exported with KerrGeoPy 0.9.3 and PyGRO 1.0.3. They contain four cases:

.. list-table:: Reference cases
   :header-rows: 1
   :widths: 30 44 26

   * - Case
     - Domain
     - Primary reference
   * - ``schwarzschild_circular_r10``
     - ``a = 0``, circular, ``r = 10 M``, one period
     - Exact circular state
   * - ``schwarzschild_isco_r6``
     - ``a = 0``, circular, ``r = 6 M``, one period
     - Exact circular state
   * - ``schwarzschild_eccentric``
     - ``a = 0``, ``p = 10 M``, ``e = 0.2``, one radial period
     - KerrGeoPy
   * - ``kerr_stable_bound``
     - ``a = 0.5 M``, ``p = 8 M``, ``e = 0.2``, ``x = 0.8``, one radial
       period
     - KerrGeoPy

KerrGeoPy's stable-orbit constructor rejects the marginally stable
``r = 6 M`` orbit, so that case uses the exact circular state. Each case is
integrated through the public API with ``rtol=1e-11`` and ``atol=1e-13``.
Endpoint checks use the raw integrated endpoint. Trajectory checks limit the
step to ``max_step = 0.2 T0`` and evaluate :meth:`~relatipy.geodesic.Solution.at`
on the reference mesh, so they include post-integration interpolation.

.. list-table:: Componentwise error budgets
   :header-rows: 1
   :widths: 16 10 10 10 10 11 11 11 11

   * - Budget
     - ``t``
     - ``r``
     - ``theta``
     - ``phi``
     - ``u^t``
     - ``u^r``
     - ``u^theta``
     - ``u^phi``
   * - Endpoint
     - ``1e-7``
     - ``1e-8``
     - ``1e-9``
     - ``1e-8``
     - ``1e-9``
     - ``1e-9``
     - ``1e-10``
     - ``1e-10``
   * - Trajectory
     - ``5e-6``
     - ``5e-8``
     - ``5e-8``
     - ``5e-8``
     - ``5e-9``
     - ``5e-9``
     - ``5e-10``
     - ``5e-10``

The recorded report ``tests/fixtures/orbit_peer_reference_results.json``
lists all twelve endpoint and twelve trajectory results (four cases, three
methods) within these budgets and has status ``pass``. It was generated with
RelatiPy 1.0.0, Python 3.11, NumPy 2.4.6, SciPy 1.17.1, and Astropy 8.0.1 on
Linux x86-64. The test module also checks endpoint convergence when the
tolerances are tightened and trajectory convergence when ``max_step`` is
reduced.

Invariant drift
---------------

``tests/fixtures/orbit_invariant_drift_results.json`` records two cases: an
inclined eccentric Kerr orbit (``a = 0.5 M``, ``p = 8 M``, ``e = 0.2``,
``x = 0.8``) over five radial periods, and a circular Schwarzschild orbit at
``r = 10 M`` over five revolutions. Each runs with ``radau``, ``dop853``, and
``dp45`` at a coarse setting (``rtol=1e-6``, ``atol=1e-8``) and a fine
setting (``rtol=1e-10``, ``atol=1e-12``). With the fine setting, the scaled
drift of every invariant and ``|u.u + 1|`` must stay below ``2e-8``; for the
Kerr case the fine-to-coarse drift ratio must be at most ``0.05``. The
recorded status is ``pass``. Conservation alone does not bound phase or
coordinate error.

Domain limits
-------------

``tests/fixtures/orbit_domain_extremes_results.json`` records near-axis,
near-horizon, and radial-infall cases at ``rtol=1e-9`` and ``atol=1e-12``.
Successful cases must keep the norm error below ``2e-9`` and each invariant
drift below ``2e-9 * max(1, |initial value|)``. Exact Schwarzschild infall
must match the exact radius within ``2e-11``. The endpoints of different
methods must agree within ``2e-8``; because every method uses the same native
right-hand side, this checks numerical consistency only. An integration stage
that crosses the polar axis must end as a numerical failure that retains the
initial state, and a case that fails from horizon conditioning must return a
finite exterior partial solution; no accuracy is claimed for either.

Limitations
-----------

* All cases are timelike and, apart from the domain-limit cases, bound,
  non-polar, and non-extremal. No result is asserted outside that domain.
* Trajectory errors include interpolation of the stored samples.
* The ``projection_radau`` method, the native-integrator comparison, and the
  cross-method checks share the RelatiPy right-hand side; they do not replace
  the independent KerrGeoPy and PyGRO comparisons.
* The frozen reports describe the environment in which they were generated.
  Rerun the tests to check a different build or platform.

Running the checks
------------------

From the repository root, with the package built in the environment:

.. code-block:: console

   $ uv run --with pytest pytest tests/reference

The public-orbit, invariant, projection, and domain tests need only RelatiPy
and its runtime dependencies. The oracle tests skip when KerrGeoPy or PyGRO
is missing; when both are installed, KerrGeoPy must be version 0.9.3 and
PyGRO at least 1.0.3. The native tests skip without a C compiler.
