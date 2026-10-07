.. _integration-roadmap:

Orbit integration: scope and limits
===================================

Current implementation
----------------------

:class:`~relatipy.metrics.Kerr` and
:class:`~relatipy.geodesic.Orbit` provide the first scalar,
timelike orbit interface. :meth:`~relatipy.geodesic.Orbit.integrate`
updates one orbit's current point;
:meth:`~relatipy.geodesic.Orbit.solve` starts from its saved initial
conditions and returns a read-only :class:`~relatipy.geodesic.Solution`.
The specialized :class:`~relatipy.observables.KerrMcmcModel` remains a separate
observable evaluator with its own interface.

The Python frontend validates inputs and uses Astropy for physical units.
Initial-condition conversion, Kerr geometry, and geodesic integration run in
C through a Cython binding. The integration loop makes no Python callbacks.
The backend uses normalized geometric units ``G = c = M = 1``; Python uses
the black-hole mass to convert length by ``r_g = G M / c²`` and time by
``T0 = G M / c³``. Its eight-component Boyer--Lindquist state is
``(t/T0, R/r_g, Theta, Phi, u^t, u^R, u^Theta, u^Phi)``, with proper time
normalized by ``T0``. Public position and velocity remain physical
quantities. The six-value Cartesian ``state_vector`` is distinct from this
native integration state.

The general methods are:

- ``"radau"`` (default);
- ``"dop853"``;
- ``"dp45"``;
- ``"projection_radau"``.

These method names are case-sensitive. The common
controls are ``method``, ``rtol``, ``atol``,
``first_step``, and ``max_step``. Omitted tolerances are automatic:
``rtol`` becomes ``1e-10`` and ``atol`` becomes ``rtol`` times a
characteristic scale of each native component at the starting state (see
:ref:`integration-tolerances`). An ``atol`` array has shape ``(8,)`` in the
native state order.
The native controller scales local error componentwise; this is no guarantee
of global trajectory error. The native integrators and their standalone
tests are described in :doc:`/development/architecture`. The output solution holds
accepted samples or validated requested samples; no native dense output is
retained. :meth:`relatipy.geodesic.Solution.at` interpolates stored
Cartesian position and coordinate velocity in Python and asks C to
reconstruct the other physical views in one batch.

The saved solution exposes Boyer--Lindquist coordinates and the input
coordinate family immediately. Other coordinate families and osculating
elements are reconstructed in C on first access and then cached. Bound Kerr
elements specify the initial orbit; they do not constitute an additional
immediate output family. Exact ``Solution.at`` queries select saved states;
intermediate queries use Cartesian interpolation before native reconstruction.

``"projection_radau"`` shares the Radau stepper. Before storing an accepted
exterior step, it fixes the position and projects the four-velocity toward
the initial timelike norm ``-1``, energy, axial angular momentum, and Carter
constant. Projection failure is reported as numerical failure and preserves
the last valid point. Projection does not apply to intermediate stages or
``Solution.at`` interpolation. Requested ``tau_eval`` samples between
accepted steps are interpolated and need not preserve the invariants to the
same tolerance. Invariant conservation does not guarantee global trajectory
accuracy.

The spin parameter is dimensionless with ``0 <= spin <= 1``. Retrograde
motion is selected by the initial velocity or inclination. The six input
orbital elements use a right-handed Cartesian frame centered on the black
hole, with ``z`` along the spin axis. The first orbit must start outside the
outer Kerr horizon. See :doc:`/user-guide/usage` and :doc:`/reference/index` for signatures and
examples.

Osculating elements and degeneracies
------------------------------------

``Orbit.orbital_elements()``, ``State.orbital_elements()``, and
``Solution.orbital_elements()`` return classical Kepler elements computed
from the instantaneous Cartesian position and coordinate velocity. A radial
or parabolic osculating conic describes that instantaneous representation;
it does not classify the complete Kerr trajectory.

Exactly zero computed Kepler energy gives ``a = +inf`` in the usual length
unit. No near-parabolic tolerance is introduced: floating-point roundoff can
give a finite ``a`` even for an analytically parabolic input. For zero or
numerically negligible angular momentum, the existing native threshold is
retained and the undefined ``inc``, ``Omega``, ``omega``, and ``f`` are
``NaN`` angular quantities. The conventions for ordinary circular and
equatorial conics remain unchanged.

The ``defined`` property of :class:`~relatipy.coordinates.OrbitalElements`
is a read-only Boolean mask in field order ``(a, e, inc, Omega, omega, f)``.
Its shape is ``(6,)`` for one state and ``(n, 6)`` for a state series.
``+inf`` semi-major axis is a defined extended value and has mask value
``True``; an undefined ``NaN`` field has mask value ``False``. Python
derives this mask from the stored elements; the native 29-column state
layout is unchanged.

These auxiliary values do not reject an otherwise valid physical state.
Construction, integration, copying, solution selection, and
``Solution.at`` preserve the representation. Physical coordinates,
coordinate velocities, and four-velocities must remain finite. Each
interpolated state has its own reconstructed elements; interpolation can
change the instantaneous conic classification. Classical Kepler elements
cannot be supplied through a semi-latus rectum instead of ``a``, so
parabolic initial conditions must use a supported position/velocity family.
The bound Kerr ``p`` input of :class:`~relatipy.coordinates.KerrOrbitalElements`
is a different parameterization for stable bound geodesics.

Plotting scope
--------------

``Orbit.preview`` uses the current osculating Kepler elements to draw a
classical conic. It marks the current point and can draw the outer horizon
and both equatorial ISCO radii as reference circles in the oblate Cartesian
``xy`` plane. These circles are not the full Kerr surfaces, and the preview
is not an integrated trajectory. ``Solution.plot`` draws the actual stored
Cartesian samples returned by the numerical integration; it does not
resample or interpolate them. Both operations accept:

- ``"views"``, the default three-panel static figure;
- ``"3d"``;
- ``"xy"``, ``"xz"``, and ``"yz"``.

Static figures use Matplotlib; ``"3d"`` defaults to an interactive Plotly
figure. Both libraries are installed with RelatiPy.
See :doc:`/reference/plotting` for the return forms.
``Orbit.preview`` raises ``ValueError`` if required element angles are
undefined or the osculating conic is parabolic. ``Solution.plot`` can still
draw the finite stored Cartesian trajectory.

Termination and failure
-----------------------

``Solution.status`` follows the ``-1, 0, 1`` convention:

- ``-1``: numerical failure;
- ``0``: requested endpoint reached;
- ``1``: internal terminal event.

A positive terminal
status has structured ``Solution.termination``. For a confirmed crossing
of the outer horizon by an accepted native step, the reported termination
state is the last valid point outside the horizon, and
``termination.tau == termination.state.tau``. It is not the exact crossing
time. If an intermediate method stage fails before an accepted crossing can
confirm the event, the solver reports numerical failure instead. Detection
tests ``R <= r_+`` after each accepted native step; crossing localization
and an additional event tolerance are not part of this contract.

``Orbit.integrate`` raises
:class:`~relatipy.geodesic.IntegrationTerminated` on a terminal event
and :class:`~relatipy.geodesic.IntegrationError` on numerical failure,
retaining its last valid point. ``Orbit.solve`` returns a partial solution
with ``status == -1`` on an integration failure when its stored states can
be reconstructed. Invalid arguments and units raise Python validation errors
before native integration.

Both operations issue :class:`~relatipy.geodesic.IntegrationWarning` before
integration in these cases:

- an explicit ``rtol`` or ``atol`` is unfit for the starting state; omitted
  tolerances never trigger it;
- for ``Orbit.solve`` only, adjacent requested proper-time samples are
  farther apart than one tenth of the initial osculating Kepler period,
  provided that the initial elements describe a finite bound conic.

This period is a Newtonian sampling estimate, not an error bound for a Kerr
trajectory. Warnings leave the return status and exceptions unchanged.

Current domain limits
---------------------

* The present Boyer--Lindquist chart has a polar-axis singularity. Its domain
  policy treats ``Theta`` outside ``(0, pi)`` as invalid and requires an error
  when ``abs(sin(Theta)) <= 64 * DBL_EPSILON`` (about ``1.42e-14`` for IEEE
  754 double precision). This numerical guard does not establish a physical
  axis crossing. An internal-stage error before an accepted step is numerical
  failure, not a terminal event. No alternate chart or axis continuation is
  approved. Full validation of this guard across paths remains pending;
  An integration stage may fail before the threshold is reached.
* Horizon detection deliberately does not localize a crossing. A failed stage
  can leave the final state outside the horizon with a numerical failure status.
* The complete physical-domain and invariant validation of the corrected
  timelike ``(x, u)`` right-hand side remains pending. Short native tests
  cover particular trajectories, not all chart domains or long-time drift.
* General ``Orbit`` calls have no supported batch constructor, Python event
  callbacks, dense-output object, method-specific tuning keywords, or
  observation/fitting interface.

Implementation and distribution decisions
-----------------------------------------

The native integrators are reentrant at the C level and use independent
configuration, state, and counters. The general binding still requires
separate validation before promising Python-level concurrent calls with the
GIL released. C output buffers are copied once to NumPy-owned storage at the
binding and freed; public state arrays are read-only. The precise native object layout and
lifetime design beyond these operations remains an open decision.

The repository currently builds its private Cython extension with
``setuptools.build_meta``. The following remain proposals or pending
decisions; they are not distribution guarantees:

- a migration to scikit-build-core and CMake;
- C11 as a required minimum;
- supported platforms and wheel matrix;
- generated analytic-expression strategy.

The scientific reference checks and their
limits are described in :doc:`/development/peer-validation`.
