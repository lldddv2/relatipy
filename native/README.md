# Experimental native Kerr reference

This directory contains C references for Kerr geometry, legacy-compatible and
corrected null-geodesic right-hand sides, and native DOP853 and Radau IIA(5)
integrators. A dedicated timelike MCMC evaluator uses these native sources
through the private `relatipy._mcmc_core` Cython extension, built with
`setuptools.build_meta`, to support `KerrMcmcModel`. The general `Kerr` and
`Orbit` API remains experimental and unavailable. These internal C entry
points are not an installed C API or a stable ABI; their use by the specialized
extension does not settle the general architecture or its build choices.

The dedicated MCMC evaluator accepts dimensionless observed arrival times and
observation terms converted to the internal geometric convention. It
integrates, matches arrival times with the Rømer delay, and writes three values
per date: `(alpha_rad, delta_rad, v_los_over_c)`. Its C computation includes
astrometric offsets and linear drifts and the velocity offset. Python prepares
the input conversions and restores public arcsec and km/s on return; it does
not project a stored trajectory afterward.

## Local convention

- Boyer--Lindquist coordinates ordered as `(t, r, theta, phi)`;
- geometrized units `G = c = 1`;
- `M` and `a` are geometric lengths, with `M > 0` and `|a| <= M`;
- metric signature `(-,+,+,+)`;
- caller-owned fixed-size output buffers, with no allocation or retained
  pointers;
- input and output buffers may overlap where each operation explicitly permits
  it because inputs are copied before an output is cleared;
- the rotation axis and roots of `Delta` are reported as coordinate
  singularities; the Kerr ring (`Sigma = 0`) is reported separately as a
  physical singularity.

Every successful call guarantees finite output. Finite inputs whose
intermediate or final values exceed the representable `double` range return
`RP_KERR_STATUS_NUMERICAL_RANGE` and leave a validated output buffer zeroed.

## Experimental interface

The internal headers `include/relatipy/kerr_geometry.h` and
`include/relatipy/kerr_geodesic.h` declare seven reference operations:

| Function | Caller-owned output |
| --- | --- |
| `rp_kerr_metric` | covariant metric `g_(mu nu)` |
| `rp_kerr_inverse_metric` | contravariant metric `g^(mu nu)` |
| `rp_kerr_christoffel` | legacy-compatible tensor `Gamma^lambda_(mu nu)` |
| `rp_kerr_four_velocity` | future-directed four-velocity from `(dr/dt, dtheta/dt, dphi/dt)` |
| `rp_kerr_geodesic_rhs` | legacy-compatible affine derivatives `(dx^mu/dlambda, du^mu/dlambda)` |
| `rp_kerr_null_tangent` | null tangent from `(E, L_z, Q)` and radial/polar signs |
| `rp_kerr_null_geodesic_rhs` | corrected null affine derivatives `(dx^mu/dlambda, dk^mu/dlambda)` |

`RP_KERR_STATUS_OK` is the only success status. Null pointers, non-finite
inputs, invalid parameters or buffer relationships, coordinate singularities,
the physical Kerr ring,
non-timelike coordinate velocities, and numerical range failures have distinct
status values. In particular,
`RP_KERR_STATUS_PHYSICAL_SINGULARITY` reports `Sigma = 0`, while
`RP_KERR_STATUS_NUMERICAL_RANGE` reports an overflow or non-finite derived
value from otherwise finite inputs, and
`RP_KERR_STATUS_NON_TIMELIKE_VELOCITY` reports a non-positive normalization
radicand, while `RP_KERR_STATUS_NO_REAL_NULL_TANGENT` reports a negative radial
or polar separated potential for the requested null constants and point.
Whenever the output pointer itself is valid, every failure leaves the
entire output buffer zero. Input values are copied first, so an input and its
output storage may overlap. The interface does not expose metric derivatives.

These functions evaluate local geometry, perform isolated timelike and null
initial-tangent construction, and evaluate one geodesic derivative. They
do not construct a complete initial state, locate integration events, allocate
solution storage, expose Python objects, or select a production normalization
or canonical integration state. The private integrator module described below
can advance this RHS for validation without changing those limitations.

The internal utilities under `src/utils/` provide rank-independent operations
on contiguous `double` tensor storage and shared scalar predicates. In
particular, `tensor.c` clears tensor buffers and checks all requested
components for finiteness, while `numeric.c` provides the scaled
machine-epsilon check used near singular surfaces. Metric implementations can
reuse them without copying rank-specific loops. They are not installed APIs,
allocate no memory, and leave shape and ownership with the caller.

The Kerr metric implementation is split by responsibility. The principal
`src/metric/kerr.c` file validates and copies caller inputs, enforces output
and error contracts, and dispatches each operation. Physical quantities and
formulas live in `src/metric/physic/kerr.c`: prepared point quantities
(`sin(theta)`, `cos(theta)`, `Sigma`, and `Delta`), the metric, inverse metric,
legacy-compatible Christoffel tensor, and coordinate-velocity normalization.
The adjacent `physic/kerr.h` is a private implementation boundary and is not
an installed header.

The geodesic implementation mirrors the same boundary. The principal
`src/geodesic/kerr.c` file owns validation, input copies, output initialization,
aliasing rules, and all three geodesic operations. The prepared affine state,
separated null-tangent formulas, legacy-compatible optimized `dx`/`du`
formulas, and corrected Levi--Civita contraction live in
`src/geodesic/physic/kerr.c`, behind its own private `kerr.h`.

The physical and geodesic evaluators compute `sin(theta)` and `cos(theta)`
once while preparing a point. Multiple-angle quantities are then obtained
algebraically: `sin(2 theta) = 2 sin(theta) cos(theta)`,
`cos(2 theta) = cos(theta)^2 - sin(theta)^2`, and
`1 - cos(4 theta) = 8 sin(theta)^2 cos(theta)^2`. Likewise,
`cot(theta) = cos(theta) / sin(theta)`, so the hot evaluators make no
additional `sin`, `cos`, or `tan` calls.

The implementation is a readable reference, not the eventual optimized hot
kernel. The metric and inverse metric use the corrected Kerr geometry. The
Christoffel evaluator is a separate compatibility reference: it directly
translates the 20 independent algebraic expressions and 12 lower-index
symmetry copies in the frozen legacy method. The resulting tensor has 32
explicitly populated entries; every component absent from that method remains
zero. It neither evaluates metric derivatives nor reconstructs the tensor from
the metric exposed by this interface.

## Legacy provenance and sign mapping

The formulas were audited against the frozen snapshot at
`../../relatipy-legacy-local-2026-09-01/`, in particular:

- `src/relatipy/numeric/metrics/kerr_metric.py`;
- `src/relatipy/symbolic/metrics/kerr_metric.py`;
- `notebooks/utils/symbolic_to_numeric.py`;
- `tests/numeric/metrics/test_kerr_metric.py`.

The legacy constructor receives a dimensionless spin and converts it to its
internal spin length with `a_internal = a_dimensionless R_s / 2`; its formulas
use `R_s = 2 M`. The effective numeric and symbolic matrices in that snapshot
have signature `(+,-,-,-)`, despite a `(-,+,+,+)` line element in the numeric
module docstring. Therefore the regression test compares the new metric to the
**negative** of the legacy matrix after mapping `R_s = 2 M`. Multiplying a
metric globally by `-1` leaves its Levi-Civita connection unchanged.

`rp_kerr_christoffel` is a literal translation of
`Kerr._get_christoffel_symbols` in that snapshot. It maps `R_s = 2 M`; the C
argument `spin` is the legacy method's internal, length-like `a`, not the
dimensionless constructor argument. The translation preserves all 20
independent expressions and all 12 symmetry assignments exactly.

`rp_kerr_four_velocity` evaluates the semantic C99 candidate derived in
`thesis/docs/v1/004-geodesica_optimizada.ipynb` from the frozen legacy metric.
For an input `v^i = dx^i/dt`, it selects the positive `u^t` branch and returns
`u^i = u^t v^i`. In geometric units, `dr/dt` is dimensionless while
`dtheta/dt` and `dphi/dt` have inverse-length units. The legacy derivation has
signature `(+,-,-,-)` and norm
`+1`; therefore the same vector has norm `-1` under this reference's globally
negated `(-,+,+,+)` metric. The operation requires a positive normalization
radicand and otherwise reports a typed non-timelike-velocity failure.

`rp_kerr_geodesic_rhs` evaluates the optimized semantic candidate derived in
the same notebook. It writes `dx^mu/dlambda = u^mu` and contracts the frozen
legacy-compatible connection directly for
`du^mu/dlambda = -Gamma^mu_(alpha beta) u^alpha u^beta`, without allocating or
materializing the Christoffel tensor. This is a local right-hand-side
reference; it does not adopt `(x, u)` as the production canonical state or
provide an integration method by itself.

`rp_kerr_null_tangent` constructs an affinely scaled tangent from the null
Kerr first integrals. With
`P = E (r^2 + a^2) - a L_z`, its radial and polar potentials are
`R = P^2 - Delta ((L_z - a E)^2 + Q)` and
`Theta = Q - cos(theta)^2 (L_z^2 / sin(theta)^2 - a^2 E^2)`.
The caller supplies the signs of the square roots. Negative potentials return
`RP_KERR_STATUS_NO_REAL_NULL_TANGENT`; roundoff-sized negative values at a
turning point are clamped to zero with the shared scaled-epsilon predicate.
The constants retain the affine scale explicitly: `(E, L_z, Q)` transforms as
`(c E, c L_z, c^2 Q)` under a positive tangent rescaling by `c`.

`rp_kerr_null_geodesic_rhs` deliberately does **not** use the frozen legacy
connection. It computes analytic first derivatives of the corrected metric,
contracts its Levi--Civita connection without materializing a rank-three
Christoffel tensor, and advances `(x^mu, k^mu)`. Despite its historical name,
this contraction also accepts a timelike tangent. Its timelike use is tested
against a connection derived independently from metric derivatives. Null
callers normally construct their initial tangent with `rp_kerr_null_tangent`.
Adaptive Runge--Kutta and collocation stages need not remain exactly on their
initial norm constraint.

## Experimental integrator module

`src/geodesic/integrators/` implements a private registry and a shared C-only
endpoint contract. It registers four adaptive method names:

- `dop853`, the explicit Dormand--Prince 8(5,3) formula;
- `dp45`, the explicit Dormand--Prince 5(4) formula;
- `radau`, the implicit three-stage, fifth-order Radau IIA formula;
- `projection_radau`, the same Radau formula with a native projection callback.

All methods accept the same native RHS callback and opaque C context. The
registry has no default method: callers must select one explicitly. A future
method can be added with a new implementation and registry entry without
changing Kerr physics. The local timelike Kerr adapter now calls the corrected
metric contraction with its experimental `(x,u)` layout; the null adapter
uses the same contraction with `(x,k)`. The approved production state uses
normalized `(x,u)`, but this endpoint prototype still uses `G = c = 1`
with an explicit mass and does not implement that production conversion.

The prototype returns one endpoint, allocates no dynamic memory, retains no
pointer, and uses no mutable global data. Dense output, event localization,
trajectory storage, public options and binding remain outside its contract.
These are current implementation facts, not a stable ABI or concurrency
promise.

The DOP853 coefficients translated from SciPy 1.17.1 live under
`vendor/scipy-dop853/`, together with complete SciPy/DOP license notices,
upstream hashes, local changes and vendored checksums. The Radau source is an
original implementation of the standard published collocation equations; its
numerical reference and the licensing boundary are recorded in the module's
`README.md`. DP45 is original RelatiPy C code using coefficients checked
against SciPy's RK45 source; the source and license boundary are recorded
in `dp45.c` and the module's `README.md`.

This compatibility choice also preserves a known historical discrepancy. The
legacy symbolic/generated Christoffel path contains a spin-sign inconsistency
relative to the corrected numeric metric. Consequently, the tensor returned by
`rp_kerr_christoffel` must not be presented as an independently corrected
Levi-Civita connection of `rp_kerr_metric`. Tests compare all 64 stored tensor
components against frozen outputs from the legacy method at three regular
points. PyGRO is an independent oracle for the complete metric only; it is not
used as a Christoffel oracle. The metric/inverse tests additionally verify
their algebraic identities and the Schwarzschild limit.

## Direct validation

From the `relatipy/` repository root:

```console
mkdir -p /tmp/relatipy-native-tests
cc -std=c99 -Wall -Wextra -Wpedantic -Werror \
  -Inative/include native/src/utils/numeric.c native/src/utils/tensor.c \
  native/src/metric/physic/kerr.c native/src/metric/kerr.c \
  native/tests/unit/metric/test_kerr_geometry.c -lm \
  -o /tmp/relatipy-native-tests/test-kerr-c99
/tmp/relatipy-native-tests/test-kerr-c99

cc -std=c11 -Wall -Wextra -Wpedantic -Werror \
  -Inative/include native/src/utils/numeric.c native/src/utils/tensor.c \
  native/src/metric/physic/kerr.c native/src/metric/kerr.c \
  native/tests/unit/metric/test_kerr_geometry.c -lm \
  -o /tmp/relatipy-native-tests/test-kerr-c11
/tmp/relatipy-native-tests/test-kerr-c11

cc -std=c11 -Wall -Wextra -Wpedantic -Werror -Inative/include -Inative/src \
  native/src/metric/kerr_properties.c \
  native/tests/unit/metric/test_kerr_properties.c -lm \
  -o /tmp/relatipy-native-tests/test-kerr-properties
/tmp/relatipy-native-tests/test-kerr-properties

cc -std=c99 -Wall -Wextra -Wpedantic -Werror \
  -Inative/include native/src/utils/numeric.c native/src/utils/tensor.c \
  native/src/metric/physic/kerr.c native/src/metric/kerr.c \
  native/src/geodesic/physic/kerr.c native/src/geodesic/kerr.c \
  native/tests/unit/geodesic/test_kerr_geodesic.c -lm \
  -o /tmp/relatipy-native-tests/test-kerr-geodesic-c99
/tmp/relatipy-native-tests/test-kerr-geodesic-c99

cc -std=c11 -Wall -Wextra -Wpedantic -Werror \
  -Inative/include native/src/utils/numeric.c native/src/utils/tensor.c \
  native/src/metric/physic/kerr.c native/src/metric/kerr.c \
  native/src/geodesic/physic/kerr.c native/src/geodesic/kerr.c \
  native/tests/unit/geodesic/test_kerr_geodesic.c -lm \
  -o /tmp/relatipy-native-tests/test-kerr-geodesic-c11
/tmp/relatipy-native-tests/test-kerr-geodesic-c11

cc -std=c11 -Wall -Wextra -Wpedantic -Werror \
  -Inative/src native/src/geodesic/integrators/integrator.c \
  native/src/geodesic/integrators/dop853.c \
  native/src/geodesic/integrators/dp45.c \
  native/src/geodesic/integrators/radau.c \
  native/tests/unit/geodesic/test_integrators.c -lm \
  -o /tmp/relatipy-native-tests/test-integrators-c11
/tmp/relatipy-native-tests/test-integrators-c11

cc -std=c11 -Wall -Wextra -Wpedantic -Werror \
  -Inative/include -Inative/src \
  native/src/utils/numeric.c native/src/utils/tensor.c \
  native/src/metric/physic/kerr.c native/src/metric/kerr.c \
  native/src/geodesic/physic/kerr.c native/src/geodesic/kerr.c \
  native/src/geodesic/integrators/integrator.c \
  native/src/geodesic/integrators/dop853.c \
  native/src/geodesic/integrators/dp45.c \
  native/src/geodesic/integrators/radau.c \
  native/src/geodesic/integrators/kerr.c \
  native/tests/unit/geodesic/test_kerr_integrators.c -lm \
  -o /tmp/relatipy-native-tests/test-kerr-integrators-c11
/tmp/relatipy-native-tests/test-kerr-integrators-c11

cc -std=c11 -Wall -Wextra -Wpedantic -Werror \
  -Inative/include -Inative/src \
  native/src/utils/numeric.c native/src/utils/tensor.c \
  native/src/metric/physic/kerr.c native/src/metric/kerr.c \
  native/src/geodesic/physic/kerr.c native/src/geodesic/kerr.c \
  native/src/geodesic/integrators/integrator.c \
  native/src/geodesic/integrators/dop853.c \
  native/src/geodesic/integrators/dp45.c \
  native/src/geodesic/integrators/radau.c \
  native/src/geodesic/integrators/kerr.c \
  native/tests/unit/geodesic/test_kerr_dp45.c -lm \
  -o /tmp/relatipy-native-tests/test-kerr-dp45-c11
/tmp/relatipy-native-tests/test-kerr-dp45-c11

cc -std=c11 -Wall -Wextra -Wpedantic -Werror \
  -Inative/include -Inative/src \
  native/src/utils/numeric.c native/src/utils/tensor.c \
  native/src/metric/physic/kerr.c native/src/metric/kerr.c \
  native/src/geodesic/physic/kerr.c native/src/geodesic/kerr.c \
  native/src/geodesic/integrators/kerr.c \
  native/tests/unit/geodesic/test_timelike_rhs_consistency.c -lm \
  -o /tmp/relatipy-native-tests/test-timelike-rhs-c11
/tmp/relatipy-native-tests/test-timelike-rhs-c11

cc -std=c11 -Wall -Wextra -Wpedantic -Werror \
  -Inative/include -Inative/src \
  native/src/utils/numeric.c native/src/utils/tensor.c \
  native/src/metric/physic/kerr.c native/src/metric/kerr.c \
  native/src/geodesic/physic/kerr.c native/src/geodesic/kerr.c \
  native/src/geodesic/physic/kerr_observables.c \
  native/src/geodesic/initial/convert.c \
  native/src/geodesic/integrators/integrator.c \
  native/src/geodesic/integrators/dop853.c \
  native/src/geodesic/integrators/dp45.c \
  native/src/geodesic/integrators/radau.c \
  native/src/geodesic/kerr_mcmc.c \
  native/tests/unit/geodesic/test_kerr_mcmc.c -lm \
  -o /tmp/relatipy-native-tests/test-kerr-mcmc-c11
/tmp/relatipy-native-tests/test-kerr-mcmc-c11
```

On compilers that provide AddressSanitizer and UndefinedBehaviorSanitizer:

```console
cc -std=c11 -Wall -Wextra -Wpedantic -Werror \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -Inative/include native/src/utils/numeric.c native/src/utils/tensor.c \
  native/src/metric/physic/kerr.c native/src/metric/kerr.c \
  native/tests/unit/metric/test_kerr_geometry.c -lm \
  -o /tmp/relatipy-native-tests/test-kerr-sanitized
ASAN_OPTIONS=detect_leaks=0 /tmp/relatipy-native-tests/test-kerr-sanitized

cc -std=c11 -Wall -Wextra -Wpedantic -Werror \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -Inative/include native/src/utils/numeric.c native/src/utils/tensor.c \
  native/src/metric/physic/kerr.c native/src/metric/kerr.c \
  native/src/geodesic/physic/kerr.c native/src/geodesic/kerr.c \
  native/tests/unit/geodesic/test_kerr_geodesic.c -lm \
  -o /tmp/relatipy-native-tests/test-kerr-geodesic-sanitized
ASAN_OPTIONS=detect_leaks=0 \
  /tmp/relatipy-native-tests/test-kerr-geodesic-sanitized
```

LeakSanitizer is disabled in this command because the reference performs no
dynamic allocation and LeakSanitizer cannot run in ptrace-restricted
environments. AddressSanitizer and UndefinedBehaviorSanitizer remain active.

Passing these commands establishes source compatibility for this isolated
reference under the selected compiler; it does not choose the project's final
C standard, build system, platforms, or ABI.
