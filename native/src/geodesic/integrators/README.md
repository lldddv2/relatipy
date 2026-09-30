# Native geodesic integrators

This private module registers four adaptive method names:

- **dop853**, the explicit Dormand--Prince 8(5,3) method for non-stiff
  problems;
- **dp45**, the explicit Dormand--Prince embedded 5(4) method;
- **radau**, the three-stage, fifth-order Radau IIA method for stiff problems;
- **projection_radau**, the same Radau solver with a required native
  candidate projector before each accepted step is committed.

`integrator.c` owns validation and registry dispatch.  Each method implements
the same C-only endpoint contract in `integrator.h`.  Adding a method requires
a registry entry, convergence and failure tests, and updated documentation.
The projection callback uses the same native context as the RHS and does not
copy Kerr physics into the Radau solver.

The current contract is deliberately internal.  It does not select a public
default, canonical state, event policy, or dense-output policy.  Each solver
uses scalar `absolute_tolerance` unless `absolute_tolerances` points to one
value per state component.  Error scales use that component's absolute
tolerance plus `relative_tolerance` times the state magnitude.  An optional
native observer receives each newly accepted state, excluding the initial
state.  It may stop at that state or report an error; both cases keep the state
and statistics at the last accepted point.  This observer does not define
horizon localization or polar-axis policy.  `kerr.c` adapts the experimental
packed
`(x, u)` state to the corrected Kerr contraction for native validation.
The approved production state uses normalized `(x, u)`; the current endpoint
adapter does not yet implement that normalization.

## Numerical and licensing provenance

The DOP853 tableau and embedded error vectors were translated from SciPy
1.17.1.  The coefficients, full applicable license notices, upstream hashes,
local modifications, and vendored checksums live in
`native/vendor/scipy-dop853/`.

The Radau implementation is original RelatiPy C code built from the standard
three-stage Radau IIA Butcher tableau and collocation equations.  Its numerical
reference is E. Hairer and G. Wanner, *Solving Ordinary Differential Equations
II: Stiff and Differential-Algebraic Problems*, second edition, Springer,
Section IV.8.  No third-party Radau source file is copied into this module, so
there is no additional source-code license to vendor.

DP45 is original RelatiPy C code using the Dormand--Prince 5(4) coefficients
checked against SciPy's RK45 source. The original paper citation, exact SciPy
source location, and applicable existing license notice are listed in
`dp45.c`. No SciPy implementation is copied into this module.

The prototype allocates no dynamic memory, uses no mutable global state, and
retains no callback or context pointer.  These are implementation properties,
not yet a promise of a stable public ABI.
