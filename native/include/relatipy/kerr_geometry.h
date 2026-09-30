/**
 * @file kerr_geometry.h
 * @brief Experimental, allocation-free Kerr geometry reference.
 *
 * This internal reference uses Boyer--Lindquist coordinates ordered as
 * `(t, r, theta, phi)`, geometrized units `G = c = 1`, mass and spin-length
 * parameters `M` and `a`, and metric signature `(-,+,+,+)`.
 *
 * @warning This header is not a public RelatiPy API or stable ABI.  The
 * normalization, installed C interface, build system, and binding remain
 * architectural decisions outside this reference implementation.
 *
 * @note A validated output buffer is zeroed on every error path.  Input
 * coordinates and coordinate velocities may overlap output storage because
 * they are copied before the output is cleared.  A null output pointer cannot
 * be cleared.
 *
 * @note Compared with the frozen legacy RelatiPy snapshot, this reference
 * uses the opposite metric signature.  Its metric is therefore the negative
 * of the legacy numeric metric after converting legacy spin to a geometric
 * spin length.  The Christoffel evaluator intentionally reproduces the
 * explicit generated expressions in the legacy implementation rather than
 * deriving them from the metric exposed here.
 */

#ifndef RELATIPY_KERR_GEOMETRY_H
#define RELATIPY_KERR_GEOMETRY_H

/** Number of spacetime dimensions used by the reference implementation. */
#define RP_KERR_DIM 4

/** Status returned by every Kerr reference operation. */
typedef enum rp_kerr_status {
    /** The operation completed and the output buffer contains a result. */
    RP_KERR_STATUS_OK = 0,
    /** An input or output pointer was null. */
    RP_KERR_STATUS_NULL_POINTER = 1,
    /** At least one scalar input was NaN or infinite. */
    RP_KERR_STATUS_NONFINITE_INPUT = 2,
    /** Scalar parameters or buffer relationships violate an operation contract. */
    RP_KERR_STATUS_INVALID_PARAMETER = 3,
    /** The point lies on, or is numerically indistinguishable from, a chart singularity. */
    RP_KERR_STATUS_COORDINATE_SINGULARITY = 4,
    /** The point lies on, or is numerically indistinguishable from, the Kerr ring. */
    RP_KERR_STATUS_PHYSICAL_SINGULARITY = 5,
    /** Finite inputs produced an intermediate or result outside the double range. */
    RP_KERR_STATUS_NUMERICAL_RANGE = 6,
    /** Coordinate velocities do not define a timelike four-velocity. */
    RP_KERR_STATUS_NON_TIMELIKE_VELOCITY = 7,
    /** Constants of motion do not define a real null tangent at this point. */
    RP_KERR_STATUS_NO_REAL_NULL_TANGENT = 8
} rp_kerr_status;

/**
 * Evaluate the covariant Kerr metric `g_(mu nu)`.
 *
 * @param mass Geometric mass `M`, with `M > 0`.
 * @param spin Spin length `a`, with `|a| <= M`.
 * @param coordinates Caller-owned input `(t, r, theta, phi)` in geometric
 *        units; angles are in radians.
 * @param metric Caller-owned `4 x 4` output buffer.  On failure after this
 *        pointer has been validated, every component is set to zero.
 * @return A typed status code.  `RP_KERR_STATUS_OK` is the only success code.
 *
 * @note The function allocates no memory and retains no pointer.
 * @note The Boyer--Lindquist axis and horizons (`Delta = 0`) are rejected as
 *       chart singularities.  The Kerr ring (`Sigma = 0`) is reported as a
 *       physical singularity.
 * @note Input and output storage may overlap: the four coordinates are copied
 *       before the output buffer is cleared or written.
 *
 * @par Example
 * @code
 * const double x[RP_KERR_DIM] = {0.0, 8.0, 1.1, 0.3};
 * double g[RP_KERR_DIM][RP_KERR_DIM];
 * rp_kerr_status status = rp_kerr_metric(1.0, 0.5, x, g);
 * if (status != RP_KERR_STATUS_OK) {
 *     // g is all zero because its pointer was valid.
 * }
 * @endcode
 */
rp_kerr_status rp_kerr_metric(
    double mass,
    double spin,
    const double coordinates[RP_KERR_DIM],
    double metric[RP_KERR_DIM][RP_KERR_DIM]
);

/**
 * Evaluate the contravariant Kerr metric `g^(mu nu)`.
 *
 * @param mass Geometric mass `M`, with `M > 0`.
 * @param spin Spin length `a`, with `|a| <= M`.
 * @param coordinates Caller-owned input `(t, r, theta, phi)` in geometric
 *        units; angles are in radians.
 * @param inverse_metric Caller-owned `4 x 4` output buffer.  On failure after
 *        this pointer has been validated, every component is set to zero.
 * @return A typed status code.
 *
 * @note The caller owns all input and output storage.  No allocation occurs.
 *
 * @par Example
 * @code
 * const double x[RP_KERR_DIM] = {0.0, 8.0, 1.1, 0.3};
 * double inverse[RP_KERR_DIM][RP_KERR_DIM];
 * rp_kerr_status status = rp_kerr_inverse_metric(1.0, 0.5, x, inverse);
 * @endcode
 */
rp_kerr_status rp_kerr_inverse_metric(
    double mass,
    double spin,
    const double coordinates[RP_KERR_DIM],
    double inverse_metric[RP_KERR_DIM][RP_KERR_DIM]
);

/**
 * Evaluate the legacy-compatible Christoffel tensor `Gamma^lambda_(mu nu)`.
 *
 * The implementation is a direct C translation of the explicit expressions
 * and lower-index symmetry assignments in the frozen legacy RelatiPy method
 * `Kerr._get_christoffel_symbols`.  It maps the legacy Schwarzschild radius
 * as `R_s = 2 * mass`; `spin` already is the length-like legacy parameter
 * `a`.
 *
 * @param mass Geometric mass `M`, with `M > 0`.
 * @param spin Spin length `a`, with `|a| <= M`.
 * @param coordinates Caller-owned input `(t, r, theta, phi)` in geometric
 *        units; angles are in radians.
 * @param christoffel Caller-owned output indexed as
 *        `christoffel[lambda][mu][nu]`.  On failure after this pointer has
 *        been validated, every component is zero.
 * @return A typed status code.
 *
 * @note This is deliberately a legacy-compatibility tensor.  It is not
 *       recomputed from metric derivatives and must not be presented as an
 *       independently corrected connection for the metric exposed above.
 * @note The legacy method provides 20 independent expressions and 12 explicit
 *       lower-index symmetry copies, producing 32 populated entries.
 *       Components absent from that method remain zero.
 *
 * @par Example
 * @code
 * double storage[RP_KERR_DIM * RP_KERR_DIM * RP_KERR_DIM] = {
 *     0.0, 8.0, 1.1, 0.3
 * };
 * rp_kerr_status status = rp_kerr_christoffel(
 *     1.0,
 *     0.5,
 *     storage,
 *     (double (*)[RP_KERR_DIM][RP_KERR_DIM])storage
 * );
 * // Input and output are permitted to overlap; storage now holds Gamma.
 * @endcode
 */
rp_kerr_status rp_kerr_christoffel(
    double mass,
    double spin,
    const double coordinates[RP_KERR_DIM],
    double christoffel[RP_KERR_DIM][RP_KERR_DIM][RP_KERR_DIM]
);

/**
 * Convert coordinate-time velocities to a future-directed four-velocity.
 *
 * The input is `(dr/dt, dtheta/dt, dphi/dt)` in Boyer--Lindquist
 * coordinates.  The output is ordered as `(u^t, u^r, u^theta, u^phi)` and
 * satisfies `u^i / u^t = dx^i / dt`.  For this reference's `(-,+,+,+)`
 * metric signature, successful output satisfies `g_(mu nu) u^mu u^nu = -1`.
 *
 * @param mass Geometric mass `M`, with `M > 0`.
 * @param spin Spin length `a`, with `|a| <= M`.
 * @param coordinates Caller-owned input `(t, r, theta, phi)` in geometric
 *        units; angles are in radians.
 * @param coordinate_velocity Caller-owned input `(dr/dt, dtheta/dt,
 *        dphi/dt)` in the same geometric convention.  The radial component
 *        is dimensionless; the angular components have inverse-length units.
 * @param four_velocity Caller-owned four-component output.  On failure after
 *        this pointer has been validated, every component is set to zero.
 * @return A typed status code.  Non-timelike coordinate velocities return
 *         `RP_KERR_STATUS_NON_TIMELIKE_VELOCITY`.
 *
 * @note The future-directed branch is selected by requiring `u^t > 0`.
 * @note The function allocates no memory and retains no pointer.
 * @note Either input buffer may overlap the output buffer because all inputs
 *       are copied before the output is cleared.
 *
 * @par Example
 * @code
 * const double x[RP_KERR_DIM] = {0.0, 8.0, 1.1, 0.3};
 * const double v[3] = {-0.01, 0.002, 0.02};
 * double u[RP_KERR_DIM];
 * rp_kerr_status status = rp_kerr_four_velocity(1.0, 0.5, x, v, u);
 * @endcode
 */
rp_kerr_status rp_kerr_four_velocity(
    double mass,
    double spin,
    const double coordinates[RP_KERR_DIM],
    const double coordinate_velocity[3],
    double four_velocity[RP_KERR_DIM]
);

#endif /* RELATIPY_KERR_GEOMETRY_H */
