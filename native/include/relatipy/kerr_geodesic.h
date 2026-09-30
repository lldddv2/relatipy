/**
 * @file kerr_geodesic.h
 * @brief Experimental, allocation-free Kerr geodesic operations.
 *
 * This internal reference uses Boyer--Lindquist coordinates ordered as
 * `(t, r, theta, phi)`, geometrized units `G = c = 1`, and the explicit
 * legacy-compatible Christoffel expressions exposed by `kerr_geometry.h`.
 * The null-geodesic operations instead use the Levi--Civita connection of
 * the corrected Kerr metric exposed by that header.
 *
 * @warning This header is not a public RelatiPy API or stable ABI.  In
 * particular, it does not select `(x, u)` as the production integrator's
 * canonical state.
 */

#ifndef RELATIPY_KERR_GEODESIC_H
#define RELATIPY_KERR_GEODESIC_H

#include "relatipy/kerr_geometry.h"

/**
 * Evaluate the affine Kerr geodesic right-hand side.
 *
 * For state `(x^mu, u^mu)`, the outputs are
 * `dx^mu/dlambda = u^mu` and
 * `du^mu/dlambda = -Gamma^mu_(alpha beta) u^alpha u^beta`.  The connection
 * is the explicit frozen-legacy compatibility reference used by
 * `rp_kerr_christoffel`; it is not reconstructed from metric derivatives.
 *
 * @param mass Geometric mass `M`, with `M > 0`.
 * @param spin Spin length `a`, with `|a| <= M`.
 * @param coordinates Caller-owned input `(t, r, theta, phi)` in geometric
 *        units; angles are in radians.
 * @param four_velocity Caller-owned input `(u^t, u^r, u^theta, u^phi)`.
 * @param coordinate_derivative Caller-owned four-component output `dx`.
 * @param four_velocity_derivative Caller-owned four-component output `du`.
 * @return A typed status code. `RP_KERR_STATUS_OK` is the only success code.
 *
 * @note The function allocates no memory and retains no pointer.
 * @note Input storage may overlap either output because both inputs are copied
 *       before the outputs are cleared.  The two output buffers themselves
 *       must not overlap.
 * @note A validated output buffer is zeroed on every error path.
 *
 * @par Example
 * @code
 * const double x[RP_KERR_DIM] = {0.0, 8.0, 1.1, 0.3};
 * const double u[RP_KERR_DIM] = {1.1, -0.01, 0.002, 0.02};
 * double dx[RP_KERR_DIM];
 * double du[RP_KERR_DIM];
 * rp_kerr_status status = rp_kerr_geodesic_rhs(
 *     1.0, 0.5, x, u, dx, du
 * );
 * @endcode
 */
rp_kerr_status rp_kerr_geodesic_rhs(
    double mass,
    double spin,
    const double coordinates[RP_KERR_DIM],
    const double four_velocity[RP_KERR_DIM],
    double coordinate_derivative[RP_KERR_DIM],
    double four_velocity_derivative[RP_KERR_DIM]
);

/**
 * Construct an affinely scaled Kerr null tangent from constants of motion.
 *
 * The output `k^mu = dx^mu/dlambda` satisfies the separated Kerr null
 * equations for energy `E`, axial angular momentum `L_z`, and Carter
 * constant `Q`.  `radial_direction` and `polar_direction` select the signs
 * of `k^r` and `k^theta` and must each be either `-1` or `+1`.  At an exact
 * turning point the corresponding component is zero for either sign.
 *
 * The scale of `lambda` remains explicit: replacing `(E, L_z, Q)` with
 * `(c E, c L_z, c^2 Q)` rescales the returned tangent by `c` for positive
 * `c`.  This operation therefore does not choose a production affine-
 * parameter normalization.
 *
 * @param mass Geometric mass `M`, with `M > 0`.
 * @param spin Spin length `a`, with `|a| <= M`.
 * @param coordinates Caller-owned input `(t, r, theta, phi)` in geometric
 *        units; angles are in radians.
 * @param energy Conserved energy `E` associated with stationarity.
 * @param axial_angular_momentum Conserved axial angular momentum `L_z`.
 * @param carter_constant Carter constant `Q` for a null geodesic.  `E`,
 *        `L_z` and `Q` must be finite and not all zero; otherwise the call
 *        returns `RP_KERR_STATUS_NONFINITE_INPUT` or
 *        `RP_KERR_STATUS_INVALID_PARAMETER`, respectively.
 * @param radial_direction Sign of radial motion, either `-1` or `+1`.
 * @param polar_direction Sign of polar motion, either `-1` or `+1`.
 * @param tangent Caller-owned output `(k^t, k^r, k^theta, k^phi)`.  On
 *        failure after this pointer has been validated, every component is
 *        zero.
 * @return `RP_KERR_STATUS_NO_REAL_NULL_TANGENT` when either separated
 *         potential is negative at the requested point; otherwise a typed
 *         status code with `RP_KERR_STATUS_OK` as the only success value.
 *
 * @note The function allocates no memory and retains no pointer.
 * @note `coordinates` and `tangent` may overlap because the coordinates are
 *       copied before the output is cleared.
 */
rp_kerr_status rp_kerr_null_tangent(
    double mass,
    double spin,
    const double coordinates[RP_KERR_DIM],
    double energy,
    double axial_angular_momentum,
    double carter_constant,
    int radial_direction,
    int polar_direction,
    double tangent[RP_KERR_DIM]
);

/**
 * Evaluate the corrected Kerr null-geodesic right-hand side.
 *
 * For state `(x^mu, k^mu)`, the outputs are
 * `dx^mu/dlambda = k^mu` and
 * `dk^mu/dlambda = -Gamma^mu_(alpha beta) k^alpha k^beta`, where `Gamma` is
 * the Levi--Civita connection derived from the corrected Kerr metric in
 * `rp_kerr_metric`.  A tangent produced by `rp_kerr_null_tangent` is null;
 * this function deliberately does not recheck the null constraint at every
 * ODE stage.
 *
 * @param mass Geometric mass `M`, with `M > 0`.
 * @param spin Spin length `a`, with `|a| <= M`.
 * @param coordinates Caller-owned input `(t, r, theta, phi)`.
 * @param tangent Caller-owned input `(k^t, k^r, k^theta, k^phi)`.
 * @param coordinate_derivative Caller-owned output `dx/dlambda`.
 * @param tangent_derivative Caller-owned output `dk/dlambda`.
 * @return A typed status code. `RP_KERR_STATUS_OK` is the only success code.
 *
 * @note Inputs may overlap outputs after being copied.  The two output
 *       buffers themselves must not overlap.
 * @note The caller is responsible for providing a null initial tangent.  The
 *       absence of per-stage constraint rejection allows adaptive Runge--
 *       Kutta and collocation methods to evaluate off-manifold stage states.
 */
rp_kerr_status rp_kerr_null_geodesic_rhs(
    double mass,
    double spin,
    const double coordinates[RP_KERR_DIM],
    const double tangent[RP_KERR_DIM],
    double coordinate_derivative[RP_KERR_DIM],
    double tangent_derivative[RP_KERR_DIM]
);

#endif /* RELATIPY_KERR_GEODESIC_H */
