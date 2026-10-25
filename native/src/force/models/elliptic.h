/**
 * @file elliptic.h
 * @brief Complete elliptic integrals and Carlson symmetric forms.
 *
 * Private to the native force models.  The parameter convention is that of
 * `K(m) = int_0^{pi/2} (1 - m sin^2 t)^{-1/2} dt` (so `m = k^2`).  Every
 * routine is pure: no allocation, no retained pointer and no global state.
 */

#ifndef RELATIPY_FORCE_MODELS_ELLIPTIC_H
#define RELATIPY_FORCE_MODELS_ELLIPTIC_H

typedef enum {
    RP_ELLIPTIC_STATUS_OK = 0,
    RP_ELLIPTIC_STATUS_NULL_POINTER,
    /** `m` is NaN or `m > 1`. */
    RP_ELLIPTIC_STATUS_INVALID_ARGUMENT
} rp_elliptic_status;

/**
 * Carlson's `R_F(x, y, z)` by duplication (DLMF 19.36.1), relative error
 * below about 2e-16.  Requires finite `x, y, z >= 0` with at most one zero;
 * returns NaN otherwise and `+inf` when two arguments are zero.
 */
double rp_carlson_rf(double x, double y, double z);

/**
 * Carlson's `R_D(x, y, z)` (symmetric in `x, y`) by duplication
 * (DLMF 19.36.2).  Requires finite `x, y >= 0` with at most one zero and
 * `z > 0`; returns NaN otherwise.
 */
double rp_carlson_rd(double x, double y, double z);

/**
 * Complete integrals `K(m)` and `E(m)` for `m <= 1`, accurate to a few ulp.
 * Uses the arithmetic-geometric mean for `0 <= m <= 1/2` and Carlson's forms
 * otherwise.  At `m = 1`, `K = +inf` and `E = 1`.  Either output may be NULL
 * only if the other is not.  Failure leaves outputs unchanged.
 */
rp_elliptic_status rp_elliptic_ke(double m, double *k, double *e);

/**
 * Same as rp_elliptic_ke() from the complementary parameter `m1 = 1 - m`,
 * given directly so that `m -> 1` keeps full relative accuracy.  Requires
 * `m1 >= 0`; `m1 = 0` gives `K = +inf`, `E = 1`.
 */
rp_elliptic_status rp_elliptic_ke_complementary(double m1, double *k, double *e);

/**
 * Hot-path AGM for `0 <= m < 1` with both `m` and `m1 = 1 - m` supplied
 * (each exact to rounding, as when `m = 4Rd/p^2` and `m1 = q^2/p^2`).
 * Writes `K(m)` and `K(m) - E(m)`; the latter has no cancellation for
 * small `m`.  Inputs are not validated.
 */
void rp_elliptic_agm(double m, double m1, double *k, double *k_minus_e);

#endif /* RELATIPY_FORCE_MODELS_ELLIPTIC_H */
