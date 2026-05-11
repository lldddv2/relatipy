/*
 * fujita_core.c
 *
 * Kerr geodesic integrator using Mino time lambda (Fujita parametrization), where
 *
 *   d tau / d lambda = Sigma(r, theta),   Sigma = r^2 + a^2 cos^2(theta)
 *
 * The radial and polar equations of motion decouple in Mino time:
 *
 *   (dr/dl)^2  = R(r)
 *   (dth/dl)^2 = Theta(theta)
 *
 * where R and Theta are the standard Kerr potentials.  The full system
 * integrated here is 7-dimensional:
 *
 *   y[0] = r
 *   y[1] = pr   = dr/dlambda   (NOT the covariant momentum)
 *   y[2] = theta
 *   y[3] = pth  = dtheta/dlambda
 *   y[4] = phi
 *   y[5] = t    (coordinate time)
 *   y[6] = tau  (proper time; dτ/dλ = Σ)
 *
 * Integration uses classic 4th-order Runge-Kutta in lambda.
 * The step size dλ is set so that the integration covers [tau0, tau_end]
 * in proper time; it is estimated as dλ = (tau_end-tau0) / (n_steps * Sigma0).
 * The integrator stops when tau >= tau_end.
 *
 * Output is the standard 8-state [t, r, theta, phi, u^t, u^r, u^theta, u^phi]
 * reconstructed at each stored step using u^mu = (dx^mu/dlambda) / Sigma.
 *
 * Build:
 *   cc -O3 -march=native -shared -fPIC -o fujita_core.so fujita_core.c -lm
 *
 * References
 * ----------
 * Mino, Y. (2003). Phys. Rev. D, 67, 084027.
 * Drasco, S. & Hughes, S. A. (2004). Phys. Rev. D, 69, 044708.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* -------------------------------------------------------------------------
 * Kerr potentials and their derivatives.
 * Rs = 2M is the Schwarzschild radius; a is the spin parameter.
 * -------------------------------------------------------------------------*/

static inline double kerr_Delta(double Rs, double a, double r)
{
    return r * r - Rs * r + a * a;
}

/* P(r) = E(r^2 + a^2) - a Lz */
static inline double kerr_P(double a, double E, double Lz, double r)
{
    return E * (r * r + a * a) - a * Lz;
}

/* Radial potential R(r) = P^2 - Delta * [(Lz - aE)^2 + Q + r^2] */
static inline double kerr_R(double Rs, double a, double E, double Lz,
                             double Q, double r)
{
    double Delta = kerr_Delta(Rs, a, r);
    double P     = kerr_P(a, E, Lz, r);
    double C     = (Lz - a * E) * (Lz - a * E) + Q + r * r;
    return P * P - Delta * C;
}

/* dR/dr */
static inline double kerr_dR_dr(double Rs, double a, double E, double Lz,
                                 double Q, double r)
{
    double Delta  = kerr_Delta(Rs, a, r);
    double P      = kerr_P(a, E, Lz, r);
    double C      = (Lz - a * E) * (Lz - a * E) + Q + r * r;
    double dP_dr  = 2.0 * E * r;
    double dDelta = 2.0 * r - Rs;
    double dC_dr  = 2.0 * r;
    return 2.0 * P * dP_dr - dDelta * C - Delta * dC_dr;
}

/* Polar potential Theta(theta) = Q - cos^2(th)*[a^2(1-E^2) + Lz^2/sin^2(th)] */
static inline double kerr_Theta(double a, double E, double Lz, double Q,
                                 double theta)
{
    double s  = sin(theta);
    double c  = cos(theta);
    double s2 = s * s;
    double c2 = c * c;
    if (s2 < 1e-30) s2 = 1e-30;
    return Q - c2 * (a * a * (1.0 - E * E) + Lz * Lz / s2);
}

/*
 * dTheta/dtheta
 *
 *  Theta = Q - cos^2(th) * [a^2(1-E^2) + Lz^2/sin^2(th)]
 *
 *  Let A = a^2(1-E^2), B = Lz^2
 *  Theta = Q - c^2*A - B*c^2/s^2
 *
 *  d/dth [ c^2 ] = -sin(2th)
 *  d/dth [ c^2/s^2 ] = d/dth [cot^2(th)] = -2*cot(th)/sin^2(th) = -2*cos(th)/sin^3(th)
 *
 *  dTheta/dth = sin(2th)*A + B * 2*cos(th)/sin^3(th)
 */
static inline double kerr_dTheta_dtheta(double a, double E, double Lz,
                                         double theta)
{
    double s  = sin(theta);
    double c  = cos(theta);
    double s3 = s * s * s;
    if (fabs(s3) < 1e-30) s3 = copysign(1e-30, s3);
    double A = a * a * (1.0 - E * E);
    double B = Lz * Lz;
    return sin(2.0 * theta) * A + 2.0 * B * c / s3;
}

/* -------------------------------------------------------------------------
 * RHS of the 7-ODE system in Mino time lambda.
 * y = [r, pr, theta, pth, phi, t, tau]
 * -------------------------------------------------------------------------*/
static void mino_rhs(double Rs, double a, double E, double Lz, double Q,
                     const double *y, double *f)
{
    double r     = y[0];
    double pr    = y[1];
    double theta = y[2];
    double pth   = y[3];

    double s     = sin(theta);
    double c     = cos(theta);
    double s2    = s * s;
    double c2    = c * c;
    if (s2 < 1e-30) s2 = 1e-30;

    double Delta = kerr_Delta(Rs, a, r);
    if (fabs(Delta) < 1e-30) Delta = copysign(1e-30, Delta);

    double P     = kerr_P(a, E, Lz, r);
    double Sigma = r * r + a * a * c2;

    /* dr/dl, dpr/dl */
    f[0] = pr;
    f[1] = 0.5 * kerr_dR_dr(Rs, a, E, Lz, Q, r);

    /* dtheta/dl, dpth/dl */
    f[2] = pth;
    f[3] = 0.5 * kerr_dTheta_dtheta(a, E, Lz, theta);

    /*
     * The codebase uses the (+,-,-,-) metric signature, so E = -p_t < 0 and
     * Lz = p_phi < 0 for a prograde equatorial orbit.  Under the substitution
     * E_std = -E, Lz_std = -Lz (where E_std > 0 is the standard convention),
     * the standard Mino equations transform as follows:
     *
     *   dphi/dl_std = -(a*E_std - Lz_std/sin^2) + a*P_std/Delta
     *              = +(a*E - Lz/sin^2) - a*P/Delta
     *
     *   dt/dl_std   = -a*(a*E_std*sin^2 - Lz_std) + (r^2+a^2)*P_std/Delta
     *              = +a*(a*E*sin^2 - Lz) - (r^2+a^2)*P/Delta
     */

    /* dphi/dl = +(aE - Lz/sin^2) - a*P/Delta */
    f[4] = (a * E - Lz / s2) - a * P / Delta;

    /* dt/dl = +a*(a*E*sin^2 - Lz) - (r^2+a^2)*P/Delta */
    f[5] = a * (a * E * s2 - Lz) - (r * r + a * a) * P / Delta;

    /* dtau/dl = Sigma */
    f[6] = Sigma;
}

/* -------------------------------------------------------------------------
 * Classic 4th-order Runge-Kutta step.
 * -------------------------------------------------------------------------*/
static void rk4_step(double Rs, double a, double E, double Lz, double Q,
                     const double *y, double dl, double *y_out)
{
    double k1[7], k2[7], k3[7], k4[7], ytmp[7];
    int i;

    mino_rhs(Rs, a, E, Lz, Q, y, k1);

    for (i = 0; i < 7; i++) ytmp[i] = y[i] + 0.5 * dl * k1[i];
    mino_rhs(Rs, a, E, Lz, Q, ytmp, k2);

    for (i = 0; i < 7; i++) ytmp[i] = y[i] + 0.5 * dl * k2[i];
    mino_rhs(Rs, a, E, Lz, Q, ytmp, k3);

    for (i = 0; i < 7; i++) ytmp[i] = y[i] + dl * k3[i];
    mino_rhs(Rs, a, E, Lz, Q, ytmp, k4);

    for (i = 0; i < 7; i++)
        y_out[i] = y[i] + (dl / 6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
}

/* -------------------------------------------------------------------------
 * Reconstruct the standard 8-state [t, r, theta, phi, u^t, u^r, u^th, u^ph]
 * from the Mino state y[7] and the current RHS f[7].
 *
 * u^mu = (dx^mu / dlambda) / Sigma
 * -------------------------------------------------------------------------*/
static void reconstruct_8state(const double *y, const double *f, double *out8)
{
    double r     = y[0];
    double theta = y[2];
    double c     = cos(theta);
    double c2    = c * c;
    double Sigma = r * r + f[6];   /* f[6] = Sigma stored in dtau/dl slot */
    /* Actually Sigma = r^2 + a^2*cos^2(theta); f[6] is dtau/dl = Sigma */
    Sigma = f[6];

    out8[0] = y[5];          /* t */
    out8[1] = y[0];          /* r */
    out8[2] = y[2];          /* theta */
    out8[3] = y[4];          /* phi */
    if (Sigma < 1e-30) Sigma = 1e-30;
    out8[4] = f[5] / Sigma;  /* u^t  = (dt/dl) / Sigma */
    out8[5] = f[0] / Sigma;  /* u^r  = (dr/dl) / Sigma = pr / Sigma */
    out8[6] = f[2] / Sigma;  /* u^th = (dth/dl) / Sigma = pth / Sigma */
    out8[7] = f[4] / Sigma;  /* u^ph = (dphi/dl) / Sigma */
}

/* -------------------------------------------------------------------------
 * Public integration function.
 *
 * Parameters
 * ----------
 * Rs, a       : Kerr parameters (Rs = 2M)
 * E, Lz, Q    : conserved quantities
 * state0[8]   : initial [t, r, theta, phi, u^t, u^r, u^theta, u^phi]
 * tau0        : initial proper time
 * tau_end     : final proper time
 * n_steps     : number of RK4 steps in lambda
 * stride      : store output every 'stride' steps
 * out_tau[*]  : preallocated array for output proper times (size >= n_steps/stride + 2)
 * out_y[*]    : preallocated array for output states (size >= (n_steps/stride+2)*8)
 * n_out       : set to the number of stored points
 * -------------------------------------------------------------------------*/
void kerr_mino(double Rs, double a, double E, double Lz, double Q,
               const double *state0,
               double tau0, double tau_end,
               int n_steps, int stride,
               double *out_tau, double *out_y, int *n_out)
{
    /* Convert standard 8-state to 7-component Mino state */
    double r0    = state0[1];
    double th0   = state0[2];
    double c0    = cos(th0);
    double Sigma0 = r0 * r0 + a * a * c0 * c0;
    if (Sigma0 < 1e-30) Sigma0 = 1e-30;

    double y[7];
    y[0] = r0;
    y[1] = state0[5] * Sigma0;   /* pr  = u^r * Sigma */
    y[2] = th0;
    y[3] = state0[6] * Sigma0;   /* pth = u^theta * Sigma */
    y[4] = state0[3];             /* phi */
    y[5] = state0[0];             /* t */
    y[6] = tau0;                  /* tau */

    /* Step size in Mino time: dl = (tau_end - tau0) / (n_steps * Sigma0) */
    double dtau  = tau_end - tau0;
    double dl    = dtau / ((double)n_steps * Sigma0);
    if (dl <= 0.0) { *n_out = 0; return; }

    int count = 0;

    /* Store initial state */
    {
        double f0[7];
        mino_rhs(Rs, a, E, Lz, Q, y, f0);
        double s8[8];
        reconstruct_8state(y, f0, s8);
        out_tau[count] = y[6];
        memcpy(&out_y[count * 8], s8, 8 * sizeof(double));
        count++;
    }

    double ytmp[7];
    int step;
    for (step = 0; step < n_steps; step++) {
        /* Adaptive last step: do not overshoot tau_end */
        double tau_now = y[6];
        if (tau_now >= tau_end) break;

        /* Recompute dl so we don't overshoot in proper time.
         * Since dtau/dl = Sigma, estimate remaining dl as:
         *   dl_remain = (tau_end - tau_now) / Sigma_now  (upper bound ~ dl)
         */
        double r_now    = y[0];
        double th_now   = y[2];
        double c_now    = cos(th_now);
        double Sigma_now = r_now * r_now + a * a * c_now * c_now;
        if (Sigma_now < 1e-30) Sigma_now = 1e-30;
        double dl_remain = (tau_end - tau_now) / Sigma_now;
        double dl_use    = (dl < dl_remain) ? dl : dl_remain;

        rk4_step(Rs, a, E, Lz, Q, y, dl_use, ytmp);
        memcpy(y, ytmp, 7 * sizeof(double));

        if ((step + 1) % stride == 0 || y[6] >= tau_end) {
            double f_cur[7];
            mino_rhs(Rs, a, E, Lz, Q, y, f_cur);
            double s8[8];
            reconstruct_8state(y, f_cur, s8);
            out_tau[count] = y[6];
            memcpy(&out_y[count * 8], s8, 8 * sizeof(double));
            count++;
        }
    }

    /* Ensure the last stored time equals tau_end exactly */
    if (count > 0)
        out_tau[count - 1] = tau_end;

    *n_out = count;
}
