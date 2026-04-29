/**
 * @file radau_core.c
 * @brief Radau IIA three-stage (order 5) geodesic integrator for the Kerr metric in Boyer–Lindquist coordinates.
 *
 * This module integrates the geodesic equations in **covariant phase space** \f$(q^\mu, p_\mu)\f$ so that the
 * right-hand side remains finite away from the polar axis, avoiding the \f$\cot\theta \to \infty\f$ singularity
 * that appears if one writes the equations directly in terms of coordinate velocities. The stored trajectory
 * is expressed as \f$[q^\mu, u^\mu]\f$ with \f$u^\mu = g^{\mu\nu} p_\nu\f$.
 *
 * At each stored output point, the four-velocity is projected onto the simultaneous solution of the unit
 * normalization, conserved energy \f$E\f$, axial angular momentum \f$L_z\f$, and Carter constant \f$Q\f$ using
 * a small Newton iteration.
 *
 * @par Build (example)
 * @code
 * cc -O3 -march=native -shared -fPIC -o radau_core.so radau_core.c -lm
 * @endcode
 *
 * @par References
 * Hairer, E., & Wanner, G. (1996). *Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems*.
 * Springer, Section IV.5 (Radau IIA methods).
 */

#include <math.h>
#include <string.h>


/* ===================================================================== */
/*  KERR METRIC INFRASTRUCTURE  (Boyer-Lindquist: [t, r, theta, phi])   */
/* ===================================================================== */

/**
 * @brief Evaluate Kerr metric components in Boyer–Lindquist coordinates.
 *
 * Coordinates are \f$q^\mu = (t, r, \theta, \phi)\f$. With Schwarzschild radius \f$R_s\f$ and spin \f$a\f$,
 * \f$\Sigma = r^2 + a^2\cos^2\theta\f$, \f$\Delta = r^2 - R_s r + a^2\f$, and the non-zero components are
 * \f[
 *   g_{00} = 1 - \frac{R_s r}{\Sigma}, \quad
 *   g_{11} = -\frac{\Sigma}{\Delta}, \quad
 *   g_{22} = -\Sigma,
 * \f]
 * \f[
 *   g_{33} = -\left(r^2 + a^2 + \frac{R_s r a^2}{\Sigma}\sin^2\theta\right)\sin^2\theta, \quad
 *   g_{03} = g_{30} = \frac{R_s r a}{\Sigma}\sin^2\theta.
 * \f]
 *
 * @param[in] Rs Schwarzschild radius \f$R_s\f$.
 * @param[in] a Kerr spin parameter \f$a\f$.
 * @param[in] q Pointer to \f$[t,r,\theta,\phi]\f$.
 * @param[out] g00 Component \f$g_{00}\f$.
 * @param[out] g11 Component \f$g_{11}\f$.
 * @param[out] g22 Component \f$g_{22}\f$.
 * @param[out] g33 Component \f$g_{33}\f$.
 * @param[out] g03 Component \f$g_{03}=g_{30}\f$.
 */
static void kerr_metric(double Rs, double a,
                         const double *q,
                         double *g00, double *g11, double *g22,
                         double *g33, double *g03)
{
    double r     = q[1], theta = q[2];
    double s     = sin(theta), c = cos(theta);
    double s2    = s * s, c2 = c * c;
    double r2    = r * r, a2 = a * a;
    double Sigma = r2 + a2 * c2;
    double Delta = r2 - Rs * r + a2;

    *g00 = 1.0 - Rs * r / Sigma;
    *g11 = -Sigma / Delta;
    *g22 = -Sigma;
    *g33 = -(r2 + a2 + Rs * r * a2 / Sigma * s2) * s2;
    *g03 = Rs * r * a / Sigma * s2;
}

/**
 * @brief Raise covector: \f$u^\mu = g^{\mu\nu} p_\nu\f$ using analytic \f$2\times 2\f$ inversion in the \f$(t,\phi)\f$ block.
 *
 * The \f$(r,\theta)\f$ directions decouple: \f$u^1 = p_1/g_{11}\f$, \f$u^2 = p_2/g_{22}\f$.
 * For the \f$(t,\phi)\f$ block, \f$\det(g_{00}g_{33}-g_{03}^2)\f$ is inverted explicitly.
 *
 * @param[in] Rs Schwarzschild radius \f$R_s\f$.
 * @param[in] a Kerr spin \f$a\f$.
 * @param[in] q Position \f$q^\mu\f$.
 * @param[in] p Covector \f$p_\mu\f$.
 * @param[out] u Contravariant four-velocity components \f$u^\mu\f$.
 */
static void kerr_ginv_dot_p(double Rs, double a,
                             const double *q, const double *p, double *u)
{
    double g00, g11, g22, g33, g03;
    kerr_metric(Rs, a, q, &g00, &g11, &g22, &g33, &g03);
    double det2 = g00 * g33 - g03 * g03;

    u[0] = ( g33 * p[0] - g03 * p[3]) / det2;
    u[1] = p[1] / g11;
    u[2] = p[2] / g22;
    u[3] = (-g03 * p[0] + g00 * p[3]) / det2;
}

/**
 * @brief Lower vector: \f$p_\mu = g_{\mu\nu} u^\nu\f$.
 *
 * @param[in] Rs Schwarzschild radius \f$R_s\f$.
 * @param[in] a Kerr spin \f$a\f$.
 * @param[in] q Position \f$q^\mu\f$.
 * @param[in] u Contravariant four-velocity \f$u^\mu\f$.
 * @param[out] p Covector \f$p_\mu\f$.
 */
static void kerr_g_dot_u(double Rs, double a,
                          const double *q, const double *u, double *p)
{
    double g00, g11, g22, g33, g03;
    kerr_metric(Rs, a, q, &g00, &g11, &g22, &g33, &g03);
    p[0] = g00 * u[0] + g03 * u[3];
    p[1] = g11 * u[1];
    p[2] = g22 * u[2];
    p[3] = g03 * u[0] + g33 * u[3];
}

/**
 * @brief Geodesic “kick” force for covariant momentum: \f$F_\mu = p_\alpha \Gamma^\alpha_{\mu\beta} u^\beta\f$.
 *
 * Christoffel symbols \f$\Gamma^\alpha_{\mu\beta}\f$ are implemented in closed form in Boyer–Lindquist
 * coordinates and contracted with \f$p_\alpha\f$ and \f$u^\beta\f$ as required by Hamiltonian geodesic flow
 * in \f$(q,p)\f$ space.
 *
 * @param[in] Rs Schwarzschild radius \f$R_s\f$.
 * @param[in] a Kerr spin \f$a\f$.
 * @param[in] q Position \f$q^\mu\f$.
 * @param[in] p Covector \f$p_\mu\f$.
 * @param[in] u Contravariant four-velocity \f$u^\mu\f$.
 * @param[out] F Components \f$F_\mu\f$.
 */
static void kerr_kick_force(
    double Rs, double a,
    const double *q, const double *p, const double *u, double *F)
{
    double r = q[1], theta = q[2];
    double s    = sin(theta), c = cos(theta);
    double s2   = s * s, c2 = c * c;
    double s4   = s2 * s2, c4 = c2 * c2;
    double s2th = sin(2.0 * theta);
    double c4th = cos(4.0 * theta);
    double cot  = c / s;

    double Rs2  = Rs * Rs;
    double a2   = a * a, a4 = a2 * a2, a6 = a2 * a4;
    double r2   = r * r, r3 = r2 * r, r4 = r2 * r2;
    double r5   = r2 * r3, r6 = r3 * r3, r7 = r4 * r3;

    double a2c2 = a2 * c2, a4c2 = a4 * c2;
    double a2r2 = a2 * r2, a4c4 = a4 * c4;

    double Sigma  = r2 + a2c2;
    double Sigma2 = Sigma * Sigma;
    double Sigma3 = Sigma2 * Sigma;
    double Delta  = r2 - Rs * r + a2;
    double den_DS = Delta * Sigma;

    double aux1           = -a2 * s2th / (2.0 * Sigma);
    double one_minus_c4th = 1.0 - c4th;

    double K0_01 = Rs * (a2r2*s2 + a4*s2 - a4 + r4) / (2.0*Sigma2*Delta);
    double K0_02 = 2.0 * Rs * aux1 * r / (2.0 * Sigma);
    double K0_13 = Rs*a*s2*(-a2c2*r2 - a2r2 - 3.0*r4 + a4c2) / (2.0*Sigma2*Delta);
    double K0_23 = Rs * a2 * a * c * (s * s2) * r / Sigma2;

    double K1_00 = Rs*(Rs*a2c2*r - Rs*r3 - a2c2*r2 - a2r2 + a4c2 + r4) / (2.0*Sigma3);
    double K1_03 = Rs*a*s2*(-Rs*a2c2*r + Rs*r3 - r4 + a2c2*r2 + a2r2 - a4c2)
                   / (2.0*Sigma3);
    double K1_11 = (-Rs*a2c2/2.0 + Rs*r2/2.0 + a2*r + a2c2*r) / den_DS;
    double K1_12 = aux1;
    double K1_22 = r * (Rs*r + a2 - r2) / Sigma;
    double K1_33 = s2 * (
        Rs*a2*s2*r4/2.0 + 2.0*Rs*a2c2*r4 + Rs*a4c4*r2
        - Rs*a4*s2*r2/2.0 - Rs*a4*r2*one_minus_c4th/16.0
        + Rs*a6*one_minus_c4th/16.0 + Rs*r6
        - Rs2*a2*s2*r3/2.0 + Rs2*a4*r*one_minus_c4th/16.0
        + a2*r5 - 2.0*a2c2*r5 - a4c4*r3 + 2.0*a4c2*r3 + a6*c4*r - r7
    ) / Sigma3;

    double K2_00 = 4.0 * Rs * aux1 * r / (4.0 * Sigma2);
    double K2_03 = 4.0 * Rs * a * s2th * r * (a2 + r2) / (8.0 * Sigma3);
    double K2_11 = -a2 * c * s / den_DS;
    double K2_12 = r / Sigma;
    double K2_22 = aux1;
    double K2_33 = c * s * (
        -2.0*Rs*r3*a2*s2 - 2.0*r2*a4c2
        + (-Rs*r)*a4*s4 + (-Rs*r)*a4*one_minus_c4th/4.0
        + (-r4)*a2 + 2.0*(-r4)*a2c2 - a4c4*r2 - a6*c4 - r6
    ) / Sigma3;

    double K3_01 = Rs * a * (-a2c2 + r2) / (2.0*Sigma2*Delta);
    double K3_02 = cot * (-Rs * r) * a / Sigma2;
    double K3_13 = (
        Rs*(-a2c2*r2) + Rs*(-a2r2)*s2/2.0 + Rs*(-r4)
        + Rs*a4*one_minus_c4th/16.0 + 2.0*a2c2*r3 + a4c4*r + r5
    ) / (Sigma2 * Delta);
    double K3_23 = cot * (
        Rs*a2*s2*r + 2.0*r2*a2 + 2.0*(-a2r2)*s2
        - 2.0*a4*s2 + a4*s4 + a4 + r4
    ) / Sigma2;

    double p0 = p[0], p1 = p[1], p2 = p[2], p3 = p[3];
    double u0 = u[0], u1 = u[1], u2 = u[2], u3 = u[3];

    F[0] = p0*(K0_01*u1 + K0_02*u2) + p1*(K1_00*u0 + K1_03*u3)
         + p2*(K2_00*u0 + K2_03*u3) + p3*(K3_01*u1 + K3_02*u2);

    F[1] = p0*(K0_01*u0 + K0_13*u3) + p1*(K1_11*u1 + K1_12*u2)
         + p2*(K2_11*u1 + K2_12*u2) + p3*(K3_01*u0 + K3_13*u3);

    F[2] = p0*(K0_02*u0 + K0_23*u3) + p1*(K1_12*u1 + K1_22*u2)
         + p2*(K2_12*u1 + K2_22*u2) + p3*(K3_02*u0 + K3_23*u3);

    F[3] = p0*(K0_13*u1 + K0_23*u2) + p1*(K1_03*u0 + K1_33*u3)
         + p2*(K2_03*u0 + K2_33*u3) + p3*(K3_13*u1 + K3_23*u2);
}


/* ===================================================================== */
/*  KERR CONSTRAINT PROJECTION (Newton iteration)                        */
/* ===================================================================== */

/**
 * @brief Solve a \f$4\times 4\f$ linear system \f$A x = b\f$ by Gaussian elimination with partial pivoting.
 *
 * @param[in] A Coefficient matrix \f$A\f$ (copied internally; not modified in place).
 * @param[in] b Right-hand side \f$b\f$.
 * @param[out] x Solution vector \f$x\f$.
 * @return 0 on success, \f$-1\f$ if the matrix is singular (within a tiny numerical threshold).
 */
static int solve4x4(double A[4][4], const double *b, double *x)
{
    double aug[4][5];
    int i, j, k, pivot;
    double factor, tmp, max_val;

    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++) aug[i][j] = A[i][j];
        aug[i][4] = b[i];
    }
    for (k = 0; k < 4; k++) {
        max_val = fabs(aug[k][k]); pivot = k;
        for (i = k+1; i < 4; i++)
            if (fabs(aug[i][k]) > max_val) { max_val = fabs(aug[i][k]); pivot = i; }
        if (max_val < 1e-30) return -1;
        if (pivot != k)
            for (j = 0; j <= 4; j++) { tmp=aug[k][j]; aug[k][j]=aug[pivot][j]; aug[pivot][j]=tmp; }
        for (i = k+1; i < 4; i++) {
            factor = aug[i][k] / aug[k][k];
            for (j = k; j <= 4; j++) aug[i][j] -= factor * aug[k][j];
        }
    }
    for (i = 3; i >= 0; i--) {
        x[i] = aug[i][4];
        for (j = i+1; j < 4; j++) x[i] -= aug[i][j] * x[j];
        x[i] /= aug[i][i];
    }
    return 0;
}

/**
 * @brief Project \f$u^\mu\f$ onto the simultaneous zeros of four Kerr geodesic constraints via Newton iteration.
 *
 * The residuals are
 * \f[
 *   C_1 = g_{\mu\nu} u^\mu u^\nu - 1, \quad
 *   C_2 = E - E_0, \quad C_3 = L_z - L_{z,0},
 * \f]
 * with \f$E = -p_0\f$, \f$L_z = p_3\f$, \f$p_\mu = g_{\mu\nu} u^\nu\f$, and Carter’s relation
 * \f[
 *   C_4 = p_\theta^2 + \cos^2\theta\left[a^2(1-E^2) + \frac{L_z^2}{\sin^2\theta}\right] - Q_0.
 * \f]
 * The Jacobian of \f$\mathbf{C}\f$ with respect to \f$u^\mu\f$ is assembled and a Newton step solves
 * \f$J\,\delta u = -\mathbf{C}\f$ until \f$\|\mathbf{C}\|_\infty < \texttt{tol}\f$ or \p max_iter is reached.
 *
 * @param[in] Rs Schwarzschild radius \f$R_s\f$.
 * @param[in] a Kerr spin \f$a\f$.
 * @param[in] q Position \f$q^\mu\f$ (read-only).
 * @param[in,out] u Four-velocity \f$u^\mu\f$; overwritten by the projected values.
 * @param[in] E0 Target conserved energy parameter \f$E_0\f$.
 * @param[in] Lz0 Target axial angular momentum \f$L_{z,0}\f$.
 * @param[in] Q0 Target Carter constant \f$Q_0\f$.
 * @param[in] tol Convergence tolerance on the maximum absolute residual.
 * @param[in] max_iter Maximum Newton iterations.
 */
void kerr_project_constraints(
    double Rs, double a,
    const double *q, double *u,
    double E0, double Lz0, double Q0,
    double tol, int max_iter)
{
    int iter, i;
    double g00, g11, g22, g33, g03;
    double theta, s, c, s2, c2;
    double p0, p2, p3, E, Lz;
    double C[4], dlt[4], J[4][4], dQ0, dQ3;

    kerr_metric(Rs, a, q, &g00, &g11, &g22, &g33, &g03);
    theta = q[2]; s = sin(theta); c = cos(theta); s2 = s*s; c2 = c*c;

    for (iter = 0; iter < max_iter; iter++) {
        p0 = g00*u[0] + g03*u[3];
        p2 = g22*u[2];
        p3 = g03*u[0] + g33*u[3];
        E  = -p0; Lz = p3;

        double full_norm = g00*u[0]*u[0] + 2.0*g03*u[0]*u[3]
                         + g11*u[1]*u[1] + g22*u[2]*u[2] + g33*u[3]*u[3];
        double Q = p2*p2 + c2*(a*a*(1.0 - E*E) + Lz*Lz/s2);

        C[0] = full_norm - 1.0;
        C[1] = E  - E0;
        C[2] = Lz - Lz0;
        C[3] = Q  - Q0;

        double max_c = 0.0;
        for (i = 0; i < 4; i++) { double ac = fabs(C[i]); if (ac > max_c) max_c = ac; }
        if (max_c < tol) break;

        double p_all[4] = { p0, g11*u[1], p2, p3 };
        for (i = 0; i < 4; i++) J[0][i] = 2.0 * p_all[i];
        J[1][0]=-g00; J[1][1]=0.0; J[1][2]=0.0; J[1][3]=-g03;
        J[2][0]= g03; J[2][1]=0.0; J[2][2]=0.0; J[2][3]= g33;
        dQ0 = c2*(2.0*a*a*E*g00 + 2.0*Lz*g03/s2);
        dQ3 = c2*(2.0*a*a*E*g03 + 2.0*Lz*g33/s2);
        J[3][0]=dQ0; J[3][1]=0.0; J[3][2]=2.0*g22*g22*u[2]; J[3][3]=dQ3;

        double neg_C[4];
        for (i = 0; i < 4; i++) neg_C[i] = -C[i];
        if (solve4x4(J, neg_C, dlt) != 0) break;
        for (i = 0; i < 4; i++) u[i] += dlt[i];
    }
}

/**
 * @brief Apply kerr_project_constraints() along a trajectory of \p N points.
 *
 * @param[in] Rs Schwarzschild radius \f$R_s\f$.
 * @param[in] a Kerr spin \f$a\f$.
 * @param[in] q_arr Positions, row-major \f$N\times 4\f$ array (\f$q^\mu\f$ per row).
 * @param[in,out] u_arr Four-velocities, row-major \f$N\times 4\f$; overwritten in place.
 * @param[in] E0 Target energy \f$E_0\f$.
 * @param[in] Lz0 Target \f$L_{z,0}\f$.
 * @param[in] Q0 Target Carter \f$Q_0\f$.
 * @param[in] tol Newton tolerance (passed through).
 * @param[in] max_iter Maximum Newton iterations (passed through).
 * @param[in] N Number of points along the trajectory.
 */
void kerr_project_trajectory(
    double Rs, double a,
    const double *q_arr, double *u_arr,
    double E0, double Lz0, double Q0,
    double tol, int max_iter, int N)
{
    int n;
    for (n = 0; n < N; n++)
        kerr_project_constraints(Rs, a, q_arr+n*4, u_arr+n*4,
                                  E0, Lz0, Q0, tol, max_iter);
}


/* ===================================================================== */
/*  RADAU IIA 3-STAGE, ORDER 5                                           */
/* ===================================================================== */

/** Butcher matrix \f$A\f$ for the three-stage Radau IIA method (row \f$s\f$, stage \f$r\f$). */
static const double RADAU3_A[3][3] = {
    { 0.19681547722366597, -0.06553542585019844,  0.02377097434422321 },
    { 0.39442431473908988,  0.29207341166522843, -0.04154875212599793 },
    { 0.37640306270046727,  0.51248582618842162,  0.11111111111111111 }
};
/** Stage abscissas \f$c_s\f$; the last stage has \f$c_3 = 1\f$ so the final stage value is the step end. */
static const double RADAU3_C[3] = {
    0.15505102572168219, 0.64494897427831781, 1.0
};

/**
 * @brief Right-hand side for geodesic flow in phase space \f$y = (q^\mu, p_\mu)\f$.
 *
 * With affine parameter \f$\tau\f$,
 * \f[
 *   \frac{\mathrm{d} q^\mu}{\mathrm{d}\tau} = u^\mu = g^{\mu\nu} p_\nu, \qquad
 *   \frac{\mathrm{d} p_\mu}{\mathrm{d}\tau} = F_\mu = p_\alpha \Gamma^\alpha_{\mu\beta} u^\beta.
 * \f]
 *
 * @param[in] Rs Schwarzschild radius \f$R_s\f$.
 * @param[in] a Kerr spin \f$a\f$.
 * @param[in] y8 Phase vector \f$[q^0,\ldots,q^3,p_0,\ldots,p_3]\f$.
 * @param[out] dy8 Time derivative \f$\mathrm{d}y/\mathrm{d}\tau\f$ (same layout as \p y8).
 */
static void kerr_geodesic_rhs_qp(double Rs, double a,
                                   const double *y8, double *dy8)
{
    double u[4];
    kerr_ginv_dot_p(Rs, a, y8, y8 + 4, u);        /* u^mu = g^{mu nu} p_nu */
    for (int i = 0; i < 4; i++) dy8[i] = u[i];    /* dq/dtau = u           */
    kerr_kick_force(Rs, a, y8, y8 + 4, u, dy8 + 4); /* dp/dtau = F_mu      */
}

/**
 * @brief One step of the three-stage Radau IIA method with fixed-point iteration on the implicit stages.
 *
 * Stages \f$Y_s\f$ satisfy \f$Y_s = y_n + h \sum_r A_{sr} f(Y_r)\f$. The initial guess uses the explicit
 * Euler increment along \f$f(y_n)\f$ at abscissas \f$c_s\f$, then Picard iteration refines \f$Y_s\f$ until
 * \p n_fix_iter passes. The accepted step is the third stage \f$Y_3\f$ (since \f$c_3=1\f$).
 *
 * @param[in] Rs Schwarzschild radius \f$R_s\f$.
 * @param[in] a Kerr spin \f$a\f$.
 * @param[in,out] y8 Phase vector at \f$\tau_n\f$; overwritten with \f$y_{n+1}\f$.
 * @param[in] h Step size \f$h\f$.
 * @param[in] n_fix_iter Number of fixed-point iterations per step.
 */
static void kerr_radau3_step_qp(double Rs, double a,
                                  double *y8, double h, int n_fix_iter)
{
    double Y[3][8], K[3][8], f0[8];
    int s, rr, i, iter;

    kerr_geodesic_rhs_qp(Rs, a, y8, f0);

    for (s = 0; s < 3; s++)
        for (i = 0; i < 8; i++)
            Y[s][i] = y8[i] + RADAU3_C[s] * h * f0[i];

    for (iter = 0; iter < n_fix_iter; iter++) {
        for (s = 0; s < 3; s++) kerr_geodesic_rhs_qp(Rs, a, Y[s], K[s]);
        for (s = 0; s < 3; s++)
            for (i = 0; i < 8; i++) {
                Y[s][i] = y8[i];
                for (rr = 0; rr < 3; rr++) Y[s][i] += h * RADAU3_A[s][rr] * K[rr][i];
            }
    }

    for (i = 0; i < 8; i++) y8[i] = Y[2][i];
}

/**
 * @brief Integrate a Kerr geodesic with Radau IIA (three stages, order 5) and store a subsampled trajectory.
 *
 * The routine converts the initial state from \f$[q^\mu, u^\mu]\f$ to \f$[q^\mu, p_\mu]\f$, projects
 * \f$u^\mu\f$ onto the constraint surface, advances with kerr_radau3_step_qp(), and at each stored time
 * recomputes \f$u^\mu = g^{\mu\nu} p_\nu\f$, projects again, and writes \f$[q^\mu, u^\mu]\f$ to \p out_y.
 *
 * @param[in] Rs Schwarzschild radius \f$R_s\f$.
 * @param[in] a Kerr spin \f$a\f$.
 * @param[in] state0 Initial \f$[q^0,\ldots,q^3,u^0,\ldots,u^3]\f$ (contravariant velocity).
 * @param[in] tau0 Initial proper time \f$\tau_0\f$.
 * @param[in] tau_end Final proper time.
 * @param[in] n_steps Number of Radau steps (uniform step \f$h = (\tau_{\mathrm{end}}-\tau_0)/\texttt{n\_steps}\f$).
 * @param[in] stride Store every \p stride-th completed step (the initial state is always stored).
 * @param[in] n_fix_iter Fixed-point iterations per implicit step (e.g. 5).
 * @param[in] E0 Target conserved energy \f$E_0\f$.
 * @param[in] Lz0 Target \f$L_{z,0}\f$.
 * @param[in] Q0 Target Carter constant \f$Q_0\f$.
 * @param[out] out_t Proper times at stored points, length at least \f$\texttt{n\_steps}+2\f$.
 * @param[out] out_y Row-major \f$(N_{\mathrm{out}}, 8)\f$ array \f$[q^\mu, u^\mu]\f$ per row.
 * @param[out] n_out Number of rows written to \p out_t and \p out_y.
 */
void kerr_radau2(
    double Rs, double a,
    const double *state0,
    double tau0, double tau_end,
    int n_steps, int stride, int n_fix_iter,
    double E0, double Lz0, double Q0,
    double *out_t, double *out_y, int *n_out)
{
    double yqp[8], u[4];
    int i, step, count = 0;

    /* Convert [q, u] -> [q, p] */
    for (i = 0; i < 4; i++) yqp[i] = state0[i];
    kerr_g_dot_u(Rs, a, state0, state0 + 4, yqp + 4);

    double h   = (tau_end - tau0) / (double)n_steps;
    double tau = tau0;

    /* Project initial u and store [q, u] */
    for (i = 0; i < 4; i++) u[i] = state0[4 + i];
    kerr_project_constraints(Rs, a, yqp, u, E0, Lz0, Q0, 1e-12, 20);
    kerr_g_dot_u(Rs, a, yqp, u, yqp + 4);

    out_t[0] = tau;
    for (i = 0; i < 4; i++) out_y[i]     = yqp[i];
    for (i = 0; i < 4; i++) out_y[4 + i] = u[i];
    count = 1;

    for (step = 0; step < n_steps; step++) {
        kerr_radau3_step_qp(Rs, a, yqp, h, n_fix_iter);
        tau += h;

        if ((step + 1) % stride == 0 || step + 1 == n_steps) {
            kerr_ginv_dot_p(Rs, a, yqp, yqp + 4, u);
            kerr_project_constraints(Rs, a, yqp, u, E0, Lz0, Q0, 1e-12, 20);
            kerr_g_dot_u(Rs, a, yqp, u, yqp + 4);

            out_t[count] = tau;
            for (i = 0; i < 4; i++) out_y[count*8 + i]     = yqp[i];
            for (i = 0; i < 4; i++) out_y[count*8 + 4 + i] = u[i];
            count++;
        }
    }
    *n_out = count;
}
