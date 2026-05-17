/**
 * @file radau_core.c
 * @brief Adaptive Radau IIA (3-stage, order 5) geodesic integrator for Kerr metric.
 *
 * Implements simplified Newton iteration (Hairer & Wanner §IV.8) with adaptive
 * step size, embedded error estimator, and dense output. Constraint projection
 * (normalization, E, Lz, Q) is applied after every accepted step.
 *
 * @par Build
 * @code
 * cc -std=c11 -O3 -fPIC -shared -o radau_core.so radau_core.c -lm
 * @endcode
 */

#include <math.h>
#include <string.h>
#include <complex.h>
#include <stdlib.h>

/* ===================================================================== */
/*  KERR METRIC INFRASTRUCTURE  (Boyer-Lindquist: [t, r, theta, phi])   */
/* ===================================================================== */

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
/*  CONSTRAINT PROJECTION                                                 */
/* ===================================================================== */

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
/*  GEODESIC RHS IN (q, p) COVARIANT PHASE SPACE                        */
/* ===================================================================== */

static void kerr_geodesic_rhs_qp(double Rs, double a,
                                   const double *y8, double *dy8)
{
    double u[4];
    kerr_ginv_dot_p(Rs, a, y8, y8 + 4, u);
    for (int i = 0; i < 4; i++) dy8[i] = u[i];
    kerr_kick_force(Rs, a, y8, y8 + 4, u, dy8 + 4);
}


/* ===================================================================== */
/*  ADAPTIVE RADAU IIA — CONSTANTS (from scipy integrate/_ivp/radau.py) */
/* ===================================================================== */

#define NEQ          8
#define RADAU_ORD    5
#define NMAX_ITER    6
#define MIN_FAC      0.2
#define MAX_FAC      10.0
#define SAFETY       0.9
#define MAX_STEPS    1000000

/* sqrt(6) = 2.449489742783178 */
static const double RADAU_C[3] = {
    0.15505102572168219,   /* (4 - sqrt(6)) / 10 */
    0.64494897427831781,   /* (4 + sqrt(6)) / 10 */
    1.0
};

/* Error embedding coefficients: E = [-13-7*S6, -13+7*S6, -1] / 3 */
static const double RADAU_E[3] = {
    -10.048809399827415,
     1.382142733160749,
    -0.333333333333333
};

/* Eigenvalue of A^{-1}: real */
#define MU_RE   3.6378342527444957
/* Eigenvalue of A^{-1}: complex  MU_COMPLEX = MU_CR + i*MU_CI
 * scipy: MU_COMPLEX = 2.681... - 3.050...j  (negative imaginary) */
#define MU_CR   2.6810828736277521
#define MU_CI  -3.0504301992474105

/* T, TI: transformation matrices (3x3) diagonalising A^{-1}
 * Values from scipy.integrate._ivp.radau (Hairer & Wanner §IV.8) */
static const double RAD_T[3][3] = {
    { 0.09443876248897524, -0.14125529502095421,  0.03002919410514742},
    { 0.25021312296533332,  0.20412935229379994, -0.38294211275726192},
    { 1.0,                  1.0,                  0.0}
};
static const double RAD_TI[3][3] = {
    { 4.17871859155190428,  0.32768282076106237,  0.52337644549944951},
    {-4.17871859155190428, -0.32768282076106237,  0.47662355450055044},
    { 0.50287263494578682, -2.57192694985560522,  0.59603920482822492}
};

/* Dense output polynomial coefficients P (3x3):
 * y(t+x*h) = y_old + Z^T @ P @ [x, x^2, x^3]  (Z is 3xNEQ stage-increments) */
static const double RAD_P[3][3] = {
    { 10.048809399827415, -25.629591780076972,  15.580782380249557},
    { -1.382142733160749,  10.296258446743638,  -8.914115713582889},
    {  0.333333333333333,  -2.666666666666667,   3.333333333333333}
};

#define EPS_MACH 2.2204460492503131e-16


/* ===================================================================== */
/*  8x8 REAL LU (Doolittle, partial pivoting)                            */
/* ===================================================================== */

static int lu8_decompose(double A[NEQ][NEQ], int piv[NEQ])
{
    for (int k = 0; k < NEQ; k++) {
        int max_r = k;
        double max_v = fabs(A[k][k]);
        for (int i = k+1; i < NEQ; i++)
            if (fabs(A[i][k]) > max_v) { max_v = fabs(A[i][k]); max_r = i; }
        if (max_v < 1e-30) return -1;
        piv[k] = max_r;
        if (max_r != k) {
            for (int j = 0; j < NEQ; j++) {
                double tmp = A[k][j]; A[k][j] = A[max_r][j]; A[max_r][j] = tmp;
            }
        }
        double inv = 1.0 / A[k][k];
        for (int i = k+1; i < NEQ; i++) {
            A[i][k] *= inv;
            for (int j = k+1; j < NEQ; j++) A[i][j] -= A[i][k] * A[k][j];
        }
    }
    return 0;
}

static void lu8_solve(const double LU[NEQ][NEQ], const int piv[NEQ],
                       const double *b, double *x)
{
    double tmp[NEQ];
    memcpy(tmp, b, NEQ * sizeof(double));
    for (int k = 0; k < NEQ; k++) {
        int p = piv[k];
        if (p != k) { double t = tmp[k]; tmp[k] = tmp[p]; tmp[p] = t; }
    }
    for (int i = 0; i < NEQ; i++)
        for (int j = 0; j < i; j++) tmp[i] -= LU[i][j] * tmp[j];
    for (int i = NEQ-1; i >= 0; i--) {
        for (int j = i+1; j < NEQ; j++) tmp[i] -= LU[i][j] * tmp[j];
        tmp[i] /= LU[i][i];
    }
    memcpy(x, tmp, NEQ * sizeof(double));
}


/* ===================================================================== */
/*  8x8 COMPLEX LU                                                       */
/* ===================================================================== */

static int lu8c_decompose(double _Complex A[NEQ][NEQ], int piv[NEQ])
{
    for (int k = 0; k < NEQ; k++) {
        int max_r = k;
        double max_v = cabs(A[k][k]);
        for (int i = k+1; i < NEQ; i++)
            if (cabs(A[i][k]) > max_v) { max_v = cabs(A[i][k]); max_r = i; }
        if (max_v < 1e-30) return -1;
        piv[k] = max_r;
        if (max_r != k) {
            for (int j = 0; j < NEQ; j++) {
                double _Complex tmp = A[k][j]; A[k][j] = A[max_r][j]; A[max_r][j] = tmp;
            }
        }
        double _Complex inv = 1.0 / A[k][k];
        for (int i = k+1; i < NEQ; i++) {
            A[i][k] *= inv;
            for (int j = k+1; j < NEQ; j++) A[i][j] -= A[i][k] * A[k][j];
        }
    }
    return 0;
}

static void lu8c_solve(const double _Complex LU[NEQ][NEQ], const int piv[NEQ],
                        const double _Complex *b, double _Complex *x)
{
    double _Complex tmp[NEQ];
    memcpy(tmp, b, NEQ * sizeof(double _Complex));
    for (int k = 0; k < NEQ; k++) {
        int p = piv[k];
        if (p != k) { double _Complex t = tmp[k]; tmp[k] = tmp[p]; tmp[p] = t; }
    }
    for (int i = 0; i < NEQ; i++)
        for (int j = 0; j < i; j++) tmp[i] -= LU[i][j] * tmp[j];
    for (int i = NEQ-1; i >= 0; i--) {
        for (int j = i+1; j < NEQ; j++) tmp[i] -= LU[i][j] * tmp[j];
        tmp[i] /= LU[i][i];
    }
    memcpy(x, tmp, NEQ * sizeof(double _Complex));
}


/* ===================================================================== */
/*  FINITE-DIFFERENCE JACOBIAN AND NEWTON MATRICES                       */
/* ===================================================================== */

static void compute_jac_fd(double Rs, double a,
                             const double *y, const double *f0,
                             double J[NEQ][NEQ])
{
    double yp[NEQ], fp[NEQ];
    double sq_eps = sqrt(EPS_MACH);
    memcpy(yp, y, NEQ * sizeof(double));
    for (int i = 0; i < NEQ; i++) {
        double h = sq_eps * fmax(1.0, fabs(y[i]));
        yp[i] = y[i] + h;
        kerr_geodesic_rhs_qp(Rs, a, yp, fp);
        for (int k = 0; k < NEQ; k++) J[k][i] = (fp[k] - f0[k]) / h;
        yp[i] = y[i];
    }
}

static void build_newton_matrices(const double J[NEQ][NEQ], double h,
                                   double Mr[NEQ][NEQ],
                                   double _Complex Mc[NEQ][NEQ])
{
    double gamma = MU_RE / h;
    double alpha = MU_CR / h;
    double beta  = MU_CI / h;
    for (int i = 0; i < NEQ; i++) {
        for (int j = 0; j < NEQ; j++) {
            Mr[i][j] = (i == j ? gamma : 0.0) - J[i][j];
            Mc[i][j] = (i == j ? (alpha + beta * _Complex_I) : 0.0) - (double _Complex)J[i][j];
        }
    }
}


/* ===================================================================== */
/*  SIMPLIFIED NEWTON ON COLLOCATION SYSTEM                              */
/* ===================================================================== */

/* Returns 1 if converged, 0 otherwise. Z[3][NEQ] updated in place. */
static int solve_collocation_system(
    double Rs, double a,
    const double *y, double h,
    double Z[3][NEQ],
    const double *scale,
    double newton_tol,
    const double LU_real[NEQ][NEQ], const int piv_real[NEQ],
    const double _Complex LU_cplx[NEQ][NEQ], const int piv_cplx[NEQ])
{
    double W[3][NEQ], F[3][NEQ];
    double fr[NEQ], dW0[NEQ];
    double _Complex fc[NEQ], dWc[NEQ];

    /* W = TI @ Z */
    for (int i = 0; i < 3; i++)
        for (int k = 0; k < NEQ; k++) {
            W[i][k] = 0.0;
            for (int s = 0; s < 3; s++) W[i][k] += RAD_TI[i][s] * Z[s][k];
        }

    double alpha = MU_CR / h, beta = MU_CI / h, gamma = MU_RE / h;
    double nrm_old = -1.0;
    int converged = 0;

    for (int iter = 0; iter < NMAX_ITER; iter++) {
        /* Evaluate RHS at each stage Y[s] = y + Z[s] */
        double y_stage[NEQ];
        for (int s = 0; s < 3; s++) {
            for (int k = 0; k < NEQ; k++) y_stage[k] = y[k] + Z[s][k];
            kerr_geodesic_rhs_qp(Rs, a, y_stage, F[s]);
        }

        /* Real RHS: TI[0] @ F - gamma * W[0] */
        for (int k = 0; k < NEQ; k++) {
            fr[k] = -gamma * W[0][k];
            for (int s = 0; s < 3; s++) fr[k] += RAD_TI[0][s] * F[s][k];
        }

        /* Complex RHS (see scipy solve_collocation_system) */
        for (int k = 0; k < NEQ; k++) {
            double ti1f = 0.0, ti2f = 0.0;
            for (int s = 0; s < 3; s++) {
                ti1f += RAD_TI[1][s] * F[s][k];
                ti2f += RAD_TI[2][s] * F[s][k];
            }
            double fcr = ti1f - alpha * W[1][k] + beta  * W[2][k];
            double fci = ti2f - alpha * W[2][k] - beta  * W[1][k];
            fc[k] = fcr + fci * _Complex_I;
        }

        /* Solve Newton systems */
        lu8_solve(LU_real, piv_real, fr, dW0);
        lu8c_solve(LU_cplx, piv_cplx, fc, dWc);

        /* Update W */
        for (int k = 0; k < NEQ; k++) {
            W[0][k] += dW0[k];
            W[1][k] += creal(dWc[k]);
            W[2][k] += cimag(dWc[k]);
        }

        /* Z = T @ W */
        for (int s = 0; s < 3; s++)
            for (int k = 0; k < NEQ; k++) {
                Z[s][k] = 0.0;
                for (int i = 0; i < 3; i++) Z[s][k] += RAD_T[s][i] * W[i][k];
            }

        /* Convergence: RMS norm of dW over scale */
        double nrm = 0.0;
        for (int k = 0; k < NEQ; k++) {
            double sc = scale[k];
            double d0 = dW0[k] / sc;
            double dr = creal(dWc[k]) / sc;
            double di = cimag(dWc[k]) / sc;
            nrm += d0*d0 + dr*dr + di*di;
        }
        nrm = sqrt(nrm / (3.0 * NEQ));

        if (nrm_old >= 0.0) {
            double rate = nrm / nrm_old;
            if (rate >= 1.0) break;
            if (rate / (1.0 - rate) * nrm < newton_tol) { converged = 1; break; }
        }
        if (nrm < newton_tol) { converged = 1; break; }
        nrm_old = nrm;
    }
    return converged;
}


/* ===================================================================== */
/*  ERROR ESTIMATION                                                      */
/* ===================================================================== */

static double estimate_error_norm(
    const double LU_real[NEQ][NEQ], const int piv_real[NEQ],
    const double Z[3][NEQ],
    const double *f0, double h,
    const double *y, const double *y_new,
    double rtol, double atol)
{
    /* scipy: error = solve(LU_real, f + E.dot(Z)/h)
     * E[s]*Z[s] has units of y; dividing by h gives y/time = same units as f */
    double rhs[NEQ], err[NEQ];
    double inv_h = 1.0 / h;
    for (int k = 0; k < NEQ; k++) {
        rhs[k] = f0[k];
        for (int s = 0; s < 3; s++) rhs[k] += RADAU_E[s] * Z[s][k] * inv_h;
    }
    lu8_solve(LU_real, piv_real, rhs, err);

    double nrm = 0.0;
    for (int k = 0; k < NEQ; k++) {
        double sc = atol + rtol * fmax(fabs(y[k]), fabs(y_new[k]));
        double e = err[k] / sc;
        nrm += e * e;
    }
    return sqrt(nrm / NEQ);
}


/* ===================================================================== */
/*  INITIAL STEP SELECTION (scipy select_initial_step)                   */
/* ===================================================================== */

static double select_initial_h(double Rs, double a,
                                 const double *y0, const double *f0,
                                 double t0, double t_bound,
                                 double rtol, double atol)
{
    double d0 = 0.0, d1 = 0.0;
    for (int k = 0; k < NEQ; k++) {
        double sc = atol + fabs(y0[k]) * rtol;
        d0 += (y0[k]/sc) * (y0[k]/sc);
        d1 += (f0[k]/sc) * (f0[k]/sc);
    }
    d0 = sqrt(d0 / NEQ); d1 = sqrt(d1 / NEQ);

    double h0 = (d0 < 1e-5 || d1 < 1e-5) ? 1e-6 : 0.01 * d0 / d1;
    h0 = fmin(h0, t_bound - t0);
    if (h0 < 1e-12) h0 = 1e-12;

    double y1[NEQ], f1[NEQ];
    for (int k = 0; k < NEQ; k++) y1[k] = y0[k] + h0 * f0[k];
    kerr_geodesic_rhs_qp(Rs, a, y1, f1);

    double d2 = 0.0;
    for (int k = 0; k < NEQ; k++) {
        double sc = atol + fabs(y0[k]) * rtol;
        double df = (f1[k] - f0[k]) / h0 / sc;
        d2 += df * df;
    }
    d2 = sqrt(d2 / NEQ);

    double d12 = fmax(d1, d2);
    double h1 = (d12 <= 1e-15) ? fmax(1e-6, h0 * 1e-3)
                                : pow(0.01 / d12, 1.0 / (RADAU_ORD + 1));
    return fmin(100.0 * h0, fmin(h1, t_bound - t0));
}


/* ===================================================================== */
/*  DENSE INTERPOLATION                                                   */
/* ===================================================================== */

/* Evaluate y(t_old + x*h) from stage increments Z and P polynomial. */
static void dense_eval(const double y_old[NEQ], const double Z[3][NEQ],
                        double x, double *y_out)
{
    double x2 = x * x, x3 = x2 * x;
    for (int k = 0; k < NEQ; k++) {
        double c0 = 0.0, c1 = 0.0, c2 = 0.0;
        for (int s = 0; s < 3; s++) {
            c0 += Z[s][k] * RAD_P[s][0];
            c1 += Z[s][k] * RAD_P[s][1];
            c2 += Z[s][k] * RAD_P[s][2];
        }
        y_out[k] = y_old[k] + c0*x + c1*x2 + c2*x3;
    }
}


/* ===================================================================== */
/*  ADAPTIVE DRIVER                                                       */
/* ===================================================================== */

/**
 * Adaptive Radau IIA (order 5) with simplified Newton and constraint projection.
 *
 * @param Rs        Schwarzschild radius
 * @param a         Kerr spin
 * @param y0        Initial state [q, u] (8 doubles, contravariant velocity)
 * @param t0        Initial proper time
 * @param t_bound   Final proper time
 * @param rtol      Relative tolerance
 * @param atol      Absolute tolerance
 * @param max_step  Maximum step size (<=0 means unbounded)
 * @param E0        Conserved energy
 * @param Lz0       Conserved angular momentum
 * @param Q0        Carter constant
 * @param tol_proj  Newton tolerance for projection
 * @param mip       Max iterations for projection
 * @param t_eval    If non-NULL: array of times at which to sample output (size n_eval)
 * @param n_eval    Length of t_eval
 * @param t_out     Output time array (preallocated)
 * @param y_out     Output state array, row-major (N_out x 8), [q, u] per row
 * @param step_cap  Max output rows when t_eval==NULL
 * @param n_out     Written: actual number of output rows
 * @param nst       Written: accepted steps
 * @param njac      Written: Jacobian evaluations
 * @param nlu       Written: LU factorizations
 *
 * @return 0 OK | -1 too many steps | -2 Newton/LU failure | -4 step underflow
 */
int kerr_radau2_adaptive(
    double Rs, double a,
    const double *y0,
    double t0, double t_bound,
    double rtol, double atol,
    double max_step,
    double E0, double Lz0, double Q0,
    double tol_proj, int mip,
    const double *t_eval, int n_eval,
    double *t_out, double *y_out,
    int step_cap,
    int *n_out,
    int *nst, int *njac, int *nlu)
{
    /* Working state in (q, p) covariant phase space */
    double yqp[NEQ];
    double f0[NEQ];
    double J[NEQ][NEQ];
    double LU_real_arr[NEQ][NEQ];
    double _Complex LU_cplx_arr[NEQ][NEQ];
    int piv_real[NEQ], piv_cplx[NEQ];
    double Z[3][NEQ];

    if (nst)  *nst  = 0;
    if (njac) *njac = 0;
    if (nlu)  *nlu  = 0;
    *n_out = 0;

    /* Convert [q, u] -> [q, p] */
    for (int k = 0; k < 4; k++) yqp[k] = y0[k];
    kerr_g_dot_u(Rs, a, y0, y0 + 4, yqp + 4);

    /* Project initial state */
    {
        double u_init[4];
        for (int k = 0; k < 4; k++) u_init[k] = y0[4 + k];
        kerr_project_constraints(Rs, a, yqp, u_init, E0, Lz0, Q0, tol_proj, mip);
        kerr_g_dot_u(Rs, a, yqp, u_init, yqp + 4);
    }

    int out_idx = 0;
    int eval_idx = 0;

    /* Advance eval_idx past any t_eval < t0 */
    if (t_eval)
        while (eval_idx < n_eval && t_eval[eval_idx] < t0) eval_idx++;

    /* Helper: store (q, u) at a given yqp state */
#define STORE_QU(ti, yqp_ptr) do {                                   \
        double _u_s[4];                                               \
        kerr_ginv_dot_p(Rs, a, (yqp_ptr), (yqp_ptr)+4, _u_s);        \
        kerr_project_constraints(Rs, a, (yqp_ptr), _u_s,             \
                                  E0, Lz0, Q0, tol_proj, mip);       \
        t_out[out_idx] = (ti);                                        \
        for (int _k = 0; _k < 4; _k++) y_out[out_idx*NEQ+_k] = (yqp_ptr)[_k]; \
        for (int _k = 0; _k < 4; _k++) y_out[out_idx*NEQ+4+_k] = _u_s[_k]; \
        out_idx++;                                                    \
    } while (0)

    /* Store initial point if in dump mode or t_eval starts at t0 */
    if (!t_eval) {
        if (out_idx < step_cap) {
            double u_s[4];
            kerr_ginv_dot_p(Rs, a, yqp, yqp + 4, u_s);
            t_out[out_idx] = t0;
            for (int k = 0; k < 4; k++) y_out[out_idx*NEQ+k]   = yqp[k];
            for (int k = 0; k < 4; k++) y_out[out_idx*NEQ+4+k] = u_s[k];
            out_idx++;
        }
    } else if (eval_idx < n_eval && fabs(t_eval[eval_idx] - t0) < 1e-12 * (t_bound - t0 + 1.0)) {
        STORE_QU(t_eval[eval_idx], yqp);
        eval_idx++;
    }

    /* Initial RHS */
    kerr_geodesic_rhs_qp(Rs, a, yqp, f0);

    /* Initial step size */
    double h = select_initial_h(Rs, a, yqp, f0, t0, t_bound, rtol, atol);
    if (max_step > 0.0) h = fmin(h, max_step);

    /* Initial Jacobian + LU */
    compute_jac_fd(Rs, a, yqp, f0, J);
    if (njac) (*njac)++;

    double h_lu = h;
    {
        double Mr[NEQ][NEQ];
        double _Complex Mc[NEQ][NEQ];
        build_newton_matrices(J, h, Mr, Mc);
        memcpy(LU_real_arr, Mr, sizeof(Mr));
        memcpy(LU_cplx_arr, Mc, sizeof(Mc));
        if (lu8_decompose(LU_real_arr, piv_real) != 0) return -2;
        if (lu8c_decompose(LU_cplx_arr, piv_cplx) != 0) return -2;
        if (nlu) (*nlu)++;
    }

    memset(Z, 0, sizeof(Z));

    double t = t0;
    int step_count = 0;
    int jac_age = 0;   /* steps since last J recompute */

    while (t < t_bound - 1e-13 * (t_bound - t0)) {
        if (step_count >= MAX_STEPS) return -1;

        /* Clip to not overshoot */
        if (t + h > t_bound) h = t_bound - t;
        if (h < 1e-14 * (1.0 + fabs(t))) return -4;

        /* Refactorize if h changed by >20% */
        double h_ratio = h / h_lu;
        if (fabs(h_ratio - 1.0) > 0.2) {
            double Mr[NEQ][NEQ];
            double _Complex Mc[NEQ][NEQ];
            build_newton_matrices(J, h, Mr, Mc);
            memcpy(LU_real_arr, Mr, sizeof(Mr));
            memcpy(LU_cplx_arr, Mc, sizeof(Mc));
            if (lu8_decompose(LU_real_arr, piv_real) != 0 ||
                lu8c_decompose(LU_cplx_arr, piv_cplx) != 0) {
                /* Singular: recompute J */
                kerr_geodesic_rhs_qp(Rs, a, yqp, f0);
                compute_jac_fd(Rs, a, yqp, f0, J);
                if (njac) (*njac)++;
                jac_age = 0;
                build_newton_matrices(J, h, Mr, Mc);
                memcpy(LU_real_arr, Mr, sizeof(Mr));
                memcpy(LU_cplx_arr, Mc, sizeof(Mc));
                if (lu8_decompose(LU_real_arr, piv_real) != 0) return -2;
                if (lu8c_decompose(LU_cplx_arr, piv_cplx) != 0) return -2;
            }
            h_lu = h;
            if (nlu) (*nlu)++;
        }

        /* Scale vector for this step */
        double scale[NEQ];
        for (int k = 0; k < NEQ; k++) scale[k] = atol + rtol * fabs(yqp[k]);

        /* Newton tolerance */
        double ntol = fmax(10.0 * EPS_MACH / rtol,
                           fmin(0.03, sqrt(rtol)));

        /* Newton iteration */
        int conv = solve_collocation_system(
            Rs, a, yqp, h, Z, scale, ntol,
            LU_real_arr, piv_real, LU_cplx_arr, piv_cplx);

        if (!conv) {
            /* Retry with fresh Jacobian if age > 0 */
            if (jac_age > 0) {
                kerr_geodesic_rhs_qp(Rs, a, yqp, f0);
                compute_jac_fd(Rs, a, yqp, f0, J);
                if (njac) (*njac)++;
                jac_age = 0;
                double Mr[NEQ][NEQ];
                double _Complex Mc[NEQ][NEQ];
                build_newton_matrices(J, h, Mr, Mc);
                memcpy(LU_real_arr, Mr, sizeof(Mr));
                memcpy(LU_cplx_arr, Mc, sizeof(Mc));
                if (lu8_decompose(LU_real_arr, piv_real) != 0) return -2;
                if (lu8c_decompose(LU_cplx_arr, piv_cplx) != 0) return -2;
                h_lu = h;
                if (nlu) (*nlu)++;
                memset(Z, 0, sizeof(Z));
                conv = solve_collocation_system(
                    Rs, a, yqp, h, Z, scale, ntol,
                    LU_real_arr, piv_real, LU_cplx_arr, piv_cplx);
            }
            if (!conv) {
                /* Shrink h, cold restart */
                h *= 0.25;
                memset(Z, 0, sizeof(Z));
                continue;
            }
        }

        /* y_new = y + Z[2] (last stage, c_3=1) */
        double y_new[NEQ];
        for (int k = 0; k < NEQ; k++) y_new[k] = yqp[k] + Z[2][k];

        /* Error estimate */
        double err_norm = estimate_error_norm(
            LU_real_arr, piv_real, Z, f0, h, yqp, y_new, rtol, atol);

        /* Step size factor: embedded estimator is 3rd order → exponent 1/4
         * (scipy: predict_factor uses error_norm**-0.25) */
        double factor = SAFETY * pow(err_norm + 1e-10, -0.25);
        factor = fmax(MIN_FAC, fmin(MAX_FAC, factor));
        double h_new = h * factor;
        if (max_step > 0.0) h_new = fmin(h_new, max_step);

        if (err_norm <= 1.0) {
            /* ---- ACCEPTED ---- */

            /* Dense output: fill t_eval points in (t, t+h] */
            if (t_eval) {
                while (eval_idx < n_eval && t_eval[eval_idx] <= t + h * (1.0 + 1e-10)) {
                    double te = fmin(t_eval[eval_idx], t + h);
                    double x  = (te - t) / h;
                    double y_te[NEQ];
                    dense_eval(yqp, Z, x, y_te);
                    /* Output as [q, u] with projection */
                    double u_te[4];
                    kerr_ginv_dot_p(Rs, a, y_te, y_te + 4, u_te);
                    kerr_project_constraints(Rs, a, y_te, u_te,
                                             E0, Lz0, Q0, tol_proj, mip);
                    t_out[out_idx] = te;
                    for (int k = 0; k < 4; k++) y_out[out_idx*NEQ+k]   = y_te[k];
                    for (int k = 0; k < 4; k++) y_out[out_idx*NEQ+4+k] = u_te[k];
                    out_idx++;
                    eval_idx++;
                }
            }

            /* Project y_new and update p */
            {
                double u_end[4];
                kerr_ginv_dot_p(Rs, a, y_new, y_new + 4, u_end);
                kerr_project_constraints(Rs, a, y_new, u_end,
                                         E0, Lz0, Q0, tol_proj, mip);
                kerr_g_dot_u(Rs, a, y_new, u_end, y_new + 4);

                /* Dump mode: store accepted step */
                if (!t_eval && out_idx < step_cap) {
                    t_out[out_idx] = t + h;
                    for (int k = 0; k < 4; k++) y_out[out_idx*NEQ+k]   = y_new[k];
                    for (int k = 0; k < 4; k++) y_out[out_idx*NEQ+4+k] = u_end[k];
                    out_idx++;
                }
            }

            /* Advance */
            memcpy(yqp, y_new, NEQ * sizeof(double));
            kerr_geodesic_rhs_qp(Rs, a, yqp, f0);
            t += h;
            step_count++;
            if (nst) (*nst)++;
            jac_age++;

            /* Recompute Jacobian every 20 accepted steps */
            if (jac_age >= 20) {
                compute_jac_fd(Rs, a, yqp, f0, J);
                if (njac) (*njac)++;
                jac_age = 0;
            }

            /* Warm-start Z for next step */
            double scale_z = h_new / h;
            for (int s = 0; s < 3; s++)
                for (int k = 0; k < NEQ; k++) Z[s][k] *= scale_z;

            h = h_new;

        } else {
            /* ---- REJECTED ---- */
            h = h_new;
            memset(Z, 0, sizeof(Z));
            jac_age = 0;
        }
    }

    *n_out = out_idx;
    return 0;

#undef STORE_QU
}


/* ===================================================================== */
/*  COMPATIBILITY SHIM (deprecated)                                       */
/* ===================================================================== */

/**
 * Legacy fixed-step Radau2 — now forwards to kerr_radau2_adaptive.
 * n_steps and n_fix_iter are ignored; rtol=1e-9, atol=1e-12 used instead.
 */
void kerr_radau2(
    double Rs, double a,
    const double *state0,
    double tau0, double tau_end,
    int n_steps, int stride, int n_fix_iter,
    double E0, double Lz0, double Q0,
    double *out_t, double *out_y, int *n_out)
{
    int step_cap = n_steps + 2;
    if (step_cap < 10000) step_cap = 10000;
    kerr_radau2_adaptive(
        Rs, a, state0, tau0, tau_end,
        1e-9, 1e-12, 0.0,
        E0, Lz0, Q0, 1e-12, 20,
        NULL, 0,
        out_t, out_y, step_cap, n_out,
        NULL, NULL, NULL);
}
