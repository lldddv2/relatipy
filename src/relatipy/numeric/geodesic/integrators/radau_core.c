/*
 * radau_core.c
 *
 * Radau IIA 3-stage (order 5) geodesic integrator for the Kerr metric
 * in Boyer-Lindquist coordinates.
 *
 * Compile (macOS/Linux):
 *   cc -O3 -march=native -shared -fPIC -o radau_core.so radau_core.c -lm
 *
 * Integration is performed in covariant phase space (q^mu, p_mu) to avoid
 * the Boyer-Lindquist coordinate singularity at theta=0 (same strategy as
 * the Yoshida integrator). Output is stored as [q^mu, u^mu].
 *
 * Constraint projection (normalization, E, Lz, Carter Q) is applied via
 * Newton iteration at every stored output point.
 *
 * Reference: Hairer & Wanner, "Solving ODEs II", §IV.5 (Radau IIA).
 */

#include <math.h>
#include <string.h>


/* ===================================================================== */
/*  KERR METRIC INFRASTRUCTURE  (Boyer-Lindquist: [t, r, theta, phi])   */
/* ===================================================================== */
/*
 * Metric (a = spin parameter):
 *   g00 = 1 - Rs*r/Sigma,  g11 = -Sigma/Delta,  g22 = -Sigma,
 *   g33 = -(r^2+a^2 + Rs*r*a^2/Sigma*sin^2)*sin^2,
 *   g03 = g30 = Rs*r*a/Sigma * sin^2
 *   Sigma = r^2 + a^2 cos^2(theta),  Delta = r^2 - Rs*r + a^2
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

/* u^mu = g^{mu nu} p_nu  (analytic 2×2 block inversion) */
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

/* p_mu = g_{mu nu} u^nu */
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

/*
 * Kick force for Kerr: F[m] = sum_{a,r} p[a] * Gamma^a_{r m} * u[r]
 * (analytical Christoffel symbols in Boyer-Lindquist coordinates)
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
/*
 * Projects u^mu onto the 4-constraint surface:
 *   C1 = g_{mu nu} u^mu u^nu - 1 = 0       (normalization)
 *   C2 = -(g_{0mu} u^mu) - E0   = 0        (energy)
 *   C3 =   g_{3mu} u^mu  - Lz0  = 0        (angular momentum)
 *   C4 = p_theta^2 + cos^2(th)[a^2(1-E^2) + Lz^2/sin^2(th)] - Q0 = 0
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
/*  RADAU IIA 3-STAGE, ORDER 5                                           */
/* ===================================================================== */
/*
 * Fixed-point (Picard) iteration for the implicit stage equations.
 *
 * Butcher tableau (Hairer & Wanner, SODII, §IV.5):
 *   c = [(4-√6)/10, (4+√6)/10, 1]
 *   b = A[2]  (c_3 = 1  =>  Y_3 = y_{n+1} after convergence)
 *
 * Integration in phase space [q^mu, p_mu] avoids the cot(theta) -> inf
 * singularity at the polar axis in Boyer-Lindquist coordinates.
 */
static const double RADAU3_A[3][3] = {
    { 0.19681547722366597, -0.06553542585019844,  0.02377097434422321 },
    { 0.39442431473908988,  0.29207341166522843, -0.04154875212599793 },
    { 0.37640306270046727,  0.51248582618842162,  0.11111111111111111 }
};
static const double RADAU3_C[3] = {
    0.15505102572168219, 0.64494897427831781, 1.0
};

/* Geodesic RHS in phase space: dy8/dtau where y8 = [q^mu, p_mu] */
static void kerr_geodesic_rhs_qp(double Rs, double a,
                                   const double *y8, double *dy8)
{
    double u[4];
    kerr_ginv_dot_p(Rs, a, y8, y8 + 4, u);        /* u^mu = g^{mu nu} p_nu */
    for (int i = 0; i < 4; i++) dy8[i] = u[i];    /* dq/dtau = u           */
    kerr_kick_force(Rs, a, y8, y8 + 4, u, dy8 + 4); /* dp/dtau = F_mu      */
}

/* One Radau IIA step. y8_{n+1} = Y[2] after fixed-point convergence. */
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

/*
 * Main integration function.
 *
 * Parameters
 * ----------
 * Rs, a        : Kerr metric parameters
 * state0       : initial [q^0..q^3, u^0..u^3]  (contravariant velocity)
 * tau0,tau_end : proper-time span
 * n_steps      : total Radau steps
 * stride       : store every stride-th step (initial always stored)
 * n_fix_iter   : fixed-point iterations per step (5 recommended)
 * E0,Lz0,Q0   : conserved quantities (energy, z-angular momentum, Carter)
 * out_t        : pre-allocated double[n_steps+2]
 * out_y        : pre-allocated double[(n_steps+2)*8], row-major (N, 8)
 * n_out        : output: number of stored points
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
