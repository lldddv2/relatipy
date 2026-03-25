/*
 * yoshida6_core.c
 *
 * Yoshida 6th-order symplectic integrator for geodesic equations.
 * Supports Schwarzschild (spherical) and Kerr (Boyer-Lindquist) metrics.
 *
 * Compile (macOS/Linux):
 *   cc -O3 -march=native -shared -fPIC -o yoshida6_core.so yoshida6_core.c -lm
 *
 * Phase-space variables:
 *   q^mu  – contravariant coordinates
 *   p_mu  – covariant momenta  p_mu = g_{mu nu} u^nu
 *
 * Kick force:
 *   F_m = sum_{a,r} p[a] * Gamma^a_{r m} * u[r]
 *       = sum_{a,r} p[a] * Gamma[a][r][m] * u[r]
 *
 * This is equivalent to the einsum 'sa,arm,s,r->m' (g_{sa} Gamma^a_{rm} u^s u^r)
 * because p_a = g_{as} u^s, so sum_s g_{sa} u^s = p_a (g symmetric).
 */

#include <math.h>
#include <string.h>

/* ---- Yoshida 6th-order weights (7-stage palindromic) ---- */
static const double YOSHIDA_W[7] = {
     0.784513610477560,
     0.235573213359357,
    -1.17767998417887,
     1.31518632068391,
    -1.17767998417887,
     0.235573213359357,
     0.784513610477560
};

/* ===================================================================== */
/*  SCHWARZSCHILD (spherical coordinates: [t, r, theta, phi])            */
/* ===================================================================== */
/*
 * Metric: g = diag(A, -1/A, -r^2, -r^2 sin^2 theta),  A = 1 - Rs/r
 * Inverse: g^-1 = diag(1/A, -A, -1/r^2, -1/(r^2 sin^2 theta))
 *
 * Non-zero Christoffels (Gamma^a_{bc}):
 *   G[0][0][1] = G[0][1][0] = Rs / (2r(r-Rs))
 *   G[1][0][0]              = Rs(r-Rs) / (2r^3)
 *   G[1][1][1]              = Rs(Rs-r) / (2r^3 A^2)
 *   G[1][2][2]              = Rs - r
 *   G[1][3][3]              = (Rs-r) sin^2(theta)
 *   G[2][1][2] = G[2][2][1] = 1/r
 *   G[2][3][3]              = -sin(theta) cos(theta)
 *   G[3][1][3] = G[3][3][1] = 1/r
 *   G[3][2][3] = G[3][3][2] = cos(theta)/sin(theta)
 */

/* u^mu = g^{mu nu} p_nu  (diagonal inverse metric) */
static void sch_ginv_dot_p(double Rs, const double *q, const double *p, double *u)
{
    double r = q[1], theta = q[2];
    double A  = 1.0 - Rs / r;
    double s  = sin(theta);
    double r2 = r * r;
    u[0] =  p[0] / A;
    u[1] = -A * p[1];
    u[2] = -p[2] / r2;
    u[3] = -p[3] / (r2 * s * s);
}

/* p_mu = g_{mu nu} u^nu  (diagonal covariant metric) */
static void sch_g_dot_u(double Rs, const double *q, const double *u, double *p)
{
    double r = q[1], theta = q[2];
    double A  = 1.0 - Rs / r;
    double s  = sin(theta);
    double r2 = r * r;
    p[0] =  A * u[0];
    p[1] = -u[1] / A;
    p[2] = -r2 * u[2];
    p[3] = -r2 * s * s * u[3];
}

/*
 * Kick force for Schwarzschild:
 *   F[m] = sum_{a,r} p[a] * Gamma^a_{r m} * u[r]
 *
 * Expanded (only non-zero Christoffel contributions):
 *   F[0] = p[0]*G0_01*u[1]  + p[1]*G1_00*u[0]
 *   F[1] = p[0]*G0_01*u[0]  + p[1]*G1_11*u[1] + p[2]*G2_12*u[2] + p[3]*G3_13*u[3]
 *   F[2] = p[1]*G1_22*u[2]  + p[2]*G2_12*u[1] + p[3]*G3_23*u[3]
 *   F[3] = p[1]*G1_33*u[3]  + p[2]*G2_33*u[3] + p[3]*(G3_13*u[1] + G3_23*u[2])
 */
static void sch_kick_force(
    double Rs, const double *q, const double *p, const double *u, double *F)
{
    double r = q[1], theta = q[2];
    double s  = sin(theta), c = cos(theta);
    double s2 = s * s;
    double A  = 1.0 - Rs / r;
    double r2 = r * r;
    double r3 = r2 * r;

    double G0_01 =  Rs / (2.0 * r * (r - Rs));
    double G1_00 =  Rs * (r - Rs) / (2.0 * r3);
    double G1_11 =  Rs * (Rs - r) / (2.0 * r3 * A * A);
    double G1_22 =  Rs - r;
    double G1_33 = (Rs - r) * s2;
    double G2_12 =  1.0 / r;
    double G2_33 = -s * c;
    double G3_13 =  1.0 / r;
    double G3_23 =  c / s;

    F[0] = p[0] * G0_01 * u[1] + p[1] * G1_00 * u[0];
    F[1] = p[0] * G0_01 * u[0] + p[1] * G1_11 * u[1]
         + p[2] * G2_12 * u[2] + p[3] * G3_13 * u[3];
    F[2] = p[1] * G1_22 * u[2] + p[2] * G2_12 * u[1]
         + p[3] * G3_23 * u[3];
    F[3] = p[1] * G1_33 * u[3] + p[2] * G2_33 * u[3]
         + p[3] * (G3_13 * u[1] + G3_23 * u[2]);
}

/* One Störmer-Verlet DKD step (Schwarzschild) */
static void sch_verlet(double Rs, double h, double *q, double *p)
{
    double u[4], q_half[4], p_new[4], F[4];
    int i;

    /* half-drift: q_half = q + (h/2) * g_inv @ p */
    sch_ginv_dot_p(Rs, q, p, u);
    for (i = 0; i < 4; i++) q_half[i] = q[i] + 0.5 * h * u[i];

    /* kick at q_half: p_new = p + h * F(q_half, p) */
    sch_ginv_dot_p(Rs, q_half, p, u);
    sch_kick_force(Rs, q_half, p, u, F);
    for (i = 0; i < 4; i++) p_new[i] = p[i] + h * F[i];

    /* second half-drift: q = q_half + (h/2) * g_inv(q_half) @ p_new */
    sch_ginv_dot_p(Rs, q_half, p_new, u);
    for (i = 0; i < 4; i++) q[i] = q_half[i] + 0.5 * h * u[i];
    for (i = 0; i < 4; i++) p[i] = p_new[i];
}

/* One Yoshida 6th-order step (Schwarzschild) */
static void sch_yoshida(double Rs, double h, double *q, double *p)
{
    int k;
    for (k = 0; k < 7; k++)
        sch_verlet(Rs, YOSHIDA_W[k] * h, q, p);
}

/*
 * Main integration function for Schwarzschild.
 *
 * Output stored as (N_out, 8) row-major: out_y[n*8 + component].
 * Caller must pre-allocate out_t and out_y with at least (n_steps+1) entries.
 */
void yoshida6_schwarzschild(
    double Rs,
    const double *state0,    /* [q^0..q^3, u^0..u^3] */
    double tau0,
    double tau_end,
    int n_steps,
    int stride,
    double *out_t,
    double *out_y,           /* (N_out, 8) row-major */
    int *n_out)
{
    double q[4], p[4], u[4];
    int i, step, count = 0;

    for (i = 0; i < 4; i++) q[i] = state0[i];
    sch_g_dot_u(Rs, q, state0 + 4, p);     /* u -> p */

    double h   = (tau_end - tau0) / n_steps;
    double tau = tau0;

    /* store initial */
    out_t[0] = tau;
    sch_ginv_dot_p(Rs, q, p, u);
    for (i = 0; i < 4; i++) out_y[i]     = q[i];
    for (i = 0; i < 4; i++) out_y[4 + i] = u[i];
    count = 1;

    for (step = 0; step < n_steps; step++) {
        sch_yoshida(Rs, h, q, p);
        tau += h;
        if ((step + 1) % stride == 0 || step + 1 == n_steps) {
            out_t[count] = tau;
            sch_ginv_dot_p(Rs, q, p, u);
            for (i = 0; i < 4; i++) out_y[count * 8 + i]     = q[i];
            for (i = 0; i < 4; i++) out_y[count * 8 + 4 + i] = u[i];
            count++;
        }
    }
    *n_out = count;
}


/* ===================================================================== */
/*  KERR (Boyer-Lindquist coordinates: [t, r, theta, phi])               */
/* ===================================================================== */
/*
 * Metric (a = self.a = spin * Rs/2):
 *   g00 = 1 - Rs*r/Sigma,  g11 = -Sigma/Delta,  g22 = -Sigma,
 *   g33 = -(r^2+a^2 + Rs*r*a^2/Sigma*sin^2)*sin^2,
 *   g03 = g30 = Rs*r*a/Sigma * sin^2
 *   where Sigma = r^2 + a^2 cos^2(theta),  Delta = r^2 - Rs*r + a^2
 *
 * Inverse (2x2 block for 00/33, diagonal for 11/22):
 *   det2  = g00*g33 - g03^2
 *   g^00  = g33/det2,   g^33 = g00/det2,   g^03 = g^30 = -g03/det2
 *   g^11  = 1/g11 = -Delta/Sigma,   g^22 = 1/g22 = -1/Sigma
 */

/* Compute Kerr metric components at q */
static void kerr_metric(double Rs, double a,
                         const double *q,
                         double *g00, double *g11, double *g22,
                         double *g33, double *g03)
{
    double r     = q[1], theta = q[2];
    double s     = sin(theta);
    double c     = cos(theta);
    double s2    = s * s;
    double c2    = c * c;
    double r2    = r * r;
    double a2    = a * a;
    double Sigma = r2 + a2 * c2;
    double Delta = r2 - Rs * r + a2;

    *g00 = 1.0 - Rs * r / Sigma;
    *g11 = -Sigma / Delta;
    *g22 = -Sigma;
    *g33 = -(r2 + a2 + Rs * r * a2 / Sigma * s2) * s2;
    *g03 = Rs * r * a / Sigma * s2;
}

/* u^mu = g^{mu nu} p_nu  (Kerr: analytic 2x2 block inversion) */
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

/* p_mu = g_{mu nu} u^nu  (Kerr: includes off-diagonal g03) */
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
 *
 * Christoffel components (translated from kerr_metric.py):
 *   K0_01  Gamma[0][0][1] = Gamma[0][1][0]
 *   K0_02  Gamma[0][0][2] = Gamma[0][2][0]
 *   K0_13  Gamma[0][1][3] = Gamma[0][3][1]
 *   K0_23  Gamma[0][2][3] = Gamma[0][3][2]
 *   K1_00  Gamma[1][0][0]
 *   K1_03  Gamma[1][0][3] = Gamma[1][3][0]
 *   K1_11  Gamma[1][1][1]
 *   K1_12  Gamma[1][1][2] = Gamma[1][2][1]
 *   K1_22  Gamma[1][2][2]
 *   K1_33  Gamma[1][3][3]
 *   K2_00  Gamma[2][0][0]
 *   K2_03  Gamma[2][0][3] = Gamma[2][3][0]
 *   K2_11  Gamma[2][1][1]
 *   K2_12  Gamma[2][1][2] = Gamma[2][2][1]
 *   K2_22  Gamma[2][2][2]
 *   K2_33  Gamma[2][3][3]
 *   K3_01  Gamma[3][0][1] = Gamma[3][1][0]
 *   K3_02  Gamma[3][0][2] = Gamma[3][2][0]
 *   K3_13  Gamma[3][1][3] = Gamma[3][3][1]
 *   K3_23  Gamma[3][2][3] = Gamma[3][3][2]
 *
 * Expanded kick components:
 *   F[0] = p0*(K0_01*u1 + K0_02*u2) + p1*(K1_00*u0 + K1_03*u3)
 *          + p2*(K2_00*u0 + K2_03*u3) + p3*(K3_01*u1 + K3_02*u2)
 *   F[1] = p0*(K0_01*u0 + K0_13*u3) + p1*(K1_11*u1 + K1_12*u2)
 *          + p2*(K2_11*u1 + K2_12*u2) + p3*(K3_01*u0 + K3_13*u3)
 *   F[2] = p0*(K0_02*u0 + K0_23*u3) + p1*(K1_12*u1 + K1_22*u2)
 *          + p2*(K2_12*u1 + K2_22*u2) + p3*(K3_02*u0 + K3_23*u3)
 *   F[3] = p0*(K0_13*u1 + K0_23*u2) + p1*(K1_03*u0 + K1_33*u3)
 *          + p2*(K2_03*u0 + K2_33*u3) + p3*(K3_13*u1 + K3_23*u2)
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
    double a2   = a * a;
    double a4   = a2 * a2;
    double a6   = a2 * a4;
    double r2   = r * r;
    double r3   = r2 * r;
    double r4   = r2 * r2;
    double r5   = r2 * r3;
    double r6   = r3 * r3;
    double r7   = r4 * r3;

    double a2c2  = a2 * c2;            /* a^2 cos^2(theta) */
    double a4c2  = a4 * c2;
    double a2r2  = a2 * r2;
    double a4c4  = a4 * c4;

    double Sigma  = r2 + a2c2;
    double Sigma2 = Sigma * Sigma;
    double Sigma3 = Sigma2 * Sigma;
    double Delta  = r2 - Rs * r + a2;

    /* aux1 = -a^2 sin(2th) / (2*Sigma)  [denom_23 = 2*Sigma] */
    double aux1 = -a2 * s2th / (2.0 * Sigma);

    double one_minus_c4th = 1.0 - c4th;

    /* ---- mu = 0 ---- */
    double K0_01 = Rs * (a2r2 * s2 + a4 * s2 - a4 + r4)
                   / (2.0 * Sigma2 * Delta);

    double K0_02 = 2.0 * Rs * aux1 * r / (2.0 * Sigma);
    /* = Rs * aux1 * r / Sigma = -Rs*a^2*sin(2th)*r / (2*Sigma^2) */

    double K0_13 = Rs * a * s2 * (-a2c2 * r2 - a2r2 - 3.0 * r4 + a4c2)
                   / (2.0 * Sigma2 * Delta);

    double K0_23 = Rs * a2 * a * c * (s * s2) * r / Sigma2;
    /* = Rs*a^3*cos*sin^3*r / Sigma^2 */

    /* ---- mu = 1 ---- */
    double K1_00 = Rs * (Rs * a2c2 * r - Rs * r3 - a2c2 * r2 - a2r2 + a4c2 + r4)
                   / (2.0 * Sigma3);

    double K1_03 = Rs * a * s2
                   * (-Rs * a2c2 * r + Rs * r3 - r4 + a2c2 * r2 + a2r2 - a4c2)
                   / (2.0 * Sigma3);

    /* den_DS = Rs*a2c2*r + Rs*r3 - a2c2*r2 - r4 + a2r2 + a4c2 = Delta*Sigma */
    double den_DS = Delta * Sigma;

    double K1_11 = (-Rs * a2c2 / 2.0 + Rs * r2 / 2.0 + a2 * r + a2c2 * r)
                   / den_DS;

    double K1_12 = aux1;

    double K1_22 = r * (Rs * r + a2 - r2) / Sigma;

    double K1_33 = s2 * (
        Rs * a2 * s2 * r4 / 2.0
        + 2.0 * Rs * a2c2 * r4
        + Rs * a4c4 * r2
        - Rs * a4 * s2 * r2 / 2.0
        - Rs * a4 * r2 * one_minus_c4th / 16.0
        + Rs * a6 * one_minus_c4th / 16.0
        + Rs * r6
        - Rs2 * a2 * s2 * r3 / 2.0
        + Rs2 * a4 * r * one_minus_c4th / 16.0
        + a2 * r5
        - 2.0 * a2c2 * r5
        - a4c4 * r3
        + 2.0 * a4c2 * r3
        + a6 * c4 * r
        - r7
    ) / Sigma3;

    /* ---- mu = 2 ---- */
    double K2_00 = 4.0 * Rs * aux1 * r / (4.0 * Sigma2);
    /* = Rs * aux1 * r / Sigma^2 = -Rs*a^2*sin(2th)*r / (2*Sigma^3) */

    double K2_03 = 4.0 * Rs * a * s2th * r * (a2 + r2)
                   / (8.0 * Sigma3);
    /* denom_23^3 = (2*Sigma)^3 = 8*Sigma^3 */

    double K2_11 = -a2 * c * s / den_DS;

    double K2_12 = r / Sigma;

    double K2_22 = aux1;

    double K2_33 = c * s * (
        -2.0 * Rs * r3 * a2 * s2
        - 2.0 * r2 * a4c2
        + (-Rs * r) * a4 * s4
        + (-Rs * r) * a4 * one_minus_c4th / 4.0
        + (-r4) * a2
        + 2.0 * (-r4) * a2c2
        - a4c4 * r2
        - a6 * c4
        - r6
    ) / Sigma3;

    /* ---- mu = 3 ---- */
    double K3_01 = Rs * a * (-a2c2 + r2) / (2.0 * Sigma2 * Delta);

    double K3_02 = cot * (-Rs * r) * a / Sigma2;

    double K3_13 = (
        Rs * (-a2c2 * r2)
        + Rs * (-a2r2) * s2 / 2.0
        + Rs * (-r4)
        + Rs * a4 * one_minus_c4th / 16.0
        + 2.0 * a2c2 * r3
        + a4c4 * r
        + r5
    ) / (Sigma2 * Delta);

    double K3_23 = cot * (
        Rs * a2 * s2 * r
        + 2.0 * r2 * a2
        + 2.0 * (-a2r2) * s2
        - 2.0 * a4 * s2
        + a4 * s4
        + a4
        + r4
    ) / Sigma2;

    /* ---- Kick components ---- */
    double p0 = p[0], p1 = p[1], p2 = p[2], p3 = p[3];
    double u0 = u[0], u1 = u[1], u2 = u[2], u3 = u[3];

    F[0] = p0 * (K0_01 * u1 + K0_02 * u2)
         + p1 * (K1_00 * u0 + K1_03 * u3)
         + p2 * (K2_00 * u0 + K2_03 * u3)
         + p3 * (K3_01 * u1 + K3_02 * u2);

    F[1] = p0 * (K0_01 * u0 + K0_13 * u3)
         + p1 * (K1_11 * u1 + K1_12 * u2)
         + p2 * (K2_11 * u1 + K2_12 * u2)
         + p3 * (K3_01 * u0 + K3_13 * u3);

    F[2] = p0 * (K0_02 * u0 + K0_23 * u3)
         + p1 * (K1_12 * u1 + K1_22 * u2)
         + p2 * (K2_12 * u1 + K2_22 * u2)
         + p3 * (K3_02 * u0 + K3_23 * u3);

    F[3] = p0 * (K0_13 * u1 + K0_23 * u2)
         + p1 * (K1_03 * u0 + K1_33 * u3)
         + p2 * (K2_03 * u0 + K2_33 * u3)
         + p3 * (K3_13 * u1 + K3_23 * u2);
}

/* One Störmer-Verlet DKD step (Kerr) */
static void kerr_verlet(double Rs, double a, double h, double *q, double *p)
{
    double u[4], q_half[4], p_new[4], F[4];
    int i;

    kerr_ginv_dot_p(Rs, a, q, p, u);
    for (i = 0; i < 4; i++) q_half[i] = q[i] + 0.5 * h * u[i];

    kerr_ginv_dot_p(Rs, a, q_half, p, u);
    kerr_kick_force(Rs, a, q_half, p, u, F);
    for (i = 0; i < 4; i++) p_new[i] = p[i] + h * F[i];

    kerr_ginv_dot_p(Rs, a, q_half, p_new, u);
    for (i = 0; i < 4; i++) q[i] = q_half[i] + 0.5 * h * u[i];
    for (i = 0; i < 4; i++) p[i] = p_new[i];
}

/* One Yoshida 6th-order step (Kerr) */
static void kerr_yoshida(double Rs, double a, double h, double *q, double *p)
{
    int k;
    for (k = 0; k < 7; k++)
        kerr_verlet(Rs, a, YOSHIDA_W[k] * h, q, p);
}

/*
 * Main integration function for Kerr.
 * Output: (N_out, 8) row-major, same layout as Schwarzschild.
 */
void yoshida6_kerr(
    double Rs,
    double a,                /* self.a = spin_param * Rs/2 */
    const double *state0,
    double tau0,
    double tau_end,
    int n_steps,
    int stride,
    double *out_t,
    double *out_y,
    int *n_out)
{
    double q[4], p[4], u[4];
    int i, step, count = 0;

    for (i = 0; i < 4; i++) q[i] = state0[i];
    kerr_g_dot_u(Rs, a, q, state0 + 4, p);

    double h   = (tau_end - tau0) / n_steps;
    double tau = tau0;

    out_t[0] = tau;
    kerr_ginv_dot_p(Rs, a, q, p, u);
    for (i = 0; i < 4; i++) out_y[i]     = q[i];
    for (i = 0; i < 4; i++) out_y[4 + i] = u[i];
    count = 1;

    for (step = 0; step < n_steps; step++) {
        kerr_yoshida(Rs, a, h, q, p);
        tau += h;
        if ((step + 1) % stride == 0 || step + 1 == n_steps) {
            out_t[count] = tau;
            kerr_ginv_dot_p(Rs, a, q, p, u);
            for (i = 0; i < 4; i++) out_y[count * 8 + i]     = q[i];
            for (i = 0; i < 4; i++) out_y[count * 8 + 4 + i] = u[i];
            count++;
        }
    }
    *n_out = count;
}


/* ===================================================================== */
/*  KERR CONSTRAINT PROJECTION (Newton iteration)                        */
/* ===================================================================== */
/*
 * Projects u^mu onto the 4-constraint surface:
 *   C1 = g_{mu nu} u^mu u^nu - 1 = 0          (normalization)
 *   C2 = -(g_{0mu} u^mu) - E0   = 0           (energy)
 *   C3 =   g_{3mu} u^mu  - Lz0  = 0           (angular momentum)
 *   C4 = p_theta^2 + cos^2(th)[a^2(1-E^2) + Lz^2/sin^2(th)] - Q0 = 0
 *
 * Jacobian J[i][nu] = dC_i/du^nu (4x4, solved by Gaussian elim):
 *   J[0] = 2 * p         (p_nu = g_{nu mu} u^mu)
 *   J[1] = [-g00, 0, 0, -g03]
 *   J[2] = [ g30, 0, 0,  g33]
 *   J[3] = [dQ/du0, 0, 2*g22^2*u2, dQ/du3]
 *          dQ/du0 = cos^2(th)*(2*a^2*E*g00 + 2*Lz*g30/sin^2(th))
 *          dQ/du3 = cos^2(th)*(2*a^2*E*g03 + 2*Lz*g33/sin^2(th))
 */

/* Gaussian elimination with partial pivoting, 4x4 system: A x = b → x */
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
        /* partial pivoting */
        max_val = fabs(aug[k][k]);
        pivot = k;
        for (i = k + 1; i < 4; i++) {
            if (fabs(aug[i][k]) > max_val) {
                max_val = fabs(aug[i][k]);
                pivot = i;
            }
        }
        if (max_val < 1e-30) return -1;  /* singular */
        if (pivot != k) {
            for (j = 0; j <= 4; j++) {
                tmp = aug[k][j]; aug[k][j] = aug[pivot][j]; aug[pivot][j] = tmp;
            }
        }
        for (i = k + 1; i < 4; i++) {
            factor = aug[i][k] / aug[k][k];
            for (j = k; j <= 4; j++) aug[i][j] -= factor * aug[k][j];
        }
    }
    /* back-substitution */
    for (i = 3; i >= 0; i--) {
        x[i] = aug[i][4];
        for (j = i + 1; j < 4; j++) x[i] -= aug[i][j] * x[j];
        x[i] /= aug[i][i];
    }
    return 0;
}

/*
 * Single-point Newton projection.
 * q[4]: fixed position (Boyer-Lindquist)
 * u[4]: velocity, corrected in-place
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
    double C[4], dlt[4];
    double J[4][4];
    double dQ0, dQ3;

    kerr_metric(Rs, a, q, &g00, &g11, &g22, &g33, &g03);
    theta = q[2];
    s = sin(theta); c = cos(theta);
    s2 = s * s; c2 = c * c;

    for (iter = 0; iter < max_iter; iter++) {
        p0 = g00*u[0] + g03*u[3];
        p2 = g22*u[2];
        p3 = g03*u[0] + g33*u[3];  /* g30 = g03 */
        E  = -p0;
        Lz =  p3;

        double full_norm = g00*u[0]*u[0] + 2.0*g03*u[0]*u[3]
                         + g11*u[1]*u[1] + g22*u[2]*u[2] + g33*u[3]*u[3];
        double Q = p2*p2 + c2*(a*a*(1.0 - E*E) + Lz*Lz/s2);

        C[0] = full_norm - 1.0;
        C[1] = E  - E0;
        C[2] = Lz - Lz0;
        C[3] = Q  - Q0;

        /* convergence check */
        double max_c = 0.0;
        for (i = 0; i < 4; i++) {
            double ac = fabs(C[i]);
            if (ac > max_c) max_c = ac;
        }
        if (max_c < tol) break;

        /* Jacobian */
        double p_all[4];
        p_all[0] = p0;
        p_all[1] = g11 * u[1];
        p_all[2] = p2;
        p_all[3] = p3;

        /* Row 0: 2*p */
        for (i = 0; i < 4; i++) J[0][i] = 2.0 * p_all[i];
        /* Row 1: dC2/du = -g_{0nu} */
        J[1][0] = -g00; J[1][1] = 0.0; J[1][2] = 0.0; J[1][3] = -g03;
        /* Row 2: dC3/du = g_{3nu} */
        J[2][0] =  g03; J[2][1] = 0.0; J[2][2] = 0.0; J[2][3] =  g33;
        /* Row 3: dC4/du */
        dQ0 = c2 * (2.0*a*a*E*g00 + 2.0*Lz*g03/s2);
        dQ3 = c2 * (2.0*a*a*E*g03 + 2.0*Lz*g33/s2);
        J[3][0] = dQ0;
        J[3][1] = 0.0;
        J[3][2] = 2.0 * g22 * g22 * u[2];
        J[3][3] = dQ3;

        /* Newton step: J * delta = -C */
        double neg_C[4];
        for (i = 0; i < 4; i++) neg_C[i] = -C[i];
        if (solve4x4(J, neg_C, dlt) != 0) break;

        for (i = 0; i < 4; i++) u[i] += dlt[i];
    }
}

/*
 * Batch projection over N output points.
 * q_arr: (N, 4) row-major  — positions, fixed
 * u_arr: (N, 4) row-major  — velocities, corrected in-place
 */
void kerr_project_trajectory(
    double Rs, double a,
    const double *q_arr, double *u_arr,
    double E0, double Lz0, double Q0,
    double tol, int max_iter, int N)
{
    int n;
    for (n = 0; n < N; n++) {
        kerr_project_constraints(
            Rs, a,
            q_arr + n * 4,
            u_arr + n * 4,
            E0, Lz0, Q0,
            tol, max_iter);
    }
}


/* ===================================================================== */
/*  RADAU IIA 3-STAGE, ORDER 5 — coordinate space (q^mu, u^mu)          */
/* ===================================================================== */
/*
 * Fixed-point (Picard) iteration for the implicit stage equations.
 *
 * Butcher tableau (Hairer & Wanner, SODII, §IV.5):
 *   c = [(4-√6)/10,  (4+√6)/10,  1]
 *   A[i][j] as tabulated below
 *   b = A[2]  (since c_3 = 1,  Y_3 = y_{n+1} after convergence)
 *
 * Fixed-point loop:
 *   Y_i^{(k+1)} = y_n + h * sum_j A[i][j] * f(Y_j^{(k)})
 *
 * Convergence rate per iteration: h * ||A|| * L_f  where L_f is the
 * Lipschitz constant of the RHS. For orbital geodesics with
 * steps_per_period >= 20,  n_fix_iter = 5 achieves full 5th-order
 * accuracy at double precision.
 */
static const double RADAU3_A[3][3] = {
    { 0.19681547722366597, -0.06553542585019844,  0.02377097434422321 },
    { 0.39442431473908988,  0.29207341166522843, -0.04154875212599793 },
    { 0.37640306270046727,  0.51248582618842162,  0.11111111111111111 }
};
static const double RADAU3_C[3] = {
    0.15505102572168219, 0.64494897427831781, 1.0
};

/* ---- Kerr geodesic RHS in PHASE SPACE (q^mu, p_mu) ------------------
 *
 * Works with covariant momenta p_mu = g_{mu nu} u^nu, exactly like the
 * Yoshida integrator. This avoids the Boyer-Lindquist coordinate
 * singularity at theta=0 (where cot(theta) -> inf in (q, u) space).
 *
 * y8 = [q^0..q^3, p_0..p_3]
 * dy8 = [dq^mu/dtau, dp_mu/dtau]
 *      = [u^mu, F_mu]  where u^mu = g^{mu nu} p_nu, F_mu = kick force
 */
static void kerr_geodesic_rhs_qp(double Rs, double a,
                                   const double *y8, double *dy8)
{
    double u[4];
    kerr_ginv_dot_p(Rs, a, y8, y8 + 4, u);       /* u^mu = g^{mu nu} p_nu */
    for (int i = 0; i < 4; i++) dy8[i] = u[i];   /* dq/dtau = u           */
    kerr_kick_force(Rs, a, y8, y8 + 4, u, dy8 + 4); /* dp/dtau = F_mu     */
}

/* One Radau IIA step in phase space. y8_{n+1} = Y[2] after convergence. */
static void kerr_radau3_step_qp(double Rs, double a,
                                  double *y8, double h, int n_fix_iter)
{
    double Y[3][8], K[3][8], f0[8];
    int s, rr, i, iter;

    kerr_geodesic_rhs_qp(Rs, a, y8, f0);

    /* Predictor: explicit Euler at each Radau node */
    for (s = 0; s < 3; s++)
        for (i = 0; i < 8; i++)
            Y[s][i] = y8[i] + RADAU3_C[s] * h * f0[i];

    /* Fixed-point iterations */
    for (iter = 0; iter < n_fix_iter; iter++) {
        for (s = 0; s < 3; s++) kerr_geodesic_rhs_qp(Rs, a, Y[s], K[s]);
        for (s = 0; s < 3; s++)
            for (i = 0; i < 8; i++) {
                Y[s][i] = y8[i];
                for (rr = 0; rr < 3; rr++) Y[s][i] += h * RADAU3_A[s][rr] * K[rr][i];
            }
    }

    /* y_{n+1} = last stage (c_3 = 1 for Radau IIA) */
    for (i = 0; i < 8; i++) y8[i] = Y[2][i];
}

/*
 * Kerr Radau2 integration with built-in constraint projection.
 *
 * Parameters
 * ----------
 * Rs, a        : Kerr metric parameters
 * state0       : initial [q^0..q^3, u^0..u^3]
 * tau0,tau_end : proper-time span
 * n_steps      : total number of Radau steps
 * stride       : store output every stride steps (initial always stored)
 * n_fix_iter   : fixed-point iterations per step (5 recommended)
 * E0, Lz0, Q0 : conserved quantities (computed from initial state in Python)
 * out_t        : pre-allocated double[n_steps+2]
 * out_y        : pre-allocated double[(n_steps+2)*8], row-major (N, 8)
 * n_out        : written with actual number of stored points
 */
/*
 * Kerr Radau2 integration with built-in constraint projection.
 *
 * Internally works in phase space [q^mu, p_mu] to avoid the
 * Boyer-Lindquist coordinate singularity at theta=0 (same strategy
 * as the Yoshida integrator). Output is stored as [q^mu, u^mu].
 *
 * Parameters
 * ----------
 * Rs, a        : Kerr metric parameters
 * state0       : initial [q^0..q^3, u^0..u^3]   (contravariant velocity)
 * tau0,tau_end : proper-time span
 * n_steps      : total number of Radau steps
 * stride       : store output every stride steps (initial always stored)
 * n_fix_iter   : fixed-point iterations per step (5 recommended)
 * E0, Lz0, Q0 : conserved quantities (computed from initial state in Python)
 * out_t        : pre-allocated double[n_steps+2]
 * out_y        : pre-allocated double[(n_steps+2)*8], row-major (N, 8)
 * n_out        : written with actual number of stored points
 */
void kerr_radau2(
    double Rs, double a,
    const double *state0,
    double tau0, double tau_end,
    int n_steps, int stride, int n_fix_iter,
    double E0, double Lz0, double Q0,
    double *out_t, double *out_y, int *n_out)
{
    double yqp[8];  /* internal state [q^mu, p_mu] */
    double u[4];    /* temporary u^mu for projection / output */
    int i, step, count = 0;

    /* Convert initial state0 = [q, u] -> [q, p] */
    for (i = 0; i < 4; i++) yqp[i] = state0[i];
    kerr_g_dot_u(Rs, a, state0, state0 + 4, yqp + 4);   /* p = g @ u */

    double h   = (tau_end - tau0) / (double)n_steps;
    double tau = tau0;

    /* Project initial u and store [q, u] */
    for (i = 0; i < 4; i++) u[i] = state0[4 + i];
    kerr_project_constraints(Rs, a, yqp, u, E0, Lz0, Q0, 1e-12, 20);
    kerr_g_dot_u(Rs, a, yqp, u, yqp + 4);               /* recompute p */

    out_t[0] = tau;
    for (i = 0; i < 4; i++) out_y[i]     = yqp[i];
    for (i = 0; i < 4; i++) out_y[4 + i] = u[i];
    count = 1;

    for (step = 0; step < n_steps; step++) {
        kerr_radau3_step_qp(Rs, a, yqp, h, n_fix_iter);   /* step in (q,p) */
        tau += h;

        if ((step + 1) % stride == 0 || step + 1 == n_steps) {
            kerr_ginv_dot_p(Rs, a, yqp, yqp + 4, u);     /* u = g_inv @ p */
            kerr_project_constraints(Rs, a, yqp, u, E0, Lz0, Q0, 1e-12, 20);
            kerr_g_dot_u(Rs, a, yqp, u, yqp + 4);        /* update p */

            out_t[count] = tau;
            for (i = 0; i < 4; i++) out_y[count * 8 + i]     = yqp[i];
            for (i = 0; i < 4; i++) out_y[count * 8 + 4 + i] = u[i];
            count++;
        }
    }
    *n_out = count;
}
