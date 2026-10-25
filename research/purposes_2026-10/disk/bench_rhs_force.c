/* Build: gcc -std=c11 -O3 -Inative/src -Inative/include bench_rhs_force.c native/src/geodesic/physic/kerr.c native/src/utils/numeric.c native/src/utils/tensor.c -lm
   Result (2026-10-25, 4-core container under load, min of 5 runs, pinned): prepare+RHS 55 ns; +frame/projection only 68 ns;
   +Miyamoto-Nagai 95 ns; +bicubic (d,z) table 123-131 ns. */
#define _POSIX_C_SOURCE 199309L
#define M_PI 3.14159265358979323846
/* Microbenchmark: library timelike RHS vs RHS + projected external (disk) force. */
#include "geodesic/physic/kerr.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define ND 512
#define NZ 512
typedef struct {
    double R[3][3];             /* rotation: disk frame xi = R x */
    double d0, idd, z0, idz;    /* table origin and inverse spacing */
    double *ad, *az;            /* tabulated a_d(d,z), a_z(d,z), row-major [iz][id] */
    double mscale;              /* disk mass factor */
    double *tab; int mode;
} disk_t;

static inline double cr(double p0, double p1, double p2, double p3, double t) {
    /* Catmull-Rom cubic (C1) */
    return p1 + 0.5 * t * (p2 - p0 + t * (2.0*p0 - 5.0*p1 + 4.0*p2 - p3 + t * (3.0*(p1 - p2) + p3 - p0)));
}
static inline double interp(const double *T, double fx, double fz) {
    int ix = (int)fx, iz = (int)fz;
    if (ix < 1) ix = 1; if (ix > ND - 3) ix = ND - 3;
    if (iz < 1) iz = 1; if (iz > NZ - 3) iz = NZ - 3;
    double tx = fx - ix, tz = fz - iz, row[4];
    for (int k = 0; k < 4; ++k) {
        const double *p = T + (size_t)(iz - 1 + k) * ND + (ix - 1);
        row[k] = cr(p[0], p[1], p[2], p[3], tx);
    }
    return cr(row[0], row[1], row[2], row[3], tz);
}

static inline void interp2(const double *T, double fx, double fz, double *o0, double *o1) {
    int ix = (int)fx, iz = (int)fz;
    if (ix < 1) ix = 1; if (ix > ND - 3) ix = ND - 3;
    if (iz < 1) iz = 1; if (iz > NZ - 3) iz = NZ - 3;
    double tx = fx - ix, tz = fz - iz, r0[4], r1[4];
    for (int k = 0; k < 4; ++k) {
        const double *p = T + 2*((size_t)(iz - 1 + k) * ND + (ix - 1));
        r0[k] = cr(p[0], p[2], p[4], p[6], tx);
        r1[k] = cr(p[1], p[3], p[5], p[7], tx);
    }
    *o0 = cr(r0[0], r0[1], r0[2], r0[3], tz);
    *o1 = cr(r1[0], r1[1], r1[2], r1[3], tz);
}
/* Miyamoto-Nagai analytic disk (for comparison): Phi = -M/sqrt(R^2+(A+sqrt(z^2+B^2))^2) */
static inline void mn_accel(const disk_t *D, const double x[3], double a[3]) {
    const double A = 1e4, B = 1e2;
    double xi0 = D->R[0][0]*x[0] + D->R[0][1]*x[1] + D->R[0][2]*x[2];
    double xi1 = D->R[1][0]*x[0] + D->R[1][1]*x[1] + D->R[1][2]*x[2];
    double xi2 = D->R[2][0]*x[0] + D->R[2][1]*x[1] + D->R[2][2]*x[2];
    double zb = sqrt(xi2*xi2 + B*B), q = A + zb;
    double den = xi0*xi0 + xi1*xi1 + q*q, f = -D->mscale / (den * sqrt(den));
    double b0 = f*xi0, b1 = f*xi1, b2 = f*xi2*q/zb;
    a[0] = D->R[0][0]*b0 + D->R[1][0]*b1 + D->R[2][0]*b2;
    a[1] = D->R[0][1]*b0 + D->R[1][1]*b1 + D->R[2][1]*b2;
    a[2] = D->R[0][2]*b0 + D->R[1][2]*b1 + D->R[2][2]*b2;
}

static inline void disk_accel(const disk_t *D, const double x[3], double a[3]) {
    double xi0 = D->R[0][0]*x[0] + D->R[0][1]*x[1] + D->R[0][2]*x[2];
    double xi1 = D->R[1][0]*x[0] + D->R[1][1]*x[1] + D->R[1][2]*x[2];
    double xi2 = D->R[2][0]*x[0] + D->R[2][1]*x[1] + D->R[2][2]*x[2];
    double d = sqrt(xi0*xi0 + xi1*xi1);
    double fx = (d - D->d0) * D->idd, fz = (fabs(xi2) - D->z0) * D->idz;
    double ad, az;
    if (D->mode == 1) { interp2(D->tab, fx, fz, &ad, &az); ad *= D->mscale; az = D->mscale * copysign(az, xi2); }
    else { ad = D->mscale * interp(D->ad, fx, fz); az = D->mscale * copysign(interp(D->az, fx, fz), xi2); }
    double inv_d = 1.0 / d;
    double b0 = ad * xi0 * inv_d, b1 = ad * xi1 * inv_d, b2 = az;
    a[0] = D->R[0][0]*b0 + D->R[1][0]*b1 + D->R[2][0]*b2;
    a[1] = D->R[0][1]*b0 + D->R[1][1]*b1 + D->R[2][1]*b2;
    a[2] = D->R[0][2]*b0 + D->R[1][2]*b1 + D->R[2][2]*b2;
}
/* f^mu = -(g^{mu nu} + u^mu u^nu) d_nu Phi, with d_nu Phi = -a . dx/dx^nu (BL, G=c=1) */
static inline void add_external(double M, double s, const struct rp_kerr_geodesic_state *st,
                                const disk_t *D, double du[4]) {
    const double r = st->x[1], sth = st->sin_theta, cth = st->cos_theta;
    const double sph = sin(st->x[3]), cph = cos(st->x[3]);
    const double rho = sqrt(r*r + s*s), r_rho = r / rho;
    const double x[3] = {rho*sth*cph, rho*sth*sph, r*cth};
    double a[3];
    if (D->mode == 3) { a[0] = 1e-12*x[0]; a[1] = 1e-12*x[1]; a[2] = 1e-12*x[2]; }
    else if (D->mode == 2) mn_accel(D, x, a); else disk_accel(D, x, a);
    const double dPr = -(a[0]*r_rho*sth*cph + a[1]*r_rho*sth*sph + a[2]*cth);
    const double dPt = -(a[0]*rho*cth*cph + a[1]*rho*cth*sph - a[2]*r*sth);
    const double dPp = -(-a[0]*rho*sth*sph + a[1]*rho*sth*cph);
    const double sig = st->sigma, del = st->delta, s2th = sth*sth;
    const double inv_sd = 1.0 / (sig * del);
    const double gtp = -2.0*M*s*r*inv_sd;
    const double grr = del / sig, gthth = 1.0 / sig;
    const double gpp = (del - s*s*s2th) * inv_sd / s2th;
    const double udP = st->u[1]*dPr + st->u[2]*dPt + st->u[3]*dPp;
    du[0] -= gtp*dPp + st->u[0]*udP;
    du[1] -= grr*dPr + st->u[1]*udP;
    du[2] -= gthth*dPt + st->u[2]*udP;
    du[3] -= gpp*dPp + st->u[3]*udP;
}
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9*t.tv_nsec; }

int main(void) {
    const double M = 1.0, s = 0.7;
    const int NS = 4096; const long NREP = 2000;
    struct rp_kerr_geodesic_state *S = malloc(NS * sizeof *S);
    for (int i = 0; i < NS; ++i) {           /* smooth orbit-like path */
        double ph = 2*M_PI*i/NS;
        S[i].x[0] = 0; S[i].x[1] = 300.0 + 3e4*(0.5 - 0.5*cos(ph)); S[i].x[2] = 1.2 + 0.2*sin(ph); S[i].x[3] = ph;
        S[i].u[0] = 1.0; S[i].u[1] = 1e-3*sin(ph); S[i].u[2] = 1e-6; S[i].u[3] = 1e-5;
    }
    disk_t D; double c = cos(0.6), sn = sin(0.6);
    double Rm[3][3] = {{1,0,0},{0,c,sn},{0,-sn,c}};
    for (int i=0;i<3;++i) for (int j=0;j<3;++j) D.R[i][j]=Rm[i][j];
    D.d0 = 0; D.idd = (ND-1)/4e4; D.z0 = 0; D.idz = (NZ-1)/4e4; D.mscale = 2.3e-4;
    D.ad = malloc(sizeof(double)*ND*NZ); D.az = malloc(sizeof(double)*ND*NZ);
    for (int iz=0; iz<NZ; ++iz) for (int id=0; id<ND; ++id) {
        double d = (id+1)/D.idd, z = (iz+0.5)/D.idz, r2 = d*d+z*z;
        D.ad[(size_t)iz*ND+id] = -d/pow(r2,1.5); D.az[(size_t)iz*ND+id] = -z/pow(r2,1.5);
    }
    D.mode = 0; D.tab = malloc(sizeof(double)*2*ND*NZ);
    for (size_t k = 0; k < (size_t)ND*NZ; ++k) { D.tab[2*k] = D.ad[k]; D.tab[2*k+1] = D.az[k]; }
    double dx[4], du[4], sink = 0, t0, tb, tp, te;
    /* A: optimized RHS only (sin/cos theta precomputed elsewhere) */
    for (int i = 0; i < NS; ++i) rp_kerr_geodesic_physic_prepare_state(M, s, &S[i]);
    t0 = now();
    for (long k = 0; k < NREP; ++k) for (int i = 0; i < NS; ++i) {
        rp_kerr_geodesic_physic_evaluate_rhs(M, s, &S[i], dx, du); sink += du[1];
    }
    tb = (now() - t0) / (NREP*(double)NS) * 1e9;
    /* B: prepare_state (sin/cos theta, checks) + RHS, as in the public RHS */
    t0 = now();
    for (long k = 0; k < NREP; ++k) for (int i = 0; i < NS; ++i) {
        rp_kerr_geodesic_physic_prepare_state(M, s, &S[i]);
        rp_kerr_geodesic_physic_evaluate_rhs(M, s, &S[i], dx, du); sink += du[1];
    }
    tp = (now() - t0) / (NREP*(double)NS) * 1e9;
    /* C: B + projected disk force */
    t0 = now();
    for (long k = 0; k < NREP; ++k) for (int i = 0; i < NS; ++i) {
        rp_kerr_geodesic_physic_prepare_state(M, s, &S[i]);
        rp_kerr_geodesic_physic_evaluate_rhs(M, s, &S[i], dx, du);
        add_external(M, s, &S[i], &D, du); sink += du[1];
    }
    te = (now() - t0) / (NREP*(double)NS) * 1e9;
    printf("RHS only            %6.1f ns/call\n", tb);
    printf("prepare + RHS       %6.1f ns/call\n", tp);
    printf("prepare + RHS + disk%6.1f ns/call  (+%.0f%%; disk term %.1f ns)\n", te, 100*(te/tp-1), te-tp);

    for (int mode = 1; mode <= 3; ++mode) {
        D.mode = mode; t0 = now();
        for (long k = 0; k < NREP; ++k) for (int i = 0; i < NS; ++i) {
            rp_kerr_geodesic_physic_prepare_state(M, s, &S[i]);
            rp_kerr_geodesic_physic_evaluate_rhs(M, s, &S[i], dx, du);
            add_external(M, s, &S[i], &D, du); sink += du[1];
        }
        double tm = (now() - t0) / (NREP*(double)NS) * 1e9;
        printf("%s %6.1f ns/call  (+%.0f%%)\n", mode == 1 ? "  + disk table joint " : mode == 2 ? "  + Miyamoto-Nagai   " : "  + frame/proj only  ", tm, 100*(tm/tp-1));
    }
    printf("sink %g\n", sink);
    return 0;
}
