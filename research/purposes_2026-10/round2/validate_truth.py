"""Validate the truth generator (first-order T: harmonic position + 2M tan(psi/2) lensing, Romer+Shapiro,
1+z = u^t(1+v.gradT)) against EXACT Schwarzschild ray tracing at S2 epochs and S301-like pericentres (a*=0)."""
import sys, numpy as np
sys.path.insert(0, "/home/user/relatipy/research/purposes_2026-10/round2")
import core as C, rays
R = np.pi/180
def run(lab, p, tobs_yr):
    a_g, ty, us = C.geo(p); tg = (np.asarray(tobs_yr) - p[C.IX['t0']])/ty
    T, X, V, Rb, s = C.exact_traj(p, tg.min() - 2*a_g - 100, tg.max() + 2*a_g + 100)
    XE, VE = C.emit(T, X, V, tg, 1.0)
    rows = []
    for k in range(len(tg)):
        x, v = XE[k], VE[k]; ut = C.ut_kerr(x, v, np.eye(3), 0.0)
        ie, ze, Te, ng = rays.exact_obs(x, v); i1, z1, T1 = rays.first_order(x, v, ut)
        rr = np.linalg.norm(x); psi = np.degrees(np.arccos(-x[2]/rr))
        # what the two other models would have: straight harmonic (model i) / straight BL (model ii) image
        istr = np.array([x[1], x[0]]); ibl = istr*(rr + 1)/rr
        rows.append((tobs_yr[k], rr + 1, psi, np.linalg.norm(i1 - ie)*us, np.linalg.norm(istr - ie)*us,
                     np.linalg.norm(ibl - ie)*us, (z1 - ze)*C.C_KMS, Te - T1))
    rows = np.array(rows); rows[:, 7] -= rows[0, 7]
    print(f"\n{lab}:  t_obs     r_BL[M]  psi[deg] | img err [uas]: 1st-order  straight-harm(i)  straight-BL(ii) | RV err 1st-order [km/s] | dT(exact-1st) rel. [M]")
    for r_ in rows:
        print(f"   {r_[0]:10.4f} {r_[1]:9.1f} {r_[2]:7.1f} | {r_[3]:9.4f} {r_[4]:12.2f} {r_[5]:14.2f} | {r_[6]:+9.4f} | {r_[7]:+.2e}")
    return rows
p = np.zeros(18); p[:8] = [4.297, 8.277, 125.1, 0.8843, 134.567*R, 228.171*R, 66.263*R, 2010.354]
p[C.IX['fsp']] = 1; p[C.IX['f']] = 1
run("S2 (a*=0 truth)", p, np.array([2010.35, 2017.5, 2018.2, 2018.33, 2018.37, 2018.39, 2018.42, 2018.5, 2019.0, 2025.0]))
# S301-like: P=8.7yr, rp~280; orientations chosen to put the pericentre at psi = 30..150 deg
for inc_deg, om_deg in ((90, 60), (90, 30), (90, 0), (90, -30), (90, -60), (60, -90), (90, -90)):
    q = p.copy(); q[2] = 83.07; q[3] = 0.98273; q[4] = inc_deg*R; q[5] = 0.7; q[6] = om_deg*R; q[7] = 2035.75
    a_g, ty, us = C.geo(q); Pyr = 2*np.pi*a_g**1.5*ty
    tp = q[7] - Pyr/2
    run(f"S301-like inc={inc_deg} om={om_deg}", q, tp + np.array([-3, -0.3, -0.03, 0.0, 0.03, 0.3, 3])/365.25)
