"""Unit checks of core.py: mirror trick, harmonic map, u^t, Hermite accuracy, PN vs exact (no fit)."""
import sys, time, numpy as np
sys.path.insert(0, "/home/user/relatipy/research/purposes_2026-10/round2")
import core as C
from astropy import units as u
from astropy.constants import c
R = np.pi/180
def herm_arr(T, X, V, tt):
    o = [C._herm(T, X, V, t) for t in tt]; return np.array([a for a, b in o]), np.array([b for a, b in o])
for lab, a_mas, e, span in (("S2", 125.1, 0.8843, (-18.5, 15.5)), ("S301", 83.07, 0.98273, (-9.8, 5.0))):
  for spin in (0.0, 0.7):
    p = np.zeros(18); p[:8] = [4.297, 8.277, a_mas, e, 134.567*R, 228.171*R, 66.263*R, 2010.354]
    p[C.IX['fsp']] = 1; p[C.IX['f']] = 1; p[C.IX['ast']] = spin; p[C.IX['ths']] = 1.1; p[C.IX['phs']] = 2.0
    a_g, ty, us = C.geo(p); tlo, thi = span[0]/ty, span[1]/ty
    t0 = time.perf_counter(); T, X, V, Rb, s, UT = C.exact_traj(p, tlo, thi, return_ut=True); t1 = time.perf_counter()
    T2, X2, V2, _, _ = C.exact_traj(p, tlo, thi, hw=0.0025)
    dut = max(abs(C.ut_kerr(X[i], V[i], Rb, s) - UT[i]) for i in range(0, len(T), 10))
    tt = np.sort(np.random.default_rng(1).uniform(T[1], T[-2], 4000))
    # add points clustered at pericentres
    rr = np.linalg.norm(X, axis=1); ip = [i for i in range(1, len(T)-1) if rr[i] < rr[i-1] and rr[i] < rr[i+1]]
    for i in ip: tt = np.concatenate([tt, T[i] + np.linspace(-30, 30, 400)*rr[i]**1.5])
    tt = tt[(tt > T[1]) & (tt < T[-2])]
    xa, va = herm_arr(T, X, V, tt); xb, vb = herm_arr(T2, X2, V2, tt)
    ex = np.linalg.norm(xa - xb, axis=1).max(); ev = np.linalg.norm(va - vb, axis=1).max()*C.C_KMS
    # mirror closure: forward integration from the earliest (backward-generated) state to t~0
    bh = C._kerr(p[0], s); rg = bh.r_g; tg = rg/c
    xb0 = Rb.T @ X[0]; vb0 = Rb.T @ V[0]; r = np.linalg.norm(xb0); th = np.arccos(xb0[2]/r); ph = np.arctan2(xb0[1], xb0[0])
    rd = xb0@vb0/r; thd = -(vb0[2]-np.cos(th)*rd)/(r*np.sin(th)); phd = (xb0[0]*vb0[1]-xb0[1]*vb0[0])/(xb0[0]**2+xb0[1]**2)
    o = bh.orbit(R=(r+1)*rg, Theta=th*u.rad, Phi=ph*u.rad, vR=rd*c, vTheta=thd*u.rad/tg, vPhi=phd*u.rad/tg)
    sol = o.solve(tau_span=(0*u.s, (-T[0]*0.999)*tg), method="dop853", rtol=1e-13, max_step=2000*tg)
    tf = (sol.t[-1]/tg).decompose().value + T[0]
    Rr = sol.R[-1].to_value(rg); Tt = sol.Theta[-1].to_value(u.rad); Pp = sol.Phi[-1].to_value(u.rad)
    xr = Rb @ ((Rr-1)*np.array([np.sin(Tt)*np.cos(Pp), np.sin(Tt)*np.sin(Pp), np.cos(Tt)]))
    xm, _ = C._herm(T, X, V, tf)
    # 1PN harmonic (fsp=1, no spin) vs exact harmonic, same initial state, no fit
    x0, v0 = C.apo_state(a_g, p[3], p[4], p[5], p[6])
    Tp, Xp, Vp = C.pn_traj(x0, v0, tlo, thi, 1.0, np.zeros(3), False, 2*np.pi/3000, 400000)
    t3 = np.linspace(tlo*0.98, thi*0.98, 4000)
    dx = np.linalg.norm(herm_arr(Tp, Xp, Vp, t3)[0] - herm_arr(T, X, V, t3)[0], axis=1)
    print(f"{lab} a*={spin}: nsteps {len(T)} ({1e3*(t1-t0):.0f} ms) |ut_metric-ut_rel|max {dut:.1e} | interp(hw .01 vs .0025): "
          f"x {ex:.1e} rg, v {ev:.1e} km/s | mirror closure {np.linalg.norm(xr-xm):.1e} rg | 1PN-exact(no fit) max {dx.max()*us:.2f} uas", flush=True)
