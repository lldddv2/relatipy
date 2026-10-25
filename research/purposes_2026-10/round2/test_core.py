"""Unit checks of core.py: mirror trick, harmonic map, ut, Hermite accuracy, PN vs exact."""
import sys, time, numpy as np
sys.path.insert(0, "/home/user/relatipy/research/purposes_2026-10/round2")
import core as C
from astropy import units as u
from astropy.constants import c
R = np.pi/180
p = np.zeros(18); p[:8] = [4.297, 8.277, 125.1, 0.8843, 134.567*R, 228.171*R, 66.263*R, 2010.354]
p[C.IX['fsp']] = 1; p[C.IX['f']] = 1
for spin in (0.0, 0.7):
    p[C.IX['ast']] = spin; p[C.IX['ths']] = 1.1; p[C.IX['phs']] = 2.0
    a_g, ty, us = C.geo(p)
    tlo, thi = -18.5/ty, 15.5/ty
    t0 = time.perf_counter()
    T, X, V, Rb, s, UT = C.exact_traj(p, tlo, thi, npo=4000, return_ut=True)
    t1 = time.perf_counter()
    T2, X2, V2, _, _ = C.exact_traj(p, tlo, thi, npo=16000)
    # u^t from metric vs relatipy samples
    ut_m = np.array([C.ut_kerr(X[i], V[i], Rb, s) for i in range(0, len(T), 50)])
    dut = np.abs(ut_m - UT[::50]).max()
    # hermite accuracy: evaluate coarse interpolant at the fine samples
    idx = np.arange(1, len(T2)-1, 7)
    err_x = max(np.linalg.norm(C._herm(T, X, V, T2[i])[0] - X2[i]) for i in idx)
    err_v = max(np.linalg.norm(C._herm(T, X, V, T2[i])[1] - V2[i]) for i in idx)*C.C_KMS
    # mirror check: restart forward from earliest backward sample, integrate to t=0
    from relatipy import Kerr
    bh = C._kerr(p[0], s); rg = bh.r_g; tg = rg/c
    xb = Rb.T @ X[0]; vb = Rb.T @ V[0]; r = np.linalg.norm(xb); th = np.arccos(xb[2]/r); ph = np.arctan2(xb[1], xb[0])
    rd = xb@vb/r; thd = -(vb[2]-np.cos(th)*rd)/(r*np.sin(th)); phd = (xb[0]*vb[1]-xb[1]*vb[0])/(xb[0]**2+xb[1]**2)
    o = bh.orbit(R=(r+1)*rg, Theta=th*u.rad, Phi=ph*u.rad, vR=rd*c, vTheta=thd*u.rad/tg, vPhi=phd*u.rad/tg)
    taus = np.linspace(0, 1.0, 2)  # placeholder
    sol = o.solve(tau_span=(0*u.s, (-T[0]*1.2)*tg), method="dop853", rtol=1e-12)
    tt = (sol.t/tg).decompose().value; k = np.searchsorted(tt, -T[0])
    st = sol.at(t=(-T[0])*tg)
    Rr = st.R.to_value(rg); Tt = st.Theta.to_value(u.rad); Pp = st.Phi.to_value(u.rad)
    xr = (Rr-1)*np.array([np.sin(Tt)*np.cos(Pp), np.sin(Tt)*np.sin(Pp), np.cos(Tt)])
    xr = Rb @ xr
    i0 = np.argmin(np.abs(T)); 
    print(f"spin {spin}: n={len(T)} exact_traj {1e3*(t1-t0):.0f} ms | max|ut_metric-ut_relatipy|={dut:.1e} | "
          f"Hermite(npo4000) err x {err_x:.1e} rg, v {err_v:.1e} km/s | mirror closure |dx|={np.linalg.norm(xr-X[i0]):.2e} rg (t={T[i0]:.1e})")
    # PN (harmonic 1PN, fsp=1, no spin) vs exact harmonic, same initial state, no fit
    x0, v0 = C.apo_state(a_g, p[3], p[4], p[5], p[6])
    Tp, Xp, Vp = C.pn_traj(x0, v0, tlo, thi, 1.0, np.zeros(3), False, 2*np.pi/3000, 200000)
    tt = np.linspace(tlo*0.98, thi*0.98, 3000)
    dx = np.array([np.linalg.norm(C._herm(Tp, Xp, Vp, t)[0] - C._herm(T, X, V, t)[0]) for t in tt])
    xb2 = np.array([C._herm(Tp, Xp, Vp, t)[0] for t in tt])
    # same with BL-like mapping (r -> r+1) to show the gauge offset is what matters
    print(f"   1PN-harmonic vs exact-harmonic (no fit): max |dx| = {dx.max():.3f} rg = {dx.max()*us:.3f} uas; rms {np.sqrt((dx**2).mean())*us:.3f} uas")
