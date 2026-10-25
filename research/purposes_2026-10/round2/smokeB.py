import sys, time, numpy as np
sys.path.insert(0, "/home/user/relatipy/research/purposes_2026-10/round2")
import core as C, taskB as B
pt = B.p_s301(0); d, per = B.s301_data(pt, 100.0)
print("peri (t_obs, r_BL, psi):", [(round(a,4), round(b,1), round(c,1)) for a,b,c in per], "n_ep", d.na)
for lab, f in (("PN", lambda p, dd: C.model_PN(p, dd, spin_on=True)), ("R", C.model_R), ("T", C.model_T)):
    f(pt, d); t0 = time.perf_counter(); o = f(pt, d); dt = time.perf_counter() - t0
    print(lab, f"{1e3*dt:.0f} ms", np.round(o[2][:3], 2))
yT = C.model_T(pt, d); 
for lab, f in (("PN spo3000", lambda p, dd: C.model_PN(p, dd, spin_on=True)), ("R", C.model_R)):
    o = f(pt, d); w = np.abs(d.t_rv - per[0][0]) < 1/365.25
    print(lab, "raw diff vs truth: ast max %.1f uas, RV max %.1f km/s, RV max(+-1d) %.1f" % (max(np.abs(o[0]-yT[0]).max(), np.abs(o[1]-yT[1]).max()), np.abs(o[2]-yT[2]).max(), np.abs(o[2]-yT[2])[w].max()))
a = C.model_PN(pt, d, spin_on=True, spo=3000); b = C.model_PN(pt, d, spin_on=True, spo=12000)
print("PN conv spo 3000 vs 12000: ast %.2e uas rv %.2e km/s" % (max(np.abs(a[0]-b[0]).max(), np.abs(a[1]-b[1]).max()), np.abs(a[2]-b[2]).max()))
a = C.model_R(pt, d, rtol=1e-11); b = C.model_R(pt, d, rtol=1e-13)
print("R conv rtol 1e-11 vs 1e-13: ast %.2e uas rv %.2e km/s" % (max(np.abs(a[0]-b[0]).max(), np.abs(a[1]-b[1]).max()), np.abs(a[2]-b[2]).max()))
a = C.model_T(pt, d, hw=0.01); b = C.model_T(pt, d, hw=0.003)
print("T conv hw .01 vs .003: ast %.2e uas rv %.2e km/s" % (max(np.abs(a[0]-b[0]).max(), np.abs(a[1]-b[1]).max()), np.abs(a[2]-b[2]).max()))
