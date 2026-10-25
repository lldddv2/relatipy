"""Q0: reproduce the GRAVITY (Nature 2026) S301 spin forecast with the Fisher matrix, no disk.
chi = 1 along S301's orbital pole (and 20 deg off), 100 uas + 1 km/s, 10 epochs/yr 2026-2035
+ 10 around pericentre; nuisance: M, R0, 6 elements; spin as Cartesian (Sx, Sy, Sz)."""
import numpy as np, time, sys
import fc, stars as st
D = fc.DEG
el = np.array(st.EL['S301A'])
dirs, (N, M_, L) = st.spin_dirs(el)

def run(label, half=0.1, past=False, rv=True, s_ast=100., s_rv=1., beta=0.0, t_end=2035.0, redshift=True):
    ta = np.concatenate([st.uniform(2026.0, t_end, 10), st.campaign(el, 2026.0, t_end, 10, half)])
    sa = np.full(ta.shape, s_ast)
    if past:
        ta = np.concatenate([st.S301_PAST, ta]); sa = np.concatenate([np.full(19, 207.), sa])
    tr = ta[ta >= 2026.0] if rv else None
    star = fc.Star('S301A', el, ta, sa, tr, None if tr is None else np.full(tr.shape, s_rv), nstep=3000)
    s_hat = np.cos(beta * D) * L + np.sin(beta * D) * N
    g = st.G0.copy(); g[2:5] = 1.0 * s_hat
    J, names = fc.jacobian([star], g, [el], gfree=[0, 1, 2, 3, 4])
    F = J.T @ J; C = np.linalg.inv(F)
    CS = C[2:5, 2:5]
    sig_chi = np.sqrt(s_hat @ CS @ s_hat)
    P = np.eye(3) - np.outer(s_hat, s_hat)
    e1 = N - (N @ s_hat) * s_hat; e1 /= np.linalg.norm(e1); e2 = np.cross(s_hat, e1)
    s1, s2 = np.sqrt(e1 @ CS @ e1), np.sqrt(e2 @ CS @ e2)
    ang = np.degrees(np.sqrt(s1 ** 2 + s2 ** 2))   # chi = 1: total angular error (rad -> deg)
    sM = np.sqrt(C[0, 0]) / g[0]; sR = np.sqrt(C[1, 1]) / g[1]
    print(f"{label:<44} sig(chi)={sig_chi:.3f}  sig_perp=({s1:.3f},{s2:.3f}) -> orient {ang:.0f} deg "
          f"(per-axis {np.degrees(s1):.0f},{np.degrees(s2):.0f})  sM/M={sM:.1e} sR0/R0={sR:.1e}", flush=True)
    return sig_chi, ang

t0 = time.time()
print("nuisance: M, R0, a, e, i, Omega, omega, tp (no offsets, no f_SP, no 2PN); spin Cartesian")
run("baseline: peri +-0.1yr, no past data")
run("peri campaign +-0.3 yr", half=0.3)
run("+ 19 past GRAVITY epochs (207 uas)", past=True)
run("spin 20 deg off the orbital pole", beta=20)
run("astrometry only", rv=False)
run("RV 5 km/s", s_rv=5.)
run("30 uas + 1 km/s", s_ast=30.)
run("baseline, data to 2041 (2 peri)", t_end=2041.0)
print(f"[{time.time()-t0:.0f}s]")
