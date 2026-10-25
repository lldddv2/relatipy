"""Validation (ii-c): full-orbit nodal precession from the disk vs Piran+ 2026 Table 1.
Newton + disk only (alpha=-1, rin=0, rout=31740 r_g, H=0.01 R', M_disk=1e3, M_BH=4.3e6), orbit set
up in the disk frame with inclination i and argument of pericentre omega; dOmega per radial
period (apocentre to apocentre, 3 orbits averaged), compared with LT 9.6e-4/sin i (S301), 3.3e-5/sin i (S2)."""
import numpy as np, diskdyn as dd, fc
D = fc.DEG
md = 1e3 / 4.3e6; rout = 31740.0; nd = np.array([0, 0, 1.0])
stars = {'S301': (16010.0, 0.9825, 9.6e-4), 'S2': (24100.0, 0.8844, 3.3e-5)}
piran = {(10, 0): (0.41, 13), (10, 90): (2.5, 19), (30, 0): (0.36, 8.9), (30, 90): (1.4, 14), (60, 0): (0.19, 4.2), (60, 90): (0.61, 6.6)}

def dOmega(a, e, inc, om, md_, S=np.zeros(3), pn=0, norb=3, nstep=4000):
    x0, v0 = fc.kepler_state(a, e, inc, 0.0, om, np.pi)
    P = 2 * np.pi * a ** 1.5
    T, X, V, A = dd.integrate(x0, v0, 0.0, (norb + 0.02) * P, fc.s_orbit(e) / nstep, nstep * (norb + 2), S, md_, nd, rout, fc.TR, fc.TZ, fc.TP, pn)
    r = np.linalg.norm(X, axis=1)
    idx = [k for k in range(1, len(r) - 1) if r[k] >= r[k - 1] and r[k] >= r[k + 1]]
    # refine apocentre states by dense output on a fine grid around each sample
    Om = []; incs = []
    for k in [0] + idx:
        ts = np.linspace(T[max(k - 1, 0)], T[min(k + 1, len(T) - 1)], 201)
        xs, vs = dd.states_at(T, X, V, A, ts); j = np.argmax(np.linalg.norm(xs, axis=1))
        L = np.cross(xs[j], vs[j]); L /= np.linalg.norm(L)
        Om.append(np.arctan2(L[0], -L[1])); incs.append(np.arccos(L[2]))
    Om = np.unwrap(Om)
    return np.diff(Om).mean(), np.diff(incs).mean()

print(f"{'star':5s} {'i':>3s} {'w':>3s}  dOmega_disk[rad/orb]  ratio(ours)  ratio(Piran T1)   di_disk")
for nm, (a, e, lt) in stars.items():
    for (i, w), (r301, r2) in piran.items():
        dO, di = dOmega(a, e, i * D, w * D, md)
        ratio = abs(dO) / (lt / np.sin(i * D))
        print(f"{nm:5s} {i:3d} {w:3d}  {dO:+.3e}            {ratio:6.2f}       {(r301 if nm=='S301' else r2):6.2f}       {di:+.1e}")
# LT check with the same machinery: chi=1, spin along the disk normal (node precession about nd)
for nm, (a, e, lt) in stars.items():
    dO, di = dOmega(a, e, 30 * D, 90 * D, 0.0, S=np.array([0, 0, 1.0]), pn=3)
    print(f"LT check {nm}: dOmega(chi=1, spin || disk normal) = {dO:.4e}  vs 4 pi/p^1.5 = {4*np.pi/(a*(1-e*e))**1.5:.4e}")
