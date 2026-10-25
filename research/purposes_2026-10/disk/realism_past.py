"""Realism check: element errors from the mock past data (GRAVITY 2017-2025, 50 uas, 10 km/s, 4/yr;
S2 2018 campaign; S301: 19 epochs at 207 uas, no RV) vs published 1-sigma errors.  M, R0, spin, disk fixed."""
import numpy as np, fc, stars as st
D = fc.DEG
g = st.G0.copy()
for name in ['S301A', 'S2', 'S29', 'S38', 'S55']:
    if name == 'S301A':
        el = st.EL[name]; ta = st.S301_PAST
        star = fc.Star(name, el, ta, 207.0, None, None, nstep=3000)
    else:
        s = st.friend(name, t_end=2025.99)
        star = fc.Star(name, s.el, s.t_ast[s.t_ast < 2026], 50.0, s.t_rv[s.t_rv < 2026], 10.0, nstep=s.nstep)
    J, nm = fc.jacobian([star], g, [np.array(star.el)], gfree=[])
    C = np.linalg.inv(J.T @ J); s_m = np.sqrt(np.diag(C))
    pub = np.array(st.EL_ERR[name])
    sc = np.array([1, 1, 1 / D, 1 / D, 1 / D, 1])
    print(f"{name:6s} mock: " + " ".join(f"{x:.2g}" for x in s_m * sc) + "  | published: " + " ".join(f"{x:.2g}" for x in pub * sc))
print("columns: a[mas] e i[deg] Omega[deg] omega[deg] tp[yr]")
