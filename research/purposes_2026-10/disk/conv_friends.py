"""Integrator convergence (max |dy| vs 4x steps) and timing per star, disk on (1e3 Msun, CWD), chi=0.7."""
import numpy as np, time, fc, stars as st
g = st.G0.copy(); g[2:5] = [0.3, -0.4, 0.5]; g[5] = 1.0
for name, ns in [('S301A', 3000), ('S2', 1500), ('S29', 2000), ('S38', 1200), ('S55', 1200)]:
    mk = (lambda n: st.s301('A', nstep=n)) if name.startswith('S301') else (lambda n, nm=name: st.friend(nm, nstep=n))
    s1, s4 = mk(ns), mk(4 * ns)
    y1 = fc.predict(s1, g, s1.el); t0 = time.time(); y1 = fc.predict(s1, g, s1.el); dt = time.time() - t0
    y4 = fc.predict(s4, g, s4.el); na = len(s1.t_ast)
    print(f"{name:6s} nstep/orbit={ns:5d}  {dt*1e3:6.1f} ms  max|dpos|={np.abs(y1-y4)[:2*na].max():.3f} uas  max|dRV|={np.abs(y1-y4)[2*na:].max():.4f} km/s  n_ast={na} n_rv={len(s1.t_rv)}")
