"""Robustness: refit S2 (a*=0 full truth) from 4 random starts (+-1 sigma per parameter) for models (i), (i'), (ii)."""
import sys, numpy as np
sys.path.insert(0, "/home/user/relatipy/research/purposes_2026-10/round2")
import core as C, taskA as A
d = A.s2_data(); pt = A.p_s2(0.0); y = C.stack(C.model_T(pt, d), d)
rng = np.random.default_rng(5)
for lab, fun, free in (("(i)", C.model_PN, A.BASE + ['fsp']), ("(i')", C.model_PN, A.BASE + ['fsp', 'f']), ("(ii)", C.model_R, A.BASE)):
    sig = C.fisher_sigma(C.jac(fun, pt, free, d))[0]
    for s in range(4):
        p0 = pt.copy(); idx = [C.IX[n] for n in free]; p0[idx] += sig*rng.normal(size=len(free))
        pb, r = C.fit(fun, p0, free, d, y, max_nfev=60)
        print(f"{lab} start {s}: chi2 {r.fun@r.fun:.4f} nfev {r.nfev} | M {(pb[0]-pt[0])/sig[0]:+.3f}s R0 {(pb[1]-pt[1])/sig[1]:+.3f}s"
              + (f" fsp {(pb[13]-1)/sig[free.index('fsp')]:+.3f}s" if 'fsp' in free else "") + (f" f {(pb[14]-1)/sig[free.index('f')]:+.3f}s" if 'f' in free else ""), flush=True)
