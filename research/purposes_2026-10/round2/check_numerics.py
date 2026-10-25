"""Numerical convergence of the three models on the S2 and S301 data sets (max |diff| between settings)."""
import sys, time, numpy as np
sys.path.insert(0, "/home/user/relatipy/research/purposes_2026-10/round2")
import core as C, taskA
d = taskA.s2_data(); p = taskA.p_s2(0.7)
def cmp(lab, f1, f2):
    t0 = time.perf_counter(); a = f1(); t1 = time.perf_counter(); b = f2()
    da = max(np.abs(a[0]-b[0]).max(), np.abs(a[1]-b[1]).max()); dv = np.abs(a[2]-b[2]).max()
    print(f"{lab:<44} max d_ast {da:.2e} uas  max d_rv {dv:.2e} km/s  [{1e3*(t1-t0):.1f} ms/eval]", flush=True)
for _ in range(2): C.model_PN(p, d)
cmp("S2 model(i) spo 3000 vs 12000", lambda: C.model_PN(p, d, spo=3000), lambda: C.model_PN(p, d, spo=12000))
cmp("S2 model(ii) rtol 1e-11 vs 1e-13", lambda: C.model_R(p, d, rtol=1e-11), lambda: C.model_R(p, d, rtol=1e-13))
cmp("S2 model(iii) hw .01 vs .003", lambda: C.model_T(p, d, hw=0.01), lambda: C.model_T(p, d, hw=0.003))
