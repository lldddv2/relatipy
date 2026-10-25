"""Interpolation accuracy of the disk table vs direct quadrature; multipole continuity; symmetry."""
import numpy as np, diskdyn as dd, disktab, sys
dw = float(sys.argv[1]) if len(sys.argv) > 1 else 0.02
dv = float(sys.argv[2]) if len(sys.argv) > 2 else 0.04
TR, TZ, tp = disktab.make_table(dw=dw, dv=dv)
xg, wg = np.polynomial.legendre.leggauss(20)
rng = np.random.default_rng(1)
def err_set(r, beta):
    e = np.empty(len(r))
    for k in range(len(r)):
        R = r[k] * np.cos(beta[k]); z = r[k] * np.sin(beta[k])
        gd = np.array(dd.disk_g_direct(R, z, 0.0, 1.0, -1.0, 0.01, xg, wg))
        sR, sZ = dd.disk_T(r[k], z, R, TR, TZ, tp); gt = np.array([sR, sZ]) / r[k]
        e[k] = np.linalg.norm(gt - gd) / np.linalg.norm(gd)
    return e
N = 4000
sets = {
 "global": (np.exp(rng.uniform(np.log(100/31740), np.log(7.9), N)), rng.uniform(-np.pi/2, np.pi/2, N)),
 "near-plane |b|<0.05": (np.exp(rng.uniform(np.log(100/31740), np.log(7.9), N)), rng.uniform(-0.05, 0.05, N)),
 "edge r=0.9-1.1,|b|<0.05": (rng.uniform(0.9, 1.1, N), rng.uniform(-0.05, 0.05, N)),
}
for k, (r, b) in sets.items():
    e = err_set(r, b)
    print(f"{k:<26} rel.err median {np.median(e):.1e}  p99 {np.percentile(e,99):.1e}  max {e.max():.1e}")
# multipole continuity at rtmax and far-field monopole limit
for rt in [7.9, 8.1, 20.0, 100.0]:
    for beta in [0.0, 0.4, 1.2]:
        R = rt*np.cos(beta); z = rt*np.sin(beta)
        gd = np.array(dd.disk_g_direct(R, z, 0.0, 1.0, -1.0, 0.01, xg, wg))
        sR, sZ = dd.disk_T(rt, z, R, TR, TZ, tp); gt = np.array([sR, sZ]) / rt
        mono = -np.array([np.cos(beta), np.sin(beta)]) / rt**2
        print(f"r={rt:6.1f} beta={beta:.1f}: |tab-direct|/|g|={np.linalg.norm(gt-gd)/np.linalg.norm(gd):.1e}  "
              f"|direct-monopole|/|g|={np.linalg.norm(gd-mono)/np.linalg.norm(gd):.1e}")
# symmetry z -> -z (table and full 3D with random disk normal; n_d -> -n_d invariance)
s_err = 0.0
for k in range(200):
    rt = np.exp(rng.uniform(np.log(0.003), np.log(7.9))); beta = rng.uniform(0, np.pi/2)
    R = rt*np.cos(beta); z = rt*np.sin(beta)
    a = np.array(dd.disk_T(rt, z, R, TR, TZ, tp)); b = np.array(dd.disk_T(rt, -z, R, TR, TZ, tp)); b[1] *= -1
    s_err = max(s_err, np.linalg.norm(a-b)/np.linalg.norm(a))
nd = np.array([0.3, -0.5, 0.8]); nd /= np.linalg.norm(nd); x = np.array([4000., -2500., 9000.])
a1 = dd.disk_acc(x, nd, 1e-4, 31740., TR, TZ, tp); a2 = dd.disk_acc(x, -nd, 1e-4, 31740., TR, TZ, tp)
print(f"symmetry: max rel |g(z)-Pg(-z)| = {s_err:.1e};  |a(nd)-a(-nd)|/|a| = {np.linalg.norm(a1-a2)/np.linalg.norm(a1):.1e}")
