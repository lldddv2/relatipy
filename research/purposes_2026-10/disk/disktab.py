"""Build / cache the unit-disk acceleration table (rout = 1, M_disk = 1, G = 1)."""
import os, time
import numpy as np
import diskdyn as dd

HERE = os.path.dirname(os.path.abspath(__file__))


def moments(rin, alpha, n):
    """<b^n> for the unit-mass disk, rout = 1."""
    if abs(alpha + 2) < 1e-12:
        C = 1 / np.log(1 / rin)
    else:
        C = (alpha + 2) / (1 - rin ** (alpha + 2))
    return C * (1 - rin ** (n + alpha + 2)) / (n + alpha + 2)


def make_table(rin=0.0, alpha=-1.0, soft=0.01, q0=0.003, dw=0.02, b0=0.002, dv=0.04,
               rtmin=5e-4, rtmax=8.0, ngl=20, cache=True, verbose=True):
    tag = f"tab2_a{alpha:+.2f}_rin{rin:.4g}_s{soft:g}_q{q0:g}_dw{dw:g}_b{b0:g}_dv{dv:g}_n{ngl}.npz"
    fn = os.path.join(HERE, "tables", tag)
    wlo = np.arcsinh(np.log(rtmin) / q0); whi = np.arcsinh(np.log(rtmax) / q0)
    nw = int(np.ceil((whi - wlo) / dw)) + 5
    w0 = wlo - 2 * dw
    vmax = np.arcsinh((np.pi / 2) / b0)
    nv = int(np.ceil(2 * vmax / dv)) + 5
    v0 = -vmax - 2 * dv
    if cache and os.path.exists(fn):
        d = np.load(fn); TR, TZ = d["TR"], d["TZ"]
    else:
        xg, wg = np.polynomial.legendre.leggauss(ngl)
        t0 = time.time()
        TR, TZ = dd.build_table(nw, w0, dw, nv, v0, dv, q0, b0, rin, alpha, soft, xg, wg)
        if verbose:
            print(f"[table {nw}x{nv} built in {time.time()-t0:.1f}s]", flush=True)
        os.makedirs(os.path.dirname(fn), exist_ok=True)
        np.savez(fn, TR=TR, TZ=TZ)
    J2 = -0.5 * moments(rin, alpha, 2); J4 = 0.375 * moments(rin, alpha, 4)
    H2 = soft ** 2 * moments(rin, alpha, 2)
    tp = np.array([w0, dw, v0, dv, q0, b0, rtmax, J2, J4, H2, rtmin])
    return TR, TZ, tp
