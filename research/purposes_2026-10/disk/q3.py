"""Q3: bias of the spin when the truth contains a disk but the model has none.
Noise-free truth (disk 1e3 / 3e2 Msun; CWD, two random orientations R1, R2), fitted by
Levenberg-Marquardt with the no-disk model (M, R0, S, elements; same priors, centred on truth).
Bias reported in units of the no-disk Fisher sigma; chi2_min = misfit (detectability)."""
import numpy as np, json, time, sys
from multiprocessing import Pool
import fc, stars as st, forecast as F
D = fc.DEG
el301 = np.array(st.EL['S301A']); dirs, (N, M_, L) = st.spin_dirs(el301)
Pv = fc.rotmat(el301[2], el301[3], el301[4])[:, 0]; Qv = np.cross(L, Pv)
SPIN = dirs['b45l45']

def fit(args):
    case, md, orient, s_rv, t_end = args
    stars, g_true, els_true, _, _ = F.build('c' if case == 'c' else 'b', s_rv=s_rv, t_end=t_end, md=md, orient=orient, spin=SPIN)
    _, g0, _, gfree, pri = F.build('c0' if case == 'c' else 'a', s_rv=s_rv, t_end=t_end, spin=SPIN)   # no-disk model layout
    ytrue = [fc.predict(s, g_true, e) for s, e in zip(stars, els_true)]
    sig = [s.sig() for s in stars]
    theta0 = np.concatenate([g0[gfree]] + [np.array(e) for e in els_true])
    def unpack(th):
        g = g0.copy(); g[gfree] = th[:len(gfree)]
        return g, [th[len(gfree) + 6 * k: len(gfree) + 6 * k + 6] for k in range(len(stars))]
    fin = np.isfinite(pri); w = np.where(fin, 1 / np.where(fin, pri, 1), 0)
    def resid(th):
        g, els = unpack(th)
        r = np.concatenate([(fc.predict(s, g, e) - y) / sg for s, e, y, sg in zip(stars, els, ytrue, sig)])
        return np.concatenate([r, w[fin] * (th - theta0)[fin]])
    th = theta0.copy(); lam = 1e-3
    r = resid(th); chi2 = r @ r; chi2_start = chi2
    for it in range(30):
        g, els = unpack(th)
        J, _ = fc.jacobian(stars, g, els, gfree)
        Jp = np.vstack([J, np.diag(w)[fin]])
        A = Jp.T @ Jp; b = -Jp.T @ r
        d = np.sqrt(np.diag(A))
        while True:
            step = np.linalg.solve(A / np.outer(d, d) + lam * np.eye(len(d)), b / d) / d
            rn = resid(th + step); c2 = rn @ rn
            if c2 < chi2: th = th + step; r = rn; dchi = chi2 - c2; chi2 = c2; lam = max(lam / 10, 1e-9); break
            lam *= 10
            if lam > 1e6: dchi = 0; break
        if dchi < 1e-4: break
    g, els = unpack(th)
    J, _ = fc.jacobian(stars, g, els, gfree); C, _ = F.cov_from(J, pri)
    CS = C[2:5, 2:5]; dS = g[2:5] - g_true[2:5]
    out = dict(case=case, md=md, orient=orient, s_rv=s_rv, t_end=t_end, iters=it + 1, chi2_start=float(chi2_start), chi2_min=float(chi2),
               ndata=int(sum(s.nobs() for s in stars)))
    for lab, u in [('P', Pv), ('Q', Qv), ('N', N), ('M', M_), ('L', L)]:
        out['b' + lab] = float(u @ dS); out['bs' + lab] = float(u @ dS / np.sqrt(u @ CS @ u))
    out['bchi'] = float(np.linalg.norm(g[2:5]) - 0.7)
    out['bM6_sig'] = float((g[0] - g_true[0]) / np.sqrt(C[0, 0]))
    return out

if __name__ == '__main__':
    t0 = time.time()
    jobs = [(case, md, orient, s_rv, t_end) for case in ('a', 'c') for md in (1.0, 0.3) for orient in ('CWD', 'R1', 'R2')
            for (s_rv, t_end) in ((1.0, 2041.0), (None, 2041.0), (1.0, 2035.0))]
    with Pool(4) as p:
        R = p.map(fit, jobs, chunksize=1)
    json.dump(R, open('q3.json', 'w'), indent=1)
    for r in R:
        print(f"{r['case']} Md={r['md']*1e3:5.0f} {r['orient']:3s} RV={r['s_rv'] or 0:g} t_end={r['t_end']:.0f}: bias/sigma P={r['bsP']:+6.2f} Q={r['bsQ']:+6.2f} L={r['bsL']:+6.2f} "
              f"(abs P={r['bP']:+.2f} Q={r['bQ']:+.2f} L={r['bL']:+.2f}) dchi={r['bchi']:+.2f} M6/sig={r['bM6_sig']:+.1f} chi2_min={r['chi2_min']:.1f} (start {r['chi2_start']:.0f}, N={r['ndata']}) it={r['iters']}", flush=True)
    print(f"[{time.time()-t0:.0f}s]")
