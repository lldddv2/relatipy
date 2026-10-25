"""Profile-likelihood check of the Fisher sigmas (strong non-linearity found along some directions,
driven by S301 RVs within ~2 weeks of pericentre).  For spin component u (P, Q, L): fix u.S at
truth + k*sigma_F (stiff penalty), minimise chi2 over ALL other parameters (Levenberg-Marquardt,
started from the Fisher conditional mean), report Delta chi2_min (Gaussian expectation k^2)."""
import numpy as np, json, sys, time
from multiprocessing import Pool
import fc, stars as st, forecast as F
el301 = np.array(st.EL['S301A']); dirs, (N, M_, L) = st.spin_dirs(el301)
Pv = fc.rotmat(el301[2], el301[3], el301[4])[:, 0]; Qv = np.cross(L, Pv)
U = {'P': Pv, 'Q': Qv, 'L': L}

def prof(args):
    case, s_rv, lab, k = args
    md = 1.0 if case in ('b', 'c') else 0.0
    stars, g0, els0, gfree, pri = F.build(case, s_rv=s_rv, t_end=2041.0, md=md, orient='CWD', spin=dirs['b45l45'])
    yt = [fc.predict(s, g0, e) for s, e in zip(stars, els0)]; sig = [s.sig() for s in stars]
    th0 = np.concatenate([g0[gfree]] + [np.array(e) for e in els0]); ng = len(gfree)
    J, _ = fc.jacobian(stars, g0, els0, gfree); C, _ = F.cov_from(J, pri)
    u = np.zeros(len(th0)); u[2:5] = U[lab]
    sF = np.sqrt(u @ C @ u)
    fin = np.isfinite(pri); w = np.where(fin, 1 / np.where(fin, pri, 1), 0)
    target = u @ th0 + k * sF; pen = 1e4 / sF
    def unpack(th):
        g = g0.copy(); g[gfree] = th[:ng]
        return g, [th[ng + 6 * q: ng + 6 * q + 6] for q in range(len(stars))]
    def resid(th):
        g, els = unpack(th)
        r = np.concatenate([(fc.predict(s, g, e) - y) / sg for s, e, y, sg in zip(stars, els, yt, sig)])
        return np.concatenate([r, w[fin] * (th - th0)[fin], [pen * (u @ th - target)]])
    th = th0 + C @ u / (u @ C @ u) * k * sF          # Fisher conditional mean (linear prediction)
    r = resid(th); chi2_lin = (r[:-1] @ r[:-1]); chi2 = r @ r; lam = 1e-3
    for it in range(80):
        g, els = unpack(th)
        Jd, _ = fc.jacobian(stars, g, els, gfree)
        Jp = np.vstack([Jd, np.diag(w)[fin], pen * u[None, :]])
        A = Jp.T @ Jp; b = -Jp.T @ r; d = np.sqrt(np.diag(A))
        ok = False
        while lam < 1e8:
            step = np.linalg.solve(A / np.outer(d, d) + lam * np.eye(len(d)), b / d) / d
            rn = resid(th + step); c2 = rn @ rn
            if c2 < chi2: dchi = chi2 - c2; th = th + step; r = rn; chi2 = c2; lam = max(lam / 10, 1e-12); ok = True; break
            lam *= 10
        if not ok or dchi < 1e-5 * max(chi2, 1): break
    return dict(case=case, s_rv=s_rv or 0, comp=lab, k=k, sigF=float(sF), chi2_linear=float(chi2_lin), chi2_profile=float(r[:-1] @ r[:-1]), iters=it + 1)

if __name__ == '__main__':
    t0 = time.time()
    jobs = [(case, s_rv, lab, k) for case, s_rv in (('c', 1.0), ('a', 1.0), ('c', 5.0)) for lab in ('P', 'Q', 'L') for k in (1.0,)]
    with Pool(4) as p: R = p.map(prof, jobs, chunksize=1)
    json.dump(R, open('profile_check.json', 'w'), indent=1)
    for r in R:
        print(f"case {r['case']} RV {r['s_rv']:g} comp {r['comp']} at +{r['k']:.0f} sigma_F ({r['sigF']:.3f}): chi2 along Fisher line {r['chi2_linear']:.1f}, profile min {r['chi2_profile']:.3f} (Gaussian: {r['k']**2:.1f}) it={r['iters']}", flush=True)
    print(f"[{time.time()-t0:.0f}s]")
