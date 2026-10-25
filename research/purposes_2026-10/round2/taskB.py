"""TASK B: S301 (P~8.7 yr, e~0.983, r_p~280 r_g), a*=0.7, 20 isotropic orbit orientations + spin axes.
Data 2026-2040.75: 8 epochs/yr (season) + dense campaign at each pericentre (nightly +-15 d + 10 within +-1 d);
astrometry 100 or 30 uas, RV 1 km/s, same epochs.  Asimov data from the first-order-T truth.
Models: (i) 1PN harm + 1.5PN SO (spin vector free), f_SP=1 [and (i') f_SP free]; (ii) relatipy, spin free;
(iii) exact + 1st-order T: Fisher at truth (RV-only / astrometry-only / both) + pipeline fit (k=0).
Item 6: post-fit RV residuals at the 2031 pericentre; first-order-T truth vs exact Schwarzschild rays there.
Usage: python taskB.py k1 k2 ...   (orientation indices)"""
import sys, time, json, numpy as np
sys.path.insert(0, "/home/user/relatipy/research/purposes_2026-10/round2")
import core as C, rays
R = np.pi/180; YR_D = 365.25
OUT = "/home/user/relatipy/research/purposes_2026-10/round2/taskB/"
BASE = ['M', 'R0', 'a', 'e', 'inc', 'Om', 'om', 't0', 'x0', 'y0', 'vx0', 'vy0', 'vz0']
SPIN = ['ast', 'ths', 'phs']

def p_s301(k):
    rng = np.random.default_rng(1000 + k)
    cosi = rng.uniform(-1, 1); Om = rng.uniform(0, 2*np.pi); om = rng.uniform(0, 2*np.pi)
    v = rng.normal(size=3); v /= np.linalg.norm(v)
    p = np.zeros(18)
    p[:8] = [4.297, 8.277, 83.07, 0.98273, np.arccos(cosi), Om, om, 2035.75]
    p[C.IX['fsp']] = 1.0; p[C.IX['f']] = 1.0
    p[C.IX['ast']] = 0.7; p[C.IX['ths']] = np.arccos(v[2]); p[C.IX['phs']] = np.arctan2(v[1], v[0])
    return p

def peri_obs_times(p):
    a_g, ty, us = C.geo(p)
    T, X, V, Rb, s = C.exact_traj(p, (2025.5 - p[7])/ty, (2041.5 - p[7])/ty)
    r = np.linalg.norm(X, axis=1); out = []
    for i in range(1, len(T) - 1):
        if r[i] < r[i-1] and r[i] <= r[i+1]:
            # refine with the Hermite interpolant (golden section on r(t))
            lo, hi = T[i-1], T[i+1]
            for _ in range(80):
                m1 = lo + 0.382*(hi - lo); m2 = lo + 0.618*(hi - lo)
                if np.linalg.norm(C._herm(T, X, V, m1)[0]) < np.linalg.norm(C._herm(T, X, V, m2)[0]): hi = m2
                else: lo = m1
            te = 0.5*(lo + hi); x, v = C._herm(T, X, V, te)
            q = np.linalg.norm(x) - x[2]
            tob = te + x[2] - 2*np.log(q/C.RREF)
            out.append((p[7] + tob*ty, np.linalg.norm(x) + 1.0, np.degrees(np.arccos(-x[2]/np.linalg.norm(x)))))
    return out

def s301_data(p, sig_ast):
    per = peri_obs_times(p)
    t = [y + f for y in range(2026, 2041) for f in np.linspace(0.15, 0.75, 8)]
    for tp, _, _ in per:
        if 2026.0 < tp < 2040.8:
            t += list(tp + np.arange(-15, 16)/YR_D)
            t += list(tp + np.array([-0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9])/YR_D)
    t = np.unique(np.round(np.array(t), 9)); t = t[(t > 2026.0) & (t < 2040.8)]
    return C.Data(t, sig_ast, t, 1.0, 2033.0), per

def summary(lab, fun, pt, pb, free, d, res, tp1):
    J = C.jac(fun, pb, free, d); sig, _ = C.fisher_sigma(J)
    na = d.na; ra = np.sqrt(np.mean((res[:2*na]*np.concatenate([d.s_ast, d.s_ast]))**2))
    rvres = res[2*na:]*d.s_rv; rr = np.sqrt(np.mean(rvres**2))
    w = np.abs(d.t_rv - tp1) < 1.0/YR_D
    out = dict(model=lab, chi2=float(res @ res), rms_ast=ra, rms_rv=rr, rv_peri_max=float(np.abs(rvres[w]).max()),
               rv_peri_rms=float(np.sqrt(np.mean(rvres[w]**2))), x=pb.tolist(), sig=dict(zip(free, sig.tolist())))
    for n in ('M', 'R0', 'ast', 'fsp'):
        if n in free:
            kk = free.index(n); b = pb[C.IX[n]] - pt[C.IX[n]]
            out['bias_' + n] = b; out['bias_sig_' + n] = b/sig[kk]; out['sig_' + n] = sig[kk]
    # spin axis error (deg) if spin free
    if 'ths' in free:
        a1 = C.axis(pb[C.IX['ths']], pb[C.IX['phs']]); a0 = C.axis(pt[C.IX['ths']], pt[C.IX['phs']])
        out['axis_err_deg'] = float(np.degrees(np.arccos(np.clip(a1 @ a0, -1, 1))))
    return out

def run_fit(lab, fun, pt, free, d, y, tp1, bounds=None, max_nfev=60):
    t0 = time.perf_counter()
    pb, r = C.fit(fun, pt.copy(), free, d, y, max_nfev=max_nfev, bounds=bounds)
    o = summary(lab, fun, pt, pb, free, d, r.fun, tp1); o['nfev'] = r.nfev; o['sec'] = time.perf_counter() - t0
    return o

def fisher_iii(pt, d):
    """sigma(a*) etc. from the exact model at truth for RV-only / astrometry-only / both, 100 & 30 uas."""
    d1 = C.Data(d.t_ast, 1.0, d.t_rv, 1.0, d.t_ref)
    free = BASE + SPIN
    J = C.jac(C.model_T, pt, free, d1)              # raw (unwhitened) Jacobian
    na = d.na; out = {}
    for sa in (100.0, 30.0):
        w = np.concatenate([np.full(2*na, 1/sa), np.full(len(d.t_rv), 1.0)])
        Jw = J*w[:, None]
        for lab, rows, cols in (("both", slice(None), free),
                                ("ast", slice(0, 2*na), [n for n in free if n != 'vz0']),
                                ("rv", slice(2*na, None), [n for n in free if n not in ('R0', 'x0', 'y0', 'vx0', 'vy0', 'Om')])):
            ci = [free.index(n) for n in cols]
            sig, Cv = C.fisher_sigma(Jw[rows][:, ci])
            o = dict(zip(cols, sig.tolist()))
            # spin-axis 1-sigma cone (deg) from (ths, phs) covariance
            ith, iph = cols.index('ths'), cols.index('phs')
            ths = pt[C.IX['ths']]
            cc = Cv[np.ix_([ith, iph], [ith, iph])]*np.array([[1, np.sin(ths)], [np.sin(ths), np.sin(ths)**2]])
            o['axis_deg'] = float(np.degrees(np.sqrt(max(np.trace(cc), 0))))
            out[f"{lab}_{int(sa)}"] = o
    return out

def exact_ray_floor(pt, d, tp1):
    """RV / astrometry error of the first-order-T truth vs exact Schwarzschild rays at the 2031 campaign."""
    a_g, ty, us = C.geo(pt)
    sel = np.abs(d.t_rv - tp1) < 15.5/YR_D; tt = d.t_rv[sel]
    tg = (tt - pt[7])/ty
    T, X, V, Rb, s = C.exact_traj(pt, tg.min() - 3*a_g, tg.max() + 3*a_g)
    XE, VE = C.emit(T, X, V, tg, 1.0)
    dz = []; di = []
    for k in range(len(tt)):
        x, v = XE[k], VE[k]; ut = C.ut_kerr(x, v, Rb, s)
        try:
            ie, ze, _, _ = rays.exact_obs(x, v); i1, z1, _ = rays.first_order(x, v, ut)
            # exact-ray 1+z uses the Schwarzschild static frame; use the Kerr u^t ratio to keep dynamics exact
            ut_s = C.ut_kerr(x, v, np.eye(3), 0.0)
            dz.append((z1 - ze*ut/ut_s)*C.C_KMS); di.append(np.linalg.norm(i1 - ie)*us)
        except Exception:
            dz.append(np.nan); di.append(np.nan)
    dz = np.array(dz); di = np.array(di)
    w1 = np.abs(tt - tp1) < 1.0/YR_D
    return dict(rv_max_1d=float(np.nanmax(np.abs(dz[w1]))), rv_max_15d=float(np.nanmax(np.abs(dz))),
                ast_max=float(np.nanmax(di)), n_nan=int(np.isnan(dz).sum()))

if __name__ == "__main__":
    import os; os.makedirs(OUT, exist_ok=True)
    for k in [int(a) for a in sys.argv[1:]]:
        t0 = time.perf_counter()
        pt = p_s301(k)
        res = dict(k=k, p=pt.tolist())
        d100, per = s301_data(pt, 100.0); d30 = C.Data(d100.t_ast, 30.0, d100.t_rv, 1.0, d100.t_ref)
        tp1 = [x for x in per if 2030 < x[0] < 2033][0][0]
        res['peri'] = per; res['n_ep'] = int(d100.na)
        y = {}
        for lab, d in (("100", d100), ("30", d30)):
            y[lab] = C.stack(C.model_T(pt, d), d)
        res['fisher_iii'] = fisher_iii(pt, d100)
        res['floor_iii'] = exact_ray_floor(pt, d100, tp1)
        fits = []
        bnd_R = (np.r_[np.full(13, -np.inf), 0.0, -np.inf, -np.inf], np.r_[np.full(13, np.inf), 0.9999, np.inf, np.inf])
        for lab, d in (("100", d100), ("30", d30)):
            f = run_fit("(i) 1PN+SO fSP=1", lambda p, dd: C.model_PN(p, dd, spin_on=True), pt, BASE + SPIN, d, y[lab], tp1)
            f['sig_ast'] = lab; fits.append(f)
            f = run_fit("(i') 1PN+SO fSP free", lambda p, dd: C.model_PN(p, dd, spin_on=True), pt, BASE + SPIN + ['fsp'], d, y[lab], tp1)
            f['sig_ast'] = lab; fits.append(f)
            f = run_fit("(ii) relatipy spin free", C.model_R, pt, BASE + SPIN, d, y[lab], tp1, bounds=bnd_R)
            f['sig_ast'] = lab; fits.append(f)
            if k == 0:
                f = run_fit("(iii) exact+1st-order T", C.model_T, pt, BASE + SPIN, d, y[lab], tp1, max_nfev=6)
                f['sig_ast'] = lab; fits.append(f)
        res['fits'] = fits; res['sec'] = time.perf_counter() - t0
        json.dump(res, open(OUT + f"k{k:02d}.json", "w"), indent=1)
        line = f"k={k:2d} i={np.degrees(pt[4]):5.1f} psi_p1={per[0][2] if per else 0:5.1f} n_ep={d100.na} floor(iii) RV {res['floor_iii']['rv_max_1d']:.2f} km/s | "
        for f in fits:
            line += f"{f['model'][:5]}/{f['sig_ast']}: a* {f.get('bias_sig_ast', np.nan):+.1f}s M {f.get('bias_sig_M', np.nan):+.1f}s rvrms {f['rms_rv']:.2f} | "
        print(line + f"[{res['sec']:.0f}s]", flush=True)
