"""TASK A: S2 with present-day-like data.  Truth = exact Kerr (relatipy) + first-order time transfer
(harmonic position + lensing, Romer+Shapiro, RV with bending).  Asimov data.  Fit models
(i) GRAVITY-like 1PN harmonic x f_SP (+f), (ii) relatipy KerrMcmcModel as is, (iii) exact + first-order T.
Usage: python taskA.py <spin> [attrib]"""
import sys, time, json, numpy as np
sys.path.insert(0, "/home/user/relatipy/research/purposes_2026-10/round2")
import core as C
R = np.pi/180
OUT = "/home/user/relatipy/research/purposes_2026-10/round2/"

def s2_data():
    na = [y + f for y in range(1992, 2017) for f in (0.35, 0.55)]                      # NACO 0.7 mas
    gr = [y + f for y in range(2017, 2026) for f in np.linspace(0.2, 0.7, 10)]         # GRAVITY 50 uas
    gr += list(np.linspace(2018.22, 2018.68, 30))                                         # + 30 in 2018
    rv = [y + f for y in range(2000, 2026) for f in np.linspace(0.25, 0.65, 5)]        # RV 10 km/s
    rv += list(np.linspace(2018.22, 2018.68, 20))
    ta = np.array(na + gr); sa = np.array([700.0]*len(na) + [50.0]*len(gr))
    o = np.argsort(ta); tr = np.sort(np.array(rv))
    return C.Data(ta[o], sa[o], tr, 10.0, 2010.0)

def p_s2(spin, seed=11):
    p = np.zeros(18)
    P = 16.05; p[:8] = [4.297, 8.277, 125.1, 0.8843, 134.567*R, 228.171*R, 66.263*R, 2018.379 - P/2]
    p[C.IX['fsp']] = 1.0; p[C.IX['f']] = 1.0
    if spin > 0:
        rng = np.random.default_rng(seed); v = rng.normal(size=3); v /= np.linalg.norm(v)
        p[C.IX['ast']] = spin; p[C.IX['ths']] = np.arccos(v[2]); p[C.IX['phs']] = np.arctan2(v[1], v[0])
    return p

BASE = ['M', 'R0', 'a', 'e', 'inc', 'Om', 'om', 't0', 'x0', 'y0', 'vx0', 'vy0', 'vz0']
PUB = dict(M=0.012, R0=0.009, fsp=0.12)

def report(lab, fun, pt, free, d, y, p0=None, max_nfev=40):
    t0 = time.perf_counter()
    pb, r = C.fit(fun, pt.copy() if p0 is None else p0.copy(), free, d, y, max_nfev=max_nfev)
    J = C.jac(fun, pb, free, d); sig, _ = C.fisher_sigma(J)
    res = r.fun; na = d.na
    ra = np.sqrt(np.mean((res[:2*na]*np.concatenate([d.s_ast, d.s_ast]))**2))
    gm = np.concatenate([d.s_ast, d.s_ast]) < 100
    rg = np.sqrt(np.mean((res[:2*na][gm]*np.concatenate([d.s_ast, d.s_ast])[gm])**2))
    rr = np.sqrt(np.mean((res[2*na:]*d.s_rv)**2)); chi2 = float(res @ res)
    out = dict(model=lab, chi2=chi2, rms_grav_uas=rg, rms_ast_uas=ra, rms_rv=rr, nfev=r.nfev,
               x=pb.tolist(), sig=dict(zip(free, sig.tolist())))
    s = f"  {lab:<34} chi2={chi2:8.3f} rms GRAV {rg:6.2f} uas, RV {rr:5.2f} km/s |"
    for n in ('M', 'R0', 'fsp', 'f'):
        if n in free:
            k = free.index(n); b = pb[C.IX[n]] - pt[C.IX[n]]
            out['bias_' + n] = b; out['bias_sig_' + n] = b/sig[k]
            s += f" {n}: {b:+.2e} ({b/sig[k]:+.2f}s, s={sig[k]:.1e}" + (f", {b/PUB[n]:+.2f}pub)" if n in PUB else ")")
    print(s + f" [{time.perf_counter() - t0:.0f}s nfev {r.nfev}]", flush=True)
    return out

if __name__ == "__main__":
    spin = float(sys.argv[1]); attrib = len(sys.argv) > 2
    d = s2_data(); pt = p_s2(spin)
    print(f"S2 truth a*={spin} axis(th,ph)=({pt[C.IX['ths']]/R:.1f},{pt[C.IX['phs']]/R:.1f}) deg; n_ast={d.na} (NACO {np.sum(d.s_ast>100)}, GRAV {np.sum(d.s_ast<100)}), n_rv={len(d.t_rv)}", flush=True)
    res = []
    truths = [("full truth", 1.0, 1.0)]
    if attrib: truths += [("truth w/o lensing", 0.0, 1.0), ("truth w/o Shapiro+RV-bending", 1.0, 0.0), ("truth w/o lensing,Shapiro,bending", 0.0, 0.0)]
    for tl, lens, shap in truths:
        y = C.stack(C.model_T(pt, d, lens=lens, shap=shap), d)
        print(f"-- {tl}", flush=True)
        r1 = report("(i) 1PN harm x fSP, f=1", C.model_PN, pt, BASE + ['fsp'], d, y); r1['truth'] = tl; res.append(r1)
        r2 = report("(i') 1PN harm x fSP, f free", C.model_PN, pt, BASE + ['fsp', 'f'], d, y); r2['truth'] = tl; res.append(r2)
        if tl == "full truth":
            r3 = report("(ii) relatipy KerrMcmcModel a*=0", C.model_R, pt, BASE, d, y); r3['truth'] = tl; res.append(r3)
            pz = pt.copy(); pz[C.IX['ast']] = 0.0
            if spin > 0:
                r4 = report("(iii) exact a*=0 + 1st-order T", C.model_T, pt, BASE, d, y, p0=pz, max_nfev=12)
            else:
                r4 = report("(iii) exact + 1st-order T", C.model_T, pt, BASE, d, y, max_nfev=12)
            r4['truth'] = tl; res.append(r4)
    json.dump(res, open(OUT + f"taskA_spin{spin}{'_attrib' if attrib else ''}.json", "w"), indent=1)
