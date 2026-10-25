"""Q2: (1) within-orbit time profile of the orbital-plane change for S301 (LT vs disk);
(2) full-orbit Fisher vs secular-only Fisher for chi_perp.
Secular-only: data -> per-star uniform plane-drift rates (di/dt, dOmega/dt [, domega/dt]) via a
drift model (1PN orbit rigidly rotated: R_Z(Odot tau) R_N(idot tau) R_L(wdot tau)) fitted to the SAME
data with the same nuisances/priors; rates -> (S, Md, i_d, O_d) via orbit-averaged rates computed
from full integrations (apocentre to apocentre, 2 orbits)."""
import numpy as np, json, sys, time
import fc, stars as st, diskdyn as dd, forecast as F
D = fc.DEG
el = np.array(st.EL['S301A']); dirs, (N, M_, L) = st.spin_dirs(el)

# ------------------------------------------------------------------ (1) time profile
def osc_pole(T, X, V, A, ts):
    xs, vs = dd.states_at(T, X, V, A, ts)
    Lh = np.cross(xs, vs); return Lh / np.linalg.norm(Lh, axis=1)[:, None]

def profile(g, pn, label):
    s = fc.Star('x', el, [2027.0, 2037.0], 1.0, nstep=6000)
    P = st.period_yr(el); tp = st.peris(el, 2027.0, 2037.0)[0]
    t_a0, t_a1 = tp - P / 2, tp + P / 2
    T, X, V, A, rgau, tgy = fc.integrate_star(s, g, el, pn=pn, t_ic=t_a0 - 0.01, t_end=t_a1 + 0.01)
    c = lambda t: (np.atleast_1d(t) - fc.T_REF) / tgy
    L0, L1 = osc_pole(T, X, V, A, c([t_a0, t_a1]))
    dL = L1 - L0; u = dL / np.linalg.norm(dL)
    res = {}
    for w in (0.01, 0.05, 0.25):
        La, Lb = osc_pole(T, X, V, A, c([tp - w * P, tp + w * P]))
        res[w] = (Lb - La) @ u / np.linalg.norm(dL)
    # sky-frame node/inclination fractions
    def iO(Lv): return np.arccos(Lv[2]), np.arctan2(Lv[0], -Lv[1])
    i0, O0 = iO(L0); i1, O1 = iO(L1); La, Lb = osc_pole(T, X, V, A, c([tp - 0.05 * P, tp + 0.05 * P]))
    ia, Oa = iO(La); ib, Ob = iO(Lb)
    fO = (Ob - Oa) / (O1 - O0); fi = (ib - ia) / (i1 - i0)
    print(f"{label:<34} |dL|/orbit={np.linalg.norm(dL):.2e} rad  dOmega_sky={O1-O0:+.2e} di_sky={i1-i0:+.2e} | frac of dL within +-1%P {res[0.01]:.3f}, +-5%P {res[0.05]:.3f}, +-25%P {res[0.25]:.3f} | frac dOmega_sky(+-5%P) {fO:.3f}, di_sky {fi:.3f}", flush=True)
    return dict(label=label, dL=float(np.linalg.norm(dL)), dO=float(O1 - O0), di=float(i1 - i0), f01=res[0.01], f05=res[0.05], f25=res[0.25], fO05=float(fO), fi05=float(fi))

# ------------------------------------------------------------------ (2) secular Fisher
def drift_rate_cov(stars, els, pri_full, nrate):
    """Covariance of per-star drift rates [rad/yr] from the drift model (1PN, no spin/disk)."""
    g = st.G0.copy()
    J0, names = fc.jacobian(stars, g, els, [0, 1], pn=1)
    nobs = [s.nobs() for s in stars]; off = np.concatenate([[0], np.cumsum(nobs)])
    Jr = np.zeros((off[-1], nrate * len(stars)))
    for k, s in enumerate(stars):
        T, X, V, A, rgau, tgy = fc.integrate_star(s, g, els[k], pn=1)
        tob = (np.concatenate([s.t_ast, s.t_rv]) - fc.T_REF) / tgy
        out, te = dd.observe(T, X, V, A, tob)
        xs, vs = dd.states_at(T, X, V, A, te)
        tau = te * tgy + fc.T_REF - 0.5 * (s.t_ast.min() + s.t_ast.max())
        ax = [np.array([np.cos(els[k][3]), np.sin(els[k][3]), 0.0]), np.array([0, 0, 1.0]), fc.pole(els[k][2], els[k][3])][:nrate]
        kk = fc.uas_per_rg(g[0], g[1]); na = len(s.t_ast)
        for j, u in enumerate(ax):
            dx = tau[:, None] * np.cross(u, xs)
            dv = tau[:, None] * np.cross(u, vs) + tgy * np.cross(u, xs)
            col = np.concatenate([dx[:na, 0] * kk, dx[:na, 1] * kk, dv[na:, 2] * fc.C_KMS])
            Jr[off[k]:off[k + 1], nrate * k + j] = col / s.sig()
    J = np.hstack([J0, Jr]); pri = np.concatenate([pri_full, np.full(nrate * len(stars), np.inf)])
    C, cond = F.cov_from(J, pri)
    n0 = J0.shape[1]
    return C[n0:, n0:]

def orbit_rates(name, g, nrate):
    """Orbit-averaged sky-frame rates (idot, Odot, wdot) [rad/yr] over 2 radial orbits from apocentre."""
    e0 = np.array(st.EL[name]); P = st.period_yr(e0)
    t_apo = e0[5] + P / 2
    s = fc.Star(name, e0, [t_apo, t_apo + 2 * P], 1.0, nstep=3000)
    T, X, V, A, rgau, tgy = fc.integrate_star(s, g, e0, pn=2, t_ic=t_apo, t_end=t_apo + 2.05 * P)
    # apocentres: refine on grid
    r = np.linalg.norm(X, axis=1)
    idx = [0] + [k for k in range(1, len(r) - 1) if r[k] >= r[k - 1] and r[k] >= r[k + 1]]
    vals = []; tt = []
    for k in idx[:3]:
        ts = np.linspace(T[max(k - 1, 0)], T[min(k + 1, len(T) - 1)], 401)
        xs, vs = dd.states_at(T, X, V, A, ts); j = np.argmax(np.linalg.norm(xs, axis=1))
        Lv = np.cross(xs[j], vs[j]); Lv /= np.linalg.norm(Lv)
        inc = np.arccos(Lv[2]); Om = np.arctan2(Lv[0], -Lv[1])
        ev = np.cross(vs[j], np.cross(xs[j], vs[j])) - xs[j] / np.linalg.norm(xs[j])   # GM=1
        nvec = np.array([np.cos(Om), np.sin(Om), 0.0])
        w = np.arctan2(ev @ np.cross(Lv, nvec), ev @ nvec)
        vals.append([inc, Om, w]); tt.append(ts[j] * tgy)
    vals = np.unwrap(np.array(vals), axis=0); dt = tt[2] - tt[0]
    return (vals[2] - vals[0])[:nrate] / dt

def rate_jac(names, g, free, nrate, steps):
    J = np.zeros((nrate * len(names), len(free)))
    for c, k in enumerate(free):
        gp = g.copy(); gm = g.copy(); gp[k] += steps[k]; gm[k] -= steps[k]
        J[:, c] = np.concatenate([(orbit_rates(n, gp, nrate) - orbit_rates(n, gm, nrate)) / (2 * steps[k]) for n in names])
    return J

def secular_sigma(case, nrate, md=1.0, orient='CWD', s_rv=1.0, t_end=2041.0, s_ast=100.0):
    stars, g, els, gfree, pri = F.build(case, s_rv=s_rv, t_end=t_end, s_ast=s_ast, md=md, orient=orient, spin=dirs['perpN'])
    pri_full = np.concatenate([pri[:2], pri[len(gfree):]])     # M6, R0 + elements
    Cr = drift_rate_cov(stars, els, pri_full, nrate)
    free = [2, 3, 4] + ([5, 6, 7] if 5 in gfree else [])
    steps = {2: 0.05, 3: 0.05, 4: 0.05, 5: 0.05, 6: 1e-3, 7: 1e-3}
    Jr = rate_jac([s.name for s in stars], g, free, nrate, steps)
    Fs = Jr.T @ np.linalg.inv(Cr) @ Jr
    ev = np.linalg.eigvalsh(Fs / np.outer(np.sqrt(np.diag(Fs)), np.sqrt(np.diag(Fs))))
    rank = int(np.sum(ev > 1e-10 * ev.max()))
    if rank < len(free):
        return dict(rank=rank, npar=len(free), sN=np.inf, sM=np.inf, sL=np.inf)
    C = np.linalg.inv(Fs); CS = C[:3, :3]
    out = dict(rank=rank, npar=len(free), sN=float(np.sqrt(N @ CS @ N)), sM=float(np.sqrt(M_ @ CS @ M_)), sL=float(np.sqrt(L @ CS @ L)))
    if 5 in free: out['sMd'] = float(1e3 * np.sqrt(C[3, 3]))
    out['sig_rates_S301_urad_per_yr'] = (1e6 * np.sqrt(np.diag(Cr)[:nrate])).tolist()
    return out

if __name__ == '__main__':
    t0 = time.time(); R = {'profile': [], 'secular': []}
    print("== (1) within-orbit profile, S301 A, one orbit apocentre->apocentre (Newton + LT only / Newton + disk only)")
    for lab in ['perpN', 'perpM']:
        g = st.G0.copy(); g[2:5] = 0.7 * dirs[lab]; R['profile'].append(profile(g, 3, f"LT chi=0.7 spin {lab}"))
    for orient in ['CWD', 'P30']:
        g = st.G0.copy(); g[5] = 1.0; g[6], g[7] = F.disk_pole(orient, el); R['profile'].append(profile(g, 0, f"disk 1e3 Msun {orient}"))
    print("== (2) secular-only vs full Fisher (S301 100 uas, RV 1 km/s from 2030, 2026-2041; spin perp)")
    for case, md, orient in [('a', 0.0, '-'), ('b', 1.0, 'CWD'), ('c', 1.0, 'CWD'), ('c', 1.0, 'P30'), ('c', 0.3, 'CWD')]:
        for nrate in (2, 3):
            r = secular_sigma(case, nrate, md=md, orient=orient)
            r.update(case=case, md=md, orient=orient, nrate=nrate); R['secular'].append(r)
            print(f"secular case={case} md={md} {orient} rates={'(i,O)' if nrate==2 else '(i,O,w)'}: rank {r['rank']}/{r['npar']}  sN={r['sN']:.3f} sM={r['sM']:.3f} sL={r['sL']:.3f}" + (f" sMd={r['sMd']:.0f}" if 'sMd' in r else ''), flush=True)
    json.dump(R, open('q2.json', 'w'), indent=1, default=float)
    print(f"[{time.time()-t0:.0f}s]")
