"""Q1 production Fisher forecasts.  Usage: python forecast.py [tag]
Cases: a = S301 alone, no disk; b = S301 alone, disk in truth, (Md, i_d, O_d) marginalised;
c = S301+S2+S55+S38+S29, disk marginalised; c0 = joint, no disk.
Priors (stand-ins for past data): M, R0 from GRAVITY 2024 (4.2996+-0.0118e6 Msun, 8275.9+-8.6 pc)
in all cases (also run without for case a); friends' elements = published 1-sigma (GRAVITY 2024
Table B.1, S29 tp error placeholder 5e-4 yr), applied to the osculating elements at t_ic.
S301: 19 past epochs at 207 uas + GRAVITY+ 2026-t_end (6/yr + 20 within +-0.3 yr of each
pericentre); RV from rv_start (4/yr + 15 within +-0.3 yr of each pericentre)."""
import numpy as np, json, sys, time, itertools
from multiprocessing import Pool
import fc, stars as st
D = fc.DEG
PRIOR_G = {0: 0.0118, 1: 0.0086}
PRIOR_EL = {k: np.array(v, float) for k, v in st.EL_ERR.items()}
PRIOR_EL['S29'][5] = 5e-4
FRIENDS = ['S2', 'S55', 'S38', 'S29']


def disk_pole(kind, el301):
    if kind == 'CWD':
        return 127 * D, 99 * D
    if kind == 'P30':       # Piran+ Table 1 geometry for S301: i = 30 deg, omega = 90 deg w.r.t. the disk
        L = fc.pole(el301[2], el301[3]); Pv = fc.rotmat(el301[2], el301[3], el301[4])[:, 0]
        return fc.pole_angles(np.cos(30 * D) * L + np.sin(30 * D) * Pv)
    if kind.startswith('R'):  # random orientation, seed = int(kind[1:])
        rng = np.random.default_rng(int(kind[1:])); v = rng.normal(size=3); v /= np.linalg.norm(v)
        if v[2] < 0: v = -v
        return fc.pole_angles(v)
    raise ValueError(kind)


def build(case, sol='A', t_end=2041.0, s_ast=100.0, s_rv=None, rv_start=2030.0, md=1.0, orient='CWD',
          spin=None, friends_past='prior', friends_fut=30.0, mprior=True):
    el301 = np.array(st.EL['S301' + sol])
    stars = [st.s301(sol, t_end=t_end, s_ast=s_ast, s_rv=s_rv, rv_start=rv_start)]
    if case in ('c', 'c0'):
        for nm in FRIENDS:
            if friends_past == 'prior':
                f = st.friend(nm, t_end=t_end, s_fut=friends_fut)
                keep_a = f.t_ast >= 2026.0; keep_r = f.t_rv >= 2026.0
                stars.append(fc.Star(nm, f.el, f.t_ast[keep_a], f.s_ast[keep_a], f.t_rv[keep_r], f.s_rv[keep_r], nstep=f.nstep))
            else:
                stars.append(st.friend(nm, t_end=t_end, s_fut=friends_fut))
    g = st.G0.copy()
    g[2:5] = 0.7 * spin
    disk = case in ('b', 'c') and md > 0
    if disk:
        g[5] = md; g[6], g[7] = disk_pole(orient, el301)
    gfree = [0, 1, 2, 3, 4] + ([5, 6, 7] if disk else [])
    els = [np.array(s.el) for s in stars]
    pri = [PRIOR_G.get(k, np.inf) if mprior else np.inf for k in gfree]
    for s in stars:
        pri += list(PRIOR_EL[s.name]) if (s.name in FRIENDS and friends_past == 'prior') else [np.inf] * 6
    return stars, g, els, gfree, np.array(pri)


def cov_from(J, pri):
    F = J.T @ J + np.diag(np.where(np.isfinite(pri), 1.0 / np.where(np.isfinite(pri), pri, 1.0) ** 2, 0.0))
    d = np.sqrt(np.diag(F)); Fn = F / np.outer(d, d)
    return np.linalg.inv(Fn) / np.outer(d, d), np.linalg.cond(Fn)


def summarize(C, g, gfree, el301, cond):
    ix = {k: j for j, k in enumerate(gfree)}
    CS = C[ix[2]:ix[4] + 1, ix[2]:ix[4] + 1]
    L = fc.pole(el301[2], el301[3]); N = np.array([np.cos(el301[3]), np.sin(el301[3]), 0.]); Mv = np.cross(L, N)
    S = g[2:5]; chi = np.linalg.norm(S); sh = S / chi
    Pp = np.eye(3) - np.outer(sh, sh)
    out = dict(sSx=np.sqrt(CS[0, 0]), sSy=np.sqrt(CS[1, 1]), sSz=np.sqrt(CS[2, 2]),
               sN=np.sqrt(N @ CS @ N), sM=np.sqrt(Mv @ CS @ Mv), sL=np.sqrt(L @ CS @ L),
               schi=np.sqrt(sh @ CS @ sh), sang_deg=np.degrees(np.sqrt(np.trace(Pp @ CS @ Pp)) / chi),
               sM6=np.sqrt(C[ix[0], ix[0]]), sR0=np.sqrt(C[ix[1], ix[1]]), cond=cond)
    if 5 in ix:
        out.update(sMd=1e3 * np.sqrt(C[ix[5], ix[5]]), sid=np.degrees(np.sqrt(C[ix[6], ix[6]])), sOd=np.degrees(np.sqrt(C[ix[7], ix[7]])))
    return out


def run_one(cfg):
    stars, g, els, gfree, pri = build(**cfg['b'])
    J, names = fc.jacobian(stars, g, els, gfree, scale=cfg.get('scale', 1.0))
    C, cond = cov_from(J, pri)
    res = summarize(C, g, gfree, els[0], cond)
    res.update(cfg['meta'])
    return res


def configs(sol='A', spins=None, cases=('a', 'a_noprior', 'b', 'c', 'c0'), asts=(100.0, 30.0)):
    dirs, _ = st.spin_dirs(np.array(st.EL['S301' + sol]))
    out = []
    for case in cases:
        disks = [(0.0, '-')] if case in ('a', 'a_noprior', 'c0') else [(1.0, 'CWD'), (0.3, 'CWD'), (1.0, 'P30'), (0.3, 'P30')]
        for (md, orient), t_end, s_rv, s_ast, (slab, sv) in itertools.product(disks, (2035.0, 2041.0), (None, 20.0, 5.0, 1.0), asts, dirs.items()):
            if spins is not None and slab not in spins: continue
            b = dict(case=case.replace('_noprior', ''), sol=sol, t_end=t_end, s_ast=s_ast, s_rv=s_rv, md=md, orient=orient, spin=sv,
                     mprior=(case != 'a_noprior'))
            out.append(dict(b=b, meta=dict(case=case, sol=sol, t_end=t_end, s_ast=s_ast, s_rv=s_rv or 0, md=md, orient=orient, spin=slab)))
    return out


if __name__ == '__main__':
    tag = sys.argv[1] if len(sys.argv) > 1 else 'A'
    t0 = time.time()
    if tag == 'A':
        cf = configs('A')
    elif tag == 'B':
        cf = [c for c in configs('B', asts=(100.0,)) if c['meta']['orient'] in ('-', 'CWD') and c['meta']['md'] in (0.0, 1.0)]
    elif tag == 'steps':
        base = [c for c in configs('A', spins=['perpN'], cases=('c',), asts=(100.0,)) if c['meta']['s_rv'] == 1.0 and c['meta']['t_end'] == 2041.0 and c['meta']['orient'] == 'CWD' and c['meta']['md'] == 1.0]
        cf = []
        for sc in (0.3, 1.0, 3.0):
            c = json.loads(json.dumps(base[0], default=lambda o: o.tolist())); c['b']['spin'] = np.array(c['b']['spin']); c['scale'] = sc; c['meta']['scale'] = sc; cf.append(c)
    elif tag == 'mockpast':
        cf = [c for c in configs('A', cases=('c',), asts=(100.0,)) if c['meta']['md'] == 1.0 and c['meta']['orient'] == 'CWD']
        for c in cf: c['b']['friends_past'] = 'mock'; c['meta']['friends_past'] = 'mock'
    elif tag == 'f100':
        cf = [c for c in configs('A', cases=('c',), asts=(100.0,)) if c['meta']['md'] == 1.0 and c['meta']['orient'] == 'CWD']
        for c in cf: c['b']['friends_fut'] = 100.0; c['meta']['friends_fut'] = 100.0
    print(f"{len(cf)} configs", flush=True)
    with Pool(4) as p:
        res = p.map(run_one, cf, chunksize=1)
    json.dump(res, open(f'results_{tag}.json', 'w'), indent=0, default=float)
    print(f"[{time.time()-t0:.0f}s] -> results_{tag}.json")
