"""ONE emcee run for case (c): S301 (100 uas, RV 1 km/s from 2030) + S2/S55/S38/S29, 2026-2041,
disk 1e3 Msun (CWD) in truth, (Md, i_d, O_d) free; 38 parameters; noise-free (Asimov) data;
walkers initialised from the Fisher covariance; sampling in Fisher-whitened coordinates
theta = theta0 + Lchol z; DE + DE-snooker moves.  Compares posterior widths with Fisher."""
import numpy as np, emcee, time, json, sys, os
from multiprocessing import Pool
import fc, stars as st, forecast as F
NSTEPS = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
NW = 160
el301 = np.array(st.EL['S301A']); dirs, (N, M_, L) = st.spin_dirs(el301)
Pv = fc.rotmat(el301[2], el301[3], el301[4])[:, 0]; Qv = np.cross(L, Pv)
STARS, G, ELS, GFREE, PRI = F.build('c', s_rv=1.0, t_end=2041.0, md=1.0, orient='CWD', spin=dirs['b45l45'])
YT = [fc.predict(s, G, e) for s, e in zip(STARS, ELS)]
SIG = [s.sig() for s in STARS]
TH0 = np.concatenate([G[GFREE]] + [np.array(e) for e in ELS])
FIN = np.isfinite(PRI); W = np.where(FIN, 1 / np.where(FIN, PRI, 1), 0)
J, names = fc.jacobian(STARS, G, ELS, GFREE)
CF, cond = F.cov_from(J, PRI)
LCH = np.linalg.cholesky(CF)
NG = len(GFREE)

def lnp(z):
    th = TH0 + LCH @ z
    g = G.copy(); g[GFREE] = th[:NG]
    if not (0.5 < th[1] / G[1] < 2 and 0 < th[0]): return -np.inf
    chi2 = 0.0
    try:
        for k, s in enumerate(STARS):
            e = th[NG + 6 * k: NG + 6 * k + 6]
            if not (0 < e[1] < 0.9995 and e[0] > 0): return -np.inf
            r = (fc.predict(s, g, e) - YT[k]) / SIG[k]; chi2 += r @ r
    except Exception:
        return -np.inf
    rp = W[FIN] * (th - TH0)[FIN]
    return -0.5 * (chi2 + rp @ rp)

if __name__ == '__main__':
    t0 = time.time()
    rng = np.random.default_rng(42)
    z0 = rng.normal(size=(NW, len(TH0)))
    moves = [(emcee.moves.DEMove(), 0.8), (emcee.moves.DESnookerMove(), 0.2)]
    out = '/tmp/claude-0/-home-user-relatipy/9021de26-e198-5a73-ad04-6e893f8ac1b1/scratchpad/mcmc_c_chain.h5'
    with Pool(4) as pool:
        sampler = emcee.EnsembleSampler(NW, len(TH0), lnp, moves=moves, pool=pool)
        for i, _ in enumerate(sampler.sample(z0, iterations=NSTEPS, progress=False)):
            if (i + 1) % 250 == 0:
                ch = sampler.get_chain()
                try: tau = emcee.autocorr.integrated_time(ch[(i + 1) // 3:], quiet=True)
                except Exception: tau = np.full(len(TH0), np.nan)
                zz = ch[(i + 1) // 3:].reshape(-1, len(TH0)); th = TH0 + zz @ LCH.T
                S = th[:, 2:5]
                sP, sQ, sL_ = [np.std(S @ u) for u in (Pv, Qv, L)]
                fP, fQ, fL = [np.sqrt(u @ CF[2:5, 2:5] @ u) for u in (Pv, Qv, L)]
                print(f"step {i+1} acc={np.mean(sampler.acceptance_fraction):.2f} tau_max={np.nanmax(tau):.0f} tau_med={np.nanmedian(tau):.0f} "
                      f"MCMC/Fisher sigma: P {sP/fP:.2f} Q {sQ/fQ:.2f} L {sL_/fL:.2f} Md {np.std(th[:,5])/np.sqrt(CF[5,5]):.2f} "
                      f"i_d {np.std(th[:,6])/np.sqrt(CF[6,6]):.2f} O_d {np.std(th[:,7])/np.sqrt(CF[7,7]):.2f} | z-std med {np.median(zz.std(0)):.2f} [{time.time()-t0:.0f}s]", flush=True)
                np.save(out.replace('.h5', '.npy'), ch[::5].astype(np.float32))
    ch = sampler.get_chain(); burn = NSTEPS // 3
    zz = ch[burn:].reshape(-1, len(TH0)); th = TH0 + zz @ LCH.T
    S = th[:, 2:5]
    res = dict(fisher=dict(zip(['sP', 'sQ', 'sL', 'sMd', 'si_d', 'sO_d'], [float(np.sqrt(u @ CF[2:5, 2:5] @ u)) for u in (Pv, Qv, L)] + [float(1e3*np.sqrt(CF[5, 5])), float(np.degrees(np.sqrt(CF[6, 6]))), float(np.degrees(np.sqrt(CF[7, 7])))])),
               mcmc=dict(zip(['sP', 'sQ', 'sL', 'sMd', 'si_d', 'sO_d'], [float(np.std(S @ u)) for u in (Pv, Qv, L)] + [float(1e3*np.std(th[:, 5])), float(np.degrees(np.std(th[:, 6]))), float(np.degrees(np.std(th[:, 7])))])),
               mcmc_q16_84=dict(zip(['P', 'Q', 'L', 'Md'], [np.percentile(S @ u - G[2:5] @ u, [16, 50, 84]).tolist() for u in (Pv, Qv, L)] + [(1e3*np.percentile(th[:, 5], [16, 50, 84])).tolist()])),
               acc=float(np.mean(sampler.acceptance_fraction)), nsteps=NSTEPS, nwalkers=NW, ndim=len(TH0),
               tau=emcee.autocorr.integrated_time(ch[burn:], quiet=True).tolist())
    json.dump(res, open('mcmc_c.json', 'w'), indent=1)
    np.save('mcmc_c_spin_disk_thin.npy', th[::20, :8].astype(np.float32))
    print(json.dumps({k: res[k] for k in ('fisher', 'mcmc', 'mcmc_q16_84', 'acc')}, indent=0))
    print(f"[{time.time()-t0:.0f}s]")
