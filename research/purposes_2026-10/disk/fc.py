"""Forecast core: units, stars, epochs, forward model, Fisher.
Global parameters (8): M6 [1e6 Msun], R0 [kpc], Sx, Sy, Sz [chi components, observer frame
X=N, Y=E, Z=away], Md [1e3 Msun], i_d, O_d [rad] (disk normal written as an orbital pole:
n_d = (sin O sin i, -cos O sin i, cos i)).
Per-star parameters (6): a [mas], e, i, Omega, omega [rad], tp [yr]  (osculating Kepler
elements, GM = M, propagated to the start state at t_ic = first epoch - 0.3 yr).
alpha = -1 and rout (fixed in AU; 31740 r_g of the fiducial mass) are fixed."""
import numpy as np
import diskdyn as dd
import disktab

GM_AU = 9.870628e-9          # GM_sun/c^2 [AU]
GM_S = 4.925490947e-6        # GM_sun/c^3 [s]
YR_S = 3.15576e7
C_KMS = 299792.458
T_REF = 2020.0
M6_FID, R0_FID = 4.297, 8.277
RG_AU_FID = GM_AU * 1e6 * M6_FID
ROUT_AU = 31740.0 * RG_AU_FID      # Piran+ fiducial disk edge (S301 apocentre)
GNAMES = ['M6', 'R0', 'Sx', 'Sy', 'Sz', 'Md', 'i_d', 'O_d']
ENAMES = ['a', 'e', 'i', 'O', 'w', 'tp']
NG = len(GNAMES)
DEG = np.pi / 180

TR, TZ, TP = disktab.make_table(verbose=False)


def rg_au(M6): return GM_AU * 1e6 * M6
def tg_yr(M6): return GM_S * 1e6 * M6 / YR_S
def uas_per_rg(M6, R0): return rg_au(M6) * 1e3 / R0   # rg[AU] / (R0[kpc]*1e3 pc) arcsec * 1e6 -> uas


def rotmat(inc, Om, om):
    cO, sO = np.cos(Om), np.sin(Om); ci, si = np.cos(inc), np.sin(inc); cw, sw = np.cos(om), np.sin(om)
    return np.array([[cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci, sO * si],
                     [sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci, -cO * si],
                     [sw * si, cw * si, ci]])


def pole(inc, Om):
    return np.array([np.sin(Om) * np.sin(inc), -np.cos(Om) * np.sin(inc), np.cos(inc)])


def pole_angles(n):
    n = n / np.linalg.norm(n)
    inc = np.arccos(np.clip(n[2], -1, 1)); Om = np.arctan2(n[0], -n[1])
    return inc, Om % (2 * np.pi)


def kepler_state(a, e, inc, Om, om, Mmean):
    """GM = 1 state for mean anomaly Mmean."""
    Mm = np.mod(Mmean, 2 * np.pi)
    E = Mm if e < 0.8 else np.pi
    for _ in range(100):
        dE = (E - e * np.sin(E) - Mm) / (1 - e * np.cos(E)); E -= dE
        if abs(dE) < 1e-15: break
    cE, sE = np.cos(E), np.sin(E)
    r = a * (1 - e * cE)
    xp = a * (cE - e); yp = a * np.sqrt(1 - e * e) * sE
    vfac = 1.0 / (np.sqrt(a) * (1 - e * cE))   # n a / (1 - e cos E) with n = a^-1.5
    vxp = -vfac * sE; vyp = vfac * np.sqrt(1 - e * e) * cE
    Rm = rotmat(inc, Om, om)
    return Rm @ np.array([xp, yp, 0.0]), Rm @ np.array([vxp, vyp, 0.0])


def s_orbit(e, n=20000):
    E = (np.arange(n) + 0.5) * 2 * np.pi / n
    return np.sum(1.0 / np.sqrt(1 - e * np.cos(E))) * 2 * np.pi / n


class Star:
    """Elements + data.  t_ast, s_ast [yr, uas]; t_rv, s_rv [yr, km/s]."""
    def __init__(self, name, el, t_ast, s_ast, t_rv=None, s_rv=None, nstep=4000):
        self.name = name; self.el = np.array(el, float)
        self.t_ast = np.asarray(t_ast, float); self.s_ast = np.broadcast_to(np.asarray(s_ast, float), self.t_ast.shape).copy()
        self.t_rv = np.zeros(0) if t_rv is None else np.asarray(t_rv, float)
        self.s_rv = np.zeros(0) if s_rv is None else np.broadcast_to(np.asarray(s_rv, float), self.t_rv.shape).copy()
        self.nstep = nstep
        allt = np.concatenate([self.t_ast, self.t_rv])
        self.t_ic = allt.min() - 0.3; self.t_end = allt.max() + 0.3
        self.ds = s_orbit(self.el[1]) / nstep

    def sig(self):
        return np.concatenate([self.s_ast, self.s_ast, self.s_rv])

    def nobs(self):
        return 2 * len(self.t_ast) + len(self.t_rv)


def disk_normal(g):
    return pole(g[6], g[7])


def integrate_star(star, g, el, pn=2, t_ic=None, t_end=None):
    """Returns (T, X, V, A) in code units of the model mass, plus conversion factors."""
    M6, R0 = g[0], g[1]
    S = np.array(g[2:5], float)
    md = g[5] * 1e3 / (M6 * 1e6)
    nd = disk_normal(g)
    rgau = rg_au(M6); tgy = tg_yr(M6)
    a_rg = el[0] * 1e-3 * R0 * 1e3 / rgau          # a[AU] = a["] * R0[pc]
    P_yr = 2 * np.pi * a_rg ** 1.5 * tgy
    t_ic = star.t_ic if t_ic is None else t_ic
    t_end = star.t_end if t_end is None else t_end
    Mmean = 2 * np.pi * (t_ic - el[5]) / P_yr
    x0, v0 = kepler_state(a_rg, el[1], el[2], el[3], el[4], Mmean)
    t0c = (t_ic - T_REF) / tgy; t1c = (t_end - T_REF) / tgy
    norb = (t_end - t_ic) / P_yr
    nmax = int(star.nstep * (norb + 1) * 1.3) + 200
    rout = ROUT_AU / rgau
    T, X, V, A = dd.integrate(x0, v0, t0c, t1c, star.ds, nmax, S, md, nd, rout, TR, TZ, TP, pn)
    if T[-1] < t1c:
        raise RuntimeError("nmax too small")
    return T, X, V, A, rgau, tgy


def predict(star, g, el, pn=2):
    """Model data vector [north uas (n_ast), east uas (n_ast), rv km/s (n_rv)]."""
    T, X, V, A, rgau, tgy = integrate_star(star, g, el, pn)
    tob = (np.concatenate([star.t_ast, star.t_rv]) - T_REF) / tgy
    out, te = dd.observe(T, X, V, A, tob)
    k = uas_per_rg(g[0], g[1]); na = len(star.t_ast)
    return np.concatenate([out[:na, 0] * k, out[:na, 1] * k, out[na:, 2] * C_KMS])


# ----------------------------------------------------------------- Fisher
GSTEP = np.array([1e-6, 1e-6, 0.02, 0.02, 0.02, 0.01, 1e-4, 1e-4])   # M6, R0 relative; others absolute
ESTEP = np.array([1e-6, 1e-7, 1e-6, 1e-6, 1e-6, 2e-6])               # a relative; others absolute
GREL = np.array([1, 1, 0, 0, 0, 0, 0, 0], bool)


def jacobian(stars, g, els, gfree, efree=True, pn=2, scale=1.0, pool=None):
    """Whitened Jacobian d(y/sigma)/dtheta for free globals (index list) + all elements of
    every star (if efree).  Returns J (nobs_tot x npar), names."""
    nst = len(stars)
    names = [GNAMES[k] for k in gfree] + ([f"{s.name}.{n}" for s in stars for n in ENAMES] if efree else [])
    nobs = [s.nobs() for s in stars]; off = np.concatenate([[0], np.cumsum(nobs)])
    J = np.zeros((off[-1], len(names)))
    tasks = []
    for col, k in enumerate(gfree):
        h = GSTEP[k] * scale * (abs(g[k]) if GREL[k] else 1.0)
        for si in range(nst):
            tasks.append((si, col, k, None, h))
    if efree:
        for si in range(nst):
            for q in range(6):
                col = len(gfree) + 6 * si + q
                h = ESTEP[q] * scale * (abs(els[si][q]) if q == 0 else 1.0)
                tasks.append((si, col, None, q, h))
    def run(t):
        si, col, k, q, h = t
        res = []
        for sgn in (+1, -1):
            gg = np.array(g, float); ee = np.array(els[si], float)
            if k is not None: gg[k] += sgn * h
            else: ee[q] += sgn * h
            res.append(predict(stars[si], gg, ee, pn))
        return (res[0] - res[1]) / (2 * h) / stars[si].sig()
    outs = pool.map(_run_task, [(stars[t[0]], g, els[t[0]], t, pn) for t in tasks]) if pool is not None else [run(t) for t in tasks]
    for t, d in zip(tasks, outs):
        si, col = t[0], t[1]
        J[off[si]:off[si + 1], col] = d
    return J, names


def _run_task(args):
    star, g, el, t, pn = args
    si, col, k, q, h = t
    res = []
    for sgn in (+1, -1):
        gg = np.array(g, float); ee = np.array(el, float)
        if k is not None: gg[k] += sgn * h
        else: ee[q] += sgn * h
        res.append(predict(star, gg, ee, pn))
    return (res[0] - res[1]) / (2 * h) / star.sig()
