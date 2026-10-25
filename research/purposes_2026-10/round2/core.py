"""round2 core: truth generator + three fit models + fitting/Fisher helpers.

Geometric units G=c=M=1 inside every model evaluation (M = that model's mass).
Observer frame (relatipy convention): X=Dec (north), Y=RA (east), Z=away from observer.
alpha (RA offset) = Y/R0, delta (Dec offset) = X/R0.  Arrival time t_obs = t_e + Z (+Shapiro).

Orbit parametrisation, identical for every model: Newtonian Kepler state (GM=1) at the
APOCENTRE (true anomaly pi) at coordinate time t0.  (At apocentre PN corrections to the
energy are ~1/r_a^2, so the 'a' of every model maps to nearly the same physical period;
a pericentre-referenced Kepler state would not for S301, where PN energy terms ~ binding energy.)
The state is interpreted as: harmonic Cartesian (model i, truth, model iii) or
Boyer-Lindquist Cartesian (model ii = relatipy KerrMcmcModel as is).

Parameter vector (full, 18 entries; fits free a subset):
 M[1e6 Msun] R0[kpc] a[mas] e inc Om om [rad] t0[yr] x0 y0[uas] vx0 vy0[uas/yr] vz0[km/s]
 fsp f ast ths phs   (spin axis polar angles in observer frame)
"""
import numpy as np, warnings
from numba import njit
from astropy import units as u, constants as const
from scipy.optimize import least_squares
warnings.filterwarnings("ignore")

C_KMS = const.c.to_value(u.km / u.s)
RG_AU = (const.G * 1e6 * const.M_sun / const.c**2).to_value(u.au)    # r_g [AU] per 1e6 Msun
TG_YR = (const.G * 1e6 * const.M_sun / const.c**3).to_value(u.yr)    # r_g/c [yr] per 1e6 Msun
AU_KPC = (1 * u.kpc).to_value(u.au)
UAS_RAD = (1 * u.rad).to_value(u.uas)
NAMES = ['M', 'R0', 'a', 'e', 'inc', 'Om', 'om', 't0', 'x0', 'y0', 'vx0', 'vy0', 'vz0',
         'fsp', 'f', 'ast', 'ths', 'phs']
IX = {n: i for i, n in enumerate(NAMES)}
# finite-difference steps (absolute) for central differences
STEP = dict(M=1e-6, R0=1e-6, a=1e-6, e=1e-7, inc=1e-6, Om=1e-6, om=1e-6, t0=1e-6,
            x0=1.0, y0=1.0, vx0=1.0, vy0=1.0, vz0=1.0, fsp=1e-4, f=1e-4, ast=1e-4, ths=1e-4, phs=1e-4)
REL = {'M', 'R0', 'a'}          # these steps are relative


def uas_per_rg(M6, R0):
    return RG_AU * M6 / (R0 * AU_KPC) * UAS_RAD


def rotmat(inc, Om, om):
    cO, sO = np.cos(Om), np.sin(Om); ci, si = np.cos(inc), np.sin(inc); cw, sw = np.cos(om), np.sin(om)
    return np.array([[cO*cw - sO*sw*ci, -cO*sw - sO*cw*ci,  sO*si],
                     [sO*cw + cO*sw*ci, -sO*sw + cO*cw*ci, -cO*si],
                     [sw*si,             cw*si,             ci]])


def apo_state(a, e, inc, Om, om):
    """Kepler (GM=1) state at apocentre, observer frame; standard convention
    (ascending node = receding, Z away)."""
    R = rotmat(inc, Om, om)
    ra = a * (1 + e); va = np.sqrt((1 - e) / (a * (1 + e)))
    return R @ np.array([-ra, 0.0, 0.0]), R @ np.array([0.0, -va, 0.0])


def axis(ths, phs):
    return np.array([np.sin(ths)*np.cos(phs), np.sin(ths)*np.sin(phs), np.cos(ths)])


def body_to_observer(vec):
    """Same construction as relatipy.observables.kerr_mcmc._body_to_observer."""
    vec = np.asarray(vec, float); n = np.linalg.norm(vec)
    if n == 0: return np.eye(3)
    ax = vec / n; ref = np.array([1.0, 0.0, 0.0])
    if abs(ax[0]) > 0.9: ref = np.array([0.0, 1.0, 0.0])
    bx = ref - ref @ ax * ax; bx /= np.linalg.norm(bx); by = np.cross(ax, bx)
    return np.ascontiguousarray(np.column_stack((bx, by, ax)))


def geo(p):
    """common derived quantities: a in r_g, yr per r_g/c, uas per r_g."""
    M6, R0 = p[0], p[1]
    return p[2] * R0 / (RG_AU * M6), TG_YR * M6, uas_per_rg(M6, R0)

# ----------------------------------------------------------------------------------------
# numba kernels: PN EOM, Sundman RK4 both directions, Hermite dense output, emission solve
# ----------------------------------------------------------------------------------------
@njit(cache=True)
def _acc(x, v, fsp, S, son):
    r = np.sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]); n = x / r
    a = -n / (r*r)
    vv = v[0]*v[0] + v[1]*v[1] + v[2]*v[2]; rd = v[0]*n[0] + v[1]*n[1] + v[2]*n[2]
    a = a + fsp / (r*r) * (n * (4.0/r - vv) + 4.0 * rd * v)
    if son:
        nS = n[0]*S[0] + n[1]*S[1] + n[2]*S[2]
        nxv = np.array([n[1]*v[2]-n[2]*v[1], n[2]*v[0]-n[0]*v[2], n[0]*v[1]-n[1]*v[0]])
        vxS = np.array([v[1]*S[2]-v[2]*S[1], v[2]*S[0]-v[0]*S[2], v[0]*S[1]-v[1]*S[0]])
        a = a + (2.0/(r*r*r)) * (3.0*nS*nxv + vxS)
    return a


@njit(cache=True)
def _f(y, fsp, S, son):
    x = y[0:3]; v = y[3:6]
    r = np.sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]); g = r**1.5
    out = np.empty(7)
    out[0:3] = g * v; out[3:6] = g * _acc(x, v, fsp, S, son); out[6] = g
    return out


@njit(cache=True)
def _run(x0, v0, tlim, fsp, S, son, ds, nmax):
    T = np.empty(nmax); X = np.empty((nmax, 3)); V = np.empty((nmax, 3))
    y = np.empty(7); y[0:3] = x0; y[3:6] = v0; y[6] = 0.0
    T[0] = 0.0; X[0] = x0; V[0] = v0; k = 0
    while abs(y[6]) < abs(tlim) and k < nmax - 1:
        k1 = _f(y, fsp, S, son); k2 = _f(y + 0.5*ds*k1, fsp, S, son)
        k3 = _f(y + 0.5*ds*k2, fsp, S, son); k4 = _f(y + ds*k3, fsp, S, son)
        y = y + ds/6.0 * (k1 + 2*k2 + 2*k3 + k4); k += 1
        T[k] = y[6]; X[k] = y[0:3]; V[k] = y[3:6]
    return T[:k+1], X[:k+1], V[:k+1]


@njit(cache=True)
def pn_traj(x0, v0, tlo, thi, fsp, S, son, ds, nmax):
    Tf, Xf, Vf = _run(x0, v0, thi, fsp, S, son, ds, nmax)
    Tb, Xb, Vb = _run(x0, v0, tlo, fsp, S, son, -ds, nmax)
    nb = len(Tb); nf = len(Tf); n = nb - 1 + nf
    T = np.empty(n); X = np.empty((n, 3)); V = np.empty((n, 3))
    for i in range(nb - 1):
        T[i] = Tb[nb-1-i]; X[i] = Xb[nb-1-i]; V[i] = Vb[nb-1-i]
    for i in range(nf):
        T[nb-1+i] = Tf[i]; X[nb-1+i] = Xf[i]; V[nb-1+i] = Vf[i]
    return T, X, V


@njit(cache=True)
def _herm(T, X, V, t):
    lo = 0; hi = len(T) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if T[mid] <= t: lo = mid
        else: hi = mid
    h = T[hi] - T[lo]; s = (t - T[lo]) / h
    h00 = 2*s**3 - 3*s**2 + 1; h10 = s**3 - 2*s**2 + s; h01 = -2*s**3 + 3*s**2; h11 = s**3 - s**2
    x = h00*X[lo] + h10*h*V[lo] + h01*X[hi] + h11*h*V[hi]
    d00 = (6*s**2 - 6*s)/h; d10 = 3*s**2 - 4*s + 1; d01 = (-6*s**2 + 6*s)/h; d11 = 3*s**2 - 2*s
    v = d00*X[lo] + d10*V[lo] + d01*X[hi] + d11*V[hi]
    return x, v

RREF = 1.0e4


@njit(cache=True)
def emit(T, X, V, tobs, shap):
    """Emission states for arrival times tobs: t_obs = t_e + Z - 2*shap*ln((r-Z)/RREF)."""
    n = len(tobs); XE = np.empty((n, 3)); VE = np.empty((n, 3))
    for i in range(n):
        t = tobs[i]
        for it in range(40):
            x, v = _herm(T, X, V, t)
            r = np.sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]); q = r - x[2]
            rd = (x[0]*v[0] + x[1]*v[1] + x[2]*v[2]) / r
            F = t + x[2] - 2.0*shap*np.log(q/RREF) - tobs[i]
            dF = 1.0 + v[2] - 2.0*shap*(rd - v[2])/q
            dt = F / dF; t -= dt
            if abs(dt) < 1e-7: break
        x, v = _herm(T, X, V, t)
        XE[i] = x; VE[i] = v
    return XE, VE


@njit(cache=True)
def ut_kerr(x, v, Rb, a):
    """u^t from the Kerr metric for a harmonic observer-frame state (r_BL = r_H + 1)."""
    xb = np.empty(3); vb = np.empty(3)
    for k in range(3):
        xb[k] = Rb[0, k]*x[0] + Rb[1, k]*x[1] + Rb[2, k]*x[2]
        vb[k] = Rb[0, k]*v[0] + Rb[1, k]*v[1] + Rb[2, k]*v[2]
    rH = np.sqrt(xb[0]**2 + xb[1]**2 + xb[2]**2); R = rH + 1.0
    rho2 = xb[0]**2 + xb[1]**2
    ct = xb[2]/rH; st2 = rho2/(rH*rH); st = np.sqrt(st2)
    rd = (xb[0]*vb[0] + xb[1]*vb[1] + xb[2]*vb[2]) / rH
    thd = -(vb[2] - ct*rd) / (rH*st)
    phd = (xb[0]*vb[1] - xb[1]*vb[0]) / rho2
    Sig = R*R + a*a*ct*ct; Del = R*R - 2.0*R + a*a
    gtt = -(1.0 - 2.0*R/Sig); gtp = -2.0*a*R*st2/Sig; grr = Sig/Del; gqq = Sig
    gpp = (R*R + a*a + 2.0*a*a*R*st2/Sig) * st2
    nrm = -(gtt + 2.0*gtp*phd + grr*rd*rd + gqq*thd*thd + gpp*phd*phd)
    return 1.0/np.sqrt(nrm)


@njit(cache=True)
def obs_T(T, X, V, tobs, Rb, a, lens, shap, rvmode):
    """Truth / model (iii): first-order time-transfer observables.
    img = x_perp (1 + 2*lens/(r-Z)); 1+z = u^t (1 + v.gradT), gradT = Zhat - 2*shap*(n-Zhat)/(r-Z)."""
    XE, VE = emit(T, X, V, tobs, shap)
    n = len(tobs); img = np.empty((n, 2)); z = np.empty(n)
    for i in range(n):
        x = XE[i]; v = VE[i]
        r = np.sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]); q = r - x[2]
        fac = 1.0 + 2.0*lens/q
        img[i, 0] = x[1]*fac; img[i, 1] = x[0]*fac
        rd = (x[0]*v[0] + x[1]*v[1] + x[2]*v[2]) / r
        ut = ut_kerr(x, v, Rb, a)
        z[i] = ut*(1.0 + v[2] - 2.0*shap*(rd - v[2])/q) - 1.0
    return img, z, XE, VE


@njit(cache=True)
def obs_pn(T, X, V, tobs, f, mult):
    """Model (i): Romer only, straight projection, z = v_Z + f(v^2/2+1/r) [mult: (1+v_Z)(1+f(..))-1]."""
    XE, VE = emit(T, X, V, tobs, 0.0)
    n = len(tobs); img = np.empty((n, 2)); z = np.empty(n)
    for i in range(n):
        x = XE[i]; v = VE[i]
        r = np.sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2])
        img[i, 0] = x[1]; img[i, 1] = x[0]
        g = 0.5*(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) + 1.0/r
        if mult: z[i] = (1.0 + v[2])*(1.0 + f*g) - 1.0
        else: z[i] = v[2] + f*g
    return img, z

# ----------------------------------------------------------------------------------------
# exact Kerr trajectories with relatipy (forward + mirrored backward), harmonic mapping
# ----------------------------------------------------------------------------------------
_KERR = {}


def _kerr(M6, spin):
    from relatipy import Kerr
    key = (M6, spin)
    if key not in _KERR:
        if len(_KERR) > 64: _KERR.clear()
        _KERR[key] = Kerr(mass=M6*1e6*u.M_sun, spin=spin)
    return _KERR[key]


def tau_grid(a, e, tmax, npo):
    """proper-time samples uniform in Kepler eccentric anomaly from apocentre (1PN tau(t))."""
    dE = 2*np.pi/npo; E = [np.pi]; sa3 = np.sqrt(a**3)
    while True:
        En = E[-1] + dE
        E.append(En)
        if sa3*(En - e*np.sin(En) - np.pi) > tmax*1.002 + 50.0: break
    E = np.array(E); t = sa3*(E - e*np.sin(E) - np.pi)
    tau = t - 2.0*np.sqrt(a)*(E - np.pi) + t/(2*a)
    return tau


def exact_traj(p, tlo, thi, npo=4000, rtol=1e-12, return_ut=False):
    """Exact Kerr geodesic from the harmonic apocentre state of p, sampled for t in [tlo, thi]
    (geometric, relative to t0).  Returns T, X, V (harmonic, observer frame), Rb, spin."""
    from astropy.constants import c
    M6 = p[0]; a_g, _, _ = geo(p); e = p[3]
    spin = float(p[IX['ast']]); vec = axis(p[IX['ths']], p[IX['phs']])
    Rb = body_to_observer(vec) if spin > 0 else np.eye(3)
    x0, v0 = apo_state(a_g, e, p[4], p[5], p[6])
    xb = Rb.T @ x0; vb = Rb.T @ v0
    r = np.linalg.norm(xb); th = np.arccos(xb[2]/r); ph = np.arctan2(xb[1], xb[0])
    rd = xb @ vb / r; thd = -(vb[2] - np.cos(th)*rd)/(r*np.sin(th))
    phd = (xb[0]*vb[1] - xb[1]*vb[0])/(xb[0]**2 + xb[1]**2)
    bh = _kerr(M6, spin); rg = bh.r_g; tg = rg / c
    out = []
    for sgn, tm in ((1.0, thi), (-1.0, -tlo)):
        o = bh.orbit(R=(r + 1.0)*rg, Theta=th*u.rad, Phi=(sgn*ph)*u.rad, vR=(sgn*rd)*c,
                     vTheta=(sgn*thd)*u.rad/tg, vPhi=phd*u.rad/tg)
        tau = tau_grid(a_g, e, max(tm, 1.0), npo)
        s = o.solve(tau_eval=tau*tg, method="dop853", rtol=rtol)
        if not s.success or len(s.t) != len(tau):
            raise RuntimeError("relatipy exact integration failed")
        T = (s.t/tg).decompose().value; R = (s.R/rg).decompose().value
        Th = s.Theta.to_value(u.rad); Ph = s.Phi.to_value(u.rad)
        vR = (s.vR/c).decompose().value; vTh = (s.vTheta*tg).to_value(u.rad); vPh = (s.vPhi*tg).to_value(u.rad)
        ut = np.asarray(s.ut, float)
        if sgn < 0: T = -T; Ph = -Ph; vR = -vR; vTh = -vTh
        rH = R - 1.0; sT, cT, sP, cP = np.sin(Th), np.cos(Th), np.sin(Ph), np.cos(Ph)
        n = np.stack([sT*cP, sT*sP, cT], 1); eth = np.stack([cT*cP, cT*sP, -sT], 1)
        eph = np.stack([-sP, cP, np.zeros_like(sP)], 1)
        Xb = rH[:, None]*n; Vb = vR[:, None]*n + rH[:, None]*(vTh[:, None]*eth + (sT*vPh)[:, None]*eph)
        out.append((T, Xb @ Rb.T, Vb @ Rb.T, ut))
    (Tf, Xf, Vf, uf), (Tb, Xbk, Vbk, ub) = out
    T = np.concatenate([Tb[::-1][:-1], Tf]); X = np.concatenate([Xbk[::-1][:-1], Xf])
    V = np.concatenate([Vbk[::-1][:-1], Vf]); UT = np.concatenate([ub[::-1][:-1], uf])
    T = np.ascontiguousarray(T); X = np.ascontiguousarray(X); V = np.ascontiguousarray(V)
    if return_ut: return T, X, V, Rb, spin, UT
    return T, X, V, Rb, spin

# ----------------------------------------------------------------------------------------
# data container and the three models
# ----------------------------------------------------------------------------------------
class Data:
    def __init__(self, t_ast, s_ast, t_rv, s_rv, t_ref):
        self.t_ast = np.asarray(t_ast, float); self.t_rv = np.asarray(t_rv, float)
        self.s_ast = np.broadcast_to(np.asarray(s_ast, float), self.t_ast.shape).copy()
        self.s_rv = np.broadcast_to(np.asarray(s_rv, float), self.t_rv.shape).copy()
        self.t_ref = t_ref
        self.tall = np.concatenate([self.t_ast, self.t_rv])
        self.na = len(self.t_ast)


def _offsets(p, d, al, de, rv):
    dt = d.t_ast - d.t_ref
    return (al + p[IX['x0']] + p[IX['vx0']]*dt, de + p[IX['y0']] + p[IX['vy0']]*dt, rv + p[IX['vz0']])


def model_T(p, d, npo=4000, lens=1.0, shap=1.0):
    """Truth and model (iii): exact Kerr dynamics + first-order time transfer function."""
    a_g, ty, us = geo(p)
    tg = (d.tall - p[IX['t0']]) / ty
    T, X, V, Rb, spin = exact_traj(p, tg.min() - 1.2*a_g*(1 + p[3]) - 100, tg.max() + 1.2*a_g*(1 + p[3]) + 100, npo)
    img, z, _, _ = obs_T(T, X, V, tg, Rb, spin, lens, shap, 0)
    na = d.na
    return _offsets(p, d, img[:na, 0]*us, img[:na, 1]*us, z[na:]*C_KMS)


def model_PN(p, d, spo=3000, spin_on=False, mult=False):
    """Model (i): harmonic 1PN x f_SP (+1.5PN spin-orbit) + Romer + straight projection."""
    a_g, ty, us = geo(p)
    tg = (d.tall - p[IX['t0']]) / ty
    x0, v0 = apo_state(a_g, p[3], p[4], p[5], p[6])
    S = p[IX['ast']]*axis(p[IX['ths']], p[IX['phs']]) if spin_on else np.zeros(3)
    pad = 1.2*a_g*(1 + p[3]) + 100
    tlo, thi = tg.min() - pad, tg.max() + pad
    P = 2*np.pi*a_g**1.5
    nmax = int(((thi - tlo)/P + 2) * spo * 3) + 1000
    T, X, V = pn_traj(x0, v0, tlo, thi, p[IX['fsp']], S, spin_on, 2*np.pi/spo, nmax)
    img, z = obs_pn(T, X, V, tg, p[IX['f']], mult)
    na = d.na
    return _offsets(p, d, img[:na, 0]*us, img[:na, 1]*us, z[na:]*C_KMS)


_KM = None


def model_R(p, d, rtol=1e-11):
    """Model (ii): relatipy KerrMcmcModel as is (BL-Cartesian state, straight line, Romer)."""
    global _KM
    from relatipy import KerrMcmcModel
    from astropy.constants import c
    if _KM is None: _KM = KerrMcmcModel()
    a_g, ty, us = geo(p)
    x0, v0 = apo_state(a_g, p[3], p[4], p[5], p[6])
    rg = RG_AU*p[0]*u.au
    spin = float(p[IX['ast']]); vec = tuple(axis(p[IX['ths']], p[IX['phs']])) if spin > 0 else (0.0, 0.0, 1.0)
    orb = dict(x=x0[0]*rg, y=x0[1]*rg, z=x0[2]*rg, vx=v0[0]*c, vy=v0[1]*c, vz=v0[2]*c, t=p[IX['t0']]*u.yr)
    al, de, vr = _KM.get_ra_dec_vr(kerr=dict(mass=p[0]*1e6*u.M_sun, spin=spin, vec=vec), orbit=orb,
                                   sol={"t_obs": d.tall*u.yr, "method": "dop853", "rtol": rtol},
                                   distance=p[1]*u.kpc)
    na = d.na
    return _offsets(p, d, np.asarray(al[:na])*1e6, np.asarray(de[:na])*1e6, np.asarray(vr[na:]))

# ----------------------------------------------------------------------------------------
# fitting and Fisher
# ----------------------------------------------------------------------------------------
def stack(obs, d):
    al, de, rv = obs
    return np.concatenate([al/d.s_ast, de/d.s_ast, rv/d.s_rv])


def steps(p, free):
    return np.array([STEP[n]*(abs(p[IX[n]]) if n in REL else 1.0) for n in free])


def jac(fun, p, free, d, h=None):
    """central-difference Jacobian of whitened observables wrt free params."""
    h = steps(p, free) if h is None else h
    cols = []
    for k, n in enumerate(free):
        pp = p.copy(); pm = p.copy(); pp[IX[n]] += h[k]; pm[IX[n]] -= h[k]
        cols.append((stack(fun(pp, d), d) - stack(fun(pm, d), d)) / (2*h[k]))
    return np.array(cols).T


def fisher_sigma(J):
    F = J.T @ J
    C = np.linalg.pinv(F, rcond=1e-15)
    return np.sqrt(np.clip(np.diag(C), 0, None)), C


def fit(fun, p0, free, d, ydata, max_nfev=60, verbose=0, bounds=None, x_scale=None):
    """least squares of model fun to noiseless data ydata (whitened stack) over free params."""
    idx = [IX[n] for n in free]
    def res(q):
        p = p0.copy(); p[idx] = q
        try:
            return stack(fun(p, d), d) - ydata
        except Exception:
            return np.full(ydata.shape, 1e6)
    def jj(q):
        p = p0.copy(); p[idx] = q
        return jac(fun, p, free, d)
    q0 = p0[idx].copy()
    if x_scale is None:
        x_scale = fisher_sigma(jj(q0))[0]
        x_scale = np.where(x_scale > 0, x_scale, 1.0)
    kw = dict(jac=jj, x_scale=x_scale, method='trf', ftol=1e-15, xtol=1e-15, gtol=1e-15,
              max_nfev=max_nfev, verbose=verbose)
    if bounds is not None: kw['bounds'] = bounds
    r = least_squares(res, q0, **kw)
    p = p0.copy(); p[idx] = r.x
    return p, r
