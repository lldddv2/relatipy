"""Fast PN test-particle model (numba).  G=c=M=1, lengths r_g, times r_g/c.
Sundman time dt/ds = r^1.5, classical RK4 in s, cubic-Hermite dense output,
Newton-solved Romer delay.  1PN terms optionally scaled by f_SP (GRAVITY style),
1.5PN spin-orbit term a=(2/r^3)[3(n.S)(n x v)+(v x S)] (sign validated vs exact Kerr)."""
import numpy as np
from numba import njit

@njit(cache=True)
def _acc(x, v, order, fsp, S):
    r = np.sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]); n = x/r
    a = -n/(r*r)
    if order >= 1:
        vv = v[0]*v[0]+v[1]*v[1]+v[2]*v[2]; rd = v[0]*n[0]+v[1]*n[1]+v[2]*n[2]
        a = a + fsp*(1.0/(r*r))*(n*(4.0/r - vv) + 4.0*rd*v)
    if order >= 2:
        nS = n[0]*S[0]+n[1]*S[1]+n[2]*S[2]
        nxv = np.array([n[1]*v[2]-n[2]*v[1], n[2]*v[0]-n[0]*v[2], n[0]*v[1]-n[1]*v[0]])
        vxS = np.array([v[1]*S[2]-v[2]*S[1], v[2]*S[0]-v[0]*S[2], v[0]*S[1]-v[1]*S[0]])
        a = a + (2.0/(r*r*r))*(3.0*nS*nxv + vxS)
    return a

@njit(cache=True)
def _f(y, order, fsp, S):
    x = y[0:3]; v = y[3:6]
    r = np.sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]); g = r**1.5
    out = np.empty(7)
    out[0:3] = g*v; out[3:6] = g*_acc(x, v, order, fsp, S); out[6] = g
    return out

@njit(cache=True)
def integrate(x0, v0, t_end, order, fsp, S, ds):
    nmax = int(1e7)
    y = np.empty(7); y[0:3] = x0; y[3:6] = v0; y[6] = 0.0
    T = [0.0]; X = [x0.copy()]; V = [v0.copy()]
    k = 0
    while y[6] < t_end and k < nmax:
        k1 = _f(y, order, fsp, S); k2 = _f(y+0.5*ds*k1, order, fsp, S)
        k3 = _f(y+0.5*ds*k2, order, fsp, S); k4 = _f(y+ds*k3, order, fsp, S)
        y = y + ds/6.0*(k1+2*k2+2*k3+k4)
        T.append(y[6]); X.append(y[0:3].copy()); V.append(y[3:6].copy()); k += 1
    n = len(T); Ta = np.empty(n); Xa = np.empty((n,3)); Va = np.empty((n,3))
    for i in range(n):
        Ta[i] = T[i]; Xa[i] = X[i]; Va[i] = V[i]
    return Ta, Xa, Va

@njit(cache=True)
def _hermite(T, X, V, t):
    lo = 0; hi = len(T)-1
    while hi-lo > 1:
        mid = (lo+hi)//2
        if T[mid] <= t: lo = mid
        else: hi = mid
    h = T[hi]-T[lo]; s = (t-T[lo])/h
    h00 = 2*s**3-3*s**2+1; h10 = s**3-2*s**2+s; h01 = -2*s**3+3*s**2; h11 = s**3-s**2
    x = h00*X[lo]+h10*h*V[lo]+h01*X[hi]+h11*h*V[hi]
    dh00 = (6*s**2-6*s)/h; dh10 = 3*s**2-4*s+1; dh01 = (-6*s**2+6*s)/h; dh11 = 3*s**2-2*s
    v = dh00*X[lo]+dh10*V[lo]+dh01*X[hi]+dh11*V[hi]
    return x, v

@njit(cache=True)
def sky(T, X, V, ghat, e1, e2, t_obs):
    """Return (n,2) sky positions in r_g and emission times; Romer via Newton."""
    n = len(t_obs); out = np.empty((n,2)); te = np.empty(n)
    for i in range(n):
        t = t_obs[i]
        for it in range(30):
            x, v = _hermite(T, X, V, t)
            F = t + (x[0]*ghat[0]+x[1]*ghat[1]+x[2]*ghat[2]) - t_obs[i]
            dF = 1.0 + (v[0]*ghat[0]+v[1]*ghat[1]+v[2]*ghat[2])
            dt = F/dF; t -= dt
            if abs(dt) < 1e-9: break
        x, v = _hermite(T, X, V, t)
        out[i,0] = x[0]*e1[0]+x[1]*e1[1]+x[2]*e1[2]
        out[i,1] = x[0]*e2[0]+x[1]*e2[1]+x[2]*e2[2]
        te[i] = t
    return out, te

@njit(cache=True)
def integrate_pre(x0, v0, t_end, order, fsp, S, ds, nmax):
    """Preallocated variant of integrate()."""
    T = np.empty(nmax); X = np.empty((nmax,3)); V = np.empty((nmax,3))
    y = np.empty(7); y[0:3] = x0; y[3:6] = v0; y[6] = 0.0
    T[0] = 0.0; X[0] = x0; V[0] = v0; k = 0
    while y[6] < t_end and k < nmax-1:
        k1 = _f(y, order, fsp, S); k2 = _f(y+0.5*ds*k1, order, fsp, S)
        k3 = _f(y+0.5*ds*k2, order, fsp, S); k4 = _f(y+ds*k3, order, fsp, S)
        y = y + ds/6.0*(k1+2*k2+2*k3+k4); k += 1
        T[k] = y[6]; X[k] = y[0:3]; V[k] = y[3:6]
    return T[:k+1], X[:k+1], V[:k+1]
