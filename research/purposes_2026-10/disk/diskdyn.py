"""diskdyn: numba kernels for the S301 + friends full-orbit spin/disk forecast.

Units: G = c = M_BH = 1 (lengths r_g, times r_g/c).  Observer frame: X = north (Dec),
Y = east (RA), Z = away from the observer (receding RV > 0); same convention as the
GRAVITY orbital elements (Thiele-Innes).
Dynamics (test particle, harmonic gauge): Newton + 1PN + 1.5PN spin-orbit
  a_SO = (2/r^3)[3(n.S)(n x v) + (v x S)]   (validated vs exact Kerr, see validate_kerr.py)
+ Newtonian razor-thin power-law disk (Sigma ~ R^alpha, rin..rout, ring softening H = soft*R').
2PN (and the spin quadrupole) are omitted CONSISTENTLY in truth and model.
Integrator: Sundman time dt/ds = r^1.5, classical RK4 in s (fixed ds => smooth in parameters),
quintic-Hermite dense output (x, v, a stored), Newton-solved Romer delay t_obs = t_e + Z(t_e).
Disk force: 2D table of r*g_R, r*g_z (unit disk mass, rout = 1) on a stretched grid in
(w, v): ln(r/rout) = q0 sinh(w), beta = b0 sinh(v) (beta = latitude above the disk plane),
Catmull-Rom bicubic interpolation; multipole (l = 0, 2, 4) beyond r = rtmax*rout.
"""
import math
import numpy as np
from numba import njit, prange

PI = math.pi


# ----------------------------------------------------------------- elliptic integrals / rings
@njit(cache=True)
def ellipKE(m):
    """Complete elliptic integrals K(m), E(m), m = k^2 < 1, by the AGM."""
    a = 1.0; b = math.sqrt(1.0 - m); s = 0.5 * m; pw = 0.5
    for _ in range(60):
        c = 0.5 * (a - b)
        an = 0.5 * (a + b); b = math.sqrt(a * b); a = an
        pw *= 2.0; s += pw * c * c
        if c < 1e-17 * a:
            break
    K = PI / (2.0 * a)
    return K, K * (1.0 - s)


@njit(cache=True)
def ring_g(R, z, b, H):
    """Acceleration (g_R, g_z) at (R, z) of a unit-mass ring of radius b in z = 0,
    softened |x - x'|^2 -> |x - x'|^2 + H^2 (G = 1)."""
    zeta2 = z * z + H * H
    D2 = (R + b) ** 2 + zeta2
    D = math.sqrt(D2)
    m = 4.0 * R * b / D2
    den = (R - b) ** 2 + zeta2
    K, E = ellipKE(m)
    gz = -2.0 * z * E / (PI * D * den)
    if R > 1e-12 * b:
        gR = (E * (b * b - R * R + zeta2) / den - K) / (PI * R * D)
    else:
        gR = 0.0
    return gR, gz


@njit(cache=True)
def disk_g_direct(R, z, rin, rout, alpha, soft, xg, wg):
    """Unit-mass disk acceleration by piecewise Gauss-Legendre over ring radius b, with
    breakpoints clustered geometrically around b = R on the softening scale.
    R < 0 (table ghost cells beyond the pole) = odd continuation of g_R, even of g_z."""
    if R < 0.0:
        gR, gz = disk_g_direct(-R, z, rin, rout, alpha, soft, xg, wg)
        return -gR, gz
    if abs(alpha + 2.0) < 1e-12:
        C = 1.0 / math.log(rout / rin)
    else:
        C = (alpha + 2.0) / (rout ** (alpha + 2.0) - rin ** (alpha + 2.0))
    h0 = math.sqrt(z * z + (soft * R) ** 2)
    if h0 < 1e-7 * rout:
        h0 = 1e-7 * rout
    bp = np.empty(256); nb = 0
    bp[0] = rin; bp[1] = rout; nb = 2
    if R > rin and R < rout:
        bp[nb] = R; nb += 1
    d = 0.25 * h0
    span = max(abs(R - rin), abs(rout - R))
    while d < span and nb < 250:
        x1 = R - d
        if x1 > rin and x1 < rout:
            bp[nb] = x1; nb += 1
        x2 = R + d
        if x2 > rin and x2 < rout:
            bp[nb] = x2; nb += 1
        d *= 2.0
    bs = np.sort(bp[:nb])
    gR = 0.0; gz = 0.0
    for k in range(nb - 1):
        lo = bs[k]; hi = bs[k + 1]
        if hi <= lo:
            continue
        half = 0.5 * (hi - lo); mid = 0.5 * (hi + lo)
        for j in range(xg.shape[0]):
            b = mid + half * xg[j]
            w = wg[j] * half * C * b ** (alpha + 1.0)
            r1, r2 = ring_g(R, z, b, soft * b)
            gR += w * r1; gz += w * r2
    return gR, gz


@njit(parallel=True, cache=True)
def build_table(nw, w0, dw, nv, v0, dv, q0, b0, rin, alpha, soft, xg, wg):
    """Table of r~ g_R, r~ g_z for a unit-mass disk with rout = 1 (r~ = r/rout)."""
    TR = np.empty((nw, nv)); TZ = np.empty((nw, nv))
    for i in prange(nw):
        rt = math.exp(q0 * math.sinh(w0 + i * dw))
        for j in range(nv):
            beta = b0 * math.sinh(v0 + j * dv)
            R = rt * math.cos(beta); z = rt * math.sin(beta)
            gR, gz = disk_g_direct(R, z, rin, 1.0, alpha, soft, xg, wg)
            TR[i, j] = rt * gR; TZ[i, j] = rt * gz
    return TR, TZ


# ----------------------------------------------------------------- table interpolation
@njit(cache=True)
def _crw(t):
    t2 = t * t; t3 = t2 * t
    return (-t3 + 2.0 * t2 - t) * 0.5, (3.0 * t3 - 5.0 * t2 + 2.0) * 0.5, \
           (-3.0 * t3 + 4.0 * t2 + t) * 0.5, (t3 - t2) * 0.5


@njit(cache=True)
def _row(T, i, j, b0, b1, b2, b3):
    return T[i, j - 1] * b0 + T[i, j] * b1 + T[i, j + 1] * b2 + T[i, j + 2] * b3


@njit(cache=True)
def multipole_T(rt, sb, cb, J2, J4, H2):
    """r~ g_R, r~ g_z of the unit disk (rout = 1) from l = 0, 2, 4 + softening monopole term.
    rt = r/rout, sb = sin(beta) = z/r, cb = cos(beta) = R/r."""
    mu = sb
    # Phi = -1/r - J2 P2/r^3 - J4 P4/r^5 + H2/(2 r^3);  g = -grad Phi
    P2 = 0.5 * (3 * mu * mu - 1); dP2 = 3 * mu
    P4 = (35 * mu ** 4 - 30 * mu * mu + 3) / 8.0; dP4 = (140 * mu ** 3 - 60 * mu) / 8.0
    r = rt
    dPhidr = 1.0 / r ** 2 + 3 * J2 * P2 / r ** 4 + 5 * J4 * P4 / r ** 6 - 1.5 * H2 / r ** 4
    dPhidmu = -J2 * dP2 / r ** 3 - J4 * dP4 / r ** 5
    # dmu/dR = -z R / r^3 = -sb cb / r ; dmu/dz = R^2/r^3 = cb^2 / r
    gR = -(dPhidr * cb + dPhidmu * (-sb * cb / r))
    gz = -(dPhidr * sb + dPhidmu * (cb * cb / r))
    return rt * gR, rt * gz


@njit(cache=True)
def disk_T(rt, z, R, TR, TZ, tp):
    """Interpolated (r~ g_R, r~ g_z), unit disk with rout = 1. tp = table params array:
    [w0, dw, v0, dv, q0, b0, rtmax, J2, J4, H2, rtmin]"""
    w0 = tp[0]; dw = tp[1]; v0 = tp[2]; dv = tp[3]; q0 = tp[4]; b0 = tp[5]
    if rt >= tp[6]:
        r = math.sqrt(R * R + z * z)
        return multipole_T(rt, z / r, R / r, tp[7], tp[8], tp[9])
    if rt < tp[10]:
        rt = tp[10]
    w = math.asinh(math.log(rt) / q0)
    beta = math.atan2(z, R)
    v = math.asinh(beta / b0)
    nw = TR.shape[0]; nv = TR.shape[1]
    fw = (w - w0) / dw; i = int(math.floor(fw))
    if i < 1: i = 1
    if i > nw - 3: i = nw - 3
    t = fw - i
    fv = (v - v0) / dv; j = int(math.floor(fv))
    if j < 1: j = 1
    if j > nv - 3: j = nv - 3
    u = fv - j
    a0, a1, a2, a3 = _crw(t); c0, c1, c2, c3 = _crw(u)
    sR = (a0 * _row(TR, i - 1, j, c0, c1, c2, c3) + a1 * _row(TR, i, j, c0, c1, c2, c3)
          + a2 * _row(TR, i + 1, j, c0, c1, c2, c3) + a3 * _row(TR, i + 2, j, c0, c1, c2, c3))
    sZ = (a0 * _row(TZ, i - 1, j, c0, c1, c2, c3) + a1 * _row(TZ, i, j, c0, c1, c2, c3)
          + a2 * _row(TZ, i + 1, j, c0, c1, c2, c3) + a3 * _row(TZ, i + 2, j, c0, c1, c2, c3))
    return sR, sZ


@njit(cache=True)
def disk_acc(x, nd, md, rout, TR, TZ, tp):
    """3D disk acceleration at x (code units). nd = unit disk normal, md = Mdisk/M_BH,
    rout = outer radius [r_g]."""
    r = math.sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2])
    z = x[0] * nd[0] + x[1] * nd[1] + x[2] * nd[2]
    Rx0 = x[0] - z * nd[0]; Rx1 = x[1] - z * nd[1]; Rx2 = x[2] - z * nd[2]
    R = math.sqrt(Rx0 * Rx0 + Rx1 * Rx1 + Rx2 * Rx2)
    sR, sZ = disk_T(r / rout, z / rout, R / rout, TR, TZ, tp)
    fac = md / (rout * r)
    out = np.empty(3)
    if R > 0.0:
        cr = fac * sR / R
        out[0] = cr * Rx0 + fac * sZ * nd[0]
        out[1] = cr * Rx1 + fac * sZ * nd[1]
        out[2] = cr * Rx2 + fac * sZ * nd[2]
    else:
        out[0] = fac * sZ * nd[0]; out[1] = fac * sZ * nd[1]; out[2] = fac * sZ * nd[2]
    return out


# ----------------------------------------------------------------- dynamics
@njit(cache=True)
def acc(x, v, S, md, nd, rout, TR, TZ, tp, pn):
    """Total acceleration.  pn = 0: Newton only; 1: +1PN; 2: +1PN+1.5PN SO; 3: Newton + SO only."""
    r = math.sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2])
    n0 = x[0] / r; n1 = x[1] / r; n2 = x[2] / r
    ir2 = 1.0 / (r * r)
    a0 = -n0 * ir2; a1 = -n1 * ir2; a2 = -n2 * ir2
    if pn == 1 or pn == 2:
        vv = v[0] * v[0] + v[1] * v[1] + v[2] * v[2]
        rd = v[0] * n0 + v[1] * n1 + v[2] * n2
        c = 4.0 / r - vv
        a0 += ir2 * (n0 * c + 4.0 * rd * v[0])
        a1 += ir2 * (n1 * c + 4.0 * rd * v[1])
        a2 += ir2 * (n2 * c + 4.0 * rd * v[2])
    if pn >= 2:
        nS = n0 * S[0] + n1 * S[1] + n2 * S[2]
        f = 2.0 * ir2 / r
        a0 += f * (3.0 * nS * (n1 * v[2] - n2 * v[1]) + (v[1] * S[2] - v[2] * S[1]))
        a1 += f * (3.0 * nS * (n2 * v[0] - n0 * v[2]) + (v[2] * S[0] - v[0] * S[2]))
        a2 += f * (3.0 * nS * (n0 * v[1] - n1 * v[0]) + (v[0] * S[1] - v[1] * S[0]))
    out = np.empty(3)
    out[0] = a0; out[1] = a1; out[2] = a2
    if md != 0.0:
        ad = disk_acc(x, nd, md, rout, TR, TZ, tp)
        out[0] += ad[0]; out[1] += ad[1]; out[2] += ad[2]
    return out


@njit(cache=True)
def _deriv(y, S, md, nd, rout, TR, TZ, tp, pn, out, aout):
    x = y[0:3]; v = y[3:6]
    r = math.sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]); g = r * math.sqrt(r)
    a = acc(x, v, S, md, nd, rout, TR, TZ, tp, pn)
    for k in range(3):
        out[k] = g * v[k]; out[3 + k] = g * a[k]; aout[k] = a[k]
    out[6] = g


@njit(cache=True)
def integrate(x0, v0, t0, t_end, ds, nmax, S, md, nd, rout, TR, TZ, tp, pn):
    """RK4 in Sundman time from t0 to >= t_end.  Returns T, X, V, A (stored at nodes)."""
    T = np.empty(nmax); X = np.empty((nmax, 3)); V = np.empty((nmax, 3)); A = np.empty((nmax, 3))
    y = np.empty(7); y[0:3] = x0; y[3:6] = v0; y[6] = t0
    k1 = np.empty(7); k2 = np.empty(7); k3 = np.empty(7); k4 = np.empty(7)
    yt = np.empty(7); aa = np.empty(3); dum = np.empty(3)
    k = 0
    _deriv(y, S, md, nd, rout, TR, TZ, tp, pn, k1, aa)
    T[0] = t0; X[0] = x0; V[0] = v0; A[0] = aa
    while y[6] < t_end and k < nmax - 1:
        for q in range(7): yt[q] = y[q] + 0.5 * ds * k1[q]
        _deriv(yt, S, md, nd, rout, TR, TZ, tp, pn, k2, dum)
        for q in range(7): yt[q] = y[q] + 0.5 * ds * k2[q]
        _deriv(yt, S, md, nd, rout, TR, TZ, tp, pn, k3, dum)
        for q in range(7): yt[q] = y[q] + ds * k3[q]
        _deriv(yt, S, md, nd, rout, TR, TZ, tp, pn, k4, dum)
        for q in range(7): y[q] = y[q] + ds / 6.0 * (k1[q] + 2.0 * k2[q] + 2.0 * k3[q] + k4[q])
        k += 1
        _deriv(y, S, md, nd, rout, TR, TZ, tp, pn, k1, aa)
        T[k] = y[6]; X[k] = y[0:3]; V[k] = y[3:6]; A[k] = aa
    return T[:k + 1], X[:k + 1], V[:k + 1], A[:k + 1]


@njit(cache=True)
def _locate(T, t):
    lo = 0; hi = T.shape[0] - 1
    if t <= T[0]: return 0
    if t >= T[hi]: return hi - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if T[mid] <= t: lo = mid
        else: hi = mid
    return lo


@njit(cache=True)
def dense(T, X, V, A, t):
    """Quintic Hermite position, its derivative for velocity."""
    lo = _locate(T, t); hi = lo + 1
    h = T[hi] - T[lo]; s = (t - T[lo]) / h
    s2 = s * s; s3 = s2 * s; s4 = s3 * s; s5 = s4 * s
    H0 = 1 - 10 * s3 + 15 * s4 - 6 * s5; H1 = s - 6 * s3 + 8 * s4 - 3 * s5
    H2 = 0.5 * (s2 - 3 * s3 + 3 * s4 - s5); H3 = 10 * s3 - 15 * s4 + 6 * s5
    H4 = -4 * s3 + 7 * s4 - 3 * s5; H5 = 0.5 * (s3 - 2 * s4 + s5)
    D0 = -30 * s2 + 60 * s3 - 30 * s4; D1 = 1 - 18 * s2 + 32 * s3 - 15 * s4
    D2 = 0.5 * (2 * s - 9 * s2 + 12 * s3 - 5 * s4); D3 = 30 * s2 - 60 * s3 + 30 * s4
    D4 = -12 * s2 + 28 * s3 - 15 * s4; D5 = 0.5 * (3 * s2 - 8 * s3 + 5 * s4)
    x = np.empty(3); v = np.empty(3)
    hh = h * h
    for q in range(3):
        x[q] = (H0 * X[lo, q] + H1 * h * V[lo, q] + H2 * hh * A[lo, q]
                + H3 * X[hi, q] + H4 * h * V[hi, q] + H5 * hh * A[hi, q])
        v[q] = (D0 * X[lo, q] + D1 * h * V[lo, q] + D2 * hh * A[lo, q]
                + D3 * X[hi, q] + D4 * h * V[hi, q] + D5 * hh * A[hi, q]) / h
    return x, v


@njit(cache=True)
def observe(T, X, V, A, tobs):
    """For each observing time: (X_north, Y_east) [r_g], RV [c] = v_Z + v^2/2 + 1/r, and t_emit.
    Romer: t_obs = t_e + Z(t_e)."""
    n = tobs.shape[0]
    out = np.empty((n, 3)); te = np.empty(n)
    for i in range(n):
        x, v = dense(T, X, V, A, tobs[i])
        t = tobs[i] - x[2]
        for it in range(50):
            x, v = dense(T, X, V, A, t)
            F = t + x[2] - tobs[i]
            dt = F / (1.0 + v[2])
            t -= dt
            if abs(dt) < 1e-13 * (1.0 + abs(t)):
                break
        x, v = dense(T, X, V, A, t)
        r = math.sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2])
        out[i, 0] = x[0]; out[i, 1] = x[1]
        out[i, 2] = v[2] + 0.5 * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) + 1.0 / r
        te[i] = t
    return out, te


@njit(cache=True)
def states_at(T, X, V, A, ts):
    n = ts.shape[0]; xs = np.empty((n, 3)); vs = np.empty((n, 3))
    for i in range(n):
        x, v = dense(T, X, V, A, ts[i]); xs[i] = x; vs[i] = v
    return xs, vs
