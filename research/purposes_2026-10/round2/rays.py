"""Exact Schwarzschild light propagation (M=1) for validation: primary-image impact parameter,
photon emission direction, exact redshift of a moving emitter, exact time of flight."""
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq


def u_at(b, psi):
    f = lambda ph, y: [y[1], -y[0] + 3.0*y[0]**2]
    s = solve_ivp(f, [0, psi], [0.0, 1.0/b], rtol=1e-12, atol=1e-16, method="DOP853")
    return s.y[0, -1], s.y[1, -1]


def impact(r, psi):
    """b of the primary image of a point at Schwarzschild radius r, angle psi from the observer axis."""
    bf = r*np.sin(psi)
    lo = max(0.5*bf, 5.3); hi = bf + 4*np.tan(psi/2) + 30
    return brentq(lambda b: u_at(b, psi)[0] - 1.0/r, lo, hi, xtol=1e-13, rtol=1e-15)


def _I(r1, r2, b):
    g = lambda r: 1.0/((1 - 2/r)*np.sqrt(max(1 - b*b*(1 - 2/r)/(r*r), 1e-300)))
    return quad(g, r1, r2, epsabs=1e-12, epsrel=1e-13, limit=400)[0]


def _I_turn(rmin, r2, b):
    # r = rmin + s^2 removes the 1/sqrt singularity at the turning point
    def g(s):
        r = rmin + s*s
        q = 1 - b*b*(1 - 2/r)/(r*r)
        return 2*s/((1 - 2/r)*np.sqrt(max(q, 1e-300)))
    return quad(g, 0.0, np.sqrt(r2 - rmin), epsabs=1e-12, epsrel=1e-13, limit=400)[0]


def tof(r, b, outgoing, D):
    """exact coordinate time of flight from r to D along impact parameter b."""
    rc = 1e4*max(r, b)
    tail = (D - rc) + 2*np.log(D/rc) + (4 + b*b/2)*(1/rc - 1/D)
    if outgoing:
        return _I(r, rc, b) + tail
    roots = np.roots([1.0, 0.0, -b*b, 2*b*b]); rmin = max(x.real for x in roots if abs(x.imag) < 1e-9)
    return _I_turn(rmin, r, b) + _I_turn(rmin, min(rc, 10*r), b) + (_I(10*r, rc, b) if rc > 10*r else 0) + tail


def exact_obs(xH, vH, D=1e9):
    """Exact (Schwarzschild, a=0) observables for a star at harmonic observer-frame position xH,
    coordinate velocity vH (observer at -Z).  Returns image vector (2: [Y, X] = [RA, Dec]) in M,
    1+z, time of flight (to distance D) and the photon direction."""
    rH = np.linalg.norm(xH); r = rH + 1.0; n = xH/rH
    psi = np.arccos(-n[2])
    b = impact(r, psi)
    perp = np.array([xH[0], xH[1]]); bh = perp/np.linalg.norm(perp)
    ub, dub = u_at(b, psi)
    outgoing = dub > 0
    sa = b*np.sqrt(1 - 2/r)/r; ca = (1 if outgoing else -1)*np.sqrt(max(1 - sa*sa, 0))
    o = np.array([0.0, 0.0, -1.0]); et = o - (o @ n)*n; et /= np.linalg.norm(et)
    ng = ca*n + sa*et
    rd = vH @ n; ndot = (vH - rd*n)/rH
    vloc = rd/(1 - 2/r)*n + r*ndot/np.sqrt(1 - 2/r)
    g = 1/np.sqrt(1 - vloc @ vloc)
    onepz = g*(1 - vloc @ ng)/np.sqrt(1 - 2/r)
    img = b*np.array([bh[1], bh[0]])
    # plane-wave observer: integrate to the plane at harmonic distance D-1 on the axis
    return img, onepz, tof(r, b, outgoing, np.sqrt((D - 1.0)**2 + b*b) + 1.0), ng


def first_order(xH, vH, ut, D=1e9):
    """first-order time transfer observables (the truth generator's formulas), plus finite-D T1."""
    r = np.linalg.norm(xH); Z = xH[2]; q = r - Z
    img = np.array([xH[1], xH[0]])*(1 + 2/q)
    rd = vH @ (xH/r)
    onepz = ut*(1 + vH[2] - 2*(rd - vH[2])/q)
    rO = D - 1.0
    T1 = rO + Z + 2*np.log(2*rO) - 2*np.log(q)      # plane-wave limit of R + 2ln((re+ro+R)/(re+ro-R))
    return img, onepz, T1
