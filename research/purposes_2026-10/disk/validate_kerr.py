"""Validation (i): numba PN (1PN + 1.5PN SO, harmonic) vs relatipy EXACT Kerr nodal precession
for an S301 orbit (a = 16010 r_g, e = 0.9825), spin chi along z, orbit inclined 60 deg.
Spin flip realised by rotating the ICs by pi about x (relatipy needs spin >= 0).
Odd-in-chi part isolates LT (removes the chi^2 quadrupole nodal term)."""
import numpy as np, time
from astropy import units as u
from astropy.constants import c
from relatipy import Kerr
import diskdyn as dd, fc

a, e, chi = 16010.0, 0.9825, 0.7
inc, Om, om = 60 * fc.DEG, 30 * fc.DEG, 90 * fc.DEG
x0, v0 = fc.kepler_state(a, e, inc, Om, om, np.pi)        # apocentre
P = 2 * np.pi * a ** 1.5
NORB = 3
Rx = np.diag([1.0, -1.0, -1.0])

def nodes_from(X, V):
    r = np.linalg.norm(X, axis=1)
    L = np.cross(X, V); L /= np.linalg.norm(L, axis=1)[:, None]
    # apocentre samples: local maxima of r
    idx = [k for k in range(1, len(r) - 1) if r[k] >= r[k - 1] and r[k] >= r[k + 1]]
    Lh = L[idx]
    Omn = np.unwrap(np.arctan2(Lh[:, 0], -Lh[:, 1])); incn = np.arccos(Lh[:, 2])
    return Omn, incn

def exact(sign):
    M = 4.297e6 * u.M_sun
    bh = Kerr(mass=M, spin=chi)
    RG_S = (bh.r_g / c).to_value(u.s)
    xx, vv = (x0, v0) if sign > 0 else (Rx @ x0, Rx @ v0)
    orb = bh.orbit(x=xx[0] * bh.r_g, y=xx[1] * bh.r_g, z=xx[2] * bh.r_g, vx=vv[0] * c, vy=vv[1] * c, vz=vv[2] * c)
    tau = np.linspace(0, (NORB + 0.3) * P * RG_S, 6000 * NORB) * u.s
    sol = orb.solve(tau_span=(tau[0], tau[-1]), tau_eval=tau, method="dop853", rtol=1e-13)
    assert sol.success, sol.message
    X = (sol.xyz / bh.r_g).to_value(u.dimensionless_unscaled).T
    V = (sol.vxyz / c).to_value(u.dimensionless_unscaled).T
    if X.shape[0] == 3 and X.shape[1] != 3: X = X.T; V = V.T
    if sign < 0: X = X @ Rx.T; V = V @ Rx.T
    return nodes_from(X, V)

def pn(sign, pnord=2):
    S = np.array([0, 0, sign * chi])
    ds = fc.s_orbit(e) / 4000
    T, X, V, A = dd.integrate(x0, v0, 0.0, (NORB + 0.3) * P, ds, 4000 * (NORB + 2), S, 0.0, np.array([0, 0, 1.0]), 1.0, fc.TR, fc.TZ, fc.TP, pnord)
    ts = np.linspace(0, (NORB + 0.3) * P, 6000 * NORB)
    Xs, Vs = dd.states_at(T, X, V, A, ts)
    return nodes_from(Xs, Vs)

t0 = time.time()
Ep, Ip = exact(+1); Em, Im = exact(-1)
Pp, PIp = pn(+1); Pm, PIm = pn(-1)
dE_odd = 0.5 * (np.diff(Ep) - np.diff(Em)); dE_even = 0.5 * (np.diff(Ep) + np.diff(Em))
dP_odd = 0.5 * (np.diff(Pp) - np.diff(Pm))
p = a * (1 - e * e)
print(f"S301 a={a} e={e} chi={chi} inc={inc/fc.DEG:.0f}deg  ({time.time()-t0:.1f}s)")
print("per-orbit dOmega (rad): exact odd ", np.array2string(dE_odd, precision=5), " exact even", np.array2string(dE_even, precision=3))
print("                        PN odd    ", np.array2string(dP_odd, precision=5))
print(f"rel diff PN/exact - 1 (odd): {np.mean(dP_odd)/np.mean(dE_odd)-1:+.2e};  analytic 4 pi chi/p^1.5 = {4*np.pi*chi/p**1.5:.5e}")
print(f"exact even part vs -3 pi chi^2 cos(i)/p^2 = {-3*np.pi*chi**2*np.cos(inc)/p**2:.3e}")
print(f"inclination drift per orbit: exact {np.mean(np.diff(Ip)):+.1e}, PN {np.mean(np.diff(PIp)):+.1e} rad")
