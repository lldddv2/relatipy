"""Forward models for the S301 + friends mocks.

A forward model is any callable

    model(elements, t_ast, t_rv, bh) -> (ra_mas, dec_mas, v_kms)

where
    elements : dict with a_mas, e, i_deg, Omega_deg, omega_deg, t_p, element_epoch
               (the keys of stars.json; extra keys are ignored),
    t_ast    : 1-D array of astrometric epochs [decimal yr, observed time],
    t_rv     : 1-D array of spectroscopic epochs [decimal yr, observed time],
    bh       : dict with M_1e6Msun, R0_kpc, chi and spin_vec (3 components in
               observer axes (Dec, RA, away); normalised by the model).
It returns RA and Dec offsets [mas] at t_ast and the line-of-sight velocity
[km/s] at t_rv.  A Kerr + disk model will later be plugged in with the same
signature (e.g. a closure carrying the disk parameters).

The default, KerrForward, wraps relatipy.KerrMcmcModel.get_ra_dec_vr (exact
Kerr geodesic, homogeneous: no extended mass).  Its observation model is
relatipy's: straight-line light with Romer delay, no Shapiro delay or lensing;
v is the redshift-equivalent velocity (Doppler + gravitational redshift).

Initial state: the osculating Kepler elements are converted by relatipy at the
element (osculation) epoch t_el, with the true anomaly at t_el obtained from
Kepler's equation with GM = M and the published t_p.  Using t_el = t_p for
elements that are osculating at apocentre (GRAVITY 2024 Table B.1) costs ~1 mas
on S2 (validation.md).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
DEG = np.pi / 180.0


def load_stars(path: Path | str = HERE / "stars.json") -> dict:
    with open(path) as fh:
        return json.load(fh)


def period_yr(a_mas: float, M6: float, R0_kpc: float) -> float:
    """Kepler period [yr] for a [mas] at distance R0 around M [1e6 Msun]."""
    a_au = a_mas * R0_kpc  # mas * kpc = au
    return float(np.sqrt(a_au ** 3 / (M6 * 1e6)))


def true_anomaly(t: float, el: dict, M6: float, R0_kpc: float) -> float:
    """Keplerian true anomaly [rad] at epoch t from (a, e, t_p)."""
    e = el["e"]
    P = period_yr(el["a_mas"], M6, R0_kpc)
    Mm = np.mod(2 * np.pi * (t - el["t_p"]) / P + np.pi, 2 * np.pi) - np.pi
    # E - e sin E is monotonic on [-pi, pi]: bracketed root, safe at e ~ 0.98
    E = brentq(lambda x: x - e * np.sin(x) - Mm, -np.pi, np.pi, xtol=1e-15, rtol=1e-15)
    return float(2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2)))


def orbit_pole(el: dict) -> np.ndarray:
    """Unit orbital angular momentum in observer axes (Dec, RA, away)."""
    i, O = el["i_deg"] * DEG, el["Omega_deg"] * DEG
    return np.array([np.sin(O) * np.sin(i), -np.cos(O) * np.sin(i), np.cos(i)])


SPIN_DIRS = {  # (beta, lambda) [deg] w.r.t. the S301 frame (node N, L x N, L); as the prototype
    "par": (0, 0), "perpN": (90, 0), "perpM": (90, 90), "b45l45": (45, 45),
    "b120l200": (120, 200), "b70l310": (70, 310),
}


def spin_vector(name_or_vec, ref_el: dict) -> np.ndarray:
    """Spin unit vector in observer axes from a SPIN_DIRS name (relative to the orbital
    frame of ref_el, normally S301A) or an explicit 3-vector (Dec, RA, away)."""
    if isinstance(name_or_vec, str):
        b, l = (x * DEG for x in SPIN_DIRS[name_or_vec])
        L = orbit_pole(ref_el)
        O = ref_el["Omega_deg"] * DEG
        N = np.array([np.cos(O), np.sin(O), 0.0])
        Mv = np.cross(L, N)
        v = np.sin(b) * np.cos(l) * N + np.sin(b) * np.sin(l) * Mv + np.cos(b) * L
    else:
        v = np.asarray(name_or_vec, dtype=float)
    return v / np.linalg.norm(v)


class KerrForward:
    """Default forward model: relatipy KerrMcmcModel.get_ra_dec_vr (exact Kerr)."""

    def __init__(self, method: str = "dop853", rtol: float = 1e-11, atol: float = 1e-13,
                 max_steps: int = 1_000_000, use_element_epoch: bool = True):
        from relatipy import KerrMcmcModel
        self._m = KerrMcmcModel()
        self.sol_opts = dict(method=method, rtol=rtol, atol=atol, max_steps=max_steps)
        self.use_element_epoch = use_element_epoch
        self.n_retries = 0

    def __call__(self, el: dict, t_ast, t_rv, bh: dict):
        import astropy.units as u
        t_ast = np.atleast_1d(np.asarray(t_ast, float))
        t_rv = np.atleast_1d(np.asarray(t_rv, float))
        t_all = np.concatenate([t_ast, t_rv])
        if t_all.size == 0:
            return np.empty(0), np.empty(0), np.empty(0)
        M6, R0 = bh["M_1e6Msun"], bh["R0_kpc"]
        chi = float(bh.get("chi", 0.0))
        vec = bh.get("spin_vec", (0.0, 0.0, 1.0))
        t_el = el["element_epoch"] if self.use_element_epoch else el["t_p"]
        f = true_anomaly(t_el, el, M6, R0)
        call = dict(
            kerr=dict(mass=M6 * 1e6 * u.Msun, spin=chi, vec=tuple(vec)),
            orbit=dict(a=el["a_mas"] * R0 * u.au, e=el["e"], inc=el["i_deg"] * u.deg,
                       Omega=el["Omega_deg"] * u.deg, omega=el["omega_deg"] * u.deg,
                       f=f * u.rad, t=t_el * u.yr),
            distance=R0 * u.kpc,
        )
        ra, de, v = self._eval(call, t_all)
        n = t_ast.size
        return ra[:n] * 1e3, de[:n] * 1e3, v[n:]

    # relatipy 1.2.0 occasionally raises "Kerr geodesic integration failed" for isolated
    # observed epochs (step-sequence dependent, e.g. S2, chi 0.7, vec (0,0,1),
    # t_obs = 2010.9470999999608 at rtol 1e-11); the same epoch succeeds at another
    # tolerance.  Retry with alternative tolerances (all agree to < 0.1 uas, < 0.01 km/s),
    # then epoch by epoch.  n_retries counts the fallbacks.
    RETRY_TOLS = ((1e-12, 1e-14), (1e-10, 1e-12), (3e-12, 3e-14), (3e-11, 3e-13))

    def _one(self, call, t, opts):
        import astropy.units as u
        ra, de, v = self._m.get_ra_dec_vr(sol=dict(t_obs=np.asarray(t) * u.yr, **opts), **call)
        return np.array(ra, float), np.array(de, float), np.array(v, float)

    def _eval(self, call, t):
        from relatipy.geodesic import IntegrationError
        try:
            return self._one(call, t, self.sol_opts)
        except IntegrationError:
            pass
        for rt, at in self.RETRY_TOLS:
            self.n_retries += 1
            try:
                return self._one(call, t, dict(self.sol_opts, rtol=rt, atol=at))
            except IntegrationError:
                continue
        out = np.empty((3, t.size))
        for k, tk in enumerate(t):
            for opts in [self.sol_opts] + [dict(self.sol_opts, rtol=rt, atol=at) for rt, at in self.RETRY_TOLS]:
                try:
                    out[:, k] = [x[0] for x in self._one(call, [tk], opts)]
                    break
                except IntegrationError:
                    self.n_retries += 1
            else:
                raise IntegrationError(f"Kerr integration failed at t_obs = {tk!r} for all tolerances")
        return out[0], out[1], out[2]
