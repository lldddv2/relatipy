"""Work-precision of pericentre precession: Relatipy vs PyGRO.

Helpers copied verbatim from docs/validation/01-Theoretical-Validations.ipynb
(cells 4, 6, 43, 45, 57). Only the solve call is timed.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

import argparse
import csv
import json
import logging
import signal
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import astropy.units as u
import astropy.constants as const
import relatipy as rp
import pygro

logging.getLogger().setLevel(logging.WARNING)


# ---- notebook cell 4 -------------------------------------------------------
def time_scale(bh):
    return (bh.r_g / const.c).to(u.s)


def hermite_pericentres(tau, r, u_r, phi, u_phi):
    phi = np.unwrap(phi)
    phases = []
    for i in np.where((u_r[:-1] < 0) & (u_r[1:] >= 0))[0]:
        h = tau[i+1] - tau[i]
        c3 = 2*r[i] - 2*r[i+1] + h*(u_r[i] + u_r[i+1])
        c2 = -3*r[i] + 3*r[i+1] - h*(2*u_r[i] + u_r[i+1])
        roots = np.roots([3*c3, 2*c2, h*u_r[i]])
        s = next(x.real for x in roots
                 if abs(x.imag) < 1e-12 and -1e-9 <= x.real <= 1 + 1e-9 and 6*c3*x.real + 2*c2 > 0)
        phases.append(
            (2*s**3 - 3*s**2 + 1) * phi[i]
            + (s**3 - 2*s**2 + s) * h * u_phi[i]
            + (-2*s**3 + 3*s**2) * phi[i+1]
            + (s**3 - s**2) * h * u_phi[i+1]
        )
    return np.array(phases)


# ---- notebook cell 6 -------------------------------------------------------
_pygro_engines = {}


def pygro_engine(chi, integrator="rkf78"):
    """Notebook cell 6, with the integrator as an argument (one engine per spin and integrator)."""
    if (chi, integrator) not in _pygro_engines:
        metric = pygro.Metric(
            name="Kerr spacetime",
            coordinates=["t", "r", "theta", "phi"],
            line_element=(
                "-(1-2*M*r/(Sigma(r,theta)))*dt**2"
                "+Sigma(r,theta)/Delta(r)*dr**2"
                "+Sigma(r,theta)*dtheta**2"
                "+(r**2+a**2+2*M*r*a**2/Sigma(r,theta)*sin(theta)**2)*sin(theta)**2*dphi**2"
                "-4*M*r*a*sin(theta)**2/Sigma(r,theta)*dt*dphi"
            ),
            transform=[
                "t",
                "sqrt(r**2+a**2)*sin(theta)*cos(phi)",
                "sqrt(r**2+a**2)*sin(theta)*sin(phi)",
                "r*cos(theta)",
            ],
            Sigma="r**2+a**2*cos(theta)**2",
            Delta="r**2-2*M*r+a**2",
            M=1,
            a=chi,
        )
        engine = pygro.GeodesicEngine(metric, backend="autowrap", integrator=integrator)
        engine.set_stopping_criterion("r**2-2*M*r+a**2 > 1e-4", "horizon")
        _pygro_engines[chi, integrator] = engine
    return _pygro_engines[chi, integrator]


def pygro_geodesic(bh, orbit, integrator):
    """Untimed part of notebook pygro_solve: build Geodesic with initial state."""
    T0 = time_scale(bh)
    engine = pygro_engine(bh.spin, integrator)
    geo = pygro.Geodesic("time-like", engine, verbose=False)
    geo.initial_x = [0.0, (orbit.R / bh.r_g).to_value(u.one),
                     orbit.Theta.to_value(u.rad), orbit.Phi.to_value(u.rad)]
    geo.initial_u = [float(orbit.ut.value), (orbit.uR / const.c).to_value(u.one),
                     (orbit.uTheta * T0).to_value(u.rad), (orbit.uPhi * T0).to_value(u.rad)]
    return engine, geo


# ---- notebook cells 43, 45 -------------------------------------------------
def exact_precession(chi, E, Lz, p, e, n=4096):
    u_p = 1 / p
    x = Lz - chi * E
    cos_phase = np.cos((np.arange(n) + 0.5) * 2*np.pi / n)
    dphi_dchi = u_p**0.5 * (x + chi*E - 2*u_p*x*(1 + e*cos_phase)) / (
        np.sqrt(1 + u_p**2 * x**2 * (e**2 - 2*e*cos_phase - 3))
        * (1 - 2*u_p*(1 + e*cos_phase) + chi**2 * u_p**2 * (1 + e*cos_phase)**2)
    )
    return 2*np.pi * np.mean(dphi_dchi) - 2*np.pi


def mean_advance(phases):
    return (phases[-1] - phases[0]) / (phases.size - 1) - 2*np.pi


# ---- setup -----------------------------------------------------------------
TIMEOUT = 120.0
N_PERIODS = 100
REPEATS = 3


class Timeout(Exception):
    pass


def _alarm(signum, frame):
    raise Timeout()


signal.signal(signal.SIGALRM, _alarm)


def prepare(bh, a_M, e):
    """Orbit + exact advance + tau_end, as in notebook measure_precession."""
    a = a_M * bh.r_g
    p = a * (1 - e**2)
    p_M = (p / bh.r_g).to_value(u.one)
    orbit = bh.orbit(p=p, e=e, x=1.0)
    T0 = time_scale(bh)
    tau_end = N_PERIODS * (orbit.get_keplerian_period() / T0).to_value(u.one)
    exact = exact_precession(bh.spin, orbit.get_E(), orbit.get_Lz(), p_M, e)
    return orbit, tau_end, exact


def run_relatipy(bh, orbit, tau_end, method, rtol):
    T0 = time_scale(bh)
    span = (0 * T0, tau_end * T0)
    t0 = time.perf_counter()
    sol = orbit.solve(tau_span=span, method=method, rtol=rtol)
    wall = time.perf_counter() - t0
    return wall, sol


def eps_relatipy(bh, sol, exact):
    T0 = time_scale(bh)
    measured = mean_advance(hermite_pericentres(
        (sol.tau / T0).to_value(u.one), (sol.R / bh.r_g).to_value(u.one),
        (sol.uR / const.c).to_value(u.one), sol.Phi.to_value(u.rad),
        (sol.uPhi * T0).to_value(u.rad)))
    return abs(measured / exact - 1)


def run_pygro(bh, orbit, tau_end, integrator, goal):
    engine, geo = pygro_geodesic(bh, orbit, integrator)
    signal.setitimer(signal.ITIMER_REAL, TIMEOUT)
    try:
        t0 = time.perf_counter()
        engine.integrate(geo, tau_end, 1.0, accuracy_goal=goal, precision_goal=goal)
        wall = time.perf_counter() - t0
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    return wall, geo


def eps_pygro(geo, exact):
    measured = mean_advance(hermite_pericentres(geo.tau, geo.x[:, 1], geo.u[:, 1], geo.x[:, 3], geo.u[:, 3]))
    return abs(measured / exact - 1)


RELATIPY_METHODS = ("projection_radau", "radau", "dop853", "dp45")
RELATIPY_RTOLS = (1e-6, 1e-8, 1e-10, 1e-12, 1e-13)
PYGRO_INTEGRATORS = ("rkf78", "dp853", "rkf45", "dp45", "ck45")
PYGRO_GOALS = (6, 8, 10, 12)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=20, help="first n orbits of the cell-57 sample")
    ap.add_argument("--out", default=str(Path(__file__).with_name("workprec.csv")))
    ap.add_argument("--relatipy", nargs="*", default=list(RELATIPY_METHODS), help="Relatipy methods")
    ap.add_argument("--pygro", nargs="*", default=["rkf78", "dp853"], help="PyGRO integrators")
    ap.add_argument("--dump-states", help="write the initial states of the sample to this JSON file and exit"
                    " (input of work_precision_gyoto.py)")
    args = ap.parse_args()
    settings = (
        [("relatipy", m, r) for m in args.relatipy for r in RELATIPY_RTOLS]
        + [("pygro", m, g) for m in args.pygro for g in PYGRO_GOALS]
    )
    warnings.simplefilter("ignore")

    bh = rp.Kerr(mass=4.3e6*u.Msun, spin=0.9)
    # notebook cell 57: draw all 100, keep the first n
    rng = np.random.default_rng(2026)
    n_orbits = 100
    mc_a = np.exp(rng.uniform(np.log(30), np.log(1e4), n_orbits))
    mc_e = rng.uniform(0.05, 0.9, n_orbits)

    t_start = time.perf_counter()
    for integrator in args.pygro:
        t0 = time.perf_counter()
        pygro_engine(0.9, integrator)
        print(f"pygro engine compile ({integrator}): {time.perf_counter() - t0:.2f} s", flush=True)

    prepared = [prepare(bh, mc_a[i], mc_e[i]) for i in range(args.n)]
    if args.dump_states:
        T0 = time_scale(bh)

        def coordinate_end(orbit, tau_end):
            """Coordinate time [M] at tau_end (GYOTO stops at a coordinate time)."""
            sol = orbit.solve(tau_span=(0 * T0, tau_end * T0), method="dop853", rtol=1e-13)
            return (sol.t[-1] / T0).to_value(u.one)

        states = [
            dict(
                orbit_id=i,
                a=mc_a[i],
                e=mc_e[i],
                y0=[0.0, (orbit.R / bh.r_g).to_value(u.one), orbit.Theta.to_value(u.rad),
                    orbit.Phi.to_value(u.rad), float(orbit.ut.value), (orbit.uR / const.c).to_value(u.one),
                    (orbit.uTheta * T0).to_value(u.rad), (orbit.uPhi * T0).to_value(u.rad)],
                tau_end=tau_end,
                t_end=coordinate_end(orbit, tau_end),
                exact=exact,
            )
            for i, (orbit, tau_end, exact) in enumerate(prepared)
        ]
        Path(args.dump_states).write_text(json.dumps(dict(spin=bh.spin, orbits=states), indent=1))
        return

    # warm-up: one untimed solve per setting on orbit 0
    orbit0, tau0, _ = prepared[0]
    for lib, method, tol in settings:
        if lib == "relatipy":
            run_relatipy(bh, orbit0, tau0, method, tol)
        else:
            try:
                run_pygro(bh, orbit0, tau0, method, tol)
            except Timeout:
                pass
    print("warm-up done", flush=True)

    fields = ["library", "method", "tol", "orbit_id", "a", "e", "wall_s", "n_steps", "eps", "status"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for lib, method, tol in settings:
            ts = time.perf_counter()
            for i in range(args.n):
                orbit, tau_end, exact = prepared[i]
                row = dict(library=lib, method=method, tol=tol, orbit_id=i, a=mc_a[i], e=mc_e[i],
                           wall_s=np.nan, n_steps=np.nan, eps=np.nan, status="ok")
                walls, result = [], None
                try:
                    for _ in range(REPEATS):
                        if lib == "relatipy":
                            wall, result = run_relatipy(bh, orbit, tau_end, method, tol)
                        else:
                            wall, result = run_pygro(bh, orbit, tau_end, method, tol)
                        walls.append(wall)
                        if wall > TIMEOUT:
                            break
                    row["wall_s"] = min(walls)
                    if min(walls) > TIMEOUT:
                        row["status"] = "timeout"
                    if lib == "relatipy":
                        row["n_steps"] = result.integration.n_steps
                        if result.status != 0:
                            row["status"] = f"solver_status_{result.status}"
                        row["eps"] = eps_relatipy(bh, result, exact)
                    else:
                        row["n_steps"] = len(result.tau) - 1
                        if getattr(result, "exit", None) not in (None, "done"):
                            row["status"] = f"exit_{result.exit}"
                        row["eps"] = eps_pygro(result, exact)
                except Timeout:
                    row["status"] = "timeout"
                    row["wall_s"] = min(walls) if walls else TIMEOUT
                except Exception as exc:  # record and continue
                    prev = "" if row["status"] == "ok" else row["status"] + ";"
                    row["status"] = prev + f"error:{type(exc).__name__}:{str(exc)[:80]}"
                    if walls:
                        row["wall_s"] = min(walls)
                w.writerow(row)
                fh.flush()
            print(f"{lib} {method} {tol}: {time.perf_counter() - ts:.1f} s", flush=True)
    print(f"total: {time.perf_counter() - t_start:.1f} s", flush=True)


if __name__ == "__main__":
    main()
