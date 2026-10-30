"""Work-precision of pericentre precession: GYOTO.

Runs with the Python that has GYOTO (not the relatipy venv). Reads the initial states written by
`work_precision.py --dump-states` and writes rows with the columns of work-precision.csv.
Helpers copied verbatim from work_precision.py (notebook cells 4 and 45). Only xFill is timed.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

import argparse
import csv
import json
import time
import warnings

import numpy as np
import gyoto.core
import gyoto.std


# ---- notebook cell 4 -------------------------------------------------------
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


# ---- notebook cell 45 ------------------------------------------------------
def mean_advance(phases):
    return (phases[-1] - phases[0]) / (phases.size - 1) - 2*np.pi


# ---- setup -----------------------------------------------------------------
REPEATS = 3
INTEGRATORS = ("runge_kutta_fehlberg78", "runge_kutta_dopri5", "runge_kutta_cash_karp54",
               "runge_kutta_cash_karp54_classic")    # Legacy: no proper time, ignores RelTol
TOLERANCES = (1e-6, 1e-8, 1e-10, 1e-12, 1e-13)
DELTA_MAX_OVER_R = 1e10    # default 1 caps the step at r; lifted so that the tolerance sets the step


def make_star(spin, integrator, tol):
    """GYOTO Star in Kerr (Boyer-Lindquist, M = 1); RelTol = AbsTol = tol."""
    metric = gyoto.std.KerrBL()
    metric.spin(spin)
    star = gyoto.std.Star()
    star.metric(metric)
    star.integrator(integrator)
    star.set("RelTol", tol)
    star.set("AbsTol", tol)
    star.set("MaxIter", 10**8)
    star.set("DeltaMaxOverR", DELTA_MAX_OVER_R)
    return star


def run_gyoto(star, y0, t_end):
    star.initCoord(list(y0))
    t0 = time.perf_counter()
    star.xFill(t_end)
    return time.perf_counter() - t0


def nodes(star):
    """Stored states: proper time, r, u^r, phi, u^phi (dx/dtau)."""
    n = star.get_nelements()
    tau = np.zeros(n)
    star.get_tau(tau)
    t = np.zeros(n)
    star.get_t(t)
    r, th, ph, ut, ur, uth, uph = (np.zeros(n) for _ in range(7))
    star.getCoord(t, r, th, ph, ut, ur, uth, uph)
    order = np.argsort(tau)
    return tau[order], r[order], ur[order], ph[order], uph[order]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("states", help="JSON from work_precision.py --dump-states")
    ap.add_argument("--n", type=int, default=100, help="first n orbits")
    ap.add_argument("--out", required=True)
    ap.add_argument("--integrators", nargs="*", default=list(INTEGRATORS))
    args = ap.parse_args()
    warnings.simplefilter("ignore")

    sample = json.load(open(args.states))
    spin, orbits = sample["spin"], sample["orbits"][:args.n]
    settings = [(m, tol) for m in args.integrators for tol in TOLERANCES]

    t_start = time.perf_counter()
    o = orbits[0]
    for method, tol in settings:    # warm-up: one untimed solve per setting on orbit 0
        run_gyoto(make_star(spin, method, tol), o["y0"], o["t_end"])
    print("warm-up done", flush=True)

    fields = ["library", "method", "tol", "orbit_id", "a", "e", "wall_s", "n_steps", "eps", "status"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for method, tol in settings:
            ts = time.perf_counter()
            for o in orbits:
                row = dict(library="gyoto", method=method, tol=tol, orbit_id=o["orbit_id"], a=o["a"], e=o["e"],
                           wall_s=np.nan, n_steps=np.nan, eps=np.nan, status="ok")
                try:
                    star = make_star(spin, method, tol)
                    # GYOTO stops at a coordinate time: t_end, the coordinate time of tau_end
                    row["wall_s"] = min(run_gyoto(star, o["y0"], o["t_end"]) for _ in range(REPEATS))
                    tau, r, ur, ph, uph = nodes(star)
                    if tau[-1] < 0.999 * o["tau_end"]:
                        row["status"] = "short"
                    row["n_steps"] = tau.size - 1
                    measured = mean_advance(hermite_pericentres(tau, r, ur, ph, uph))
                    row["eps"] = abs(measured / o["exact"] - 1)
                except Exception as exc:  # record and continue
                    prev = "" if row["status"] == "ok" else row["status"] + ";"
                    row["status"] = prev + f"error:{type(exc).__name__}:{str(exc)[:80]}"
                w.writerow(row)
                fh.flush()
            print(f"gyoto {method} {tol}: {time.perf_counter() - ts:.1f} s", flush=True)
    print(f"total: {time.perf_counter() - t_start:.1f} s", flush=True)


if __name__ == "__main__":
    main()
