"""Wall time of the closed-form pericentre advance of KerrGeoPy on the work-precision sample.

KerrGeoPy does not integrate: it gives the fundamental frequencies in closed form, and the advance
per radial period is 2 pi Omega_phi / Omega_r - 2 pi (notebook section 3.2). Its error is not
comparable with the integrated ones, so only the time is measured, with the protocol of
work_precision.py (one core, one untimed warm-up, best of 3 repeats per orbit).
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

import argparse
import csv
import time
import warnings
from pathlib import Path

import numpy as np
import kerrgeopy as kg

CHI = 0.9
REPEATS = 3


def advance(p, e):
    """Orbit from (p, e, x=1) and its advance per radial period [rad]; both steps are timed."""
    orbit = kg.StableOrbit(CHI, p, e, 1.0)
    omega_r, _, omega_phi = orbit.fundamental_frequencies()
    return 2*np.pi * omega_phi / omega_r - 2*np.pi


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=100, help="first n orbits of the cell-57 sample")
    ap.add_argument("--out", default=str(Path(__file__).with_name("work-precision-kerrgeopy.csv")))
    args = ap.parse_args()
    warnings.simplefilter("ignore")

    # notebook cell 57: draw all 100, keep the first n
    rng = np.random.default_rng(2026)
    n_orbits = 100
    mc_a = np.exp(rng.uniform(np.log(30), np.log(1e4), n_orbits))
    mc_e = rng.uniform(0.05, 0.9, n_orbits)

    advance(mc_a[0] * (1 - mc_e[0]**2), mc_e[0])    # warm-up, untimed
    fields = ["library", "method", "tol", "orbit_id", "a", "e", "wall_s", "n_steps", "eps", "status"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for i in range(args.n):
            p = mc_a[i] * (1 - mc_e[i]**2)
            walls = []
            for _ in range(REPEATS):
                t0 = time.perf_counter()
                advance(p, mc_e[i])
                walls.append(time.perf_counter() - t0)
            w.writerow(dict(library="kerrgeopy", method="fundamental_frequencies", tol=np.nan, orbit_id=i,
                            a=mc_a[i], e=mc_e[i], wall_s=min(walls), n_steps=np.nan, eps=np.nan, status="ok"))
    print(f"wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
