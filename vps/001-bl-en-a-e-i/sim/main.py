#!/usr/bin/env python3
"""
Kerr geodesic bank in Boyer-Lindquist coordinates for the Sgr A* spin study.

Each sampled orbit is defined by the Keplerian initial conditions
``(a, e, i)`` drawn by :mod:`orbit_sampling` (Latin-Hypercube design in
``(ln a, log10(1-e), i)`` with a periapsis-stability rejection), and is
integrated in the Kerr spacetime with the C-backed ``Radau2`` integrator.
The stored trajectory is the Boyer-Lindquist state

.. math::

    x^\\mu(\\tau) = (t, r, \\theta, \\phi), \\qquad
    u^\\mu(\\tau) = \\frac{\\mathrm{d}x^\\mu}{\\mathrm{d}\\tau},

sampled on a uniform proper-time grid
``tau_j = j / (n_samples - 1) * tau_max`` with ``tau_max = periods * T_kepler``
and ``T_kepler = 2 pi a^{3/2}``. Geometric units throughout: ``G = c = M = 1``.

The initial conditions ``(a, e, i)`` are *labels* of each trajectory, not
outputs: they are stored alongside every orbit in the shard file and in
``index.csv`` so that any trajectory can be traced back to the point of the
LHS design that generated it.

Design notes for the VPS (ARM64, 4 vCPU, no GPU):

* One process per orbit through a ``multiprocessing`` pool; BLAS is pinned to
  one thread so the workers do not oversubscribe the 4 vCPU.
* Results are written in shards, so partial output is usable while the job runs
  and the local 2-minute sync always finds complete files.
* ``out/checkpoint.json`` records the last completed shard; a resumed run skips
  everything already on disk.

Run it through ``run.sh`` (which bootstraps the environment)::

    vpsjob run --dir ./sim --cpu 350% -- ./run.sh main.py --n-orbits 100000
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import signal
import sys
import time
from pathlib import Path

# Must be set before numpy/scipy import their BLAS: with 4 vCPU and 3 worker
# processes, a 4-thread BLAS per process only makes the workers fight.
for _var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
             "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_var, "1")

import multiprocessing as mp  # noqa: E402

import numpy  # noqa: E402

import orbit_sampling as osamp  # noqa: E402

try:  # available inside the VPS job, absent on a laptop
    from vpsprogress import report
except ImportError:  # pragma: no cover

    def report(*_args, **_kwargs):
        return False


SCHEMA_VERSION = 1

STOP = False


def _sigterm(_signo, _frame):
    """``vpsjob stop`` sends SIGTERM: finish the current shard, then checkpoint."""
    global STOP
    STOP = True
    print("[sim] SIGTERM received, closing after the current shard", flush=True)


# --------------------------------------------------------------------- worker

_KERR = None
_POINTS = None
_CFG = None


def _worker_init(points, cfg):
    """Build the per-process Kerr metric once and pin the worker's signals."""
    global _KERR, _POINTS, _CFG
    # The parent decides when to stop; a worker that dies mid-task would lose
    # the whole shard.
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
    signal.signal(signal.SIGINT, signal.SIG_IGN)

    from relatipy.numeric.metrics import Kerr

    _POINTS = points
    _CFG = cfg
    _KERR = Kerr(1.0, cfg["spin"])


def _constants_of_motion(kerr, ys0):
    """
    Energy, axial angular momentum and Carter constant of an initial state.

    Parameters
    ----------
    kerr : relatipy.numeric.metrics.Kerr
        Kerr metric, spin ``a`` in geometric units.
    ys0 : ndarray, shape (8,)
        Initial Boyer-Lindquist state ``[x^mu, u^mu]``.

    Returns
    -------
    tuple of float
        ``(E, Lz, Q)`` per unit test-particle mass, in the usual physical
        convention: ``E < 1`` for a bound orbit, ``Lz > 0`` prograde, and
        ``Q = L^2 sin^2 i`` in the weak-field limit.

    Notes
    -----
    Same construction the ``Radau2`` integrator uses internally to define the
    constraints it projects onto, so these values label exactly the surface the
    trajectory is held on. The signs differ from the integrator's internal
    ``E0``/``Lz0`` because relatipy carries the metric with signature
    ``(+,-,-,-)``, which flips ``p_t`` and ``p_phi``.
    """
    g0 = kerr.metric(ys0[:4])
    p0 = g0 @ ys0[4:]
    E = p0[0]
    Lz = -p0[3]
    c2 = numpy.cos(ys0[2]) ** 2
    s2 = numpy.sin(ys0[2]) ** 2
    Q = p0[2] ** 2 + c2 * (kerr.a ** 2 * (1.0 - E ** 2) + Lz ** 2 / s2)
    return float(E), float(Lz), float(Q)


def _norm_residual(kerr, sol):
    """Worst ``|g_{mu nu} u^mu u^nu - 1|`` over a few samples of the path."""
    worst = 0.0
    for j in (0, sol.shape[1] // 2, sol.shape[1] - 1):
        u = sol[4:, j]
        worst = max(worst, abs(float(u @ kerr.metric(sol[:4, j]) @ u) - 1.0))
    return worst


def _integrate_one(idx):
    """
    Integrate one orbit of the design and return its trajectory plus labels.

    Parameters
    ----------
    idx : int
        Row of the LHS design to integrate.

    Returns
    -------
    dict
        ``t`` (float64), ``xsp`` and ``u`` (float32) trajectories plus the
        scalar metadata written to ``index.csv``. On failure the arrays are
        filled with NaN and ``status`` is 1.
    """
    from relatipy.numeric.coordinates import OrbitalElements

    a, e, inc = (float(v) for v in _POINTS[idx])
    ns = _CFG["samples"]
    rec = {
        "orbit_id": int(idx),
        "a": a,
        "e": e,
        "inc_rad": inc,
        "status": 1,
        "error": "",
        "period": float("nan"),
        "tau_max": float("nan"),
        "E": float("nan"),
        "Lz": float("nan"),
        "Q": float("nan"),
        "r_min": float("nan"),
        "r_max": float("nan"),
        "norm_err": float("nan"),
        "t": None,
        "xsp": None,
        "u": None,
    }

    t0 = time.perf_counter()
    try:
        oe = OrbitalElements(
            t=0.0, a=a, e=e, inc=numpy.degrees(inc), mass=1.0
        )  # OrbitalElements takes angles in degrees
        period = float(oe._get_period())
        tau_max = _CFG["periods"] * period
        taus = numpy.linspace(0.0, tau_max, ns)

        bl = oe.convert_to("BoyerLindquist", **_KERR.kwargs)
        ys0 = _KERR.get_4state_vector(bl)
        E, Lz, Q = _constants_of_motion(_KERR, ys0)

        # Private entry point on purpose: get_path() wraps the result in a
        # coordinate object that keeps dx^i/dt but drops the four-velocity,
        # and u^mu is exactly what a redshift/astrometry model needs.
        sol = numpy.asarray(
            _KERR.geodesic._get_path_from_4state_vector(
                ys0, taus,
                integrator=_CFG["integrator"], adaptative=False,
                period=period, rtol=_CFG["rtol"], atol=_CFG["atol"],
            ),
            dtype=float,
        )
        if sol.shape != (8, ns):
            raise RuntimeError(f"unexpected solution shape {sol.shape}")
        if not numpy.all(numpy.isfinite(sol)):
            raise RuntimeError("non-finite values in the solution")

        rec.update(
            status=0,
            period=period,
            tau_max=float(tau_max),
            E=E, Lz=Lz, Q=Q,
            r_min=float(sol[1].min()),
            r_max=float(sol[1].max()),
            norm_err=_norm_residual(_KERR, sol),
            # t stays float64: it reaches ~1e8 M and float32 would quantise it
            # to a few M, i.e. tens of seconds for Sgr A*.
            t=sol[0].copy(),
            xsp=sol[1:4].astype(numpy.float32),
            u=sol[4:8].astype(numpy.float32),
        )
    except Exception as exc:  # a bad orbit must not kill the shard
        rec["error"] = f"{type(exc).__name__}: {exc}"[:200]
        rec["t"] = numpy.full(ns, numpy.nan)
        rec["xsp"] = numpy.full((3, ns), numpy.nan, dtype=numpy.float32)
        rec["u"] = numpy.full((4, ns), numpy.nan, dtype=numpy.float32)

    rec["wall_s"] = time.perf_counter() - t0
    return rec


# ----------------------------------------------------------------------- I/O

INDEX_COLUMNS = [
    "orbit_id", "shard", "row",
    "a", "e", "inc_rad", "inc_deg", "r_p", "r_a",
    "period", "tau_max", "n_samples",
    "E", "Lz", "Q", "r_min", "r_max", "norm_err",
    "bound", "wall_s", "status", "error",
]


def _write_shard(out_dir, shard_id, recs, cfg):
    """Write one shard of trajectories as a single ``.npz``."""
    path = out_dir / "shards" / f"orbits_{shard_id:05d}.npz"
    tmp = path.with_suffix(".npz.tmp")
    # a file object, not a path: numpy.savez would append a second ".npz"
    with tmp.open("wb") as fh:
        _savez_shard(fh, recs, cfg)
    os.replace(tmp, path)
    return path


def _savez_shard(fh, recs, cfg):
    numpy.savez(
        fh,
        schema_version=numpy.int32(SCHEMA_VERSION),
        orbit_id=numpy.array([r["orbit_id"] for r in recs], dtype=numpy.int32),
        ci=numpy.array([[r["a"], r["e"], r["inc_rad"]] for r in recs], dtype=numpy.float64),
        tau_max=numpy.array([r["tau_max"] for r in recs], dtype=numpy.float64),
        period=numpy.array([r["period"] for r in recs], dtype=numpy.float64),
        constants=numpy.array([[r["E"], r["Lz"], r["Q"]] for r in recs], dtype=numpy.float64),
        status=numpy.array([r["status"] for r in recs], dtype=numpy.int8),
        bound=numpy.array([r["E"] < 1.0 for r in recs], dtype=numpy.int8),
        t=numpy.stack([r["t"] for r in recs]),
        xsp=numpy.stack([r["xsp"] for r in recs]),
        u=numpy.stack([r["u"] for r in recs]),
        spin=numpy.float64(cfg["spin"]),
        n_samples=numpy.int32(cfg["samples"]),
        periods=numpy.float64(cfg["periods"]),
    )


def _index_rows(recs, shard_id, cfg):
    for row, r in enumerate(recs):
        yield {
            "orbit_id": r["orbit_id"],
            "shard": shard_id,
            "row": row,
            "a": f"{r['a']:.10e}",
            "e": f"{r['e']:.10e}",
            "inc_rad": f"{r['inc_rad']:.10e}",
            "inc_deg": f"{numpy.degrees(r['inc_rad']):.6f}",
            "r_p": f"{r['a'] * (1.0 - r['e']):.10e}",
            "r_a": f"{r['a'] * (1.0 + r['e']):.10e}",
            "period": f"{r['period']:.10e}",
            "tau_max": f"{r['tau_max']:.10e}",
            "n_samples": cfg["samples"],
            "E": f"{r['E']:.12e}",
            "Lz": f"{r['Lz']:.12e}",
            "Q": f"{r['Q']:.12e}",
            "r_min": f"{r['r_min']:.10e}",
            "r_max": f"{r['r_max']:.10e}",
            "norm_err": f"{r['norm_err']:.3e}",
            # The Keplerian (a, e, i) -> Boyer-Lindquist map is Newtonian, so a
            # design point with a tight periapsis can come out relativistically
            # unbound (E >= 1). Those orbits are integrated and stored anyway;
            # this flag is what you filter on.
            "bound": int(r["E"] < 1.0) if numpy.isfinite(r["E"]) else 0,
            "wall_s": f"{r['wall_s']:.4f}",
            "status": r["status"],
            "error": r["error"],
        }


def _truncate_index(path, n_keep):
    """Drop index rows beyond the checkpoint so a resume cannot duplicate them."""
    if not path.exists():
        return
    with path.open("r", newline="", encoding="utf-8") as fh:
        rows = [row for row in csv.DictReader(fh) if int(row["orbit_id"]) < n_keep]
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=INDEX_COLUMNS)
        w.writeheader()
        w.writerows(rows)
        fh.flush()


def _config_signature(cfg):
    keys = ("n_orbits", "spin", "seed", "samples", "periods", "integrator",
            "rtol", "atol", "rp_floor", "shard_size")
    return json.dumps({k: cfg[k] for k in keys}, sort_keys=True)


# ---------------------------------------------------------------------- main


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    p.add_argument("--n-orbits", type=int, default=100_000,
                   help="orbits in the LHS design")
    p.add_argument("--spin", type=float, default=0.9,
                   help="dimensionless Kerr spin a/M assumed for Sgr A*")
    p.add_argument("--seed", type=int, default=20260731,
                   help="seed of the Latin-Hypercube design")
    p.add_argument("--samples", type=int, default=2048,
                   help="proper-time samples stored per orbit")
    p.add_argument("--periods", type=float, default=2.0,
                   help="Keplerian periods integrated per orbit")
    p.add_argument("--rp-floor", type=float, default=10.0,
                   help="conservative periapsis floor [M] of the stability cut")
    p.add_argument("--integrator", default="Radau2")
    p.add_argument("--rtol", type=float, default=1e-12)
    p.add_argument("--atol", type=float, default=1e-13)
    p.add_argument("--workers", type=int, default=3,
                   help="worker processes (the VPS has 4 vCPU)")
    p.add_argument("--shard-size", type=int, default=500,
                   help="orbits per .npz shard; also the resume granularity")
    p.add_argument("--out", default=os.environ.get("VPSJOB_OUT", "out"))
    p.add_argument("--fresh", action="store_true",
                   help="ignore any existing checkpoint and start over")
    args = p.parse_args()

    signal.signal(signal.SIGTERM, _sigterm)
    signal.signal(signal.SIGINT, _sigterm)

    cfg = {
        "n_orbits": args.n_orbits,
        "spin": args.spin,
        "seed": args.seed,
        "samples": args.samples,
        "periods": args.periods,
        "integrator": args.integrator,
        "rtol": args.rtol,
        "atol": args.atol,
        "rp_floor": args.rp_floor,
        "shard_size": args.shard_size,
    }

    out = Path(args.out)
    (out / "shards").mkdir(parents=True, exist_ok=True)
    ck_path = out / "checkpoint.json"
    index_path = out / "index.csv"

    # ---- LHS design (deterministic: same seed -> same points) --------------
    print(f"[sim] sampling {args.n_orbits} orbits (seed={args.seed}, "
          f"spin={args.spin}, rp_floor={args.rp_floor})", flush=True)
    t_sample = time.perf_counter()
    points, sample_stats = osamp.sample_orbits_lhs(
        args.n_orbits, spin=args.spin, seed=args.seed,
        rp_floor=args.rp_floor, use_separatrix=True, return_stats=True,
    )
    print(f"[sim] design ready in {time.perf_counter() - t_sample:.1f} s: "
          f"{sample_stats}", flush=True)
    numpy.save(out / "initial_conditions.npy", points)
    numpy.savetxt(out / "initial_conditions.csv", points, delimiter=",",
                  header="a,e,inc_rad", comments="")

    # ---- resume ------------------------------------------------------------
    n_shards = -(-args.n_orbits // args.shard_size)
    start_shard = 0
    if ck_path.exists() and not args.fresh:
        ck = json.loads(ck_path.read_text(encoding="utf-8"))
        if ck.get("config") != _config_signature(cfg):
            print("[sim] checkpoint belongs to a different configuration; "
                  "rerun with --fresh to discard it", flush=True)
            return 2
        start_shard = int(ck["next_shard"])
        _truncate_index(index_path, start_shard * args.shard_size)
        print(f"[sim] resuming at shard {start_shard}/{n_shards}", flush=True)

    (out / "run.json").write_text(json.dumps({
        "schema_version": SCHEMA_VERSION,
        "job_id": os.environ.get("VPSJOB_ID"),
        "started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "units": "geometric, G = c = M = 1",
        "coordinates": "Boyer-Lindquist (t, r, theta, phi)",
        "signature": "(+,-,-,-); timelike normalisation g_uu u^u u^v = +1",
        "initial_conditions": "Keplerian (a, e, i); i=0 prograde equatorial",
        "tau_grid": "tau_j = j/(n_samples-1) * tau_max, tau_max = periods * 2*pi*a**1.5",
        "n_shards": n_shards,
        "sample_stats": sample_stats,
        **cfg,
    }, indent=2, default=float), encoding="utf-8")

    # ---- integration -------------------------------------------------------
    new_index = not index_path.exists() or index_path.stat().st_size == 0
    ctx = mp.get_context("fork")
    done = start_shard * args.shard_size
    t_start = time.perf_counter()
    n_fail = 0
    n_unbound = 0
    rc = 0

    with index_path.open("a", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=INDEX_COLUMNS)
        if new_index:
            writer.writeheader()
            fh.flush()

        # close()/join() rather than the context manager: __exit__ calls
        # terminate(), which signals SIGTERM to workers that deliberately
        # ignore it, and the parent would then wait on them forever.
        pool = ctx.Pool(processes=args.workers, initializer=_worker_init,
                        initargs=(points, cfg))
        try:
            for shard_id in range(start_shard, n_shards):
                lo = shard_id * args.shard_size
                hi = min(lo + args.shard_size, args.n_orbits)
                t_shard = time.perf_counter()

                recs = pool.map(_integrate_one, range(lo, hi), chunksize=4)

                path = _write_shard(out, shard_id, recs, cfg)
                writer.writerows(_index_rows(recs, shard_id, cfg))
                fh.flush()
                os.fsync(fh.fileno())
                # checkpoint last: a crash before this only redoes one shard
                tmp = ck_path.with_suffix(".json.tmp")
                tmp.write_text(json.dumps({
                    "next_shard": shard_id + 1,
                    "n_shards": n_shards,
                    "orbits_done": hi,
                    "config": _config_signature(cfg),
                }), encoding="utf-8")
                os.replace(tmp, ck_path)

                done = hi
                n_fail += sum(r["status"] for r in recs)
                n_unbound += sum(1 for r in recs if not r["E"] < 1.0)
                dt_shard = time.perf_counter() - t_shard
                elapsed = time.perf_counter() - t_start
                rate = (done - start_shard * args.shard_size) / max(elapsed, 1e-9)
                mb = path.stat().st_size / 1e6

                print(f"[sim] shard {shard_id + 1}/{n_shards}  "
                      f"orbits {done}/{args.n_orbits}  "
                      f"{dt_shard:.1f} s  {rate:.1f} orb/s  "
                      f"{mb:.1f} MB  fails={n_fail}  unbound={n_unbound}", flush=True)
                report(step=done, total=args.n_orbits,
                       note=f"shard {shard_id + 1}/{n_shards}",
                       orbits_per_s=round(rate, 2),
                       gb_written=round(mb * n_shards / 1000.0, 3),
                       fails=n_fail, unbound=n_unbound)

                if STOP:
                    print(f"[sim] stopped after shard {shard_id + 1}, "
                          "checkpoint saved", flush=True)
                    rc = 143
                    break
        finally:
            pool.close()
            pool.join()

    total = time.perf_counter() - t_start
    print(f"[sim] {'stopped' if rc else 'finished'}: {done}/{args.n_orbits} orbits, "
          f"{n_fail} failures, {n_unbound} unbound (E >= 1), "
          f"{total / 3600:.2f} h", flush=True)
    return rc


if __name__ == "__main__":
    sys.exit(main())
