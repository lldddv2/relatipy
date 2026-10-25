"""Measure the five-star likelihood and a few emcee steps, and project a full run.

Builds the five-star likelihood of ``bench_kerr_mcmc_multi.py`` (S2, S29,
S38, S55, S301; 100-150 epochs per star) on noise-free data, without a force
(40 parameters) and with the fiducial ThinDisk (1e3 Msun, ``r_out`` at the
S301 apocentre, 45 parameters with ``disk_mass`` and the two normal angles),
then:

1. reports the astrometric and velocity error of the chosen tolerance
   against a ``rtol = 1e-12`` reference at the truth;
2. times likelihood calls (median);
3. runs a few emcee iterations with ``2 x n_dim`` walkers, in one process
   and with a ``--processes`` pool;
4. projects the wall time of ``--walkers-per-parameter`` walkers per
   parameter for ``--steps`` steps from the measured per-evaluation times.

The table of the disk is built once per process (about 1 s; pool workers
build their own). The whole script takes about one minute. Run from the
repository root after building the extensions::

    PYTHONPATH=src python tools/benchmarks/bench_mcmc_projection.py

Measured on the shared 4-core development container (2026-10-25, dop853,
4-process spawn pool, 3 iterations of 2 x n_dim walkers):

===========  ======  =====  ===============  =========  =============  =========
rtol/atol    force   n_dim  tolerance error  call (ms)  serial eval/s  4-proc /s
===========  ======  =====  ===============  =========  =============  =========
1e-9/1e-12   none    40     0.006 uas        21.8       48             123
1e-9/1e-12   disk    43     2.6 uas          13.4       78             182
1e-10/1e-13  none    40     0.0004 uas       22.2       45             121
1e-10/1e-13  disk    43     0.07 uas         13.5       73             177
===========  ======  =====  ===============  =========  =============  =========

The forced path is faster than the homogeneous one because its active-force
right-hand side is the generated Kerr kernel. ``rtol = 1e-10`` costs almost
nothing here (the step count is set by the epoch spacing) and keeps the
disk's numerical error far below 50 uas astrometry, so it is the default.
Projection for 4 walkers per parameter (172 walkers with the disk) and
2e4 steps, 3.4e6 evaluations: about 13 h in one process and 5.4 h with four
processes (speed-up 2.4 on the shared machine); with an autocorrelation time
of 300 steps that is about 1.1e4 independent samples.
"""

from __future__ import annotations

import argparse
import multiprocessing
import statistics
import sys
import time
from pathlib import Path

import numpy as np
from astropy import units as u

sys.path.insert(0, str(Path(__file__).parent))

from bench_kerr_mcmc_multi import ELEMENTS, EPOCHS, SHARED  # noqa: E402

from relatipy import (  # noqa: E402
    KerrMcmcGaussianLikelihood, KerrMcmcMultiModel, PerturbedKerrMcmcMultiModel,
    ThinDisk,
)

DISK = ThinDisk(mass=1e3 * u.M_sun, r_out=1362.5 * u.au)
DISK_VALUES = {"force.disk_mass": 1e3, "force.normal_theta": 1.0, "force.normal_phi": 0.5}
SPIN = {"spin_dec": 0.3, "spin_ra": -0.2, "spin_away": 0.5}


def build(with_disk, rtol, atol):
    """Return the model, the truth vector and its likelihood on noise-free data."""
    astrometry = {s: np.linspace(*span, n) for s, (span, n, _, _) in EPOCHS.items()}
    spectroscopy = {s: np.linspace(*span, n) for s, (_, _, span, n) in EPOCHS.items()}
    if with_disk:
        model = PerturbedKerrMcmcMultiModel(
            astrometry, spectroscopy, reference_epoch=2000.0, force=DISK
        )
    else:
        model = KerrMcmcMultiModel(astrometry, spectroscopy, reference_epoch=2000.0)
    model.set_solver(method="dop853", rtol=rtol, atol=atol)
    values = dict(SHARED) | SPIN | (DISK_VALUES if with_disk else {})
    for star, (t_p, a, e, inc, node, peri) in ELEMENTS.items():
        values |= {
            f"{star}.t_p": t_p, f"{star}.a": a, f"{star}.e": e,
            f"{star}.inc": np.deg2rad(inc), f"{star}.Omega": np.deg2rad(node),
            f"{star}.omega": np.deg2rad(peri),
        }
    truth = model.pack(values)
    data = {}
    for star, (alpha, delta, v_los) in model(truth).items():
        data[star] = {
            "alpha": alpha, "delta": delta, "v_los": v_los,
            "alpha_err": np.full(alpha.size, 5e-5),
            "delta_err": np.full(delta.size, 5e-5),
            "v_los_err": np.full(v_los.size, 10.0),
        }
    return model, truth, KerrMcmcGaussianLikelihood(model, data)


def tolerance_error(model, truth, rtol, atol):
    """Largest astrometric (uas) and velocity (km/s) change against rtol 1e-12."""
    model.set_solver(method="dop853", rtol=1e-12, atol=1e-15)
    reference = model(truth)
    model.set_solver(method="dop853", rtol=rtol, atol=atol)
    current = model(truth)
    astrometry = max(
        float(np.max(np.abs(current[s][i] - reference[s][i]))) for s in current for i in (0, 1)
    )
    velocity = max(float(np.max(np.abs(current[s][2] - reference[s][2]))) for s in current)
    return 1e6 * astrometry, velocity


def median_call(function, argument, repeat):
    function(argument)
    samples = []
    for _ in range(repeat):
        start = time.perf_counter()
        function(argument)
        samples.append(time.perf_counter() - start)
    return statistics.median(samples)


def emcee_seconds_per_evaluation(likelihood, truth, iterations, pool=None):
    """Seconds per likelihood evaluation over a short emcee run."""
    import emcee

    dimension = truth.size
    walkers = 2 * dimension
    rng = np.random.default_rng(1)
    scale = np.maximum(np.abs(truth), 1e-3) * 1e-9
    start = truth + scale * rng.standard_normal((walkers, dimension))
    sampler = emcee.EnsembleSampler(walkers, dimension, likelihood, pool=pool)
    begin = time.perf_counter()
    sampler.run_mcmc(start, iterations, progress=False)
    elapsed = time.perf_counter() - begin
    # emcee evaluates the initial ensemble once, then every walker per step.
    return elapsed / (walkers * (iterations + 1)), walkers


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--rtol", type=float, default=1e-10)
    parser.add_argument("--atol", type=float, default=1e-13)
    parser.add_argument("--repeat", type=int, default=15)
    parser.add_argument("--iterations", type=int, default=3)
    parser.add_argument("--processes", type=int, default=4)
    parser.add_argument("--walkers-per-parameter", type=int, default=4)
    parser.add_argument("--steps", type=int, default=20000)
    parser.add_argument("--autocorrelation", type=float, default=300.0)
    args = parser.parse_args()

    print(f"dop853 rtol={args.rtol:g} atol={args.atol:g}; "
          f"{args.iterations} emcee iterations with 2 x n_dim walkers")
    for with_disk in (False, True):
        label = "ThinDisk" if with_disk else "no force"
        model, truth, likelihood = build(with_disk, args.rtol, args.atol)
        astrometric, velocity = tolerance_error(model, truth, args.rtol, args.atol)
        call = median_call(likelihood, truth, args.repeat)
        serial, walkers = emcee_seconds_per_evaluation(likelihood, truth, args.iterations)
        context = multiprocessing.get_context("spawn")
        with context.Pool(args.processes) as pool:
            pool.map(likelihood, [truth] * args.processes)  # warm the workers
            pooled, _ = emcee_seconds_per_evaluation(
                likelihood, truth, args.iterations, pool=pool
            )
        total_walkers = args.walkers_per_parameter * truth.size
        evaluations = total_walkers * args.steps
        print(f"{label}: {truth.size} parameters")
        print(f"  tolerance error vs rtol 1e-12: {astrometric:.3g} uas, "
              f"{velocity:.3g} km/s")
        print(f"  likelihood call: {1e3 * call:.2f} ms (median of {args.repeat})")
        print(f"  emcee, {walkers} walkers: {1e3 * serial:.2f} ms/eval serial "
              f"({1 / serial:.0f} eval/s), {1e3 * pooled:.2f} ms/eval with "
              f"{args.processes} processes ({1 / pooled:.0f} eval/s, "
              f"x{serial / pooled:.2f})")
        print(f"  projection, {total_walkers} walkers x {args.steps} steps = "
              f"{evaluations:.3g} evaluations: {evaluations * serial / 3600:.1f} h serial, "
              f"{evaluations * pooled / 3600:.1f} h with {args.processes} processes; "
              f"about {total_walkers * args.steps / args.autocorrelation:.3g} "
              f"independent samples for tau = {args.autocorrelation:g} steps")


if __name__ == "__main__":
    main()
