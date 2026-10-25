"""Time five-star perturbed and homogeneous multi-star MCMC calls, interleaved.

Uses the stars, elements and epochs of ``bench_kerr_mcmc_multi.py`` and
compares one call of:

- ``KerrMcmcMultiModel`` (homogeneous);
- ``PerturbedKerrMcmcMultiModel`` with the built-in fiducial thin disk
  (``r_out = 1362.5 au``) of zero mass and of ``--disk-mass`` solar masses;
- a Plummer sphere of ``--plummer-mass`` solar masses written as a run-time
  ``Potential`` (needs SymPy and a C compiler; skipped otherwise).

Per-star times come from single-star calls of the same models. Calls are
interleaved and the median of ``--repeat`` rounds is reported. Run from the
repository root after building the extensions::

    PYTHONPATH=src python tools/benchmarks/bench_perturbed_kerr_mcmc_multi.py
"""

from __future__ import annotations

import argparse
import statistics
import sys
import time
from pathlib import Path

import numpy as np
from astropy import units as u

sys.path.insert(0, str(Path(__file__).parent))

from bench_kerr_mcmc_multi import ELEMENTS, EPOCHS, SHARED  # noqa: E402

from relatipy import (  # noqa: E402
    KerrMcmcMultiModel, PerturbedKerrMcmcMultiModel, Potential, ThinDisk,
)


def values(spin, force_values=None):
    named = dict(SHARED) | dict(zip(("spin_dec", "spin_ra", "spin_away"), spin))
    for star, (t_p, a, e, inc, node, peri) in ELEMENTS.items():
        named |= {
            f"{star}.t_p": t_p, f"{star}.a": a, f"{star}.e": e,
            f"{star}.inc": np.deg2rad(inc), f"{star}.Omega": np.deg2rad(node),
            f"{star}.omega": np.deg2rad(peri),
        }
    return named | (force_values or {})


def models(stars, args):
    astrometry = {s: np.linspace(*EPOCHS[s][0], EPOCHS[s][1]) for s in stars}
    spectroscopy = {s: np.linspace(*EPOCHS[s][2], EPOCHS[s][3]) for s in stars}
    radius = args.plummer_radius * u.mpc
    disk = ThinDisk(mass=0 * u.M_sun, r_out=1362.5 * u.au)
    cases = {"homogeneous": (KerrMcmcMultiModel(
        astrometry, spectroscopy, reference_epoch=2000.0), {})}
    for label, mass in (("disk zero", 0.0), ("disk", args.disk_mass)):
        cases[label] = (PerturbedKerrMcmcMultiModel(
            astrometry, spectroscopy, reference_epoch=2000.0, force=disk,
        ), {"force.disk_mass": mass, "force.normal_theta": 0.0,
            "force.normal_phi": 0.0})
    try:
        potential = Potential(
            "-m / sqrt(x**2 + y**2 + z**2 + b**2)",
            parameters={"m": 0 * u.M_sun, "b": radius},
        )
        compiled = PerturbedKerrMcmcMultiModel(
            astrometry, spectroscopy, reference_epoch=2000.0, force=potential,
        )
        cases["compiled potential"] = (compiled, {
            "force.m": args.plummer_mass, "force.b": args.plummer_radius,
        })
    except (ImportError, RuntimeError) as exc:
        print(f"compiled potential skipped: {exc}")
    prepared = {}
    for label, (model, force_values) in cases.items():
        model.set_solver(method=args.method, rtol=args.rtol, atol=args.atol)
        named = {k: v for k, v in values(args.spin, force_values).items()
                 if k in model.parameter_names}
        vector = model.pack(named)
        model(vector)
        prepared[label] = (model, vector)
    return prepared


def interleaved(prepared, repeat):
    samples = {label: [] for label in prepared}
    for _ in range(repeat):
        for label, (model, vector) in prepared.items():
            start = time.perf_counter()
            model(vector)
            samples[label].append(time.perf_counter() - start)
    return {label: 1e3 * statistics.median(times) for label, times in samples.items()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--repeat", type=int, default=30)
    parser.add_argument("--method", default="dop853", choices=("dop853", "radau"))
    parser.add_argument("--rtol", type=float, default=1e-9)
    parser.add_argument("--atol", type=float, default=1e-12)
    parser.add_argument("--spin", type=float, nargs=3, default=(0.3, -0.2, 0.5))
    parser.add_argument("--disk-mass", type=float, default=1e3,
                        help="thin-disk mass in solar masses")
    parser.add_argument("--plummer-mass", type=float, default=4e3,
                        help="Plummer mass in solar masses")
    parser.add_argument("--plummer-radius", type=float, default=1.0,
                        help="Plummer scale radius in mpc")
    args = parser.parse_args()

    stars = tuple(ELEMENTS)
    print(f"solver {args.method} rtol={args.rtol:g} atol={args.atol:g}; "
          f"disk {args.disk_mass:g} Msun; "
          f"Plummer {args.plummer_mass:g} Msun, b = {args.plummer_radius:g} mpc")
    totals = interleaved(models(stars, args), args.repeat)
    print("five stars (ms, median):")
    for label, value in totals.items():
        ratio = value / totals["homogeneous"]
        print(f"  {label:20s} {value:8.3f}  x{ratio:.2f}")
    print("per star (ms, median):")
    header = "  ".join(f"{label:>18s}" for label in totals)
    print(f"  {'star':6s}{header}")
    for star in stars:
        row = interleaved(models((star,), args), args.repeat)
        print(f"  {star:6s}" + "  ".join(f"{row[label]:18.3f}" for label in totals))


if __name__ == "__main__":
    main()
