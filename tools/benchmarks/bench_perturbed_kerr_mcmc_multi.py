"""Time five-star perturbed and homogeneous multi-star MCMC calls, interleaved.

Stars and elements:

- S2, S29, S38, S55: GRAVITY Collaboration (2024), A&A 692 A242, Table B.1,
  with the angles used as given there (the convention of the KerrMcmcModel
  tests).
- S301: solution A (a = 83.0 mas, e = 0.9832, i = 124.09 deg,
  Omega = 73.8 deg, omega = 293.4 deg, t_P = 2023.126, P = 8.68 yr).

Shared: M = 4.297e6 Msun, R0 = 8277 pc. Epochs are evenly spaced over
realistic spans with about 100-150 epochs per star. One call of each model
is timed:

- ``KerrMcmcMultiModel`` (homogeneous);
- ``PerturbedKerrMcmcMultiModel`` with the fiducial thin disk of Piran et
  al. (2026) (``r_out = 1362.5 au``, the S301 apocentre) of zero mass and of
  ``--disk-mass`` solar masses.

Per-star times come from single-star calls of the same models. Calls are
interleaved and the median of ``--repeat`` rounds is reported. Run from the
repository root after building the extensions::

    PYTHONPATH=src python tools/benchmarks/bench_perturbed_kerr_mcmc_multi.py
"""

from __future__ import annotations

import argparse
import statistics
import time

import numpy as np
from astropy import units as u

from relatipy import KerrMcmcMultiModel, PerturbedKerrMcmcMultiModel, ThinDisk

# star: (t_p [yr], a [arcsec], e, inc [deg], Omega [deg], omega [deg])
ELEMENTS = {
    "S2": (2018.3789, 0.12502, 0.88444, 134.67, 228.21, 66.279),
    "S29": (2021.4102, 0.39025, 0.96880, 144.24, 4.9259, 203.68),
    "S38": (2022.6843, 0.14249, 0.81807, 168.69, 122.43, 40.065),
    "S55": (2009.4738, 0.10424, 0.72980, 159.59, 319.43, 327.77),
    "S301": (2023.126, 0.0830, 0.9832, 124.09, 73.8, 293.4),
}
# star: (astrometry span, astrometry count, spectroscopy span, count)
EPOCHS = {
    "S2": ((1992.2, 2024.5), 110, (2000.0, 2024.5), 35),
    "S29": ((2002.3, 2024.5), 100, (2007.6, 2024.5), 20),
    "S38": ((2004.2, 2024.5), 120, (2008.3, 2024.5), 15),
    "S55": ((2004.5, 2024.5), 100, (2010.0, 2024.5), 10),
    "S301": ((2016.0, 2025.5), 100, (2019.0, 2025.5), 20),
}
SHARED = {
    "D": 8.277, "M": 4.297, "xS0": -7.2e-4, "yS0": -1.2e-4,
    "vxS0": 7.8e-5, "vyS0": 4.1e-5, "v_LSR": -2.0,
}


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
    disk = ThinDisk(mass=0 * u.M_sun, r_out=1362.5 * u.au)
    cases = {"homogeneous": (KerrMcmcMultiModel(
        astrometry, spectroscopy, reference_epoch=2000.0), {})}
    for label, mass in (("disk zero", 0.0), ("disk", args.disk_mass)):
        cases[label] = (PerturbedKerrMcmcMultiModel(
            astrometry, spectroscopy, reference_epoch=2000.0, force=disk,
        ), {"force.disk_mass": mass, "force.normal_theta": 0.0,
            "force.normal_phi": 0.0})
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
    parser.add_argument("--rtol", type=float, default=1e-10)
    parser.add_argument("--atol", type=float, default=1e-13)
    parser.add_argument("--spin", type=float, nargs=3, default=(0.3, -0.2, 0.5))
    parser.add_argument("--disk-mass", type=float, default=1e3,
                        help="thin-disk mass in solar masses")
    args = parser.parse_args()

    stars = tuple(ELEMENTS)
    print(f"solver {args.method} rtol={args.rtol:g} atol={args.atol:g}; "
          f"disk {args.disk_mass:g} Msun")
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
