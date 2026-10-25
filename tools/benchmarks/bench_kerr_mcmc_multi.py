"""Time one five-star KerrMcmcMultiModel call and split native from Python time.

Stars and elements:

- S2, S29, S38, S55: GRAVITY Collaboration (2024), A&A 692 A242, Table B.1,
  with the angles used as given there (the convention of the KerrMcmcModel
  tests).
- S301: solution A (a = 83.0 mas, e = 0.9832, i = 124.09 deg,
  Omega = 73.8 deg, omega = 293.4 deg, t_P = 2023.126, P = 8.68 yr).

Shared: M = 4.297e6 Msun, R0 = 8277 pc. Epochs are evenly spaced over
realistic spans with about 100-150 epochs per star.

Run from the repository root after building the extensions::

    PYTHONPATH=src python tools/benchmarks/bench_kerr_mcmc_multi.py

The native time is measured by wrapping ``relatipy._mcmc_core.evaluate``.
The Python overhead is measured directly by replacing it with a stub that
returns precomputed arrays, and indirectly as total minus native time.
"""

from __future__ import annotations

import argparse
import statistics
import time

import numpy as np

from relatipy import KerrMcmcGaussianLikelihood, KerrMcmcMultiModel
from relatipy import _mcmc_core


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


def build(spin, method, rtol, atol):
    """Return the model, the truth vector and a likelihood on noise-free data."""
    astrometry = {s: np.linspace(*span, n) for s, (span, n, _, _) in EPOCHS.items()}
    spectroscopy = {s: np.linspace(*span, n) for s, (_, _, span, n) in EPOCHS.items()}
    model = KerrMcmcMultiModel(astrometry, spectroscopy, reference_epoch=2000.0)
    model.set_solver(method=method, rtol=rtol, atol=atol)
    values = dict(SHARED)
    values |= dict(zip(("spin_dec", "spin_ra", "spin_away"), spin))
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


def timed(function, argument, repeat):
    """Return per-call wall times in seconds."""
    function(argument)
    samples = []
    for _ in range(repeat):
        start = time.perf_counter()
        function(argument)
        samples.append(time.perf_counter() - start)
    return samples


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--repeat", type=int, default=50)
    parser.add_argument("--method", default="dop853", choices=("dop853", "radau"))
    parser.add_argument("--rtol", type=float, default=1e-9)
    parser.add_argument("--atol", type=float, default=1e-12)
    parser.add_argument("--spin", type=float, nargs=3, default=(0.3, -0.2, 0.5))
    args = parser.parse_args()

    model, truth, likelihood = build(args.spin, args.method, args.rtol, args.atol)
    counts = {
        star: (model.astrometry_times[star].size, model.spectroscopy_times[star].size)
        for star in model.stars
    }
    print(f"stars and (astrometry, spectroscopy) epochs: {counts}")
    print(f"parameters: {model.n_parameters}; solver {args.method} "
          f"rtol={args.rtol:g} atol={args.atol:g}; spin={tuple(args.spin)}")

    original = _mcmc_core.evaluate
    native = []
    captured = []

    def measured(*arguments):
        start = time.perf_counter()
        result = original(*arguments)
        native.append(time.perf_counter() - start)
        captured.append(result)
        return result

    try:
        _mcmc_core.evaluate = measured
        model(truth)
        per_star = dict(zip(model.stars, native))
        native.clear()
        totals = timed(model, truth, args.repeat)
        native_per_call = [
            sum(native[i * len(model.stars):(i + 1) * len(model.stars)])
            for i in range(args.repeat + 1)
        ][1:]
        outputs = captured[: len(model.stars)]
        cursor = iter(range(10**9))

        def stub(*arguments):
            return outputs[next(cursor) % len(outputs)].copy()

        _mcmc_core.evaluate = stub
        python_model = timed(model, truth, args.repeat * 20)
        python_likelihood = timed(likelihood, truth, args.repeat * 20)
    finally:
        _mcmc_core.evaluate = original
    likelihoods = timed(likelihood, truth, args.repeat)

    def ms(samples):
        return 1e3 * statistics.median(samples)

    total = ms(totals)
    native_total = ms(native_per_call)
    print(f"model call: {total:.3f} ms median ({1e3 * min(totals):.3f} ms min)")
    print(f"  native (wrapped evaluate): {native_total:.3f} ms "
          f"= {100 * native_total / total:.1f}%")
    print(f"  python (total - native):   {total - native_total:.3f} ms")
    print(f"  python (stubbed native):   {1e3 * ms(python_model):.1f} us")
    print(f"likelihood call: {ms(likelihoods):.3f} ms median; "
          f"python with stubbed native: {1e3 * ms(python_likelihood):.1f} us")
    print("native time per star at the truth (first call, ms): " + ", ".join(
        f"{star} {1e3 * seconds:.2f}" for star, seconds in per_star.items()
    ))


if __name__ == "__main__":
    main()
