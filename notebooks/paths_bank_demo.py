"""
Minimal PathsBank demo.

Run::

    python notebooks/paths_bank_demo.py

First execution integrates; second returns from cache. The named group
``demo`` bundles the trajectory.
"""

import time

import astropy.units as u
import numpy as np

from relatipy.numeric.coordinates import OrbitalElements
from relatipy.numeric.metrics import Schwarzschild
from relatipy.numeric.utils.banks import PathsBank, init_bank


def main() -> None:
    bank = init_bank(PathsBank("notebooks/paths_bank"))
    metric = Schwarzschild(1.989e30 * u.kg)
    ic = OrbitalElements(a=50.0, e=0.2, mass=1.989e30 * u.kg)
    taus = np.linspace(0.0, 200.0, 400)
    kw = dict(integrator="Yoshida6", adaptative=False, rtol=1e-9, atol=1e-12)

    t0 = time.perf_counter()
    p1 = bank.get_path(metric, ic, taus, **kw)
    t_miss = time.perf_counter() - t0

    t0 = time.perf_counter()
    p2 = bank.get_path(metric, ic, taus, **kw)
    t_hit = time.perf_counter() - t0

    print(f"miss: {t_miss*1e3:.1f} ms   hit: {t_hit*1e3:.1f} ms")
    print(f"hash: {p1._bank_hash}")
    print(f"trajectory shape: xs={p1.xs.shape}, vs={p1.vs.shape}")
    print(f"identical arrays: {np.array_equal(p1.xs, p2.xs)}")

    if "demo" not in bank.list_groups():
        bank.save_group("demo", p1, notes="paths_bank smoke demo")
    bank.show_groups()


if __name__ == "__main__":
    main()
