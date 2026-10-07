"""Small public-API builders for frozen INF-030 null references.

Fixtures are local to RelatiPy; the brain path is provenance only. Geometric
length and coordinate time are expressed as GM/c² and GM/c³, respectively.
Private canonical states are used only for numerical observations in tests.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.constants import c

from relatipy import Kerr


def load_fixture(name: str) -> dict:
    """Read a frozen reference without accessing the thesis knowledge base."""
    path = Path(__file__).resolve().parents[1] / "fixtures" / "null" / name
    return json.loads(path.read_text(encoding="utf-8"))


def make_metric(spin: float) -> Kerr:
    """Use one solar mass for presentation and normalize by its own scales."""
    return Kerr(mass=1 * u.Msun, spin=spin)


def scales(metric: Kerr):
    """Return the geometric length (metres) and coordinate time (seconds)."""
    length = metric.r_g.to(u.m)
    return length, (length / c).to(u.s)


def constants_ray(metric, r, b, eta, theta=np.pi / 2,
                  radial_sign=-1, polar_sign=1):
    """Build a future null ray from normalized Boyer–Lindquist constants."""
    length, _ = scales(metric)
    return metric.null(
        R=r * length, Theta=theta * u.rad, Phi=0 * u.rad,
        b=b * length, eta=eta * length**2,
        radial_sign=radial_sign, polar_sign=polar_sign,
    )


def solve_ray(ray, t_final, method="radau", rtol=1e-10, atol=1e-12,
              r_escape=None):
    """Solve without t_eval; the last row alone may be a time-cut interpolant."""
    length, time = scales(ray._metric)
    options = {} if r_escape is None else {"r_escape": r_escape * length}
    return ray.solve(
        t_span=(ray.t, t_final * time), method=method,
        rtol=rtol, atol=atol, **options,
    )


def canonical(solution) -> np.ndarray:
    """Observe the documented private normalized (x, k) state in tests only."""
    return solution._state._canonical
