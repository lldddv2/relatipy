"""Schwarzschild scattering against frozen INF-030:D003 quadrature values.

INF-030 manual review remains pending. Geometric units are G = c = M = 1.
The finite initial and escape radii are corrected with independent analytic
tail integrals, rather than compared directly with an asymptotic deflection.
"""

import numpy as np
import pytest
from scipy.integrate import quad

from .null_reference import (
    canonical,
    constants_ray,
    load_fixture,
    make_metric,
    solve_ray,
)


REFERENCE = load_fixture("deflection.json")
METHODS = ("radau", "dop853", "dp45")


def _asymptotic_tail(b: float, radius: float) -> float:
    """Evaluate the missing C006 azimuth from a finite radius to infinity.

    With w=b/r, C006 gives integral_0^(b/r) dw/sqrt(1-w²+2w³/b).
    This independent reference formula is used only by the test oracle.
    """
    value, _ = quad(
        lambda w: 1 / np.sqrt(1 - w**2 + 2 * w**3 / b),
        0, b / radius, epsabs=1e-14, epsrel=1e-14,
    )
    return value


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("case", REFERENCE["cases"], ids=lambda case: case["id"])
def test_schwarzschild_deflection_matches_asymptotic_quadrature(case, method):
    """INF-030:C006/C008, D003: finite-radius-corrected weak-field scattering."""
    row = case["source_row"]
    b = float(row["b_over_M"])
    initial_radius = case["initial_radius_over_M"]
    escape_radius = case["escape_radius_over_M"]
    metric = make_metric(0)
    ray = constants_ray(metric, initial_radius, b, 0)
    solution = solve_ray(
        ray, case["final_time_over_M"], method=method,
        rtol=REFERENCE["solver"]["rtol"],
        atol=np.array(REFERENCE["solver"]["atol"]),
        r_escape=escape_radius,
    )
    assert solution.success
    assert solution.status == 1
    assert solution.termination.reason == "escape"

    # Escape is detected on an accepted step. Use its actual radius: the
    # event overshoots r_escape and neither endpoint is a t_eval interpolant.
    states = canonical(solution)
    assert states[0, 5] < 0 < states[-1, 5]
    assert states[-1, 1] >= escape_radius
    tails = (_asymptotic_tail(b, states[0, 1])
             + _asymptotic_tail(b, states[-1, 1]))
    alpha = states[-1, 3] - states[0, 3] + tails - np.pi
    error = abs(alpha - float(row["alpha_exact_quad_rad"]))
    # Per-case bounds round up 10× the measured absolute error; calibration
    # values and the exact numerical controls are frozen in deflection.json.
    assert error <= case["measurements"][method]["tolerance_rad"]
