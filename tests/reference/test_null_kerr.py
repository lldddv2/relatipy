"""Kerr null reference gates from INF-030, whose manual review is pending.

The frozen fixture copies D005 values into this repository. CI never reads the
knowledge base. These short integrations constrain the selected unstable
spherical trajectories and an inclined infall, not arbitrary Kerr photons.
"""

from __future__ import annotations

import numpy as np
import pytest

from relatipy import _core

from .null_reference import canonical, constants_ray, load_fixture, make_metric, solve_ray


REFERENCE = load_fixture("kerr.json")
METHODS = tuple(REFERENCE["measurement"]["methods"])


def _polar_bounds(spin: float, b: float, eta: float) -> tuple[float, float]:
    """Derive the two polar roots from INF-030:C011 (S005, Eq. 8)."""
    # Set z=cos(theta)^2 in Theta(theta)=0. The rationalized positive
    # quadratic root avoids subtracting two nearly equal positive numbers.
    coefficient = spin**2 - eta - b**2
    z = 2 * eta / (np.sqrt(coefficient**2 + 4 * spin**2 * eta) - coefficient)
    lower = float(np.arccos(np.sqrt(z)))
    return lower, float(np.pi - lower)


@pytest.mark.parametrize("case", REFERENCE["spherical_cases"], ids=lambda case: case["id"])
@pytest.mark.parametrize("method", METHODS)
def test_null_kerr_spherical_radius_and_polar_turns(case: dict, method: str) -> None:
    """INF-030:C013/C011: D005 double-root rays stay spherical and turn in theta."""
    row = case["reference_row"]
    spin = float(row["a_over_M"])
    radius = float(row["r_over_M"])
    b = float(row["lambda_over_M"])
    eta = float(row["eta_over_M2"])
    photon = constants_ray(
        make_metric(spin), radius, b, eta, theta=case["initial_theta_rad"],
        radial_sign=case["radial_sign"], polar_sign=case["polar_sign"],
    )
    solution = solve_ray(photon, case["t_final_over_M"], method=method)
    assert solution.status == 0 and solution.success
    assert solution.termination is None
    # The last row is a time-cut interpolant, even when t_eval is absent.
    accepted = canonical(solution)[:-1]
    assert len(accepted) > 2
    tolerance = case["measurements"][method]["tolerances"]
    # These orbits are unstable. Restrict to 40M and use ten times the
    # measured radial error; initial effective k^r is recorded in the fixture.
    radius_error = float(np.max(np.abs(accepted[:, 1] - radius)))
    assert radius_error <= tolerance["radius_absolute"]

    lower, upper = _polar_bounds(spin, b, eta)
    theta = accepted[:, 2]
    # Containment allows only ten times the polar-root float64 roundoff,
    # measured independently against mpmath at 80 decimal digits. The much
    # larger finite-sampling extrema gap is never used for containment.
    bound_tolerance = tolerance["theta_bounds_absolute"]
    assert np.min(theta) >= lower - bound_tolerance
    assert np.max(theta) <= upper + bound_tolerance
    extrema_gap = max(float(np.min(theta) - lower), float(upper - np.max(theta)))
    assert extrema_gap <= tolerance["theta_extrema_gap_absolute"]
    # Both limits are approached and at least two real turning points occur.
    # A stationary or monotonic polar trajectory must not satisfy this gate.
    polar_signs = np.sign(accepted[:, 6])
    polar_signs = polar_signs[polar_signs != 0]
    assert np.count_nonzero(np.diff(polar_signs)) >= 2
    assert np.min(theta) < np.pi / 2 < np.max(theta)


@pytest.mark.parametrize("method", METHODS)
def test_null_kerr_inclined_infall_conserves_invariants(method: str) -> None:
    """INF-030:C010/C011: inclined a=0.9 rays conserve E, Lz, Q and zero norm."""
    case = REFERENCE["invariant_case"]
    photon = constants_ray(
        make_metric(case["spin"]), case["r_over_M"], case["b_over_M"],
        case["eta_over_M2"], theta=case["theta_rad"],
        radial_sign=case["radial_sign"], polar_sign=case["polar_sign"],
    )
    solution = solve_ray(photon, case["t_final_over_M"], method=method)
    assert solution.status == 0 and solution.success
    assert solution.termination is None
    # Observe the documented private state only. _core computes invariants
    # on accepted native rows; interpolated momenta can dominate their errors.
    accepted = canonical(solution)[:-1]
    assert len(accepted) > 2
    values = []
    for state in accepted:
        invariants, status = _core.null_invariants(case["spin"], state)
        assert status == 0
        assert all(np.isfinite(value) for value in invariants.values())
        values.append(invariants)
    initial = values[0]
    measurement = case["measurements"][method]
    for name, expected in (("impact_parameter", case["b_over_M"]),
                           ("eta", case["eta_over_M2"])):
        tolerance_name = "b" if name == "impact_parameter" else "eta"
        tolerance = measurement["initial_constants_tolerance"][tolerance_name]
        assert abs(initial[name] - expected) <= tolerance
    # Relative drift scales each nonzero conserved quantity by its initial
    # value. Fixture bounds are ten times the measured errors, not rtol.
    for name in ("energy", "axial_angular_momentum", "carter_constant"):
        assert initial[name] != 0
        drift = max(abs(value[name] - initial[name]) / abs(initial[name]) for value in values)
        assert drift <= measurement["tolerances"][name]
    # The relative norm uses sum |g_mu_nu k^mu k^nu|, as architecture 6.13
    # requires; also bound the raw absolute null norm for this fixed case.
    norm = max(abs(value["relative_norm"]) for value in values)
    raw_norm = max(abs(value["norm"]) for value in values)
    assert norm <= measurement["tolerances"]["relative_norm"]
    assert raw_norm <= measurement["tolerances"]["absolute_norm"]
    assert np.ptp(accepted[:, 2]) > 0
    assert np.min(accepted[:, 1]) < case["r_over_M"]
