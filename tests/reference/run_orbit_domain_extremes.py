"""Reproduce public BL domain and extreme-spin checks without a chart change.

Independent metric invariants and exact Schwarzschild radial infall supply
physical checks. Cross-method endpoints supply numerical checks only, because
all three methods use the same native RHS. Scientific manual review is pending.
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from importlib import metadata
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.constants import c

from relatipy.geodesic.exceptions import IntegrationError
from relatipy.metrics.kerr import Kerr

from .timelike_invariants import (
    kerr_metric,
    public_canonical_state,
    timelike_invariants,
)

SPINS = (0.0, 0.5, 0.99, 1.0)
METHODS = ("radau", "dop853", "dp45")
POLAR_DISTANCES = (1e-4, 1e-8)
GUARD = 64 * np.finfo(float).eps


def _solve(metric, orbit, method, end, first, maximum):
    return orbit.solve(
        tau_span=(0 * metric._time_scale, end * metric._time_scale),
        method=method, rtol=1e-9, atol=1e-12,
        first_step=first * metric._time_scale,
        max_step=maximum * metric._time_scale,
    )


def _summary(metric, solution):
    """Measure accepted samples, excluding interpolation as an error source."""
    canonical = public_canonical_state(solution, metric.mass)
    values = timelike_invariants(metric.spin, canonical)
    metric_values = kerr_metric(metric.spin, canonical[:, :4])
    drift = {
        name: float(np.max(np.abs(series - series[0])))
        for name, series in values.items()
    }
    norm_error = float(np.max(np.abs(values["norm"] + 1)))
    terms = metric_values * canonical[:, 4:, None] * canonical[:, None, 4:]
    # Coordinate matrix conditioning and cancellation are distinct. The norm
    # sum has natural target magnitude one, so this sum diagnoses cancellation.
    cancellation = float(np.max(np.sum(np.abs(terms), axis=(1, 2))))
    return {
        "status": int(solution.status),
        "steps": solution.integration.n_steps,
        "nfev": solution.integration.nfev,
        "max_norm_error": norm_error,
        "invariant_max_abs_drift": drift,
        "initial_invariants": {name: float(series[0]) for name, series in values.items()},
        "max_metric_condition_number": float(np.max(np.linalg.cond(metric_values))),
        "max_norm_absolute_term_sum": cancellation,
        "min_horizon_distance": float(np.min(canonical[:, 1])
            - (1 + np.sqrt((1 - metric.spin) * (1 + metric.spin)))),
        "min_polar_distance": float(np.min(np.minimum(canonical[:, 2], np.pi - canonical[:, 2]))),
        "endpoint_normalized": canonical[-1].tolist(),
    }


def _check_physical_summary(result, norm_tolerance=2e-9):
    assert result["status"] == 0
    assert result["max_norm_error"] < norm_tolerance
    assert result["min_horizon_distance"] > 0
    assert result["min_polar_distance"] > GUARD
    for name, drift in result["invariant_max_abs_drift"].items():
        if name == "norm":
            continue  # The fixed norm target has its own conditioning-aware bound.
        assert drift < 2e-9 * max(1, abs(result["initial_invariants"][name]))


def near_polar_case(spin, pole, distance, method):
    """Move away from either axis with nonzero radial and angular velocity."""
    metric = Kerr(mass=1 * u.Msun, spin=spin)
    direction = 1 if pole == "north" else -1
    theta = distance if pole == "north" else np.pi - distance
    orbit = metric.orbit(
        R=8 * metric.r_g, Theta=theta * u.rad, Phi=0.2 * u.rad,
        vR=-0.005 * c,
        vTheta=direction * distance * 0.1 * u.rad / metric._time_scale,
        vPhi=0.02 * u.rad / metric._time_scale,
    )
    initial = public_canonical_state(orbit.initial, metric.mass)
    solution = _solve(metric, orbit, method, 0.1, 0.002, 0.005)
    result = _summary(metric, solution)
    _check_physical_summary(result)
    actual_change = result["endpoint_normalized"][2] - initial[2]
    assert direction * actual_change > 0
    assert orbit.tau == orbit.initial.tau
    result.update(spin=spin, pole=pole, initial_distance=distance,
                  method=method, polar_angle_change=actual_change)
    return result


def near_horizon_case(spin, method):
    """Use timelike corotation just outside each outer horizon, including a=1."""
    metric = Kerr(mass=1 * u.Msun, spin=spin)
    radius = 1 + np.sqrt((1 - spin) * (1 + spin)) + 0.03
    theta = 0.7
    point = np.array([0, radius, theta, 0])
    independent_metric = kerr_metric(spin, point)
    angular_velocity = -independent_metric[0, 3] / independent_metric[3, 3]
    delta = radius**2 - 2 * radius + spin**2
    orbit = metric.orbit(
        R=radius * metric.r_g, Theta=theta * u.rad, Phi=0 * u.rad,
        vR=(-0.1 * delta / (radius**2 + spin**2)) * c,
        vTheta=1e-5 * u.rad / metric._time_scale,
        vPhi=angular_velocity * u.rad / metric._time_scale,
    )
    solution = _solve(metric, orbit, method, 0.003, 0.0001, 0.0005)
    result = _summary(metric, solution)
    _check_physical_summary(result)
    result.update(spin=spin, method=method, initial_horizon_distance=0.03)
    return result


def radial_infall_case(method):
    """Check E=1, Lz=Q=0 Schwarzschild infall with an exact proper-time radius.

    Integrating dr/dtau=-sqrt(2/r) gives
    r(tau)=[r(0)**(3/2)-3 sqrt(2) tau/2]**(2/3). This is a test-only
    derivation of the separated radial equation at a=0 and fixed rest mass;
    the interval ends outside the horizon. ut=1/(1-2/r) follows E=1.
    """
    metric = Kerr(mass=1 * u.Msun, spin=0)
    radius = 2.05
    orbit = metric.orbit(
        R=radius * metric.r_g, Theta=1.1 * u.rad, Phi=0.2 * u.rad,
        vR=(-(1 - 2 / radius) * np.sqrt(2 / radius)) * c,
    )
    # A smaller public maximum step resolves the rapidly growing BL ut well
    # enough to keep the independent norm error below 2e-9 for all methods.
    solution = _solve(metric, orbit, method, 0.02, 0.00025, 0.0005)
    canonical = public_canonical_state(solution, metric.mass)
    tau = (solution.tau / metric._time_scale).to_value(u.one)
    expected_radius = (radius**1.5 - 1.5 * np.sqrt(2) * tau)**(2 / 3)
    errors = {
        "radius": float(np.max(np.abs(canonical[:, 1] - expected_radius))),
        "radial_tangent": float(np.max(np.abs(canonical[:, 5] + np.sqrt(2 / expected_radius)))),
        "time_tangent": float(np.max(np.abs(canonical[:, 4] - 1 / (1 - 2 / expected_radius)))),
    }
    assert errors["radius"] < 2e-11
    assert errors["radial_tangent"] < 2e-10
    assert errors["time_tangent"] < 2e-8
    result = _summary(metric, solution)
    _check_physical_summary(result)
    result.update(method=method, exact_max_abs_error=errors)
    return result


def stage_failure_case(method):
    """A true polar stage crossing is numerical failure, not a horizon event."""
    metric = Kerr(mass=1 * u.Msun, spin=0)
    orbit = metric.orbit(
        R=8 * metric.r_g, Theta=0.01 * u.rad, Phi=0 * u.rad,
        vTheta=(-1 / np.sqrt(65 / 0.75)) * u.rad / metric._time_scale,
    )
    initial = public_canonical_state(orbit.initial, metric.mass)
    solution = _solve(metric, orbit, method, 0.02, 0.02, 0.02)
    assert solution.status == -1 and not solution.success
    assert solution.termination is None and solution.integration.n_steps == 0
    np.testing.assert_array_equal(public_canonical_state(solution, metric.mass)[0], initial)
    try:
        orbit.integrate(
            0.02 * metric._time_scale, method=method,
            rtol=1e-9, atol=1e-12,
            first_step=0.02 * metric._time_scale,
            max_step=0.02 * metric._time_scale,
        )
    except IntegrationError:
        pass
    else:
        raise AssertionError("stage failure must raise IntegrationError")
    np.testing.assert_array_equal(public_canonical_state(orbit, metric.mass), initial)
    return {"method": method, "status": int(solution.status),
            "n_steps": solution.integration.n_steps, "nfev": solution.integration.nfev,
            "termination": None, "retained_initial_state": True}


def horizon_failure_case(method):
    """A BL-conditioning failure must preserve a reconstructible exterior state.

    This case does not require horizon termination. The chart may fail before
    any crossing step is accepted; that outcome is numerical failure. A real
    accepted crossing, if observed, must preserve its preceding exterior state.
    """
    metric = Kerr(mass=1 * u.Msun, spin=0)
    radius = 2.05
    orbit = metric.orbit(
        R=radius * metric.r_g, Theta=1.1 * u.rad, Phi=0 * u.rad,
        vR=(-(1 - 2 / radius) * np.sqrt(2 / radius)) * c,
    )
    solution = _solve(metric, orbit, method, 0.1, 0.001, 0.001)
    assert solution.status in (-1, 1)
    assert np.all(np.isfinite(public_canonical_state(solution, metric.mass)))
    assert np.all((solution.R / metric.r_g).to_value(u.one) > 2)
    if solution.status == 1:
        assert solution.termination.reason == "outer_horizon"
        assert solution.termination.tau == solution.termination.state.tau
        assert solution.termination.state.R > metric.horizons.event
        np.testing.assert_array_equal(
            public_canonical_state(solution.termination.state, metric.mass),
            public_canonical_state(solution, metric.mass)[-1],
        )
    else:
        assert solution.termination is None
    result = _summary(metric, solution)
    result.update(method=method, final_tau=float((solution.tau[-1] / metric._time_scale).value),
                  termination=solution.termination.reason if solution.termination else None,
                  physical_accuracy="not_claimed_at_the_chart_conditioning_floor")
    return result


def rejected_initial_conditions():
    """Exercise the public chart, horizon, and spin bounds without integrating."""
    results = []
    angles = (0, GUARD, -0.01, np.pi, np.pi - GUARD, np.pi + 0.01, np.nan, np.inf)
    for spin in SPINS:
        metric = Kerr(mass=1 * u.Msun, spin=spin)
        discriminant = np.sqrt((1 - spin) * (1 + spin))
        np.testing.assert_allclose(
            (metric.horizons.event / metric.r_g).to_value(u.one),
            1 + discriminant, rtol=0, atol=2e-15,
        )
        np.testing.assert_allclose(
            (metric.horizons.cauchy / metric.r_g).to_value(u.one),
            1 - discriminant, rtol=0, atol=2e-15,
        )
        for angle in angles:
            try:
                metric.orbit(R=8 * metric.r_g, Theta=angle * u.rad, Phi=0 * u.rad)
            except ValueError:
                results.append({"kind": "polar", "spin": spin, "angle": str(angle), "rejected": True})
            else:
                raise AssertionError(f"invalid polar angle accepted: {spin}, {angle}")
        for radius in (metric.horizons.event, 0.99 * metric.horizons.event):
            try:
                metric.orbit(R=radius, Theta=0.7 * u.rad, Phi=0 * u.rad)
            except ValueError:
                results.append({"kind": "horizon", "spin": spin,
                                "radius_over_M": float((radius / metric.r_g).value), "rejected": True})
            else:
                raise AssertionError("non-exterior initial condition accepted")
        for theta in (np.nextafter(GUARD, np.inf), np.nextafter(np.pi - GUARD, -np.inf)):
            orbit = metric.orbit(R=8 * metric.r_g, Theta=theta * u.rad, Phi=0 * u.rad)
            assert np.isfinite(public_canonical_state(orbit.initial, metric.mass)).all()
    for spin in (-0.01, 1.01, np.nan, np.inf, -np.inf):
        try:
            Kerr(mass=1 * u.Msun, spin=spin)
        except ValueError:
            results.append({"kind": "spin", "spin": str(spin), "rejected": True})
        else:
            raise AssertionError(f"invalid spin accepted: {spin}")
    return results


def run():
    start = time.perf_counter()
    polar = [near_polar_case(spin, pole, distance, method)
             for spin in SPINS for pole in ("north", "south")
             for distance in POLAR_DISTANCES for method in METHODS]
    horizons = [near_horizon_case(spin, method) for spin in SPINS for method in METHODS]
    infall = [radial_infall_case(method) for method in METHODS]
    stage = [stage_failure_case(method) for method in METHODS]
    failures = []
    for method in METHODS:
        try:
            failures.append(horizon_failure_case(method))
        except ValueError as error:
            failures.append({"method": method, "status": "public_reconstruction_failure",
                             "error": str(error)})
    rejection = rejected_initial_conditions()
    differences = []
    for group, stride in ((polar, 3), (horizons, 3)):
        for index in range(0, len(group), stride):
            endpoints = np.array([item["endpoint_normalized"] for item in group[index:index + stride]])
            difference = float(np.max(np.abs(endpoints - endpoints[0])))
            assert difference < 2e-8
            differences.append(difference)
    gaps = [item for item in failures if item["status"] == "public_reconstruction_failure"]
    return {
        "schema_version": "1.0", "status": "known_failure" if gaps else "pass",
        "scientific_manual_review": "pending", "elapsed_seconds": time.perf_counter() - start,
        "python": sys.version, "platform": platform.platform(),
        "versions": {name: metadata.version(name) for name in ("numpy", "scipy", "astropy", "relatipy")},
        "units": "G=c=M=mu=1; BL (x,u) order; angles rad",
        "numerical_options": {
            "common": {"rtol": 1e-9, "atol": 1e-12},
            "near_polar": {"tau_end": 0.1, "first_step": 0.002, "max_step": 0.005},
            "near_horizon": {"tau_end": 0.003, "first_step": 0.0001, "max_step": 0.0005},
            "exact_schwarzschild_infall": {"tau_end": 0.02, "first_step": 0.00025, "max_step": 0.0005},
            "polar_stage_failure": {"tau_end": 0.02, "first_step": 0.02, "max_step": 0.02},
            "horizon_conditioning_failure": {"tau_end": 0.1, "first_step": 0.001, "max_step": 0.001},
        },
        "acceptance_limits": {
            "successful_case_norm_error": 2e-9,
            "successful_case_invariant_drift": "2e-9 * max(1, abs(initial invariant)); norm uses its fixed target",
            "exact_infall_radius_error": 2e-11,
            "exact_infall_radial_tangent_error": 2e-10,
            "exact_infall_time_tangent_error": 2e-8,
            "cross_method_endpoint_difference": 2e-8,
            "horizon_conditioning_failure": "finite exterior partial state and structured outcome; no invariant accuracy guarantee",
        },
        "source_ids": ["SCHMIDT_2002", "INF012_C001", "INF021_C001"],
        "source_locators": {
            "SCHMIDT_2002": "timelike_invariants.py; Schmidt (2002), Eqs. (3),(4),(6),(8),(17),(18); arXiv:gr-qc/0202090",
            "INF012_C001": "brain/012 record.json C001, S001 section II Eq. (2.6), p.349",
            "INF021_C001": "brain/021 record.json C001, S001 section II.1 Eq. (1) and following paragraphs",
        },
        "claim_locators": {
            "DOM001": {"kind": "derived", "source_ids": ["SCHMIDT_2002"],
                       "locator": "cases.near_polar; moving near-axis invariant conservation"},
            "DOM002": {"kind": "derived", "source_ids": ["SCHMIDT_2002", "INF012_C001"],
                       "locator": "cases.near_horizon; short exterior corotating invariant checks"},
            "DOM003": {"kind": "derived", "source_ids": ["SCHMIDT_2002"],
                       "locator": "cases.exact_schwarzschild_infall; integrated a=0 radial first integral"},
            "DOM004": {"kind": "derived", "source_ids": ["LOCAL_TEST"],
                       "locator": "cases.polar_stage_failure and cases.horizon_conditioning_failure; runtime contracts"},
            "DOM005": {"kind": "derived", "source_ids": ["LOCAL_TEST", "INF012_C001", "INF021_C001"],
                       "locator": "cases.initial_rejections; approved public bounds and chart guards"},
        },
        "local_test_source": "tests/reference/run_orbit_domain_extremes.py; output is generated numerical evidence, manual review pending",
        "accepted_guard_adjacent_initial_count": 8,
        "cases": {"near_polar": polar, "near_horizon": horizons,
                  "exact_schwarzschild_infall": infall, "polar_stage_failure": stage,
                  "horizon_conditioning_failure": failures, "initial_rejections": rejection},
        "cross_method_endpoint_max_abs_difference": max(differences),
        "reference_limits": [
            "Cross-method differences use the same production RHS and only test numerical consistency.",
            "Independent metric invariants constrain conservation, not uniqueness of a full trajectory.",
            "Large metric condition numbers near the axis do not alone measure norm cancellation.",
            "Short valid near-horizon intervals do not establish successful horizon crossing.",
            "BL chart conditioning can cause numerical failure before an accepted crossing; no alternate chart is used.",
            "Local tolerance controls do not guarantee normalization or global accuracy near a singular chart.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(run(), indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
        print(args.output)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
