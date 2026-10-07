"""Reproduce the frozen INF-030 null measurements and their regression budgets.

From relatipy/: PYTHONPATH=tests uv run python -m reference.run_null_calibration
--case all --check --output /tmp/null-calibration.json

Output contains regenerated fixtures plus every original Kerr candidate report.
Inputs are frozen decimal source rows, not frozen measurements. No thesis brain
or external scratch directory is needed. Scientific manual review stays pending.
--check compares float64 bits, integer counts, classifications and budgets;
elapsed seconds/runtime are reported separately and never equality gates.
The command does not update fixtures: --output writes a separate JSON report.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import platform
from importlib.metadata import version
from pathlib import Path
import time

import numpy as np
from astropy import units as u
from astropy.constants import c
from scipy.integrate import quad
from mpmath import mp

from relatipy import Kerr, _core
from .null_reference import make_metric, constants_ray, solve_ray, canonical, scales, load_fixture

METHODS = ("radau", "dop853", "dp45")
SCRIPT = "tests/reference/run_null_calibration.py"
LINEARIZATION_RULE = (
    "max(10 * perturbed.max_linear_error, 100 * roundoff_growth_floor), "
    "rounded upward to two significant figures. The floor prevents a "
    "roundoff-dominated regression threshold; only dop853 changes."
)
BUDGET_RATIONALE = (
    "Ten times each measured nonzero error, rounded upward to two significant "
    "figures. Zero-error circular radius uses ten times the largest perturbed "
    "linearization discrepancy, plus the separately stated roundoff growth "
    "floor. Circular linearization uses max(10 times each method's measured "
    "linearization error, 100 times roundoff_growth_floor), rounded upward to "
    "two significant figures: radau and dp45 retain their budgets; dop853 "
    "gains margin above roundoff. These are local regression budgets, not "
    "accuracy guarantees."
)

# All nine E1 spherical candidates, copied as decimal strings from INF-030
# D005, physical CSV lines (header is line 1). Source SHA-256:
# 1d3cc56fa27e880fdcd99e0c1568f66d5caa0c10771897527f19b2fb2a4e24fd
# Claims C011/C013 and exact source evidence are retained in kerr.json.
KERR_CANDIDATES = [
    {
        "csv_line": 3,
        "a_over_M": "0.5",
        "r_over_M": "2.4953954216968726",
        "lambda_over_M": "3.0337434076525634",
        "eta_over_M2": "10.134263804565824",
        "R_residual": "1.2829270608442539e-49",
        "dRdr_residual": "1.9062695218548791e-49"
    },
    {
        "csv_line": 4,
        "a_over_M": "0.5",
        "r_over_M": "2.6434944880598845",
        "lambda_over_M": "1.92323053761101",
        "eta_over_M2": "18.165240408922266",
        "R_residual": "6.8422776578360209e-49",
        "dRdr_residual": "9.1068047523812475e-49"
    },
    {
        "csv_line": 5,
        "a_over_M": "0.5",
        "r_over_M": "2.7915935544228965",
        "lambda_over_M": "0.75487200666177653",
        "eta_over_M2": "23.823409183604733",
        "R_residual": "3.0790249460262094e-48",
        "dRdr_residual": "4.5338069607926888e-48"
    },
    {
        "csv_line": 7,
        "a_over_M": "0.5",
        "r_over_M": "3.0877916871489203",
        "lambda_over_M": "-1.7808232914115635",
        "eta_over_M2": "26.373574073638304",
        "R_residual": "1.0263416486754031e-48",
        "dRdr_residual": "1.1284779937637317e-48"
    },
    {
        "csv_line": 8,
        "a_over_M": "0.5",
        "r_over_M": "3.2358907535119322",
        "lambda_over_M": "-3.1566677939593133",
        "eta_over_M2": "22.229147140954433",
        "R_residual": "2.3947971802426073e-48",
        "dRdr_residual": "2.3922713102362925e-48"
    },
    {
        "csv_line": 22,
        "a_over_M": "0.9",
        "r_over_M": "2.1459579553432963",
        "lambda_over_M": "1.3426451811475596",
        "eta_over_M2": "15.559475929325161",
        "R_residual": "1.2829270608442539e-49",
        "dRdr_residual": "1.6910689365457761e-49"
    },
    {
        "csv_line": 24,
        "a_over_M": "0.9",
        "r_over_M": "2.7340612832632098",
        "lambda_over_M": "-0.66425532980129748",
        "eta_over_M2": "25.564168473517438",
        "R_residual": "3.4211388289180104e-49",
        "dRdr_residual": "6.3574863722583401e-49"
    },
    {
        "csv_line": 25,
        "a_over_M": "0.9",
        "r_over_M": "3.0281129472231665",
        "lambda_over_M": "-1.9287509289054148",
        "eta_over_M2": "26.981806015203976",
        "R_residual": "3.0790249460262094e-48",
        "dRdr_residual": "2.6782495085631064e-48"
    },
    {
        "csv_line": 26,
        "a_over_M": "0.9",
        "r_over_M": "3.3221646111831232",
        "lambda_over_M": "-3.3764533978979278",
        "eta_over_M2": "24.303557625160723",
        "R_residual": "6.8422776578360209e-49",
        "dRdr_residual": "2.0865830445542649e-49"
    }
]


def measure_schwarzschild() -> dict:
    """Measure radial infall, both critical-impact sides and circular release."""
    metric = make_metric(0)
    length, _ = scales(metric)
    results = {"rtol": 1e-10, "atol": 1e-12, "radial": {}, "capture": {}, "circle": {}}
    for method in METHODS:
        started = time.perf_counter()
        photon = constants_ray(metric, 50, 0, 0)
        solution = solve_ray(photon, 60.247572207803873, method)
        states = canonical(solution)
        radius = states[:-1, 1].astype(np.longdouble)
        exact = (50 - radius) + 2 * np.log(48 / (radius - 2))
        error = float(np.max(abs(states[:-1, 0].astype(np.longdouble) - exact)))
        results["radial"][method] = {
            "status": solution.status, "accepted_t_error": error,
            "endpoint_r_error": float(abs(states[-1, 1] - 2.1)),
            "rows": len(states), "seconds": time.perf_counter() - started,
            "initial": photon._initial_y.tolist(),
        }
        for sign in (-1, 1):
            started = time.perf_counter()
            critical_b = 3 * np.sqrt(3)
            photon = constants_ray(metric, 100, critical_b * (1 + sign * 1e-3), 0)
            solution = solve_ray(photon, 400, method, r_escape=110)
            states = canonical(solution)
            results["capture"][method + str(sign)] = {
                "status": solution.status,
                "reason": None if solution.termination is None else solution.termination.reason,
                "r_final": float(states[-1, 1]), "r_min": float(np.min(states[:, 1])),
                "t_final": float(states[-1, 0]), "rows": len(states),
                "seconds": time.perf_counter() - started,
            }
        for dr in (0, 1e-7):
            started = time.perf_counter()
            photon = metric.null(R=(3 + dr) * length, Theta=np.pi / 2 * u.rad,
                                 Phi=0 * u.rad, vPhi=1 * u.rad / u.s)
            solution = solve_ray(photon, 10, method)
            states = canonical(solution)
            lyapunov = 1 / (3 * np.sqrt(3))
            initial = photon._initial_y
            f = 1 - 2 / initial[1]
            # Independent Schwarzschild connection at dr/dt=0, INF-030:C011.
            acceleration = (-f / initial[1]**2 * initial[4]**2
                            + f * initial[1] * initial[7]**2)
            effective = (abs(initial[1] - 3) + abs(initial[5] / initial[4]) / lyapunov
                         + abs(acceleration) / lyapunov**2)
            linear = 3 + dr * np.cosh(lyapunov * states[:, 0])
            results["circle"][method + str(dr)] = {
                "status": solution.status,
                "max_r_error": float(np.max(abs(states[:, 1] - 3))),
                "max_linear_error": float(np.max(abs(states[:, 1] - linear))),
                "initial_effective_perturbation": float(effective),
                "initial_acceleration": float(acceleration), "initial": initial.tolist(),
                "endpoint_r": float(states[-1, 1]), "rows": len(states),
                "seconds": time.perf_counter() - started,
            }
    return results


def ray(spin, r, b, eta, theta=np.pi / 2, radial_sign=1, polar_sign=1):
    """Preserve the original Kerr measurement's public input conversions."""
    metric = Kerr(mass=1 * u.Msun, spin=spin)
    return metric, metric.null(
        R=r * metric.r_g, Theta=theta * u.rad, Phi=0 * u.rad,
        b=b * metric.r_g, eta=eta * metric.r_g**2,
        radial_sign=radial_sign, polar_sign=polar_sign,
    )


def solve(metric, photon, t, method):
    """Solve Kerr rays with E1's explicit scalar tolerances and no t_eval."""
    return photon.solve(t_span=(0 * u.s, t * (metric.r_g / c).to(u.s)),
                        method=method, rtol=1e-10, atol=1e-12)


def theta_limits(a, b, eta):
    """INF-030:C011: solve Theta(theta)=0 in z=cos(theta)^2."""
    d = a*a - eta - b*b
    z = 2 * eta / (np.sqrt(d*d + 4*a*a*eta) - d)
    return float(np.arccos(np.sqrt(z))), float(np.pi - np.arccos(np.sqrt(z)))


def theta_root_roundoff(a, b, eta):
    """Bound the float64 polar roots against the same inputs at 80 digits."""
    with mp.workdps(80):
        am, bm, em = [mp.mpf(value) for value in (a, b, eta)]
        d = am*am - em - bm*bm
        z = 2 * em / (mp.sqrt(d*d + 4*am*am*em) - d)
        exact = mp.acos(mp.sqrt(z))
        low, high = theta_limits(a, b, eta)
        return float(max(abs(mp.mpf(low) - exact), abs(mp.mpf(high) - (mp.pi - exact))))


def drift(a, states):
    """Observe native null invariants on the original retained rows."""
    values = []
    for state in states:
        value, status = _core.null_invariants(a, state)
        if status != 0:
            raise RuntimeError(f"null_invariants failed with status {status}")
        values.append(value)
    first = values[0]
    return {
        name: float(max(abs(value[name] - first[name]) / abs(first[name]) for value in values))
        for name in ("energy", "axial_angular_momentum", "carter_constant")
    } | {
        "relative_norm": float(max(abs(value["relative_norm"]) for value in values)),
        "absolute_norm": float(max(abs(value["norm"]) for value in values)),
    }


def measure_kerr() -> list[dict]:
    """Reproduce all candidate spheres (including unused ones) and infall."""
    reports = []
    for row in KERR_CANDIDATES:
        a, r, b, eta = [float(row[key]) for key in (
            "a_over_M", "r_over_M", "lambda_over_M", "eta_over_M2")]
        try:
            metric, photon = ray(a, r, b, eta)
        except ValueError as exc:
            reports.append({"kind": "rejected", "csv_line": row["csv_line"], "reason": str(exc)})
            continue
        low, high = theta_limits(a, b, eta)
        for method in METHODS:
            started = time.perf_counter()
            solution = solve(metric, photon, 40, method)
            states = canonical(solution)[:-1]
            polar_velocity = states[:, 6]
            reports.append({
                "kind": "sphere", "csv_line": row["csv_line"], "method": method,
                "status": solution.status, "runtime": time.perf_counter() - started,
                "a": a, "r": r, "b": b, "eta": eta, "initial_kr": float(states[0, 5]),
                "max_radius_error": float(np.max(abs(states[:, 1] - r))),
                "theta_limits": [low, high],
                "theta_range": [float(np.min(states[:, 2])), float(np.max(states[:, 2]))],
                # Retain initial zero in the original sign-change count.
                "turns": int(np.count_nonzero(np.diff(np.sign(polar_velocity)))),
                "theta_bound_excess": float(max(0, low - np.min(states[:, 2]),
                                                np.max(states[:, 2]) - high)),
                "theta_root_roundoff": theta_root_roundoff(a, b, eta),
                "theta_extrema_gap": float(max(np.min(states[:, 2]) - low,
                                               high - np.max(states[:, 2]))),
                "drift": drift(a, states),
            })
    metric, photon = ray(.9, 10, 3, 5, 1.1, radial_sign=-1, polar_sign=1)
    for method in METHODS:
        started = time.perf_counter()
        solution = solve(metric, photon, 20, method)
        states = canonical(solution)[:-1]
        reports.append({
            "kind": "invariants", "method": method, "status": solution.status,
            "runtime": time.perf_counter() - started,
            "r_range": [float(states[:, 1].min()), float(states[:, 1].max())],
            "drift": drift(.9, states),
            "initial_constants_error": {
                "b": abs(_core.null_invariants(.9, states[0])[0]["impact_parameter"] - 3),
                "eta": abs(_core.null_invariants(.9, states[0])[0]["eta"] - 5),
            },
        })
    return reports


def rounded_up(value: float, digits: int = 2) -> float:
    """Keep the original float64 ceiling arithmetic for frozen budgets."""
    if value == 0:
        return 0.0
    exponent = math.floor(math.log10(value)) - digits + 1
    unit = 10 ** exponent
    return math.ceil(value / unit) * unit


def budget(error: float) -> float:
    """Round ten times a measured error upward to two significant figures."""
    return rounded_up(10 * error)


def metadata(case: str) -> dict:
    """Describe this executable and the actual numerical environment."""
    return {
        "script": SCRIPT,
        "command": "PYTHONPATH=tests uv run python -m reference.run_null_calibration "
                   f"--case {case} --check",
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "python": platform.python_version(), "platform": platform.platform(),
        **{name: version(name) for name in ("numpy", "astropy", "scipy", "mpmath")},
        "extension_sha256": hashlib.sha256(Path(_core.__file__).read_bytes()).hexdigest(),
        "timing_policy": "seconds/runtime are descriptive, remeasured, and excluded from bitwise equality.",
    }


def generate_schwarzschild(frozen: dict, measured: dict) -> dict:
    """Rebuild every Schwarzschild measurement and budget from fresh results."""
    out = copy.deepcopy(frozen)
    out["measurement"].update(metadata("schwarzschild"), budget_rationale=BUDGET_RATIONALE)
    radial, capture, circle = out["radial"], out["capture"], out["circle"]
    radial["methods"] = {
        method: {**row, "accepted_t_tolerance": budget(row["accepted_t_error"]),
                 "endpoint_r_tolerance": budget(row["endpoint_r_error"])}
        for method, row in measured["radial"].items()
    }
    errors = [abs(np.longdouble(row["t_minus_t0_over_M"]) - (
        50 - np.longdouble(row["r_over_M"]) + 2 * np.log(
            48 / (np.longdouble(row["r_over_M"]) - 2))))
        for row in out["provenance"][1]["rows"]]
    radial["csv_formula_error_measured"] = float(max(errors))
    radial["csv_formula_tolerance"] = budget(float(max(errors)))
    capture["measured"] = measured["capture"]
    capture["horizon_excess_measured"] = max(
        row["r_final"] - 2 for row in measured["capture"].values() if row["reason"] == "horizon")
    capture["horizon_excess_tolerance"] = budget(capture["horizon_excess_measured"])
    capture["classification_mismatches_measured"] = sum(
        row["status"] != 1 or row["reason"] != ("horizon" if side == -1 else "escape")
        for method in METHODS for side in (-1, 1)
        for row in [measured["capture"][method + str(side)]])
    capture["classification_mismatches_tolerance"] = 0
    circle["lyapunov"] = 1 / (3 * np.sqrt(3))
    floor = float(10 * np.finfo(float).eps * 3 * np.exp(10 / (3 * np.sqrt(3))))
    circle["roundoff_growth_floor"] = floor
    circle["linearization_tolerance_rule"] = LINEARIZATION_RULE
    linear_max = max(row["max_linear_error"] for row in measured["circle"].values())
    circle["methods"] = {}
    for method in METHODS:
        exact = measured["circle"][method + "0"]
        perturbed = measured["circle"][method + "1e-07"]
        circle["methods"][method] = {
            "exact": exact, "perturbed": perturbed,
            "exact_r_tolerance": max(budget(linear_max), floor),
            "perturbed_r_tolerance": budget(perturbed["max_r_error"]),
            "linearization_tolerance": rounded_up(max(10 * perturbed["max_linear_error"], 100 * floor)),
        }
    return out


def tail(b: float, radius: float) -> float:
    """INF-030:C006 finite-radius azimuthal tail; preserve E1 arithmetic."""
    return quad(lambda w: 1 / np.sqrt(1 - w*w + 2*w**3 / b),
                0, b / radius, epsabs=1e-14, epsrel=1e-14)[0]


def generate_deflection(frozen: dict) -> tuple[dict, list[dict]]:
    """Remeasure all nine scattering rays and regenerate their budgets."""
    out = copy.deepcopy(frozen)
    out["measurement"] = {**frozen.get("measurement", frozen["calibration"]),
                          **metadata("deflection")}
    out["measurement"]["historical_calibration"] = (
        "The calibration block is the original E1 record. measurement identifies "
        "the current reproducible generator; historical paths are not dependencies.")
    out["measurement"]["budget_rationale"] = (
        "Ten times measured error rounded upward to two significant figures; "
        "retain the original coarser one-significant-figure budget for b=1000, "
        "dop853 (7e-12). No deflection budget changes.")
    reports = []
    for case in out["cases"]:
        b = float(case["source_row"]["b_over_M"])
        case["measurements"] = {}
        for method in METHODS:
            metric = Kerr(mass=1*u.Msun, spin=0)
            photon = metric.null(R=case["initial_radius_over_M"]*metric.r_g,
                                 Theta=np.pi/2*u.rad, Phi=0*u.rad, b=b*metric.r_g,
                                 eta=0*metric.r_g**2, radial_sign=-1, polar_sign=1)
            scale = (metric.r_g/c).to(u.s)
            started = time.perf_counter()
            solution = photon.solve(
                t_span=(0*u.s, case["final_time_over_M"]*scale),
                r_escape=case["escape_radius_over_M"]*metric.r_g,
                method=method, rtol=out["solver"]["rtol"], atol=np.array(out["solver"]["atol"]))
            states = canonical(solution)
            correction = tail(b, states[0, 1]) + tail(b, states[-1, 1])
            alpha = states[-1, 3] - states[0, 3] + correction - np.pi
            error = float(abs(alpha - float(case["source_row"]["alpha_exact_quad_rad"])))
            elapsed = time.perf_counter() - started
            reason = None if solution.termination is None else solution.termination.reason
            reports.append({"b": b, "method": method, "status": solution.status,
                            "reason": reason, "r_final": float(states[-1, 1]),
                            "samples": len(solution), "correction_rad": float(correction),
                            "alpha_rad": float(alpha), "error_rad": error, "seconds": elapsed})
            case["measurements"][method] = {
                "measured_error_rad": error,
                "tolerance_rad": rounded_up(10*error, digits=1 if b == 1000 and method == "dop853" else 2),
                "measured_correction_rad": float(correction),
                "measured_final_radius_over_M": float(states[-1, 1]),
                "samples": len(solution), "seconds": elapsed,
                "status": solution.status, "termination_reason": reason,
            }
    return out, reports


def generate_kerr(frozen: dict, reports: list[dict]) -> dict:
    """Regenerate the two frozen spheres and every inclined-infall budget."""
    out = copy.deepcopy(frozen)
    out["measurement"].update(metadata("kerr"))
    for key in ("script_path", "output_path", "output_sha256"):
        out["measurement"].pop(key, None)
    out["measurement"]["output_policy"] = "--output exports all 30 candidate/invariant reports and regenerated fixtures."
    for case in out["spherical_cases"]:
        line = int(case["csv_locator"].split()[3])
        case["measurements"] = {}
        for row in reports:
            if row["kind"] != "sphere" or row["csv_line"] != line:
                continue
            case["measurements"][row["method"]] = {
                "measured": {key: row[key] for key in (
                    "max_radius_error", "theta_bound_excess", "theta_root_roundoff",
                    "theta_extrema_gap", "theta_limits", "theta_range", "turns", "initial_kr")},
                "tolerances": {"radius_absolute": 10*row["max_radius_error"],
                               "theta_bounds_absolute": 10*row["theta_root_roundoff"],
                               "theta_extrema_gap_absolute": 10*row["theta_extrema_gap"]},
            }
    case = out["invariant_case"]
    case["measurements"] = {}
    for row in reports:
        if row["kind"] != "invariants":
            continue
        case["measurements"][row["method"]] = {
            "measured": row["drift"], "tolerances": {key: 10*value for key, value in row["drift"].items()},
            "initial_constants_error": row["initial_constants_error"],
            "initial_constants_tolerance": {key: 10*value for key, value in row["initial_constants_error"].items()},
            "r_range": row["r_range"],
        }
    return out


def compare(frozen: dict, generated: dict) -> dict:
    """Report maximum absolute differences by field and enforce exact bits.

    Provenance and executable metadata are not numerical measurements. Times
    are measured anew and reported separately; every other fixture field is
    checked recursively, including counts, arrays, termination and budgets.
    """
    differences, timings, mismatches = {}, {}, []

    def visit(expected, actual, path):
        field = path[-1]
        if isinstance(expected, dict):
            if not isinstance(actual, dict) or expected.keys() != actual.keys():
                mismatches.append(".".join(path) + ": keys differ")
                return
            for key in expected:
                visit(expected[key], actual[key], (*path, key))
        elif isinstance(expected, list):
            if not isinstance(actual, list) or len(expected) != len(actual):
                mismatches.append(".".join(path) + ": lengths differ")
                return
            for index, (left, right) in enumerate(zip(expected, actual)):
                visit(left, right, (*path, str(index)))
        elif isinstance(expected, (int, float)):
            group = next((key for key in reversed(path) if not key.isdigit()), field)
            delta = abs(actual - expected)
            target = timings if field in ("seconds", "runtime") else differences
            target[group] = max(target.get(group, 0), delta)
            same = (actual.hex() == expected.hex() if isinstance(expected, float)
                    and isinstance(actual, float) else actual == expected)
            if target is differences and not same:
                mismatches.append(f"{'.'.join(path)}: {expected!r} != {actual!r} (delta={delta:.17g})")
        elif expected != actual:
            mismatches.append(f"{'.'.join(path)}: {expected!r} != {actual!r}")

    for key in frozen:
        if key in ("measurement", "calibration", "provenance"):
            continue
        expected, actual = frozen[key], generated[key]
        visit(expected, actual, (key,))
    # Enforce executable provenance without treating environment or runtimes
    # as reproducible observations across machines.
    for key in ("script", "command", "script_sha256"):
        expected = frozen.get("measurement", {}).get(key)
        actual = generated["measurement"][key]
        if expected != actual:
            mismatches.append(f"measurement.{key}: {expected!r} != {actual!r}")
    return {"status": "PASS" if not mismatches else "FAIL",
            "max_abs_difference_by_field": differences,
            "nonreproducible_timing_max_abs_difference": timings,
            "mismatches": mismatches}


def main() -> None:
    """Run selected measurements; optionally export and enforce equality."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=(*("schwarzschild", "deflection", "kerr"), "all"), default="all")
    parser.add_argument("--output", type=Path, help="Separate JSON report path; no fixture is overwritten implicitly.")
    parser.add_argument("--check", action="store_true", help="Require bitwise scientific measurements and budget equality.")
    args = parser.parse_args()
    fixture_directory = Path(__file__).resolve().parents[1] / "fixtures" / "null"
    if args.output and args.output.resolve() in fixture_directory.glob("*.json"):
        parser.error("--output must be a separate report, not an input fixture")
    cases = ("schwarzschild", "deflection", "kerr") if args.case == "all" else (args.case,)
    payload = {"manual_review": "pending", "fixtures": {}, "raw_measurements": {}, "checks": {}}
    for case in cases:
        frozen = load_fixture(case + ".json")
        if case == "schwarzschild":
            raw = measure_schwarzschild()
            generated = generate_schwarzschild(frozen, raw)
        elif case == "deflection":
            generated, raw = generate_deflection(frozen)
        else:
            raw = measure_kerr()
            generated = generate_kerr(frozen, raw)
        payload["raw_measurements"][case] = raw
        payload["fixtures"][case] = generated
        payload["checks"][case] = compare(frozen, generated)
        print(case + ": " + json.dumps(payload["checks"][case], sort_keys=True))
    payload["status"] = "PASS" if all(check["status"] == "PASS" for check in payload["checks"].values()) else "FAIL"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    if args.check and payload["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
