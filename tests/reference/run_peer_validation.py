"""Run the peer-oracle harness and emit a machine-readable result manifest."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from pathlib import Path
from typing import Any

import numpy as np

from .peer_oracles import (
    case_by_id,
    exact_schwarzschild_circular_reference,
    integrate_with_pygro,
    kerrgeopy_bound_reference,
    kerrgeopy_circular_reference,
    make_pygro_engine,
    software_versions,
)


def _recorded_relatipy_comparison() -> dict[str, Any]:
    """Link a production run, explicitly separate from this peer-only run."""
    repository = Path(__file__).resolve().parents[2]
    report = repository / "tests/fixtures/orbit_physical_validation_results.json"
    command = "PYTHONPATH=tests .venv/bin/python -m reference.run_physical_validation"
    if not report.is_file():
        return {"status": "not_run", "reason": "Production comparison is executable in its own environment.",
                "command": command}
    payload = json.loads(report.read_text(encoding="utf-8"))
    expected = dict(payload["source_sha256"])
    expected[payload["extension"]["path"]] = payload["extension"]["sha256"]
    stale = [name for name, checksum in expected.items()
             if not (repository / name).is_file()
             or hashlib.sha256((repository / name).read_bytes()).hexdigest() != checksum]
    return {
        "status": "stale" if stale else "recorded_" + payload["status"],
        "execution": "Recorded production run; this command reruns external peers only.",
        "path": str(report.relative_to(repository)),
        "report_sha256": hashlib.sha256(report.read_bytes()).hexdigest(),
        "generated_at_utc": payload["generated_at_utc"],
        "scientific_manual_review": payload["scientific_manual_review"],
        "changed_inputs": stale,
        "command": command,
    }


def _summary(result: Any) -> dict[str, Any]:
    return {
        "pygro_steps": int(result.proper_time.size - 1),
        "endpoint_error_t_r_theta_phi": result.endpoint_error.tolist(),
        "trajectory_max_abs_error_t_r_theta_phi": np.max(
            np.abs(result.trajectory_error), axis=1
        ).tolist(),
        "max_abs_trajectory_error": result.max_abs_trajectory_error,
        "max_abs_endpoint_error": result.max_abs_endpoint_error,
        "max_abs_normalization_error": result.max_abs_normalization_error,
    }


def run() -> dict[str, Any]:
    metric, engine = make_pygro_engine()
    results: dict[str, Any] = {}

    circular_case = case_by_id("schwarzschild_circular_r10")
    circular_reference = kerrgeopy_circular_reference(circular_case)
    circular = integrate_with_pygro(
        metric,
        engine,
        circular_reference,
        spin=circular_case["a_over_M"],
        goal=circular_case["pygro_goal"],
    )
    assert circular.max_abs_endpoint_error <= circular_case["tolerances"]["coordinate_max_abs"]
    assert circular.max_abs_trajectory_error <= circular_case["tolerances"]["trajectory_coordinate_max_abs"]
    assert circular.max_abs_normalization_error <= circular_case["tolerances"]["normalization_max_abs"]
    results[circular_case["id"]] = _summary(circular)

    isco_case = case_by_id("schwarzschild_isco_r6")
    isco_reference = exact_schwarzschild_circular_reference(isco_case["radius_over_M"])
    isco = integrate_with_pygro(
        metric,
        engine,
        isco_reference,
        spin=isco_case["a_over_M"],
        goal=isco_case["pygro_goal"],
    )
    isco_radius_error = float(
        np.max(np.abs(isco.coordinates[:, 1] - isco_case["radius_over_M"]))
    )
    assert isco.max_abs_endpoint_error <= isco_case["tolerances"]["coordinate_max_abs"]
    assert isco.max_abs_trajectory_error <= isco_case["tolerances"]["trajectory_coordinate_max_abs"]
    assert isco_radius_error <= isco_case["tolerances"]["radius_max_abs"]
    assert isco.max_abs_normalization_error <= isco_case["tolerances"]["normalization_max_abs"]
    results[isco_case["id"]] = {
        **_summary(isco),
        "max_abs_radius_error": isco_radius_error,
        "kerrgeopy_applicable": False,
        "kerrgeopy_exclusion": isco_case["kerrgeopy_exclusion"],
    }

    eccentric_case = case_by_id("schwarzschild_eccentric")
    eccentric_reference = kerrgeopy_bound_reference(eccentric_case)
    eccentric = integrate_with_pygro(
        metric,
        engine,
        eccentric_reference,
        spin=eccentric_case["a_over_M"],
        goal=eccentric_case["pygro_goal"],
    )
    eccentric_delta_phi = float(
        eccentric_reference.coordinates[3, -1]
        - eccentric_reference.coordinates[3, 0]
    )
    eccentric_advance = eccentric_delta_phi - 2.0 * np.pi
    assert eccentric.max_abs_endpoint_error <= eccentric_case["tolerances"]["coordinate_max_abs"]
    assert eccentric.max_abs_trajectory_error <= eccentric_case["tolerances"]["trajectory_coordinate_max_abs"]
    assert eccentric.max_abs_normalization_error <= eccentric_case["tolerances"]["normalization_max_abs"]
    assert eccentric_advance >= eccentric_case["tolerances"]["minimum_periapsis_advance_rad"]
    results[eccentric_case["id"]] = {
        **_summary(eccentric),
        "proper_time_final_M": float(eccentric_reference.proper_time[-1]),
        "delta_phi_rad": eccentric_delta_phi,
        "periapsis_advance_rad": eccentric_advance,
    }

    kerr_case = case_by_id("kerr_stable_bound")
    kerr_reference = kerrgeopy_bound_reference(kerr_case)
    coarse_goal, fine_goal = kerr_case["convergence_goals"]
    coarse = integrate_with_pygro(
        metric, engine, kerr_reference, spin=kerr_case["a_over_M"], goal=coarse_goal
    )
    fine = integrate_with_pygro(
        metric, engine, kerr_reference, spin=kerr_case["a_over_M"], goal=fine_goal
    )
    fine_to_coarse = fine.max_abs_endpoint_error / coarse.max_abs_endpoint_error
    assert fine.max_abs_endpoint_error <= kerr_case["tolerances"]["coordinate_max_abs"]
    assert fine.max_abs_trajectory_error <= kerr_case["tolerances"]["trajectory_coordinate_max_abs"]
    assert fine.max_abs_normalization_error <= kerr_case["tolerances"]["normalization_max_abs"]
    assert fine_to_coarse <= kerr_case["tolerances"]["fine_to_coarse_error_ratio_max"]
    assert fine.max_abs_normalization_error < coarse.max_abs_normalization_error
    results[kerr_case["id"]] = {
        "proper_time_final_M": float(kerr_reference.proper_time[-1]),
        "coarse_goal": coarse_goal,
        "fine_goal": fine_goal,
        "coarse": _summary(coarse),
        "fine": _summary(fine),
        "fine_to_coarse_error_ratio": fine_to_coarse,
    }

    return {
        "schema_version": "1.0",
        "status": "pass",
        "scientific_manual_review": "pending",
        "python": sys.version,
        "platform": platform.platform(),
        "versions": software_versions(),
        "results": results,
        "relatipy_comparison": _recorded_relatipy_comparison(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = run()
    rendered = json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        print(args.output)


if __name__ == "__main__":
    main()
