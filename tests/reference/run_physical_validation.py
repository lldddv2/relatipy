"""Run the selected-domain production Orbit physical validation blocks.

Independent peer references are generated separately by
``reference.export_orbit_references``. This runner needs production RelatiPy,
not the optional peer packages. It preserves numerical failure near the BL
horizon as a measured limitation, not as a successful physical trajectory.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sys
from time import perf_counter

from relatipy import _core

from . import relatipy_peer_adapter, run_orbit_domain_extremes, run_orbit_drift


REPOSITORY = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = REPOSITORY / "tests/fixtures/orbit_physical_validation_results.json"


def run() -> dict:
    """Return live component reports with source and extension checksums.

    Times are descriptive single executions, including public reconstruction.
    Passing means that the declared fixed-case regression gates pass; it does
    not establish accuracy for every exterior BL state or long integration.
    """
    start = perf_counter()
    reports = {
        "independent_references": relatipy_peer_adapter.run(),
        "invariant_drift": run_orbit_drift.run(),
        "domain_extremes": run_orbit_domain_extremes.run(),
    }
    sources = sorted({
        *REPOSITORY.glob("native/src/**/*.c"),
        *REPOSITORY.glob("native/src/**/*.h"),
        *REPOSITORY.glob("native/include/**/*.h"),
        *REPOSITORY.glob("src/relatipy/**/*.py"),
        *REPOSITORY.glob("bindings/cython/*.pyx"),
        *REPOSITORY.glob("tests/reference/*.py"),
        REPOSITORY / "tests/fixtures/orbit_peer_reference.json",
        REPOSITORY / "tests/fixtures/peer_validation_cases.json",
        REPOSITORY / "setup.py",
    })
    return {
        "schema_version": "1.0",
        "status": "pass" if all(r["status"] == "pass" for r in reports.values()) else "fail",
        "scientific_manual_review": "pending",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "elapsed_seconds": perf_counter() - start,
        "extension": {
            "path": str(Path(_core.__file__).relative_to(REPOSITORY)),
            "sha256": hashlib.sha256(Path(_core.__file__).read_bytes()).hexdigest(),
        },
        "source_sha256": {
            str(path.relative_to(REPOSITORY)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sources
        },
        "reports": reports,
        "limitations": [
            "Passing regression gates cover only the cases and tolerances in the component reports.",
            "Accepted-step invariants use integrated x/u; interpolated output has separate error budgets.",
            "Near-horizon numerical failure may contain finite exterior states with large invariant drift.",
            "A horizon event is reported only after an accepted crossing, never inferred from proximity.",
            "Single-run elapsed times are not a reproducible performance benchmark or speed guarantee.",
            "Scientific sources and computed results remain pending human review.",
        ],
    }


def main() -> None:
    """Write strict JSON and exit unsuccessfully if a component gate fails."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"{payload['status']}: {args.output}")
    if payload["status"] != "pass":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
