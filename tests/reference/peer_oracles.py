"""Adapters for independent KerrGeoPy and PyGRO reference calculations.

These helpers are test infrastructure, not a RelatiPy production backend.  All
arrays use geometric units with ``G = c = M = 1`` and Boyer--Lindquist order
``(t, r, theta, phi)``.  PyGRO receives the position and four-velocity produced
by KerrGeoPy without reconstructing either from orbital elements.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from importlib import metadata
from pathlib import Path
from typing import Any

import numpy as np


FIXTURE_PATH = (
    Path(__file__).resolve().parents[1] / "fixtures" / "peer_validation_cases.json"
)


@dataclass(frozen=True)
class ReferenceTrajectory:
    """KerrGeoPy trajectory sampled at a common Mino/proper-time mesh."""

    mino_time: np.ndarray
    proper_time: np.ndarray
    coordinates: np.ndarray
    four_velocity: np.ndarray

    @property
    def initial_state(self) -> np.ndarray:
        return np.concatenate((self.coordinates[:, 0], self.four_velocity[:, 0]))


@dataclass(frozen=True)
class PeerResult:
    """A PyGRO trajectory plus errors against a supplied reference."""

    proper_time: np.ndarray
    coordinates: np.ndarray
    four_velocity: np.ndarray
    sampled_coordinates: np.ndarray
    trajectory_error: np.ndarray
    endpoint_error: np.ndarray
    max_abs_trajectory_error: float
    max_abs_endpoint_error: float
    max_abs_normalization_error: float


def load_cases() -> dict[str, Any]:
    """Load the small, reviewed-by-humans-pending validation manifest."""

    return json.loads(FIXTURE_PATH.read_text(encoding="utf-8"))


def case_by_id(case_id: str) -> dict[str, Any]:
    """Return a named validation case from the manifest."""

    for case in load_cases()["cases"]:
        if case["id"] == case_id:
            return case
    raise KeyError(case_id)


def software_versions() -> dict[str, str]:
    """Return distribution versions without relying on module attributes."""

    return {
        "kerrgeopy": metadata.version("kerrgeopy"),
        "pygro": metadata.version("pygro"),
        "numpy": metadata.version("numpy"),
        "scipy": metadata.version("scipy"),
    }


def circular_proper_period(radius_over_m: float) -> float:
    """Proper-time period of an equatorial Schwarzschild circular orbit."""

    radius = float(radius_over_m)
    return float(2.0 * np.pi * np.sqrt(radius**3 * (1.0 - 3.0 / radius)))


def exact_schwarzschild_circular_reference(radius_over_m: float) -> ReferenceTrajectory:
    """Return the exact initial/final state for one circular proper period."""

    radius = float(radius_over_m)
    tau_final = circular_proper_period(radius)
    u_t = 1.0 / np.sqrt(1.0 - 3.0 / radius)
    u_phi = np.sqrt(1.0 / (radius**3 * (1.0 - 3.0 / radius)))
    proper_time = np.array([0.0, tau_final], dtype=float)
    coordinates = np.array(
        [
            [0.0, u_t * tau_final],
            [radius, radius],
            [np.pi / 2.0, np.pi / 2.0],
            [0.0, 2.0 * np.pi],
        ],
        dtype=float,
    )
    velocity = np.array([u_t, 0.0, 0.0, u_phi], dtype=float)
    four_velocity = np.column_stack((velocity, velocity))
    return ReferenceTrajectory(
        mino_time=np.array([0.0, tau_final / radius**2], dtype=float),
        proper_time=proper_time,
        coordinates=coordinates,
        four_velocity=four_velocity,
    )


def kerrgeopy_circular_reference(case: dict[str, Any]) -> ReferenceTrajectory:
    """Evaluate KerrGeoPy for one Schwarzschild circular proper period."""

    import kerrgeopy as kg

    radius = float(case["radius_over_M"])
    tau_final = circular_proper_period(radius)
    mino_time = np.array([0.0, tau_final / radius**2], dtype=float)
    orbit = kg.StableOrbit(a=0.0, p=radius, e=0.0, x=1.0)
    coordinates = np.vstack(
        [np.asarray(component(mino_time), dtype=float) for component in orbit.trajectory()]
    )
    four_velocity = np.vstack(
        [
            np.asarray(component(mino_time), dtype=float)
            for component in orbit.four_velocity()
        ]
    )
    # For an exactly circular equatorial orbit these derivatives vanish
    # identically.  KerrGeoPy's double-precision turning-point evaluation can
    # leave O(1e-8) radial residue; retaining it would change the initial state.
    four_velocity[1:3, :] = 0.0
    return ReferenceTrajectory(
        mino_time=mino_time,
        proper_time=np.array([0.0, tau_final], dtype=float),
        coordinates=coordinates,
        four_velocity=four_velocity,
    )


def kerrgeopy_bound_reference(case: dict[str, Any]) -> ReferenceTrajectory:
    """Evaluate a stable bound Kerr orbit and map Mino time to proper time."""

    import kerrgeopy as kg
    from scipy.integrate import solve_ivp

    spin = float(case["a_over_M"])
    orbit = kg.StableOrbit(
        spin,
        float(case["p_over_M"]),
        float(case["eccentricity"]),
        float(case["x"]),
        initial_phases=tuple(float(value) for value in case["initial_phases"]),
    )
    lambda_r = 2.0 * np.pi / orbit.upsilon_r
    mino_time = np.linspace(
        0.0,
        float(case["radial_periods"]) * lambda_r,
        int(case["samples"]),
    )
    trajectory = orbit.trajectory()
    four_velocity_functions = orbit.four_velocity()
    coordinates = np.vstack(
        [np.asarray(component(mino_time), dtype=float) for component in trajectory]
    )
    four_velocity = np.vstack(
        [
            np.asarray(component(mino_time), dtype=float)
            for component in four_velocity_functions
        ]
    )
    mapping = load_cases()["mapping"]
    proper_solution = solve_ivp(
        lambda lam, _state: [
            trajectory[1](lam) ** 2
            + spin**2 * np.cos(trajectory[2](lam)) ** 2
        ],
        (float(mino_time[0]), float(mino_time[-1])),
        [0.0],
        t_eval=mino_time,
        method=mapping["tau_mapping_method"],
        rtol=float(mapping["tau_mapping_rtol"]),
        atol=float(mapping["tau_mapping_atol"]),
    )
    if not proper_solution.success:
        raise RuntimeError(f"Mino-to-proper-time mapping failed: {proper_solution.message}")
    proper_time = np.asarray(proper_solution.y[0], dtype=float)
    if not np.all(np.diff(proper_time) > 0.0):
        raise ValueError("KerrGeoPy reference did not produce monotonic proper time")
    return ReferenceTrajectory(
        mino_time=mino_time,
        proper_time=proper_time,
        coordinates=coordinates,
        four_velocity=four_velocity,
    )


def _repair_bundled_pygro_metric(payload: dict[str, Any]) -> dict[str, Any]:
    """Adapt PyGRO 1.0.3's bundled Kerr metric to its current loader.

    The distributed metric uses legacy velocity symbols (``ut`` rather than
    ``u_t``) and only contains the time-component normalization helpers.  The
    harness directly supplies all four initial velocity components, so the
    unused spatial helpers are filled with inert expressions.
    """

    repaired = json.loads(json.dumps(payload))
    velocity = re.compile(r"\b(ut|ur|utheta|uphi)\b")

    def current_symbol(match: re.Match[str]) -> str:
        return "u_" + match.group(1)[1:]

    for key in ("eq_x", "eq_u"):
        repaired[key] = [velocity.sub(current_symbol, value) for value in repaired[key]]
    for key in ("u0_null", "u0_timelike"):
        repaired[key] = velocity.sub(current_symbol, repaired[key])
    for kind in ("null", "timelike"):
        for component in (1, 2, 3):
            repaired.setdefault(f"u{component}_{kind}", "0")
    return repaired


def make_pygro_engine() -> tuple[Any, Any]:
    """Load PyGRO's bundled Kerr-BL metric and create a DOP853 engine."""

    import pygro as pg

    metric_path = Path(pg.__file__).resolve().parent / "default_metrics" / "Kerr-BL.metric"
    payload = _repair_bundled_pygro_metric(
        json.loads(metric_path.read_text(encoding="utf-8"))
    )
    metric = pg.Metric.__new__(pg.Metric)
    metric._initialized = False
    metric._initialized_metric = False
    metric._geodesic_engine_linked = False
    metric._load_metric_from_json(payload, m=1.0, a=0.0)
    pg.Metric.instances.append(metric)
    engine = pg.GeodesicEngine(metric, backend="lambdify", integrator="dp853")
    return metric, engine


def _align_azimuth(endpoint: np.ndarray, reference: np.ndarray) -> np.ndarray:
    aligned = np.asarray(endpoint, dtype=float).copy()
    aligned[3] += 2.0 * np.pi * np.rint((reference[3] - aligned[3]) / (2.0 * np.pi))
    return aligned


def _align_azimuth_series(coordinates: np.ndarray, reference: np.ndarray) -> np.ndarray:
    aligned = np.asarray(coordinates, dtype=float).copy()
    aligned[3] += 2.0 * np.pi * np.rint(
        (reference[3] - aligned[3]) / (2.0 * np.pi)
    )
    return aligned


def integrate_with_pygro(
    metric: Any,
    engine: Any,
    reference: ReferenceTrajectory,
    *,
    spin: float,
    goal: int,
) -> PeerResult:
    """Integrate the exact reference initial state to the same proper time."""

    import pygro as pg

    metric.set_constant(m=1.0, a=float(spin))
    geodesic = pg.Geodesic("time-like", engine, verbose=False)
    geodesic.initial_x = reference.coordinates[:, 0].tolist()
    geodesic.initial_u = reference.four_velocity[:, 0].tolist()
    tau_final = float(reference.proper_time[-1])
    initial_step = min(float(load_cases()["mapping"]["pygro_initial_step_M"]), tau_final / 100.0)
    engine.integrate(
        geodesic,
        tau_final,
        initial_step=initial_step,
        verbose=False,
        accuracy_goal=int(goal),
        precision_goal=int(goal),
    )
    reference_endpoint = reference.coordinates[:, -1]
    endpoint = _align_azimuth(geodesic.x[-1], reference_endpoint)
    endpoint_error = endpoint - reference_endpoint
    sampled_coordinates = np.asarray(
        geodesic.interpolator(reference.proper_time), dtype=float
    )[:, :4].T
    sampled_coordinates = _align_azimuth_series(
        sampled_coordinates, reference.coordinates
    )
    trajectory_error = sampled_coordinates - reference.coordinates
    normalization_error = max(
        abs(float(metric.norm(position, velocity)) + 1.0)
        for position, velocity in zip(geodesic.x, geodesic.u, strict=True)
    )
    return PeerResult(
        proper_time=np.asarray(geodesic.tau, dtype=float),
        coordinates=np.asarray(geodesic.x, dtype=float),
        four_velocity=np.asarray(geodesic.u, dtype=float),
        sampled_coordinates=sampled_coordinates,
        trajectory_error=trajectory_error,
        endpoint_error=endpoint_error,
        max_abs_trajectory_error=float(np.max(np.abs(trajectory_error))),
        max_abs_endpoint_error=float(np.max(np.abs(endpoint_error))),
        max_abs_normalization_error=float(normalization_error),
    )
