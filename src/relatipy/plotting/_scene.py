"""Backend-independent description of an orbit figure.

A scene holds only what is drawn: Cartesian positions in one length unit,
the role of each element, and the metric quantities needed to annotate it.
No physics is computed here; reference radii come from the native preview.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np
from astropy import units as u

from ._common import stored_positions, trajectory_status, trajectory_title

if TYPE_CHECKING:
    from ..geodesic import Orbit, Solution
    from ..metrics import Kerr

AXES = "xyz"
PLANES = ("xy", "xz", "yz")

# Orbital sense of the equatorial ISCO about the spin axis (+z): prograde
# orbits turn counter-clockwise seen from +z, retrograde ones clockwise.
ISCO_SENSE = {"isco_prograde": 1.0, "isco_retrograde": -1.0}
_REFERENCE_ROLES = ("horizon", "isco_prograde", "isco_retrograde")
_CIRCLE_SAMPLES = 240


@dataclass(frozen=True)
class Track:
    """One drawn path: its positions, markers and look.

    ``path`` has shape ``(n, 3)``; ``points`` maps marker roles to one
    position each. ``label`` and ``color`` override the style of ``role``
    when several tracks share a figure.
    """

    role: str                          # "trajectory" or "preview"
    path: np.ndarray
    points: dict[str, np.ndarray]
    tau: np.ndarray | None = None
    status: str | None = None          # early end of integration, always shown
    label: str | None = None
    color: str | None = None

    def name(self, style) -> str:
        """Legend name: the track label, else the style name of its role."""
        return self.label if self.label is not None else style.name(self.role)

    def line(self, style) -> dict:
        """Line settings of the role, with the track colour when it has one."""
        look = dict(style.lines[self.role])
        if self.color is not None:
            look["color"] = self.color
        return look


@dataclass(frozen=True)
class Scene:
    """What an orbit figure shows, in the length unit ``unit``.

    ``tracks`` are the drawn paths; ``references`` maps circle roles to
    equatorial radii.
    """

    kind: str                          # "trajectory" or "preview"
    tracks: tuple[Track, ...]
    references: dict[str, float]
    unit: u.UnitBase
    title: str
    metric: Kerr | None = None
    tau_unit: u.UnitBase | None = None
    extent: np.ndarray = field(init=False)

    def __post_init__(self) -> None:
        """Store the Cartesian extent of paths, markers, references, and origin."""
        # Everything drawn plus the black hole at the origin.
        parts = [np.zeros((1, 3))]
        for track in self.tracks:
            parts += [track.path, *[p[None, :] for p in track.points.values()]]
        for radius in self.references.values():
            parts.append(np.array([[radius, radius, 0.0], [-radius, -radius, 0.0]]))
        object.__setattr__(self, "extent", np.vstack(parts))

    @property
    def status(self) -> str | None:
        """Early ends of integration, one line per track that ended early."""
        if len(self.tracks) == 1:
            return self.tracks[0].status
        lines = [f"{track.label}: {track.status}" for track in self.tracks
                 if track.status is not None]
        return "\n".join(lines) or None

    @property
    def r_g(self) -> u.Quantity | None:
        """Return the metric gravitational radius when metric context exists."""
        return None if self.metric is None else self.metric.r_g

    @property
    def spin(self) -> float | None:
        """Return the metric dimensionless spin when metric context exists."""
        return None if self.metric is None else self.metric.spin

    def length_scale(self, length_unit: str | None) -> tuple[float, str]:
        """Data units per displayed unit, and the displayed unit's name."""
        if length_unit == "r_g" and self.r_g is not None:
            return float(self.r_g.to_value(self.unit)), "r_g"
        return 1.0, self.unit.to_string()

    def circle(self, role: str) -> np.ndarray:
        """Equatorial reference circle of ``role`` as ``(m, 3)`` points."""
        angle = np.linspace(0.0, 2.0 * np.pi, _CIRCLE_SAMPLES)
        radius = self.references[role]
        return np.column_stack((radius * np.cos(angle), radius * np.sin(angle),
                                np.zeros_like(angle)))


def _solution_track(solution: Solution, unit: u.UnitBase, tau_unit: u.UnitBase,
                    label: str | None = None, color: str | None = None) -> Track:
    """Stored samples of a solution in ``unit``, with first and last positions."""
    points = stored_positions(solution) * float(solution.x.unit.to(unit))
    tau = np.array(solution.tau.to_value(tau_unit), dtype=float, copy=True)
    if tau.shape != (len(points),) or not np.all(np.isfinite(tau)):
        raise ValueError("solution proper times must be finite and match positions")
    return Track(role="trajectory", path=points,
                 points={"initial": points[0], "end": points[-1]}, tau=tau,
                 status=trajectory_status(solution), label=label, color=color)


def solution_scene(solution: Solution) -> Scene:
    """Stored samples of a solution, with its first and last positions."""
    unit, tau_unit = solution.x.unit, solution.tau.unit
    return Scene(kind="trajectory", tracks=(_solution_track(solution, unit, tau_unit),),
                 references={}, unit=unit, title=trajectory_title(solution),
                 metric=getattr(solution, "_metric", None), tau_unit=tau_unit)


def solutions_scene(solutions, labels, style) -> Scene:
    """Stored samples of several solutions about one black hole.

    Lengths and proper times use the units of the first solution. Each
    track takes its label from ``labels`` (default: the trajectory name
    numbered from 1) and its colour from ``style.trajectory_colors``.
    """
    from ..geodesic import Solution

    if not solutions:
        raise ValueError("plot_sols needs at least one solution")
    for solution in solutions:
        if not isinstance(solution, Solution):
            raise TypeError("plot_sols accepts only Solution objects")
    if labels is not None:
        labels = list(labels)
        if len(labels) != len(solutions) or not all(isinstance(x, str) for x in labels):
            raise ValueError("labels must give one string per solution")
    # A solution built without its metric cannot contradict the others.
    known = [m for m in (getattr(s, "_metric", None) for s in solutions) if m is not None]
    if any(metric != known[0] for metric in known[1:]):
        raise ValueError("solutions must orbit the same black hole")
    metric = known[0] if known else None
    first = solutions[0]
    unit, tau_unit = first.x.unit, first.tau.unit
    if len(solutions) == 1:
        # One solution looks as with Solution.plot, unless it is labelled.
        label = None if labels is None else labels[0]
        tracks = (_solution_track(first, unit, tau_unit, label),)
        return Scene(kind="trajectory", tracks=tracks, references={}, unit=unit,
                     title=trajectory_title(first), metric=metric, tau_unit=tau_unit)
    if labels is None:
        base = style.name("trajectory")
        labels = [f"{base} {index}" for index in range(1, len(solutions) + 1)]
    colors = style.trajectory_colors
    tracks = tuple(_solution_track(solution, unit, tau_unit, label,
                                   colors[index % len(colors)])
                   for index, (solution, label) in enumerate(zip(solutions, labels)))
    return Scene(kind="trajectory", tracks=tracks, references={}, unit=unit,
                 title="Integrated trajectories", metric=metric, tau_unit=tau_unit)


def preview_scene(orbit: Orbit, show_horizon: bool = True,
                  show_isco: bool = True) -> Scene:
    """Native osculating conic, current position and equatorial references."""
    if not isinstance(show_horizon, bool) or not isinstance(show_isco, bool):
        raise TypeError("show_horizon and show_isco must be booleans")
    path, reference_radii = orbit._osculating_preview()
    unit = path.unit
    radii = dict(zip(_REFERENCE_ROLES, reference_radii.to_value(unit)))
    shown = {"horizon": show_horizon, "isco_prograde": show_isco,
             "isco_retrograde": show_isco}
    track = Track(role="preview", path=np.array(path.value, dtype=float),
                  points={"current": np.array(orbit.xyz.to_value(unit), dtype=float)})
    return Scene(kind="preview", tracks=(track,),
                 references={role: float(radius) for role, radius in radii.items()
                             if shown[role]},
                 unit=unit, title="Osculating Kepler preview", metric=orbit._metric)


def widest_plane(extent: np.ndarray) -> tuple[str, str, str]:
    """(horizontal, vertical, third) axes of the widest projection.

    The main plane is the pair of coordinates whose bounding box has the
    largest area.
    """
    spans = dict(zip(AXES, np.ptp(extent, axis=0)))
    horizontal, vertical = max(PLANES, key=lambda pair: spans[pair[0]] * spans[pair[1]])
    third = next(axis for axis in AXES if axis not in (horizontal, vertical))
    return horizontal, vertical, third
