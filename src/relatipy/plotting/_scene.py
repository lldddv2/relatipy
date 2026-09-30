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
class Scene:
    """What an orbit figure shows, in the length unit ``unit``.

    ``path`` has shape ``(n, 3)``; ``points`` maps marker roles to one
    position each; ``references`` maps circle roles to equatorial radii.
    """

    kind: str                          # "trajectory" or "preview"
    path: np.ndarray
    points: dict[str, np.ndarray]
    references: dict[str, float]
    unit: u.UnitBase
    title: str
    metric: Kerr | None = None
    tau: np.ndarray | None = None
    tau_unit: u.UnitBase | None = None
    status: str | None = None          # early end of integration, always shown
    extent: np.ndarray = field(init=False)

    def __post_init__(self) -> None:
        """Store the Cartesian extent of paths, markers, references, and origin."""
        # Everything drawn plus the black hole at the origin.
        parts = [self.path, np.zeros((1, 3)), *[p[None, :] for p in self.points.values()]]
        for radius in self.references.values():
            parts.append(np.array([[radius, radius, 0.0], [-radius, -radius, 0.0]]))
        object.__setattr__(self, "extent", np.vstack(parts))

    @property
    def path_role(self) -> str:
        """Return the style role of the scene path."""
        return self.kind

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


def solution_scene(solution: Solution) -> Scene:
    """Stored samples of a solution, with its first and last positions."""
    points = stored_positions(solution)
    tau = np.array(solution.tau.to_value(solution.tau.unit), dtype=float, copy=True)
    if tau.shape != (len(points),) or not np.all(np.isfinite(tau)):
        raise ValueError("solution proper times must be finite and match positions")
    return Scene(kind="trajectory", path=points,
                 points={"initial": points[0], "end": points[-1]}, references={},
                 unit=solution.x.unit, title=trajectory_title(solution),
                 metric=getattr(solution, "_metric", None), tau=tau,
                 tau_unit=solution.tau.unit, status=trajectory_status(solution))


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
    return Scene(kind="preview", path=np.array(path.value, dtype=float),
                 points={"current": np.array(orbit.xyz.to_value(unit), dtype=float)},
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
