"""Internal helpers shared by the Matplotlib and Plotly renderers."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from ..geodesic import Solution

PROJECTIONS = ("3d", "xy", "xz", "yz")
_COLUMNS = {"3d": (0, 1, 2), "xy": (0, 1), "xz": (0, 2), "yz": (1, 2)}


def projection_columns(projection: str) -> tuple[int, ...]:
    """Return the Cartesian columns shown by a supported projection.

    Examples
    --------
    >>> projection_columns("xz")
    (0, 2)
    """
    if projection not in _COLUMNS:
        raise ValueError("projection must be 3d, xy, xz or yz")
    return _COLUMNS[projection]


def stored_positions(solution: Solution) -> np.ndarray:
    """Return finite stored positions of shape ``(n, 3)`` in ``solution.x.unit``."""
    points = np.array(solution.xyz.to_value(solution.x.unit), dtype=float, copy=True)
    if points.ndim != 2 or points.shape[1] != 3 or len(points) == 0:
        raise ValueError("solution positions must have nonempty shape (n, 3)")
    if not np.all(np.isfinite(points)):
        raise ValueError("solution positions must be finite")
    return points


def trajectory_status(solution: Solution) -> str | None:
    """Describe an early end of integration, or return ``None`` on success."""
    if solution.status == -1:
        return "partial: numerical failure"
    if solution.status == 1:
        return f"terminated: {solution.termination.reason.replace('_', ' ')}"
    return None


def trajectory_title(solution: Solution) -> str:
    """Describe the stored trajectory and how its integration ended."""
    status = trajectory_status(solution)
    if status is None:
        return "Integrated trajectory"
    return f"Integrated trajectory ({status})"
