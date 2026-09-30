"""Marker shapes: dots, and chevrons or arrowheads rotated to a heading."""

from __future__ import annotations

import math

import numpy as np


def heading(points: np.ndarray) -> float | None:
    """Angle in degrees of the last step between two distinct 2D points."""
    points = np.asarray(points, dtype=float)
    for index in range(len(points) - 2, -1, -1):
        step = points[-1] - points[index]
        if np.hypot(*step) > 0:
            return math.degrees(math.atan2(step[1], step[0]))
    return None


def marker_path(shape: str, angle: float | None, tip: bool = True):
    """Marker outline, rotated to ``angle`` (degrees) when directional.

    ``shape`` is ``"dot"``, ``"chevron"`` (an open ``>``) or ``"triangle"``.
    With ``tip`` the arrow tip sits on the point; otherwise it is centred.
    """
    from matplotlib.markers import MarkerStyle
    from matplotlib.path import Path
    from matplotlib.transforms import Affine2D

    if shape not in ("dot", "chevron", "triangle"):
        raise ValueError("marker shape must be 'dot', 'chevron' or 'triangle'")
    if shape == "dot" or angle is None:
        base = MarkerStyle("o")
        return base.get_path().transformed(base.get_transform())
    if shape == "chevron":
        path = Path([(-0.5, 0.5), (0.5, 0.0), (-0.5, -0.5)],
                    [Path.MOVETO, Path.LINETO, Path.LINETO])
    else:
        base = MarkerStyle(">")
        path = base.get_path().transformed(base.get_transform())
    shift = -0.5 if tip else 0.0
    return path.transformed(Affine2D().translate(shift, 0.0).rotate_deg(angle))


def style_marker(collection, settings, angle: float | None, color) -> None:
    """Apply shape, size and color to a one-point scatter collection."""
    shape = settings.get("shape", "dot")
    directional = shape in ("chevron", "triangle") and angle is not None
    collection.set_paths([marker_path(shape, angle if directional else None)])
    collection.set_sizes([settings["size"]])
    if shape == "chevron" and directional:
        collection.set_facecolor("none")
        collection.set_edgecolor(color)
        collection.set_linewidth(settings.get("linewidth", 1.0))
    else:
        collection.set_facecolor(color)
        collection.set_edgecolor(color)
        collection.set_linewidth(0.0)
