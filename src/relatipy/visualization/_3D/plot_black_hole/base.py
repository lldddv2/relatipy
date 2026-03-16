"""
Base classes and descriptors for 3D black hole plotting.
"""
import numpy as np
from abc import ABC, abstractmethod

from ..base_elements import (
    construct_path,
    construct_point,
    construct_arrow,
    construct_circular_plane_surface,
    construct_square_plane_surface,
)
from ..orbits import set_axis_equal


class OrbitPath:
    """Encapsulates a path to be plotted (stores parameters only)."""

    def __init__(
        self,
        path,
        color='red',
        opacity=0.6,
        size=5,
        show_start=True,
        show_end=True,
        label="orbit",
    ):
        self.path = path
        self.color = color
        self.opacity = opacity
        self.size = size
        self.show_start = show_start
        self.show_end = show_end
        self.label = label


class EquatorialPlane:
    """Encapsulates the equatorial plane (circular disk)."""

    def __init__(self, color='gray', opacity=0.08):
        self.color = color
        self.opacity = opacity


class SquarePlane:
    """Encapsulates a square equatorial plane."""

    def __init__(self, semi_side, color='gray', opacity=0.08):
        self.semi_side = semi_side
        self.color = color
        self.opacity = opacity


class BasePlotBlackHole(ABC):
    """Base class for 3D black hole + orbits visualization."""

    def __init__(self, metric):
        self.metric = metric
        self._paths: list[OrbitPath] = []
        self._planes: list = []

    def add_path(
        self,
        path,
        color='red',
        opacity=0.6,
        size=5,
        show_start=True,
        show_end=True,
        label="orbit",
    ) -> None:
        """Create and append an OrbitPath."""
        self._paths.append(
            OrbitPath(
                path,
                color=color,
                opacity=opacity,
                size=size,
                show_start=show_start,
                show_end=show_end,
                label=label,
            )
        )

    def add_plane(self, plane) -> None:
        """Add an EquatorialPlane or SquarePlane instance."""
        self._planes.append(plane)

    @abstractmethod
    def _build_black_hole_elements(self) -> list:
        """Return list of plotly elements for the black hole (horizon, etc.)."""
        pass

    def plot(self, show_center=True):
        """Build and return the full figure."""
        import plotly.graph_objects as go

        # 1. max_range: max of abs of spatial coords (rows 1,2,3) over all paths; fallback R_s * 10
        max_range = self.metric.R_s * 10
        for op in self._paths:
            cart = op.path.convert_to("Cartesian")
            sv = cart.state_vector
            if sv.ndim == 2:
                spatial = sv[1:4]
            else:
                spatial = np.asarray(sv[1:4]).reshape(3, -1)
            max_range = max(max_range, float(np.max(np.abs(spatial))))

        # 2. Black hole elements
        all_elements = list(self._build_black_hole_elements())

        # 3. Per OrbitPath: convert to Cartesian, construct_path, optional start point and end arrow
        for op in self._paths:
            cart = op.path.convert_to("Cartesian")
            all_elements.append(
                construct_path(cart, color=op.color, opacity=op.opacity, label=op.label)
            )
            if op.show_start:
                all_elements.append(
                    construct_point(
                        cart[1, 0],
                        cart[2, 0],
                        cart[3, 0],
                        color=op.color,
                        size=op.size,
                        opacity=min(op.opacity + 0.2, 1.0),
                        label=f"{op.label} start",
                    )
                )
            if op.show_end:
                all_elements.append(
                    construct_arrow(
                        cart[1, -1],
                        cart[2, -1],
                        cart[3, -1],
                        u=cart[4, -1],
                        v=cart[5, -1],
                        w=cart[6, -1],
                        color=op.color,
                        size=op.size,
                        opacity=min(op.opacity + 0.2, 1.0),
                        label=f"{op.label} end",
                    )
                )

        # 4. Planes
        for plane in self._planes:
            if isinstance(plane, EquatorialPlane):
                all_elements.append(
                    construct_circular_plane_surface(
                        max_range * 1.2,
                        plane.color,
                        plane.opacity,
                    )
                )
            elif isinstance(plane, SquarePlane):
                all_elements.append(
                    construct_square_plane_surface(
                        plane.semi_side,
                        plane.color,
                        plane.opacity,
                    )
                )

        # 5. Center point
        if show_center:
            all_elements.append(
                construct_point(0, 0, 0, color='black', size=10, opacity=1.0, label="center")
            )

        # 6. Build figure and set equal axes
        fig = go.Figure(data=all_elements)
        set_axis_equal(fig, max_range * 1.2)
        return fig
