"""
Base classes for 3D black-hole and orbit visualization with Plotly.

This module defines lightweight configuration objects (:class:`OrbitPath`,
:class:`EquatorialPlane`, :class:`SquarePlane`) and the abstract
:class:`BasePlotBlackHole`, which collects orbit paths and optional reference
planes, then builds a :mod:`plotly.graph_objects` figure with the event
horizon (implemented in subclasses), trajectories, markers, and equal scene
axes.

Notes
-----
Spatial extent of the scene is inferred from Cartesian samples of all added
paths; if no paths are present, a default range proportional to
:attr:`~relatipy.numeric.metrics.base.BaseMetric.R_s` is used.

Examples
--------
>>> from relatipy.visualization._3D.plot_black_hole import OrbitPath
>>> cfg = OrbitPath(None, color="red", label="orbit")
>>> cfg.label
'orbit'
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
    """
    Container for styling and display options of a single orbit path.

    Instances store only parameters; geometry is resolved when
    :meth:`BasePlotBlackHole.plot` converts each path to Cartesian coordinates
    and calls the helpers in :mod:`~relatipy.visualization._3D.base_elements`.

    Parameters
    ----------
    path : object
        Trajectory or coordinate path with a ``convert_to("Cartesian")`` API
        compatible with the plotting stack (typically a geodesic result object).
    color : str, optional
        Line and marker color for this orbit. Default is ``"red"``.
    opacity : float, optional
        Opacity for the path trace, in ``[0, 1]``. Default is ``0.6``.
    size : float, optional
        Marker size for start/end decorations. Default is ``5``.
    show_start : bool, optional
        If True, draw a point at the initial spatial position. Default is True.
    show_end : bool, optional
        If True, draw a velocity arrow at the final spatial position. Default is True.
    label : str, optional
        Legend label for the path and derived markers. Default is ``"orbit"``.

    Attributes
    ----------
    path : object
        The trajectory or path object supplied at construction.
    color : str
        Line and marker color.
    opacity : float
        Opacity for the path trace.
    size : float
        Marker size for start/end decorations.
    show_start : bool
        Whether to show the start point.
    show_end : bool
        Whether to show the end arrow.
    label : str
        Legend label.

    Examples
    --------
    >>> from relatipy.visualization._3D.plot_black_hole import OrbitPath
    >>> op = OrbitPath(None, color="blue", opacity=0.5, label="test")
    >>> op.color
    'blue'
    """

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
    """
    Configuration for a circular disk in the equatorial (x–y) plane.

    The disk radius is chosen at plot time from the spatial extent of the
    scene (see :meth:`BasePlotBlackHole.plot`).

    Parameters
    ----------
    color : str, optional
        Surface color. Default is ``"gray"``.
    opacity : float, optional
        Surface opacity. Default is ``0.08``.

    Attributes
    ----------
    color : str
        Surface color.
    opacity : float
        Surface opacity.

    Examples
    --------
    >>> from relatipy.visualization._3D.plot_black_hole import EquatorialPlane
    >>> plane = EquatorialPlane(color="gray", opacity=0.1)
    >>> plane.opacity
    0.1
    """

    def __init__(self, color='gray', opacity=0.08):
        self.color = color
        self.opacity = opacity


class SquarePlane:
    """
    Configuration for a square reference plane in the equatorial plane.

    The half-edge length ``semi_side`` sets the extent in the x and y
    directions; the plane is centered on the origin.

    Parameters
    ----------
    semi_side : float
        Half side length of the square in the same length units as the plot.
    color : str, optional
        Surface color. Default is ``"gray"``.
    opacity : float, optional
        Surface opacity. Default is ``0.08``.

    Attributes
    ----------
    semi_side : float
        Half side length of the square.
    color : str
        Surface color.
    opacity : float
        Surface opacity.

    Examples
    --------
    >>> from relatipy.visualization._3D.plot_black_hole import SquarePlane
    >>> sq = SquarePlane(10.0, color="gray")
    >>> sq.semi_side
    10.0
    """

    def __init__(self, semi_side, color='gray', opacity=0.08):
        self.semi_side = semi_side
        self.color = color
        self.opacity = opacity


class BasePlotBlackHole(ABC):
    """
    Abstract base class for a 3D Plotly figure with a black hole and orbits.

    Subclasses implement :meth:`_build_black_hole_elements` to supply the
    horizon (and related) geometry. Call :meth:`add_path` and :meth:`add_plane`
    to register content, then :meth:`plot` to assemble the figure.

    Parameters
    ----------
    metric : BaseMetric
        Spacetime metric instance (e.g. :class:`~relatipy.numeric.metrics.schwarzschild_metric.Schwarzschild`).
        Its ``R_s`` attribute scales the default scene when no paths are added.

    Attributes
    ----------
    metric : BaseMetric
        The metric backing this visualization.
    _paths : list of OrbitPath
        Registered orbit paths and their display options.
    _planes : list of EquatorialPlane or SquarePlane
        Registered reference planes.

    Examples
    --------
    >>> from relatipy.visualization._3D.plot_black_hole import PlotSchwarzschild
    >>> from relatipy.numeric.metrics import Schwarzschild
    >>> plotter = PlotSchwarzschild(Schwarzschild(1.989e30))
    >>> fig = plotter.plot()
    >>> fig is not None
    True
    """

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
        """
        Append an :class:`OrbitPath` built from the given trajectory.

        Parameters
        ----------
        path : object
            Trajectory with ``convert_to("Cartesian")`` for the plotting stack.
        color : str, optional
            Line and marker color. Default is ``"red"``.
        opacity : float, optional
            Path opacity. Default is ``0.6``.
        size : float, optional
            Marker size for start/end decorations. Default is ``5``.
        show_start : bool, optional
            Whether to draw the start point. Default is True.
        show_end : bool, optional
            Whether to draw the end arrow. Default is True.
        label : str, optional
            Legend label. Default is ``"orbit"``.

        Returns
        -------
        None

        Examples
        --------
        >>> from relatipy.visualization._3D.plot_black_hole import PlotSchwarzschild
        >>> from relatipy.numeric.metrics import Schwarzschild
        >>> p = PlotSchwarzschild(Schwarzschild(1.989e30))
        >>> p.add_path(None, label="orbit")
        >>> len(p._paths)
        1
        """
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
        """
        Register a reference plane (circular or square equatorial).

        Parameters
        ----------
        plane : EquatorialPlane or SquarePlane
            Plane configuration appended to :attr:`_planes` and drawn in
            :meth:`plot`.

        Returns
        -------
        None

        Examples
        --------
        >>> from relatipy.visualization._3D.plot_black_hole import (
        ...     PlotSchwarzschild,
        ...     SquarePlane,
        ... )
        >>> from relatipy.numeric.metrics import Schwarzschild
        >>> p = PlotSchwarzschild(Schwarzschild(1.989e30))
        >>> p.add_plane(SquarePlane(5.0))
        """
        self._planes.append(plane)

    @abstractmethod
    def _build_black_hole_elements(self) -> list:
        """
        Build Plotly traces for the black hole (e.g. event horizon).

        Returns
        -------
        list
            List of :mod:`plotly.graph_objects` traces (e.g. mesh or surface)
            rendered in addition to paths and planes.

        Examples
        --------
        >>> from relatipy.visualization._3D.plot_black_hole import PlotSchwarzschild
        >>> from relatipy.numeric.metrics import Schwarzschild
        >>> plotter = PlotSchwarzschild(Schwarzschild(1.989e30))
        >>> isinstance(plotter._build_black_hole_elements(), list)
        True
        """
        pass

    def plot(self, show_center=True):
        """
        Assemble a 3D figure with horizon, paths, planes, and optional origin.

        Cartesian samples from each path determine the scene scale; the
        half-range is at least ``max(|spatial coords|)`` over all paths, or
        ``10 * R_s`` if there are no paths. Axes are set to equal scale via
        :func:`~relatipy.visualization._3D.orbits.set_axis_equal`.

        Parameters
        ----------
        show_center : bool, optional
            If True, draw a black marker at the coordinate origin. Default is True.

        Returns
        -------
        plotly.graph_objects.Figure
            Plotly figure containing all constructed traces.

        Examples
        --------
        >>> from relatipy.visualization._3D.plot_black_hole import PlotSchwarzschild
        >>> from relatipy.numeric.metrics import Schwarzschild
        >>> plotter = PlotSchwarzschild(Schwarzschild(1.989e30))
        >>> fig = plotter.plot(show_center=True)
        >>> len(fig.data) >= 1
        True
        """
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
