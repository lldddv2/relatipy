"""
Three-dimensional orbit visualization helpers for Plotly.

This module builds composite Plotly traces for geodesic orbits in Cartesian
coordinates: polylines with optional start markers and end velocity cones,
reference planes, Schwarzschild-radius spheres, and Innermost Stable Circular
Orbit (ISCO) guides in the equatorial plane.

Notes
-----
Trajectory arrays follow the convention used elsewhere in relatipy: row index
``0`` is the time coordinate, rows ``1``--``3`` are spatial Cartesian
coordinates, and rows ``4``--``6`` store the corresponding velocity or
derivative components used for tangent visualization at the path endpoints.

Examples
--------
>>> import numpy as np
>>> from relatipy.visualization._3D.orbits import set_axis_equal
>>> import plotly.graph_objects as go
>>> fig = go.Figure()
>>> set_axis_equal(fig, 10.0)  # doctest: +SKIP
"""

from __future__ import annotations

from typing import Any, List

import numpy as np
import plotly.graph_objects as go

from .base_elements import (
    construct_arrow,
    construct_circular_plane_surface,
    construct_path,
    construct_point,
    construct_sphere_surface,
)


def construct_orbit_path(
    path: Any,
    color: str = "red",
    opacity: float = 1.0,
    size: float = 5,
    start_point: bool = True,
    end_point: bool = True,
    label: str = "orbit",
) -> List[Any]:
    """
    Build Plotly traces for a 3D orbit path with optional start and end markers.

    The trajectory is converted to Cartesian coordinates, drawn as a line, and
    optionally annotated with a start point and an end cone aligned with the
    local velocity at the final sample.

    Parameters
    ----------
    path
        Trajectory object with :meth:`convert_to` accepting ``"Cartesian"`` and
        supporting two-dimensional indexing. After conversion, row ``1``--``3``
        are :math:`x,y,z` along the parameter, and rows ``4``--``6`` (if present)
        give the components used for the end tangent.
    color : str, optional
        Color for the path and markers.
    opacity : float, optional
        Opacity of the path and markers (clamped behavior for start/end
        highlights: slightly higher opacity when ``opacity < 0.8``).
    size : float, optional
        Marker or cone size passed to :func:`~.base_elements.construct_point`
        and :func:`~.base_elements.construct_arrow`.
    start_point : bool, optional
        If True, add a scatter point at the first spatial sample.
    end_point : bool, optional
        If True, add a velocity cone at the last spatial sample using
        ``(path[4,-1], path[5,-1], path[6,-1])`` as direction components.
    label : str, optional
        Legend label for the path; start/end traces use derived labels.

    Returns
    -------
    list
        Plotly graph objects (e.g. :class:`plotly.graph_objects.Scatter3d`,
        :class:`plotly.graph_objects.Cone`) ready to pass to
        :class:`plotly.graph_objects.Figure`.

    Examples
    --------
    >>> import numpy as np
    >>> from relatipy.numeric.coordinates import Cartesian
    >>> from relatipy.visualization._3D.orbits import construct_orbit_path
    >>> n = 8
    >>> t = np.linspace(0.0, 1.0, n)
    >>> xs = np.stack([t, np.cos(2 * np.pi * t), np.sin(2 * np.pi * t), np.zeros(n)])
    >>> vs = np.zeros((3, n))
    >>> path = Cartesian(xs, vels=vs, from_dxs_dt=False)
    >>> len(construct_orbit_path(path)) >= 1
    True
    """
    path = path.convert_to("Cartesian")

    elements = [
        construct_path(path, color, opacity, label=label),
    ]
    if start_point:
        start_opacity = opacity + 0.2 if opacity < 0.8 else 1.0
        elements.append(
            construct_point(
                path[1, 0],
                path[2, 0],
                path[3, 0],
                color=color,
                size=size,
                opacity=start_opacity,
                label=f"{label} start point",
            ),
        )
    if end_point:
        end_opacity = opacity + 0.2 if opacity < 0.8 else 1.0
        elements.append(
            construct_arrow(
                path[1, -1],
                path[2, -1],
                path[3, -1],
                u=path[4, -1],
                v=path[5, -1],
                w=path[6, -1],
                color=color,
                size=size,
                opacity=end_opacity,
                label=f"{label} end point",
            ),
        )

    return elements


def construct_black_hole(
    R_s: float,
    color: str = "black",
    opacity: float = 1.0,
) -> go.Surface:
    """
    Build a spherical surface at the Schwarzschild radius for 3D plots.

    Parameters
    ----------
    R_s : float
        Schwarzschild radius (sphere radius in the same length units as the plot).
    color : str, optional
        Uniform surface color.
    opacity : float, optional
        Surface opacity in ``[0, 1]``.

    Returns
    -------
    plotly.graph_objects.Surface
        A colored sphere suitable as ``data`` for :class:`plotly.graph_objects.Figure`.

    Examples
    --------
    >>> from relatipy.visualization._3D.orbits import construct_black_hole
    >>> surf = construct_black_hole(2.0)
    >>> surf.name
    'Schwarzschild radius'
    """
    return construct_sphere_surface(R_s, color, opacity)


def set_axis_equal(fig: go.Figure, max_range: float) -> None:
    """
    Set equal axis limits on a 3D scene with cubic aspect ratio.

    Parameters
    ----------
    fig : plotly.graph_objects.Figure
        Figure whose ``scene`` layout is updated in place.
    max_range : float
        Half-width of each axis range; axes span ``[-max_range, max_range]``.

    Returns
    -------
    None

    Examples
    --------
    >>> import plotly.graph_objects as go
    >>> from relatipy.visualization._3D.orbits import set_axis_equal
    >>> fig = go.Figure()
    >>> set_axis_equal(fig, 5.0)
    """
    fig.update_layout(
        scene=dict(
            xaxis=dict(range=[-max_range, max_range]),
            yaxis=dict(range=[-max_range, max_range]),
            zaxis=dict(range=[-max_range, max_range]),
            aspectmode="cube",
        )
    )


def construct_basic_path_plot(
    R_s: float,
    path: Any,
    color_path: str = "red",
    color_black_hole: str = "black",
    opacity_path: float = 0.2,
    opacity_black_hole: float = 0.3,
    plot_plane: bool = True,
    plot_black_hole: bool = True,
    plot_center: bool = True,
) -> go.Figure:
    """
    Assemble a default 3D figure for a single orbit around a central mass.

    Optionally draws an equatorial reference disk, a sphere at ``R_s``, the
    origin, and the orbit from ``construct_orbit_path``. Axis limits are set
    from the spatial extent of the trajectory.

    Parameters
    ----------
    R_s : float
        Schwarzschild radius passed to :func:`construct_black_hole` when
        ``plot_black_hole`` is True.
    path
        Same convention as :func:`construct_orbit_path`; converted to Cartesian
        internally. :attr:`state_vector` (excluding the time component) is used
        to estimate the plot range.
    color_path : str, optional
        Orbit color.
    color_black_hole : str, optional
        Color for the horizon sphere and reference plane.
    opacity_path : float, optional
        Opacity for the orbit line and markers.
    opacity_black_hole : float, optional
        Opacity for the horizon; the reference plane uses one third of this.
    plot_plane : bool, optional
        If True, add a circular disk in the ``z = 0`` plane.
    plot_black_hole : bool, optional
        If True, add the Schwarzschild-radius sphere.
    plot_center : bool, optional
        If True, mark the origin ``(0, 0, 0)``.

    Returns
    -------
    plotly.graph_objects.Figure
        Figure containing the selected elements with equal-scaled axes.

    Examples
    --------
    >>> import numpy as np
    >>> from relatipy.numeric.coordinates import Cartesian
    >>> n = 20
    >>> t = np.linspace(0.0, 1.0, n)
    >>> xs = np.stack([t, np.cos(2 * np.pi * t), np.sin(2 * np.pi * t), np.zeros_like(t)])
    >>> path = Cartesian(xs, vels=np.zeros((3, n)), from_dxs_dt=False)
    >>> fig = construct_basic_path_plot(0.1, path, plot_black_hole=False)
    >>> fig.layout.scene.aspectmode
    'cube'
    """
    path = path.convert_to("Cartesian")
    max_range = np.max(np.abs(path.state_vector[1:]))

    elements = []
    if plot_plane:
        elements.append(
            construct_circular_plane_surface(
                max_range * 1.2,
                color_black_hole,
                opacity_black_hole / 3,
            )
        )
    if plot_black_hole:
        elements.append(construct_black_hole(R_s, color_black_hole, opacity_black_hole))
    if plot_center:
        elements.append(
            construct_point(
                0,
                0,
                0,
                color="black",
                size=10,
                opacity=1.0,
                label="center",
            )
        )

    elements.extend(construct_orbit_path(path, color_path, opacity_path))

    fig = go.Figure(data=elements)

    set_axis_equal(fig, max_range * 1.2)
    return fig


def construct_isco(
    r: float,
    prograde: bool,
    color: str,
    n_arrows: int = 8,
    line_width: int = 4,
    opacity: float = 0.9,
) -> List[go.Scatter3d]:
    """
    Build the ISCO circle in the equatorial plane with tangential rotation arrows.

    The innermost stable circular orbit is drawn as a dashed circle of radius
    ``r`` in the :math:`z = 0` plane. Several arrow glyphs indicate prograde
    (counterclockwise when viewed from :math:`+z`) or retrograde motion.

    Parameters
    ----------
    r : float
        ISCO radius in geometric (length) units consistent with the figure.
    prograde : bool
        If True, arrows indicate prograde (counterclockwise in the :math:`xy`
        plane as :math:`\\phi` increases). If False, retrograde (clockwise).
    color : str
        Color for the circle and arrow segments.
    n_arrows : int, optional
        Number of arrow glyphs distributed uniformly in azimuth.
    line_width : int, optional
        Line width in pixels for Plotly traces.
    opacity : float, optional
        Opacity of the dashed circle (arrows use full opacity).

    Returns
    -------
    list of plotly.graph_objects.Scatter3d
        The first element is the circle; the remaining elements are paired line
        segments forming arrow heads at each azimuth.

    Notes
    -----
    This is a visualization aid only; ``r`` must be supplied from the metric
    (e.g. Kerr or Schwarzschild ISCO) by the caller.

    Examples
    --------
    >>> from relatipy.visualization._3D.orbits import construct_isco
    >>> traces = construct_isco(6.0, True, "cyan", n_arrows=4)
    >>> len(traces) == 1 + 4
    True
    """
    label = "↺ prograde ISCO" if prograde else "↻ retrograde ISCO"
    sign = 1 if prograde else -1

    phis = np.linspace(0, 2 * np.pi, 120)
    circle = go.Scatter3d(
        x=r * np.cos(phis),
        y=r * np.sin(phis),
        z=np.zeros(120),
        mode="lines",
        line=dict(color=color, width=line_width, dash="dashdot"),
        opacity=opacity,
        name=label,
        legendgroup=label,
    )

    length = r * 0.3
    head_len = 0.2 * length
    spread = 0.12 * length

    arrow_phis = np.linspace(0, 2 * np.pi, n_arrows, endpoint=False)
    arrows = []
    for phi in arrow_phis:
        u_hat = sign * (-np.sin(phi))
        v_hat = sign * (np.cos(phi))

        tip_x = r * np.cos(phi)
        tip_y = r * np.sin(phi)

        bx = tip_x - head_len * u_hat
        by = tip_y - head_len * v_hat

        px, py = -v_hat, u_hat  # perpendicular

        arrows.append(
            go.Scatter3d(
                x=[tip_x, bx + spread * px, None, tip_x, bx - spread * px],
                y=[tip_y, by + spread * py, None, tip_y, by - spread * py],
                z=[0, 0, None, 0, 0],
                mode="lines",
                line=dict(color=color, width=line_width),
                opacity=1,
                hoverinfo="skip",
                showlegend=False,
                legendgroup=label,
            )
        )

    return [circle] + arrows
