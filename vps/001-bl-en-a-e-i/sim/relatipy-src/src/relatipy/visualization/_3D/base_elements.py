"""
3D Plotly geometry helpers for Relatipy visualizations.

This module provides factory functions that build ``plotly.graph_objects`` traces
(surfaces, paths, markers, and arrows) for composing interactive 3D figures.

Notes
-----
Surfaces are typically uniform-color ``go.Surface`` traces with ``showscale=False``.
Arrow helpers include a ``go.Cone``-based implementation and a ``Scatter3d``-based
alternative whose line width scales with the view.

Examples
--------
>>> import numpy as np
>>> from relatipy.visualization._3D.base_elements import construct_sphere_surface
>>> trace = construct_sphere_surface(1.0, color="gray")
>>> trace.type == "surface"
True
"""

import plotly.graph_objects as go
import numpy as np


def construct_sphere_surface(R_s, color="black", opacity=1.0):
    """
    Build a spherical surface centered at the origin.

    The sphere uses a standard spherical-coordinate parameterization on a
    ``50 x 50`` grid.

    Parameters
    ----------
    R_s : float
        Radius of the sphere (Schwarzschild radius in typical black-hole plots).
    color : str, optional
        Solid color for the surface (Plotly color string).
    opacity : float, optional
        Opacity in ``[0, 1]``.

    Returns
    -------
    plotly.graph_objs.Surface
        A ``go.Surface`` trace named ``"Schwarzschild radius"``.

    Examples
    --------
    >>> s = construct_sphere_surface(2.0, color="black", opacity=1.0)
    >>> s.name
    'Schwarzschild radius'
    """
    u = np.linspace(0, 2 * np.pi, 50)
    v = np.linspace(0, np.pi, 50)
    x_s = R_s * np.outer(np.cos(u), np.sin(v))
    y_s = R_s * np.outer(np.sin(u), np.sin(v))
    z_s = R_s * np.outer(np.ones(np.size(u)), np.cos(v))

    surface = go.Surface(
        x=x_s,
        y=y_s,
        z=z_s,
        opacity=opacity,
        showscale=False,
        name="Schwarzschild radius",
        colorscale=[[0, color], [1, color]],
    )

    return surface


def construct_circular_plane_surface(semi_side, color="black", opacity=1.0):
    """
    Build a filled disk in the plane :math:`z = 0`.

    The disk is parameterized on a ``50 x 50`` grid with radial extent set by
    ``semi_side`` (distance from the origin to the rim in the :math:`xy` plane).

    Parameters
    ----------
    semi_side : float
        Outer radius of the disk in coordinate units.
    color : str, optional
        Solid color for the surface.
    opacity : float, optional
        Opacity in ``[0, 1]``.

    Returns
    -------
    plotly.graph_objs.Surface
        A ``go.Surface`` trace for the disk.

    Examples
    --------
>>> d = construct_circular_plane_surface(1.5)
>>> d.type == "surface"
True
    """
    u = np.linspace(0, 2 * np.pi, 50)
    v = np.linspace(0, np.pi, 50)
    x_s = semi_side * np.outer(np.cos(u), np.sin(v))
    y_s = semi_side * np.outer(np.sin(u), np.sin(v))
    z_s = np.zeros_like(x_s)

    surface = go.Surface(
        x=x_s,
        y=y_s,
        z=z_s,
        opacity=opacity,
        colorscale=[[0, color], [1, color]],
        showscale=False,
    )

    return surface


def construct_square_plane_surface(semi_side, color="black", opacity=1.0):
    """
    Build a square patch in the plane :math:`z = 0`.

    The patch spans ``[-semi_side, semi_side]`` in both :math:`x` and :math:`y`
    on a ``50 x 50`` grid.

    Parameters
    ----------
    semi_side : float
        Half-width of the square along each in-plane axis.
    color : str, optional
        Solid color for the surface.
    opacity : float, optional
        Opacity in ``[0, 1]``.

    Returns
    -------
    plotly.graph_objs.Surface
        A ``go.Surface`` trace for the square patch.

    Examples
    --------
>>> p = construct_square_plane_surface(10.0, opacity=0.3)
>>> p.type == "surface"
True
    """
    x_s, y_s = np.meshgrid(
        np.linspace(-semi_side, semi_side, 50),
        np.linspace(-semi_side, semi_side, 50),
    )
    z_s = np.zeros_like(x_s)

    surface = go.Surface(
        x=x_s,
        y=y_s,
        z=z_s,
        opacity=opacity,
        colorscale=[[0, color], [1, color]],
        showscale=False,
    )

    return surface


def construct_surface(x, y, z, color="black", opacity=0.5, name="surface"):
    """
    Build a generic colored surface from precomputed coordinate meshes.

    Parameters
    ----------
    x, y, z : array_like
        2-D arrays of the same shape giving surface coordinates.
    color : str, optional
        Solid color for the surface.
    opacity : float, optional
        Opacity in ``[0, 1]``.
    name : str, optional
        Legend name for the trace.

    Returns
    -------
    plotly.graph_objs.Surface
        A ``go.Surface`` with ``showlegend=True``.

    Examples
    --------
    >>> import numpy as np
    >>> x = y = np.linspace(-1, 1, 5)
    >>> X, Y = np.meshgrid(x, y)
    >>> Z = X**2 + Y**2
    >>> t = construct_surface(X, Y, Z, name="paraboloid")
    >>> t.name
    'paraboloid'
    """
    surface = go.Surface(
        x=x,
        y=y,
        z=z,
        opacity=opacity,
        colorscale=[[0, color], [1, color]],
        showscale=False,
        name=name,
        showlegend=True,
    )
    return surface


def construct_path(path, color="red", opacity=1.0, label="line"):
    """
    Build a 3-D line trace along stored coordinates.

    Coordinate arrays are taken from ``path[1]``, ``path[2]``, and ``path[3]``
    (index ``0`` is unused), matching layouts where the first row holds an
    affine parameter or similar.

    Parameters
    ----------
    path : array_like
        Array with at least four rows; ``x = path[1]``, ``y = path[2]``,
        ``z = path[3]``.
    color : str, optional
        Line color.
    opacity : float, optional
        Opacity in ``[0, 1]``.
    label : str, optional
        Legend name for the trace.

    Returns
    -------
    plotly.graph_objs.Scatter3d
        A ``go.Scatter3d`` trace in ``lines`` mode.

    Examples
    --------
    >>> import numpy as np
    >>> t = np.linspace(0, 1, 20)
    >>> path = np.vstack([t, np.cos(t), np.sin(t), t])
    >>> tr = construct_path(path, label="helix")
    >>> tr.mode
    'lines'
    """
    return go.Scatter3d(
        x=path[1],
        y=path[2],
        z=path[3],
        mode="lines",
        line=dict(color=color, width=2),
        opacity=opacity,
        name=label,
    )


def construct_point(x, y, z, color="red", size=10, opacity=1.0, label="point"):
    """
    Build a single 3-D marker at ``(x, y, z)``.

    Parameters
    ----------
    x, y, z : float
        Coordinates of the point.
    color : str, optional
        Marker color.
    size : float, optional
        Marker size (Plotly ``marker.size``).
    opacity : float, optional
        Opacity in ``[0, 1]``.
    label : str, optional
        Legend name for the trace.

    Returns
    -------
    plotly.graph_objs.Scatter3d
        A ``go.Scatter3d`` trace in ``markers`` mode.

    Examples
    --------
    >>> pt = construct_point(0.0, 0.0, 0.0, label="origin")
    >>> pt.mode
    'markers'
    """
    return go.Scatter3d(
        x=[x],
        y=[y],
        z=[z],
        mode="markers",
        marker=dict(color=color, size=size, opacity=opacity),
        name=label,
    )


def construct_arrow(x, y, z, u=0, v=0, w=1, color="red", size=10, opacity=1.0, label="line"):
    """
    Draw a 3-D arrow using ``plotly.graph_objects.Cone``.

    The cone is anchored at ``(x, y, z)`` with direction given by ``(u, v, w)``
    (components set arrow orientation and length). Defaults point along
    :math:`+\\hat{z}`.

    Parameters
    ----------
    x, y, z : float
        Tail position of the arrow.
    u, v, w : float, optional
        Vector components defining direction and magnitude.
    color : str, optional
        Solid color for the cone.
    size : float, optional
        Scale factor; passed to ``sizeref`` as ``size / 3``.
    opacity : float, optional
        Opacity in ``[0, 1]``.
    label : str, optional
        Legend name for the trace.

    Returns
    -------
    plotly.graph_objs.Cone
        A ``go.Cone`` trace with ``anchor="tail"`` and ``sizemode="absolute"``.

    Notes
    -----
    Requires Plotly 4.x or newer (``go.Cone`` API).

    Examples
    --------
>>> c = construct_arrow(0, 0, 0, u=0, v=0, w=2, size=12)
>>> c.type == "cone"
True
    """
    return go.Cone(
        x=[x],
        y=[y],
        z=[z],
        u=[u],
        v=[v],
        w=[w],
        anchor="tail",
        showscale=False,
        sizemode="absolute",
        colorscale=[[0, color], [1, color]],
        opacity=opacity,
        sizeref=size / 3,
        name=label,
    )


def construct_arrow_scatter(
    x, y, z, u=0, v=0, w=1, color="red", line_width=4, opacity=1.0, label="arrow"
):
    """
    Draw a 3-D arrow using ``Scatter3d`` line segments (shaft + V-shaped head).

    Line width is in screen pixels, so it scales with zoom unlike ``go.Cone``.
    The arrow runs from ``(x, y, z)`` to ``(x + u, y + v, z + w)``. If the
    vector length is below ``1e-12``, only a degenerate line trace is returned.

    Parameters
    ----------
    x, y, z : float
        Tail position of the arrow.
    u, v, w : float, optional
        Vector components defining direction and magnitude.
    color : str, optional
        Line color.
    line_width : float, optional
        Line width in pixels for all segments.
    opacity : float, optional
        Opacity in ``[0, 1]``.
    label : str, optional
        Legend name; also used as ``legendgroup``.

    Returns
    -------
    list of plotly.graph_objs.Scatter3d
        One trace if the vector is (near) zero; otherwise three traces (shaft
        and two head segments).

    Examples
    --------
    >>> traces = construct_arrow_scatter(0, 0, 0, 1, 0, 0, label="x")
    >>> len(traces)
    3
    """
    tip_x, tip_y, tip_z = x + u, y + v, z + w
    length = np.sqrt(u * u + v * v + w * w)
    if length < 1e-12:
        line_trace = go.Scatter3d(
            x=[x],
            y=[y],
            z=[z],
            mode="lines",
            line=dict(color=color, width=line_width),
            opacity=opacity,
            name=label,
            legendgroup=label,
            showlegend=True,
        )
        return [line_trace]
    # Unit direction and arrowhead (V) in a plane: symmetric V with one perpendicular
    dx, dy, dz = u / length, v / length, w / length
    head_len = 0.2 * length
    spread = 0.12 * length
    back_x = tip_x - head_len * dx
    back_y = tip_y - head_len * dy
    back_z = tip_z - head_len * dz
    # One unit vector perpendicular to (dx, dy, dz); V wings = back ± spread*perp
    if abs(dz) < 0.9:
        ref = np.array([0, 0, 1])
    else:
        ref = np.array([1, 0, 0])
    perp = np.array([dy * ref[2] - dz * ref[1], dz * ref[0] - dx * ref[2], dx * ref[1] - dy * ref[0]])
    n = np.sqrt(perp[0] ** 2 + perp[1] ** 2 + perp[2] ** 2)
    if n > 1e-12:
        perp = perp / n
    else:
        perp = np.array([1, 0, 0])
    left_x = back_x + spread * perp[0]
    left_y = back_y + spread * perp[1]
    left_z = back_z + spread * perp[2]
    right_x = back_x - spread * perp[0]
    right_y = back_y - spread * perp[1]
    right_z = back_z - spread * perp[2]
    shaft = go.Scatter3d(
        x=[x, tip_x],
        y=[y, tip_y],
        z=[z, tip_z],
        mode="lines",
        line=dict(color=color, width=line_width),
        opacity=opacity,
        name=label,
        legendgroup=label,
        showlegend=True,
    )
    head_left = go.Scatter3d(
        x=[tip_x, left_x],
        y=[tip_y, left_y],
        z=[tip_z, left_z],
        mode="lines",
        line=dict(color=color, width=line_width),
        opacity=opacity,
        legendgroup=label,
        showlegend=False,
    )
    head_right = go.Scatter3d(
        x=[tip_x, right_x],
        y=[tip_y, right_y],
        z=[tip_z, right_z],
        mode="lines",
        line=dict(color=color, width=line_width),
        opacity=opacity,
        legendgroup=label,
        showlegend=False,
    )
    return [shaft, head_left, head_right]
