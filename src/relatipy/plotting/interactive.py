"""Interactive Plotly orbit figures (``interactive=True``, default in 3D).

Lengths are shown in the style's unit, also in hover labels: the data of the
traces are divided by that unit. In 3D the outer horizon and the ergosurface
are closed surfaces whose meridian profiles come from native C
(``Kerr._surface_profile``); only their revolution about the spin axis is
done here. The equatorial plane ``z = 0`` is shown translucent.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np
from astropy import units as u

from ._scene import AXES, ISCO_SENSE, Scene, preview_scene, solution_scene
from .style import Style, resolve_style

if TYPE_CHECKING:
    from ..geodesic import Orbit, Solution

_DASHES = {"-": "solid", "--": "dash", ":": "dot", "-.": "dashdot"}


def plot_solution_interactive(solution: Solution, *, projection: str = "3d",
                              style: Style | None = None, fig: Any = None) -> Any:
    """Plotly figure of stored samples, with 3D as its default projection.

    Equivalent to ``plot_solution(projection="3d", interactive=True)`` when
    its default projection is used.

    Parameters
    ----------
    solution : relatipy.geodesic.Solution
        Stored Cartesian samples and integration outcome.
    projection : {"3d", "xy", "xz", "yz"}, optional
        Cartesian view, default ``"3d"``.
    style : Style or None, optional
        Appearance.
    fig : plotly.graph_objects.Figure or None, optional
        Figure to extend, or ``None`` to create one.

    Returns
    -------
    plotly.graph_objects.Figure
        Hover labels give positions in the displayed unit and proper time.

    Raises
    ------
    ImportError
        If the optional Plotly dependency is unavailable.
    TypeError
        If ``fig`` is not a Plotly figure or ``style`` is not a
        :class:`~relatipy.plotting.Style`.
    ValueError
        If the projection, positions or proper times are invalid.

    Examples
    --------
    >>> from astropy import units as u
    >>> from astropy.constants import c
    >>> from relatipy import Kerr
    >>> from relatipy.plotting import plot_solution_interactive
    >>> bh = Kerr(mass=1 * u.Msun, spin=0.5)
    >>> solution = bh.orbit(x=12 * bh.r_g, vy=0.1 * c).solve(
    ...     tau_span=(0 * u.s, 2e-8 * u.s), method="dp45")
    >>> fig = plot_solution_interactive(solution)
    >>> type(fig).__name__
    'Figure'
    """
    return render_interactive(solution_scene(solution), projection,
                              resolve_style(style), fig)


def preview_orbit_interactive(orbit: Orbit, *, projection: str = "3d",
                              show_horizon: bool = True, show_isco: bool = True,
                              style: Style | None = None, fig: Any = None) -> Any:
    """Plotly figure of the native osculating preview and its references.

    Same as ``preview_orbit(projection="3d", interactive=True)``; see
    :func:`relatipy.plotting.preview_orbit`.

    Parameters
    ----------
    orbit : relatipy.geodesic.Orbit
        Scalar orbit whose current osculating Kepler conic is drawn.
    projection : {"3d", "xy", "xz", "yz"}, optional
        Cartesian view, default ``"3d"``.
    show_horizon, show_isco : bool, optional
        Include the outer-horizon and prograde/retrograde ISCO references,
        respectively. Both default to ``True``.
    style : Style or None, optional
        Appearance. ``None`` uses :data:`DEFAULT_STYLE`.
    fig : plotly.graph_objects.Figure or None, optional
        Figure to extend, or ``None`` to create one.

    Returns
    -------
    plotly.graph_objects.Figure
        Figure containing the instantaneous preview. It is neither an
        integrated trajectory nor an operation that advances ``orbit``.

    Raises
    ------
    ImportError
        If the optional Plotly dependency is unavailable.
    TypeError
        If ``fig`` is not a Plotly figure, ``style`` is not a
        :class:`~relatipy.plotting.Style`, or ``show_horizon`` or
        ``show_isco`` is not a bool.
    ValueError
        If the projection is invalid, the osculating angles are undefined,
        or the conic is parabolic.
    RuntimeError
        If the native preview returns invalid positions or reference radii.

    Examples
    --------
    >>> from astropy import units as u
    >>> from astropy.constants import c
    >>> from relatipy import Kerr
    >>> from relatipy.plotting import preview_orbit_interactive
    >>> bh = Kerr(mass=1 * u.Msun, spin=0.5)
    >>> orbit = bh.orbit(x=12 * bh.r_g, vy=0.1 * c)
    >>> fig = preview_orbit_interactive(orbit, projection="xy")
    >>> type(fig).__name__
    'Figure'
    """
    return render_interactive(preview_scene(orbit, show_horizon, show_isco),
                              projection, resolve_style(style), fig)


def _plotly():
    """Import and return Plotly's graph-object module with an actionable error."""
    try:
        import plotly.graph_objects as go
    except ImportError as exc:
        raise ImportError(
            "Interactive plotting requires Plotly; install relatipy[interactive] "
            "or pass interactive=False for a static figure"
        ) from exc
    return go


def render_interactive(scene: Scene, projection: str, style: Style, fig=None):
    """Draw ``scene`` with Plotly in 3D or on one plane."""
    if projection not in ("3d", "xy", "xz", "yz"):
        raise ValueError("projection must be 3d, xy, xz or yz")
    go = _plotly()
    if fig is not None and not isinstance(fig, go.Figure):
        raise TypeError("fig must be a Plotly Figure")
    if fig is None:
        fig = go.Figure()
    settings = style.interactive
    family = f"{style.text_font}, serif"
    scale, unit = scene.length_scale(style.length_unit)
    unit_html = "<i>r</i><sub>g</sub>" if unit == "r_g" else unit
    shown_axes = AXES if projection == "3d" else projection
    columns = [AXES.index(axis) for axis in shown_axes]
    keys = "xyz"[:len(columns)]
    trace = go.Scatter3d if projection == "3d" else go.Scatter

    def coordinates(points):
        """Map Cartesian points to the active Plotly coordinate fields."""
        points = np.atleast_2d(points) / scale
        return {key: points[:, column] for key, column in zip(keys, columns)}

    hover = "<br>".join(f"{name} [{unit}]: %{{{key}:.6g}}"
                        for name, key in zip(shown_axes, keys))
    if scene.tau is not None:
        hover += f"<br>tau [{scene.tau_unit.to_string()}]: %{{customdata:.6g}}"
    hover += "<extra>%{fullData.name}</extra>"

    def name(role: str) -> str:
        """Return a Plotly-safe display name for a scene role."""
        return style.name(role).replace("\n", "<br>")

    def line(role: str, width: float) -> dict:
        """Return Plotly line settings for a styled scene role."""
        look = style.lines[role]
        return {"color": look["color"], "width": width,
                "dash": _DASHES[look.get("linestyle", "-")]}

    role = scene.path_role
    fig.add_trace(trace(**coordinates(scene.path), mode="lines", name=name(role),
                        customdata=scene.tau, hovertemplate=hover,
                        line=line(role, settings.trajectory_width)))
    horizon_surface = (projection == "3d" and scene.metric is not None
                       and settings.horizon_surface is not None)
    for reference in scene.references:
        if reference == "horizon" and horizon_surface:
            continue  # drawn as a surface below
        fig.add_trace(trace(**coordinates(scene.circle(reference)), mode="lines",
                            name=name(reference), hoverinfo="skip",
                            line=line(reference, settings.reference_width)))
    for point_role, point in scene.points.items():
        marker = style.markers[point_role]
        index = 0 if point_role == "initial" else -1
        fig.add_trace(trace(
            **coordinates(point), mode="markers", name=name(point_role),
            showlegend=point_role not in style.legend_hide,
            customdata=None if scene.tau is None else [scene.tau[index]],
            hovertemplate=hover,
            marker={"symbol": "circle", "size": settings.marker_sizes.get(point_role, 3),
                    "color": marker["color"] or style.lines[role]["color"],
                    # The end is shown by an arrowhead; its marker is for hover.
                    "opacity": 0.0 if point_role == "end" and projection == "3d" else 1.0}))

    if projection == "3d":
        _decorate_3d(go, fig, scene, style, scale)
    center = style.center_marker
    if center is not None:
        fig.add_trace(trace(**coordinates(np.zeros(3)), mode="markers",
                            name=style.bh_name or "BH", showlegend=bool(style.bh_name),
                            hovertemplate="black hole<extra></extra>",
                            marker={"symbol": "x", "color": center["color"],
                                    "size": settings.bh_size if projection == "3d"
                                    else 7}))
    _layout(fig, scene, style, projection, shown_axes, unit_html, family)
    return fig


def _decorate_3d(go, fig, scene: Scene, style: Style, scale: float) -> None:
    """Equatorial plane, arrowheads, and horizon and ergosurface surfaces."""
    settings = style.interactive
    extent = scene.extent / scale
    span = float(np.max(np.ptp(extent, axis=0)))
    flat = {"ambient": 1.0, "diffuse": 0.0, "specular": 0.0}

    plane = settings.equatorial_plane
    if plane is not None:
        pad = plane["padding"] * span
        (x0, y0), (x1, y1) = extent[:, :2].min(axis=0) - pad, extent[:, :2].max(axis=0) + pad
        fig.add_trace(go.Mesh3d(x=[x0, x1, x1, x0], y=[y0, y0, y1, y1], z=[0, 0, 0, 0],
                                i=[0, 0], j=[1, 2], k=[2, 3], color=plane["color"],
                                opacity=plane["opacity"], flatshading=True,
                                showlegend=False, lighting=flat, name="Equatorial plane",
                                hovertemplate="equatorial plane, z = 0<extra></extra>"))

    def cone(point, direction, color, size, anchor):
        """Add one direction cone with size relative to the scene span."""
        # One trace per cone: cones sharing a trace are rescaled by spacing.
        fig.add_trace(go.Cone(x=[point[0]], y=[point[1]], z=[point[2]],
                              u=[direction[0]], v=[direction[1]], w=[direction[2]],
                              anchor=anchor, sizemode="absolute", sizeref=size * span,
                              colorscale=[[0, color], [1, color]], showscale=False,
                              showlegend=False, hoverinfo="skip", lighting=flat))

    if scene.kind == "trajectory":
        points = scene.path / scale
        step = next((points[-1] - points[i] for i in range(len(points) - 2, -1, -1)
                     if np.linalg.norm(points[-1] - points[i]) > 0), None)
        if step is not None:
            cone(points[-1], step / np.linalg.norm(step),
                 style.lines["trajectory"]["color"], settings.end_cone, "tip")

    arrows = style.isco_arrows
    for role, sense in ISCO_SENSE.items():
        if arrows is None or role not in scene.references:
            continue
        radius = scene.references[role] / scale
        for angle in 2 * np.pi * (np.arange(arrows["count"]) + 0.5) / arrows["count"]:
            cone(radius * np.array([np.cos(angle), np.sin(angle), 0.0]),
                 sense * np.array([-np.sin(angle), np.cos(angle), 0.0]),
                 style.lines[role]["color"], settings.isco_cone, "center")

    if scene.metric is None:
        return
    length = scene.unit * scale
    theta = np.linspace(0.0, np.pi, settings.surface_resolution)
    phi = np.linspace(0.0, 2.0 * np.pi, settings.surface_resolution)
    horizon = settings.horizon_surface
    event_rho, event_z = scene.metric._surface_profile("outer_horizon", theta * u.rad)
    event = float(np.max(event_rho.to_value(length)))
    for role, surface, look in (("ergosurface", "ergosurface", settings.ergosurface),
                                ("horizon", "outer_horizon", horizon)):
        if look is None:
            continue
        rho, z = scene.metric._surface_profile(surface, theta * u.rad)
        rho, z = rho.to_value(length), z.to_value(length)
        if role == "ergosurface" and horizon is not None:
            # The ergosurface meets the horizon at the poles; where they are
            # this close the two surfaces flicker, so the horizon shows.
            distance = np.hypot(rho, z) - np.hypot(event_rho.to_value(length),
                                                    event_z.to_value(length))
            keep = distance > settings.ergo_gap * event
            rho, z = rho[keep], z[keep]
            if rho.size < 2:
                continue
        fig.add_trace(go.Surface(
            x=rho[:, None] * np.cos(phi)[None, :], y=rho[:, None] * np.sin(phi)[None, :],
            z=np.repeat(z[:, None], phi.size, axis=1),
            colorscale=[[0, look["color"]], [1, look["color"]]], showscale=False,
            opacity=look["opacity"], name=style.name(role), showlegend=True,
            hovertemplate=f"{style.name(role).lower()}<extra></extra>",
            lighting={"ambient": 0.75, "diffuse": 0.35, "specular": 0.05, "fresnel": 0.0},
            contours={axis: {"highlight": False} for axis in "xyz"}))


def _layout(fig, scene: Scene, style: Style, projection: str, shown_axes: str,
            unit_html: str, family: str) -> None:
    """Apply axes, legend, labels, annotations, and camera to a Plotly figure."""
    settings = style.interactive
    axis_style = {
        "showgrid": True, "gridcolor": settings.grid_color, "zeroline": False,
        "showline": True, "linecolor": settings.line_color, "linewidth": 1,
        "ticks": "inside", "ticklen": 5, "tickcolor": settings.line_color,
        "nticks": settings.ticks,
        "tickfont": {"family": family, "size": settings.tick_size, "color": "black"},
    }

    def title(axis: str) -> dict:
        """Return the formatted Plotly title for one displayed Cartesian axis."""
        return {"text": f"<i>{axis}</i> [{unit_html}]",
                "font": {"family": family, "size": settings.font_size}}

    shown = sum(1 + (trace.name or "").count("<br>") for trace in fig.data
                if trace.showlegend is not False
                and trace.type not in ("cone", "mesh3d"))
    legend_on = shown >= style.legend_min_entries
    annotations = []
    from .orbits import corner_note

    text = corner_note(style, scene)
    if text:
        text = text.replace("\n", "<br>")
        # Below the legend (rows are about 1.6 font sizes high); the upper
        # right corner holds the Plotly tool bar.
        rows = shown if legend_on else 0
        annotations.append({"text": text, "x": 0.0, "y": 1.0, "xref": "paper",
                            "yref": "paper", "xanchor": "left", "yanchor": "top",
                            "yshift": -rows * 1.6 * settings.small_size - 6, "xshift": 4,
                            "showarrow": False,
                            "font": {"family": family, "size": settings.small_size}})
    layout = {
        "template": "simple_white", "height": settings.height,
        "font": {"family": family, "size": settings.font_size, "color": "black"},
        "paper_bgcolor": style.background, "plot_bgcolor": style.background,
        "margin": {"l": 0, "r": 0, "t": 10, "b": 0} if projection == "3d"
        else {"l": 10, "r": 10, "t": 10, "b": 10},
        "showlegend": legend_on, "annotations": annotations,
        "title": {"text": scene.title} if style.show_title else None,
        "legend": {"x": 0.0, "y": 1.0, "xanchor": "left", "yanchor": "top",
                   "bgcolor": "rgba(0,0,0,0)", "borderwidth": 0, "itemsizing": "trace",
                   "font": {"family": family, "size": settings.small_size}},
    }
    if projection == "3d":
        elevation = np.radians(style.static_3d.elevation)
        azimuth = np.radians(style.static_3d.azimuth)
        eye = settings.camera_distance * np.array([
            np.cos(elevation) * np.cos(azimuth), np.cos(elevation) * np.sin(azimuth),
            np.sin(elevation)])
        layout["scene"] = {
            **{f"{axis}axis": {**axis_style, "showbackground": False, "showspikes": False,
                               "title": title(axis)} for axis in shown_axes},
            "aspectmode": "data", "bgcolor": style.background,
            "camera": {"eye": dict(zip("xyz", eye.tolist())),
                       "up": {"x": 0, "y": 0, "z": 1}},
        }
    else:
        fig.update_xaxes(**axis_style, mirror="ticks", title=title(shown_axes[0]))
        fig.update_yaxes(**axis_style, mirror="ticks", title=title(shown_axes[1]),
                         scaleanchor="x", scaleratio=1)
    fig.update_layout(**layout)
