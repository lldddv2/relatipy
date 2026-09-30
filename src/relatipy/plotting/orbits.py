"""Orbit figures: stored trajectories and native osculating previews.

Stored numerical trajectories and instantaneous Kepler previews have distinct
roles and names. Preview paths and equatorial reference radii come from the
native kernel; nothing here integrates, interpolates or evaluates physics.
Matplotlib and Plotly are imported only when a figure is drawn.

``interactive`` chooses the kind of figure. The default ``"views"``
projection is static Matplotlib. Explicit 3D views are interactive Plotly
figures by default; the planes ``"xy"``, ``"xz"`` and ``"yz"`` are static
Matplotlib figures in the publication format of
:mod:`relatipy.plotting.style`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np

from ._layout import (
    annotate_iscos, corner_text, fit_height, locator, place_bh_name, place_legend,
    scaled_formatter,
)
from ._markers import heading, style_marker
from ._scene import AXES, ISCO_SENSE, Scene, preview_scene, solution_scene
from .style import MM, Style, format_axes, resolve_style, set_axis_label

if TYPE_CHECKING:
    from ..geodesic import Orbit, Solution

PROJECTIONS = ("3d", "xy", "xz", "yz", "views")


def _interactive(projection: str, interactive: bool | None) -> bool:
    """Resolve ``interactive``: ``None`` means interactive 3D, static planes."""
    if projection not in PROJECTIONS:
        raise ValueError("projection must be 3d, xy, xz, yz or views")
    if interactive is None:
        return projection == "3d"
    if not isinstance(interactive, bool):
        raise TypeError("interactive must be True, False or None")
    if projection == "views" and interactive:
        raise ValueError("projection 'views' is a static figure only")
    return interactive


def _validate_axes(projection: str, fig, ax, operation: str):
    """Check plotting arguments before creating axes or adding artists."""
    try:
        import matplotlib.pyplot as plt
        from matplotlib.axes import Axes
        from matplotlib.figure import Figure
    except ImportError as exc:
        raise ImportError(f"{operation} requires Matplotlib") from exc
    if projection == "views" and (fig is not None or ax is not None):
        raise ValueError("projection 'views' creates its own figure")
    if fig is not None and not isinstance(fig, Figure):
        raise TypeError("fig must be a Matplotlib Figure")
    if ax is not None:
        if not isinstance(ax, Axes):
            raise TypeError("ax must be a Matplotlib Axes")
        if fig is not None and ax.figure is not fig:
            raise ValueError("ax must belong to fig")
        expected = "3d" if projection == "3d" else "rectilinear"
        if ax.name != expected:
            raise ValueError("ax projection does not match projection")
    return plt


def plot_solution(
    solution: Solution, *, projection: str = "views", interactive: bool | None = None,
    style: Style | None = None, fig: Any = None, ax: Any = None,
) -> Any:
    """Plot the stored Cartesian samples of an integrated trajectory.

    Parameters
    ----------
    solution : relatipy.geodesic.Solution
        Solution whose stored positions are plotted without interpolation.
    projection : {"3d", "xy", "xz", "yz", "views"}, optional
        Cartesian view, default ``"views"``. It draws the three planes at
        one common scale (:func:`plot_views`).
    interactive : bool or None, optional
        ``True``: interactive Plotly figure; ``False``: static Matplotlib
        figure. ``None`` (default) is ``True`` for ``"3d"`` and ``False``
        for the planes and ``"views"`` (always static).
    style : Style or None, optional
        Appearance; ``None`` uses :data:`DEFAULT_STYLE` and shows lengths
        in the stored position unit.
    fig : matplotlib.figure.Figure or plotly.graph_objects.Figure, optional
        Existing figure to draw into: a Matplotlib figure for a static plot
        or a Plotly figure for an interactive plot. Default ``None`` creates
        a new figure. The ``"views"`` layout creates its own figure and does
        not accept ``fig``.
    ax : matplotlib.axes.Axes, optional
        Existing static axes belonging to ``fig``, with a 3D projection for
        ``"3d"`` and a rectilinear projection otherwise. It must be ``None``
        for interactive figures and for ``"views"``.

    Returns
    -------
    figure : plotly.graph_objects.Figure or tuple
        A Plotly figure when interactive. Otherwise a Matplotlib
        ``(fig, ax)`` pair, or ``(fig, (main, top, right))`` for
        ``"views"``.

    Raises
    ------
    ImportError
        If Plotly (interactive) or Matplotlib (static) is unavailable.
    TypeError
        If ``interactive``, ``fig``, ``ax`` or ``style`` has the wrong type.
    ValueError
        If ``projection`` is unknown, ``"views"`` is combined with
        ``interactive=True``, ``fig`` or ``ax``, ``ax`` does not match
        ``fig`` or ``projection``, or the stored positions are invalid.

    Notes
    -----
    Lines connect stored samples; they do not represent additional
    integration or dense output. All spatial axes share one physical scale.
    Figures are never shown or saved.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from astropy import units as u
    >>> from astropy.constants import c
    >>> from relatipy import Kerr
    >>> from relatipy.plotting import plot_solution
    >>> bh = Kerr(mass=1 * u.Msun, spin=0.5)
    >>> solution = bh.orbit(x=12 * bh.r_g, vy=0.1 * c).solve(
    ...     tau_span=(0 * u.s, 2e-8 * u.s), method="dp45")
    >>> fig, ax = plot_solution(solution, projection="xy")
    >>> ax.figure is fig
    True
    >>> fig, (main, top, right) = plot_solution(solution)
    >>> plt.close("all")
    """
    style = resolve_style(style)
    if _interactive(projection, interactive):
        if ax is not None:
            raise ValueError("ax is only used by static figures")
        from .interactive import render_interactive

        return render_interactive(solution_scene(solution), projection, style, fig)
    _validate_axes(projection, fig, ax, "Solution.plot")
    scene = solution_scene(solution)
    return _render_matplotlib(scene, projection, style, fig, ax)


def preview_orbit(
    orbit: Orbit, *, projection: str = "views", interactive: bool | None = None,
    show_horizon: bool = True, show_isco: bool = True,
    style: Style | None = None, fig: Any = None, ax: Any = None,
) -> Any:
    """Plot the current native osculating Kepler path and references.

    Parameters
    ----------
    orbit : relatipy.geodesic.Orbit
        Current scalar orbit whose instantaneous conic is previewed.
    projection : {"3d", "xy", "xz", "yz", "views"}, optional
        Cartesian view, default ``"views"``; see :func:`plot_solution`.
    interactive : bool or None, optional
        Rendering mode, default ``None``; see :func:`plot_solution`.
    show_horizon : bool, optional
        Show the outer-horizon reference, default ``True``.
    show_isco : bool, optional
        Show the prograde and retrograde ISCO references, default ``True``.
    style : Style or None, optional
        Appearance; ``None`` uses :data:`DEFAULT_STYLE` and shows lengths in
        the orbit's stored position unit.
    fig : matplotlib.figure.Figure or plotly.graph_objects.Figure, optional
        Existing figure; see :func:`plot_solution`.
    ax : matplotlib.axes.Axes, optional
        Existing static axes; see :func:`plot_solution`.

    Returns
    -------
    figure : plotly.graph_objects.Figure or tuple
        Same forms as :func:`plot_solution`.

    Raises
    ------
    ImportError
        If Plotly (interactive) or Matplotlib (static) is unavailable.
    TypeError
        If ``interactive``, ``fig``, ``ax`` or ``style`` has the wrong type,
        or ``show_horizon`` or ``show_isco`` is not a bool.
    ValueError
        If the plotting arguments are invalid as in :func:`plot_solution`,
        or the osculating angles are undefined or the conic is parabolic.
    RuntimeError
        If the native preview returns invalid positions or radii.

    Notes
    -----
    This instantaneous Kepler conic is not an integrated Kerr trajectory.
    Horizon and prograde/retrograde ISCO circles are equatorial references
    whose Cartesian radii, like the path, are supplied by native C. The
    orbit is never advanced.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from astropy import units as u
    >>> from astropy.constants import c
    >>> from relatipy import Kerr
    >>> from relatipy.plotting import preview_orbit
    >>> bh = Kerr(mass=1 * u.Msun, spin=0.5)
    >>> orbit = bh.orbit(x=12 * bh.r_g, vy=0.1 * c)
    >>> fig, ax = preview_orbit(orbit, projection="xy", show_isco=False)
    >>> bool(orbit.tau == orbit.initial.tau)
    True
    >>> plt.close(fig)
    """
    style = resolve_style(style)
    if _interactive(projection, interactive):
        if ax is not None:
            raise ValueError("ax is only used by static figures")
        from .interactive import render_interactive

        return render_interactive(preview_scene(orbit, show_horizon, show_isco),
                                  projection, style, fig)
    _validate_axes(projection, fig, ax, "Orbit.preview")
    scene = preview_scene(orbit, show_horizon, show_isco)
    return _render_matplotlib(scene, projection, style, fig, ax)


def _render_matplotlib(scene: Scene, projection: str, style: Style, fig, ax):
    """Render a scene with Matplotlib and return its figure and axes."""
    import matplotlib.pyplot as plt

    if projection == "views":
        from .views import render_views

        return render_views(scene, style)
    if projection == "3d":
        from ._static3d import render_static_3d

        return render_static_3d(scene, style, fig, ax)
    owns = ax is None
    if owns:
        width = style.target.width_mm * MM
        if fig is None:
            fig = plt.figure(figsize=(width, width), layout="constrained")
            margin = style.outer_margin_mm * MM
            fig.get_layout_engine().set(w_pad=margin, h_pad=margin)
        ax = fig.add_subplot(111)
    else:
        fig = ax.figure
    fig.set_facecolor(style.background)
    handles = draw_plane(scene, ax, projection[0], projection[1], style)
    finish_panel(ax, scene, style, handles, resize=owns and fig.get_layout_engine()
                 is not None)
    return fig, ax


def draw_plane(scene: Scene, ax, horizontal: str, vertical: str, style: Style,
               name_black_hole: bool = True) -> list[tuple[Any, str]]:
    """Draw the scene on one Cartesian plane; return legend candidates.

    Data stay in the scene unit; ticks show lengths in the style's unit.
    Each artist carries its role as ``gid``.
    """
    columns = (AXES.index(horizontal), AXES.index(vertical))
    candidates = []

    def project(points):
        """Select the horizontal and vertical coordinates of Cartesian points."""
        points = np.atleast_2d(points)
        return points[:, columns[0]], points[:, columns[1]]

    path_role = scene.path_role
    path_line, = ax.plot(*project(scene.path), gid=path_role, zorder=3,
                         **style.lines[path_role])
    candidates.append((path_line, path_role))
    for role in scene.references:
        line, = ax.plot(*project(scene.circle(role)), gid=role, zorder=2,
                        **style.lines[role])
        candidates.append((line, role))

    angle = heading(np.column_stack(project(scene.path)))
    for role, point in scene.points.items():
        settings = style.markers[role]
        collection = ax.scatter(*project(point), zorder=4, gid=role)
        color = settings["color"] or style.lines[path_role]["color"]
        style_marker(collection, settings, angle, color)
        candidates.append((collection, role))

    center = style.center_marker
    if center is not None:
        ax.scatter([0.0], [0.0], marker=center["marker"], color=center["color"],
                   s=center["size"], linewidths=center["linewidth"], zorder=4,
                   gid="center")
        if name_black_hole and style.bh_name:
            name = ax.annotate(style.bh_name, xy=(0.0, 0.0),
                               xytext=style.bh_name_offset, textcoords="offset points",
                               ha="left", va="bottom", fontsize=style.target.small_size,
                               fontfamily=style.text_font, color=center["color"],
                               zorder=4, gid="bh_name")
            name.set_math_fontfamily(style.math_fontset)

    scale, unit = scene.length_scale(style.length_unit)
    unit_text = r"$r_\mathrm{g}$" if unit == "r_g" else unit
    for which, axis, name in (("x", ax.xaxis, horizontal), ("y", ax.yaxis, vertical)):
        axis.set_major_locator(locator(scale))
        set_axis_label(ax, which, rf"${name}$ [{unit_text}]", style)
    format_axes(ax, style)
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_formatter(scaled_formatter(scale))
    # The frame follows the data aspect; the figure height adapts to it.
    ax.set_aspect("equal", adjustable="box")
    if style.show_title:
        ax.set_title(scene.title, loc="left", fontsize=style.target.small_size,
                     color=style.title_color, fontfamily=style.text_font)
    return candidates


def legend_candidates(candidates, style: Style, exclude=()) -> tuple[list, list]:
    """Handles and names of the roles a legend lists."""
    handles, names = [], []
    for handle, role in candidates:
        if role in style.legend_hide or role in exclude:
            continue
        handles.append(handle)
        names.append(style.name(role))
    return handles, names


def black_hole_handle(style: Style):
    """Legend handle of the black-hole mark."""
    from matplotlib.lines import Line2D

    center = style.center_marker
    return Line2D([], [], linestyle="none", marker=center["marker"], color=center["color"],
                  markersize=np.sqrt(center["size"]), markeredgewidth=center["linewidth"])


def face_on(ax) -> bool:
    """The equatorial circles are seen face-on (an xy or yx view)."""
    return {ax.get_xlabel()[1:2], ax.get_ylabel()[1:2]} == {"x", "y"}


def finish_panel(ax, scene: Scene, style: Style, candidates, resize: bool) -> None:
    """Legend, ISCO names and arrows, black-hole name and spin of one panel."""
    fig = ax.figure
    write_names = style.isco_labels and face_on(ax)
    iscos = [role for role in ISCO_SENSE if role in scene.references]
    handles, names = legend_candidates(candidates, style,
                                       exclude=iscos if write_names else ())
    names_on_mark = [text for text in ax.texts if text.get_gid() == "bh_name"]
    legend = None
    if len(handles) + bool(names_on_mark) >= style.legend_min_entries and handles:
        # With a legend, the black hole is named there, not in the plot.
        if names_on_mark:
            handles.append(black_hole_handle(style))
            names.append(names_on_mark[0].get_text())
            names_on_mark[0].remove()
        legend = place_legend(ax, handles, names, style)
    if resize:
        fit_height(fig, style.outer_margin_mm * MM)
    if face_on(ax) and (write_names or style.isco_arrows is not None):
        unlabeled = annotate_iscos(ax, style, write_names)
        if unlabeled:
            roles = [line.get_gid() for line in unlabeled]
            if legend is not None:
                legend.remove()
            place_legend(ax, handles + unlabeled, names + [style.name(r) for r in roles],
                         style)
    place_bh_name(ax)
    text = corner_note(style, scene)
    if text:
        corner_text(ax, text, style)


def corner_note(style: Style, scene: Scene) -> str | None:
    """Corner note: the Kerr spin (if shown) and, always, an early end of
    the integration (e.g. "partial: numerical failure"); ``None`` if empty."""
    from .style import number

    lines = []
    if style.spin_text is not None and scene.spin is not None:
        lines.append(style.spin_text.format(spin=number(scene.spin)))
    if scene.status is not None:
        lines.append(scene.status)
    return "\n".join(lines) or None


__all__ = ["plot_solution", "preview_orbit"]
