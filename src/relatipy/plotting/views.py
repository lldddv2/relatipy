"""Three orthogonal projections of an orbit at one exact common scale.

The widest projection is the main panel (bottom left); the panel above
shares its horizontal axis and the panel to its right its vertical axis, as
in an engineering drawing::

    top
    main  right

Panels are placed by hand: after measuring the label margins, one length
scale fills the figure width, so a metre has the same length in every panel.
With three panels a name written inside a panel would repeat, so the black
hole and the ISCOs are named once, in the legend of the free corner.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np

from ._layout import annotate_iscos, clear_joint_labels, place_bh_name
from ._scene import Scene, preview_scene, solution_scene, widest_plane
from .orbits import black_hole_handle, corner_note, draw_plane, legend_candidates
from .style import MM, Style, legend_properties, resolve_style

if TYPE_CHECKING:
    from ..geodesic import Orbit, Solution


def plot_views(source: Solution | Orbit, *, style: Style | None = None) -> tuple[Any, tuple]:
    """Draw the xy, xz and yz projections of an orbit at one common scale.

    Parameters
    ----------
    source : relatipy.geodesic.Solution or relatipy.geodesic.Orbit
        Stored trajectory, or an orbit whose native osculating preview (with
        its equatorial references) is drawn.
    style : Style or None, optional
        Appearance; the figure is ``style.views.width_factor`` times the
        target's single-figure width.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure without an implicit show or save operation.
    axes : tuple of Axes
        ``(main, top, right)``.

    Raises
    ------
    ImportError
        If Matplotlib is unavailable.
    TypeError
        If ``style`` is neither :class:`~relatipy.plotting.Style` nor
        ``None``.
    ValueError
        If the stored positions are invalid or, for an orbit, its osculating
        preview is undefined (see :func:`preview_orbit`).

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from astropy import units as u
    >>> from astropy.constants import c
    >>> from relatipy import Kerr
    >>> from relatipy.plotting import plot_views
    >>> bh = Kerr(mass=1 * u.Msun, spin=0.5)
    >>> orbit = bh.orbit(x=12 * bh.r_g, vy=0.1 * c)
    >>> fig, (main, top, right) = plot_views(orbit)
    >>> plt.close(fig)
    """
    style = resolve_style(style)
    if hasattr(source, "solve"):
        scene = preview_scene(source)
    else:
        scene = solution_scene(source)
    return render_views(scene, style)


def render_views(scene: Scene, style: Style):
    """Render three orthogonal scene projections at a shared physical scale."""
    import matplotlib.pyplot as plt

    horizontal, vertical, third = widest_plane(scene.extent)
    low, high = scene.extent.min(axis=0), scene.extent.max(axis=0)
    margin = style.views.padding * float(np.max(high - low))
    limits = {axis: (low[i] - margin, high[i] + margin) for i, axis in enumerate("xyz")}
    span = {axis: limits[axis][1] - limits[axis][0] for axis in "xyz"}

    width = style.target.width_mm * style.views.width_factor * MM
    fig = plt.figure(figsize=(width, width))
    fig.set_facecolor(style.background)
    main, top, right, corner = (fig.add_axes((0.1, 0.1, 0.1, 0.1)) for _ in range(4))
    corner.axis("off")
    panels = ((main, horizontal, vertical), (top, horizontal, third),
              (right, third, vertical))
    candidates = []
    for ax, h_axis, v_axis in panels:
        drawn = draw_plane(scene, ax, h_axis, v_axis, style, name_black_hole=False)
        if ax is main:
            candidates = drawn
        ax.set_xlim(limits[h_axis])
        ax.set_ylim(limits[v_axis])
        ax.set_aspect("auto")  # the hand layout below sets the exact scale
    top.tick_params(labelbottom=False)
    top.set_xlabel("")
    right.tick_params(labelleft=False)
    right.set_ylabel("")

    handles, names = legend_candidates(candidates, style)
    if style.center_marker is not None and style.bh_name:
        handles.append(black_hole_handle(style))
        names.append(style.bh_name)

    def corner_legend() -> None:
        """Replace the free-corner legend using current scene candidates."""
        if corner.get_legend() is not None:
            corner.get_legend().remove()
        if handles and len(handles) >= style.legend_min_entries:
            legend = corner.legend(handles, names, loc="center left",
                                   bbox_to_anchor=(0.0, 0.5), borderaxespad=0.0,
                                   **legend_properties(style))
            for text in legend.texts:
                text.set_math_fontfamily(style.math_fontset)

    corner_legend()
    gap = style.views.gap_mm * MM

    def layout(left: float, bottom: float, top_pad: float, right_pad: float) -> None:
        """Place the panels for margins in inches; the width is fixed."""
        scale = (width - left - right_pad - gap) / (span[horizontal] + span[third])
        height = bottom + scale * (span[vertical] + span[third]) + gap + top_pad
        fig.set_size_inches(width, height)
        main_w, main_h = scale * span[horizontal], scale * span[vertical]
        side = scale * span[third]
        boxes = {main: (left, bottom, main_w, main_h),
                 top: (left, bottom + main_h + gap, main_w, side),
                 right: (left + main_w + gap, bottom, side, main_h),
                 corner: (left + main_w + gap, bottom + main_h + gap, side, side)}
        for ax, (x, y, w, h) in boxes.items():
            ax.set_position((x / width, y / height, w / width, h / height))

    # Measure how far labels reach beyond the panels and make room for them.
    pad = style.outer_margin_mm * MM
    left, bottom, top_pad, right_pad = 0.6, 0.45, pad, pad
    for _ in range(4):
        layout(left, bottom, top_pad, right_pad)
        fig.canvas.draw()
        tight = fig.get_tightbbox(fig.canvas.get_renderer())
        current_width, current_height = fig.get_size_inches()
        needed = (pad - tight.x0, pad - tight.y0, tight.y1 - current_height + pad,
                  tight.x1 - current_width + pad)
        if max(abs(value) for value in needed) < 0.005:
            break
        left, bottom = left + needed[0], bottom + needed[1]
        top_pad, right_pad = top_pad + needed[2], right_pad + needed[3]

    # Touching panels: hide a secondary tick label that meets a main one.
    clear_joint_labels(top.yaxis, main.yaxis)
    clear_joint_labels(right.xaxis, main.xaxis)

    # ISCO arrows in the face-on panel; the names stay in the legend.
    for ax, h_axis, v_axis in panels:
        if {h_axis, v_axis} == {"x", "y"}:
            annotate_iscos(ax, style, write_names=False)
        place_bh_name(ax)
    text = corner_note(style, scene)
    if text:
        target = style.target
        legend = corner.get_legend()
        options = {"ha": "left", "fontsize": target.small_size,
                   "fontfamily": style.text_font}
        if legend is not None:
            fig.canvas.draw()
            box = legend.get_window_extent(fig.canvas.get_renderer())
            inner = legend.borderpad * target.small_size * fig.dpi / 72
            x, y = corner.transAxes.inverted().transform((box.x0 + inner, box.y0))
            corner.text(x, y, text, transform=corner.transAxes, va="top", **options)
        else:
            corner.text(0.0, 0.5, text, transform=corner.transAxes, va="center",
                        **options)
    return fig, (main, top, right)
