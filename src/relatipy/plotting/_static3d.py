"""Static Matplotlib 3D orbit figures, trimmed to the printed width.

Equal physical scale: the box aspect equals the axis spans. Every name goes
to the legend; text inside a 3D scene would be foreshortened.
"""

from __future__ import annotations

import numpy as np

from ._layout import locator, scaled_formatter
from ._markers import heading, marker_path, style_marker
from ._scene import ISCO_SENSE, Scene
from .orbits import black_hole_handle, corner_note, legend_candidates
from .style import MM, Style, legend_properties


def _display(ax, points: np.ndarray) -> np.ndarray:
    """Display (pixel) coordinates of 3D points for the current view."""
    from mpl_toolkits.mplot3d import proj3d

    points = np.atleast_2d(np.asarray(points, dtype=float))
    x, y, _ = proj3d.proj_transform(points[:, 0], points[:, 1], points[:, 2],
                                    ax.get_proj())
    return ax.transData.transform(np.column_stack((x, y)))


def _box_corners(ax) -> np.ndarray:
    """Return the eight Cartesian corners of the current 3D axis limits."""
    (x0, x1), (y0, y1), (z0, z1) = ax.get_xlim(), ax.get_ylim(), ax.get_zlim()
    return np.array([[x, y, z] for x in (x0, x1) for y in (y0, y1) for z in (z0, z1)])


def _blocked(ax) -> np.ndarray:
    """Display points a legend must not cover: data, box edges and labels."""
    renderer = ax.figure.canvas.get_renderer()
    parts = [_display(ax, np.column_stack(line.get_data_3d()))
             for line in ax.lines if len(line.get_data_3d()[0])]
    for collection in ax.collections:
        offsets = getattr(collection, "_offsets3d", None)
        if offsets is not None and len(offsets[0]):
            parts.append(_display(ax, np.column_stack(offsets)))
    corners = _box_corners(ax)
    fraction = np.linspace(0.0, 1.0, 25)[:, None]
    for i, start in enumerate(corners):
        for end in corners[i + 1:]:
            if np.count_nonzero(start != end) == 1:  # one box edge
                parts.append(_display(ax, start + fraction * (end - start)))
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        for text in [*axis.get_ticklabels(), axis.label]:
            if text.get_visible() and text.get_text():
                box = text.get_window_extent(renderer)
                parts.append(np.array([[box.x0, box.y0], [box.x1, box.y1],
                                       [box.x0, box.y1], [box.x1, box.y0]]))
    return np.vstack(parts)


def _content_box(ax):
    """Display box of what is drawn (Matplotlib's tight box spans the axes)."""
    from matplotlib.transforms import Bbox

    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    boxes = [line.get_window_extent(renderer) for line in ax.lines
             if line.get_visible() and len(line.get_data_3d()[0])]
    boxes += [axis.get_tightbbox(renderer) for axis in (ax.xaxis, ax.yaxis, ax.zaxis)]
    boxes += [text.get_window_extent(renderer) for text in ax.texts]
    if ax.get_legend() is not None:
        boxes.append(ax.get_legend().get_window_extent(renderer))
    points = [_display(ax, _box_corners(ax))]
    for collection in ax.collections:
        offsets = getattr(collection, "_offsets3d", None)
        if offsets is not None and len(offsets[0]):
            reach = np.sqrt(np.max(collection.get_sizes())) / 2 * fig.dpi / 72
            centers = _display(ax, np.column_stack(offsets))
            points += [centers - reach, centers + reach]
    points = np.vstack(points)
    boxes.append(Bbox([points.min(axis=0), points.max(axis=0)]))
    return Bbox.union([box for box in boxes if box is not None and box.width >= 0])


def _trim(fig, ax, width: float, pad: float):
    """Scale the figure so the content plus ``pad`` per side is ``width`` wide;
    return the region to save, in inches."""
    from matplotlib.transforms import Bbox

    for _ in range(5):
        box = _content_box(ax)
        used = box.width / fig.dpi + 2 * pad
        if abs(used - width) < 0.002:
            break
        fig.set_size_inches(fig.get_size_inches() * width / used)
    box = _content_box(ax)
    return Bbox.from_extents(box.x0 / fig.dpi - pad, box.y0 / fig.dpi - pad,
                             box.x1 / fig.dpi + pad, box.y1 / fig.dpi + pad)


def render_static_3d(scene: Scene, style: Style, fig=None, ax=None):
    """Draw ``scene`` in Matplotlib 3D; return ``(fig, ax)``.

    A figure created here is trimmed to the target width and gets
    ``fig.relatipy_bbox`` (inches) to pass as ``bbox_inches`` when saving.
    """
    import matplotlib.pyplot as plt

    settings = style.static_3d
    owns = ax is None
    width = style.target.width_mm * MM
    if owns:
        if fig is None:
            fig = plt.figure(figsize=(width, width * settings.height_factor))
        ax = fig.add_axes((0.0, 0.0, 1.0, 1.0), projection="3d")
    else:
        fig = ax.figure
    fig.set_facecolor(style.background)
    ax.set_facecolor(style.background)
    ax.view_init(elev=settings.elevation, azim=settings.azimuth)
    ax.computed_zorder = settings.depth_sorting

    frame = style.frame
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_pane_color(settings.pane_color)
        axis.pane.set_edgecolor(settings.pane_edge_color)
        axis.pane.set_linewidth(settings.pane_edge_width)
        axis.line.set_color(frame.color)
        axis.line.set_linewidth(frame.width)
        info = axis._axinfo
        info["tick"].update({"inward_factor": settings.tick_inward,
                             "outward_factor": settings.tick_outward,
                             "linewidth": {True: frame.tick_width,
                                           False: 0.8 * frame.tick_width}})
        if settings.grid is not None:
            info["grid"].update({**settings.grid, "linestyle": "-"})
    ax.grid(settings.grid is not None)
    ax.tick_params(labelsize=style.target.tick_size, labelfontfamily=style.text_font,
                   colors=frame.color, pad=0.0)

    candidates = []
    for track in scene.tracks:
        line, = ax.plot(*track.path.T, zorder=3, gid=track.role, **track.line(style))
        candidates.append((line, track.role, track.name(style)))
    for reference in scene.references:
        circle, = ax.plot(*scene.circle(reference).T, zorder=2, gid=reference,
                          **style.lines[reference])
        candidates.append((circle, reference, style.name(reference)))

    low, high = scene.extent.min(axis=0), scene.extent.max(axis=0)
    margin = settings.padding * float(np.max(high - low))
    low, high = low - margin, high + margin
    ax.set_xlim(low[0], high[0])
    ax.set_ylim(low[1], high[1])
    ax.set_zlim(low[2], high[2])
    ax.set_box_aspect(high - low)

    scale, unit = scene.length_scale(style.length_unit)
    unit_text = r"$r_\mathrm{g}$" if unit == "r_g" else unit
    for name, axis in zip("xyz", (ax.xaxis, ax.yaxis, ax.zaxis)):
        # 3D axes keep stale labels for ticks outside the view: drop them.
        axis.set_major_locator(locator(scale, nbins=settings.ticks, inset=1e-9))
        axis.set_major_formatter(scaled_formatter(scale))
        axis.set_label_text(rf"${name}$ [{unit_text}]", fontsize=style.target.label_size,
                            fontfamily=style.text_font)
        axis.label.set_math_fontfamily(style.math_fontset)
        axis.labelpad = settings.label_pad
    if style.show_title:
        ax.set_title(scene.title, loc="left", fontsize=style.target.small_size,
                     color=style.title_color, fontfamily=style.text_font)

    # The projection is final from here on: screen angles can be measured.
    fig.canvas.draw()
    for track in scene.tracks:
        angle = heading(_display(ax, track.path))
        for point_role, point in track.points.items():
            settings_marker = style.markers[point_role]
            collection = ax.scatter(*point[:, None], depthshade=False, zorder=4,
                                    gid=point_role)
            color = settings_marker["color"] or track.line(style)["color"]
            style_marker(collection, settings_marker, angle, color)
            candidates.append((collection, point_role, style.name(point_role)))
    _isco_chevrons(ax, style)
    center = style.center_marker
    if center is not None:
        ax.scatter([0.0], [0.0], [0.0], marker=center["marker"], color=center["color"],
                   s=center["size"], linewidths=center["linewidth"], depthshade=False,
                   zorder=4, gid="center")
    if owns:
        _trim(fig, ax, width, style.outer_margin_mm * MM)
    _legend(fig, ax, scene, style, candidates)
    if owns:
        fig.relatipy_bbox = _trim(fig, ax, width, style.outer_margin_mm * MM)
    return fig, ax


def _isco_chevrons(ax, style: Style) -> None:
    """Chevrons along each ISCO circle, pointing along its projected sense."""
    arrows = style.isco_arrows
    if arrows is None:
        return
    for line in list(ax.lines):
        sense = ISCO_SENSE.get(line.get_gid())
        if sense is None:
            continue
        x, y, _ = line.get_data_3d()
        radius = float(np.max(np.hypot(x, y)))
        if radius == 0:
            continue
        count = arrows["count"]
        angles = 2 * np.pi * (np.arange(count) + 0.5) / count
        positions = radius * np.column_stack((np.cos(angles), np.sin(angles),
                                              np.zeros(count)))
        tangents = sense * np.column_stack((-np.sin(angles), np.cos(angles),
                                            np.zeros(count)))
        headings = [heading(_display(ax, np.vstack((p, p + 1e-3 * radius * t))))
                    for p, t in zip(positions, tangents)]
        chevrons = ax.scatter(*positions.T, s=arrows["size"], facecolors="none",
                              edgecolors=line.get_color(), linewidths=arrows["linewidth"],
                              depthshade=False, zorder=line.get_zorder() + 0.1)
        chevrons.set_paths([marker_path("chevron", a, tip=False) for a in headings])


def _legend(fig, ax, scene: Scene, style: Style, candidates) -> None:
    """Legend in a free corner of the drawn content (else above it), with the
    black hole listed, and the spin below it."""
    handles, names = legend_candidates(candidates, style)
    if style.center_marker is not None and style.bh_name:
        handles.append(black_hole_handle(style))
        names.append(style.bh_name)
    legend = None
    if handles and len(handles) >= style.legend_min_entries:
        content = _content_box(ax)
        renderer = fig.canvas.get_renderer()
        points = _blocked(ax)
        anchor = content.transformed(fig.transFigure.inverted()).bounds
        common = {"handles": handles, "labels": names, "bbox_to_anchor": anchor,
                  "bbox_transform": fig.transFigure, "borderaxespad": 0.0,
                  **legend_properties(style)}
        for loc in ("upper left", "upper right"):
            legend = ax.legend(**common, loc=loc)
            box = legend.get_window_extent(renderer).padded(2 * fig.dpi / 72)
            if not (((points[:, 0] >= box.x0) & (points[:, 0] <= box.x1)
                     & (points[:, 1] >= box.y0) & (points[:, 1] <= box.y1)).any()):
                break
        else:
            x0, y0, w, h = anchor
            above = {**common, "bbox_to_anchor": (x0, y0 + h, w, 0.0), "loc": "lower left"}
            for columns in range(len(handles), 0, -1):
                legend = ax.legend(**above, ncols=columns, columnspacing=1.2)
                if legend.get_window_extent(renderer).width <= content.width:
                    break
        for text in legend.texts:
            text.set_math_fontfamily(style.math_fontset)
    text = corner_note(style, scene)
    if not text:
        return
    options = {"ha": "left", "va": "top", "fontsize": style.target.small_size,
               "fontfamily": style.text_font}
    if legend is not None:
        # Anchored to the legend, below its handles: follows any resize.
        inner = legend.borderpad * style.target.small_size
        ax.annotate(text, xy=(0.0, 0.0), xycoords=legend, xytext=(inner, 0.0),
                    textcoords="offset points", **options)
    else:
        content = _content_box(ax).transformed(fig.transFigure.inverted())
        ax.text2D(content.x0, content.y1, text, transform=fig.transFigure, **options)
