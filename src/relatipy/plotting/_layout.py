"""Placement helpers for 2D Matplotlib panels.

Legends, names and annotations are placed where they cover nothing drawn,
by measuring the rendered figure. Every function works on artists tagged
with a role in their ``gid`` (see ``orbits.py``).
"""

from __future__ import annotations

import numpy as np

from ._scene import ISCO_SENSE
from .style import Style, legend_properties, number


def locator(scale: float = 1.0, nbins="auto", inset: float = 0.0):
    """Round ticks in units of ``scale``, dropping those near the axis ends.

    Data stay in their unit; ticks fall on round values of ``data / scale``
    (e.g. multiples of 5 r_g). ``inset`` is the fraction of the axis kept
    free of ticks at each end (panels that touch).
    """
    from matplotlib.ticker import MaxNLocator

    class _Locator(MaxNLocator):
        """Scale tick positions and omit positions in the configured inset."""

        def tick_values(self, vmin, vmax):
            """Return scaled major tick positions inside the visible interval."""
            low, high = sorted((vmin, vmax))
            values = super().tick_values(low / scale, high / scale) * scale
            margin = inset * (high - low)
            return values[(values >= low + margin) & (values <= high - margin)]

    return _Locator(nbins=nbins, steps=[1, 2, 2.5, 5, 10],
                    min_n_ticks=1 if inset else 2)


def scaled_formatter(scale: float = 1.0):
    """Tick labels of ``value / scale`` with a typographic minus."""
    from matplotlib.ticker import FuncFormatter

    return FuncFormatter(lambda value, _: number(value / scale))


def drawn_points(ax, include_texts: bool = True, exclude=()) -> np.ndarray:
    """Display coordinates of every vertex, marker and text box in ``ax``."""
    from matplotlib.collections import LineCollection

    parts = [line.get_xydata() for line in ax.lines]
    parts += [collection.get_offsets() for collection in ax.collections
              if len(collection.get_offsets())]
    for collection in ax.collections:
        if isinstance(collection, LineCollection):
            parts += list(collection.get_segments())
    arrays = [np.asarray(part, dtype=float).reshape(-1, 2) for part in parts]
    points = ax.transData.transform(np.vstack(arrays)) if arrays else np.empty((0, 2))
    if include_texts:
        renderer = ax.figure.canvas.get_renderer()
        for text in ax.texts:
            if text in exclude:
                continue
            box = text.get_window_extent(renderer)
            points = np.vstack((points, [[box.x0, box.y0], [box.x1, box.y0],
                                         [box.x0, box.y1], [box.x1, box.y1],
                                         [(box.x0 + box.x1) / 2, (box.y0 + box.y1) / 2]]))
    return points


def _covers(box, points: np.ndarray) -> bool:
    """Return whether a display-space box contains at least one point."""
    return bool(((points[:, 0] >= box.x0) & (points[:, 0] <= box.x1)
                 & (points[:, 1] >= box.y0) & (points[:, 1] <= box.y1)).any())


_INSIDE = ("upper left", "upper right", "lower left", "lower right",
           "center left", "center right", "upper center", "lower center")


def place_legend(ax, handles, labels, style: Style):
    """First corner or side inside the panel whose legend covers nothing;
    otherwise one row above the panel."""
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    points = drawn_points(ax)
    common = {"handles": handles, "labels": labels, **legend_properties(style)}
    for loc in _INSIDE:
        legend = ax.legend(**common, loc=loc)
        if not _covers(legend.get_window_extent(renderer).padded(2 * fig.dpi / 72), points):
            return _math_font(legend, style)
    return _math_font(legend_above(ax, common), style)


def legend_above(ax, common: dict):
    """One row above the panel; fewer columns only if the row is too wide."""
    renderer = ax.figure.canvas.get_renderer()
    width = ax.get_window_extent(renderer).width
    for columns in range(len(common["handles"]), 0, -1):
        legend = ax.legend(**common, loc="lower left", bbox_to_anchor=(0.0, 1.01),
                           ncols=columns, borderaxespad=0.0, columnspacing=1.2)
        if legend.get_window_extent(renderer).width <= width:
            break
    return legend


def _math_font(legend, style: Style):
    """Apply the selected math font family to every legend label."""
    for text in legend.texts:
        text.set_math_fontfamily(style.math_fontset)
    return legend


def fit_height(fig, pad: float) -> None:
    """Resize the figure height to its content plus ``pad`` inches per side."""
    width, height = fig.get_size_inches()
    for _ in range(4):
        fig.canvas.draw()
        used = fig.get_tightbbox(fig.canvas.get_renderer()).height
        change = height - used - 2 * pad
        if abs(change) < 0.01:
            break
        height -= change
        fig.set_size_inches(width, height)


_CORNERS = ((0.97, 0.97, "right", "top"), (0.03, 0.97, "left", "top"),
            (0.97, 0.03, "right", "bottom"), (0.03, 0.03, "left", "bottom"))


def corner_text(ax, text: str, style: Style):
    """Write ``text`` in the first inner corner that covers nothing drawn."""
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    points = drawn_points(ax)
    legend = ax.get_legend()
    blocked = [legend.get_window_extent(renderer)] if legend is not None else []
    options = {"fontsize": style.target.small_size, "fontfamily": style.text_font,
               "zorder": 5}
    for x, y, ha, va in _CORNERS:
        label = ax.text(x, y, text, transform=ax.transAxes, ha=ha, va=va, **options)
        box = label.get_window_extent(renderer).padded(2 * fig.dpi / 72)
        if not _covers(box, points) and not any(box.overlaps(b) for b in blocked):
            return label
        label.remove()
    return ax.text(1.0, 1.01, text, transform=ax.transAxes, ha="right",
                   va="bottom", **options)


def clear_joint_labels(secondary, primary) -> None:
    """Blank tick labels of ``secondary`` that overlap labels of ``primary``."""
    from matplotlib.ticker import FuncFormatter

    fig = secondary.axes.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    def shown(axis):
        """Return visible, nonempty labels and locations within an axis view."""
        low, high = sorted(axis.get_view_interval())
        return [(location, label)
                for location, label in zip(axis.get_ticklocs(), axis.get_ticklabels())
                if low <= location <= high and label.get_visible() and label.get_text()]

    taken = [label.get_window_extent(renderer) for _, label in shown(primary)]
    hidden = [float(location) for location, label in shown(secondary)
              if any(label.get_window_extent(renderer).padded(fig.dpi / 72).overlaps(t)
                     for t in taken)]
    if hidden:
        formatter = secondary.get_major_formatter()
        secondary.set_major_formatter(FuncFormatter(
            lambda value, position: "" if any(np.isclose(value, h) for h in hidden)
            else formatter(value, position)))


def place_bh_name(ax) -> None:
    """Move the black-hole name to the nearest diagonal that touches no line.

    Farther steps clear small horizon circles; beyond the second step a thin
    leader keeps the name pointing at the mark.
    """
    names = [text for text in ax.texts if text.get_gid() == "bh_name"]
    if not names:
        return
    name = names[0]
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    points = drawn_points(ax, exclude=(name,))
    center = ax.transData.transform((0.0, 0.0))
    points = points[np.hypot(*(points - center).T) > 4 * fig.dpi / 72]
    dx, dy = (abs(value) for value in name.xyann)
    for step in range(1, 13):
        for sx, sy in ((1, 1), (-1, 1), (1, -1), (-1, -1)):
            name.xyann = (sx * step * dx, sy * step * dy)
            name.set_ha("left" if sx > 0 else "right")
            name.set_va("bottom" if sy > 0 else "top")
            if not _covers(name.get_window_extent(renderer).padded(fig.dpi / 72), points):
                if step > 2:
                    ax.annotate("", xy=(0.0, 0.0), xytext=name.xyann,
                                textcoords="offset points", zorder=name.get_zorder(),
                                arrowprops={"arrowstyle": "-", "color": name.get_color(),
                                            "lw": 0.5, "shrinkA": 0, "shrinkB": 4})
                return
    name.xyann = (dx, dy)
    name.set_ha("left")
    name.set_va("bottom")


def curved_text(ax, text: str, radius: float, top: bool, color, style: Style,
                fontsize: float, draw: bool = True) -> tuple[float, float]:
    """Write ``text`` along a circle about the origin, outside its line.

    Letters are spaced by their rendered widths, so the text follows the
    circle at any size. Returns the covered angular interval in radians.
    """
    from matplotlib import font_manager

    fig = ax.figure
    renderer = fig.canvas.get_renderer()
    properties = font_manager.FontProperties(family=style.text_font, size=fontsize)
    widths = [renderer.get_text_width_height_descent(char, properties, ismath=False)[0]
              for char in text]
    ascent = renderer.get_text_width_height_descent("I", properties, ismath=False)[1]
    per_unit = abs(ax.transData.transform((1.0, 0.0))[0]
                   - ax.transData.transform((0.0, 0.0))[0])
    baseline = (radius * per_unit + style.isco_label_gap * fig.dpi / 72
                + (0.0 if top else ascent))
    span = sum(widths) / baseline
    center = 0.5 * np.pi if top else 1.5 * np.pi
    direction = -1.0 if top else 1.0  # left to right along the circle
    if draw:
        position = center - direction * span / 2
        for char, width in zip(text, widths):
            middle = position + direction * width / (2 * baseline)
            ax.text(*(baseline / per_unit * np.array([np.cos(middle), np.sin(middle)])),
                    char, rotation=np.degrees(middle) + (-90.0 if top else 90.0),
                    rotation_mode="anchor", ha="center", va="baseline", color=color,
                    fontsize=fontsize, fontfamily=style.text_font, zorder=4,
                    gid="isco_name")
            position += direction * width / baseline
    return center - span / 2, center + span / 2


def annotate_iscos(ax, style: Style, write_names: bool, min_fontsize: float = 6.0,
                   max_span: float = np.pi) -> list:
    """Name each ISCO along its circle and mark its sense (face-on views).

    A name shrinks to ``min_fontsize`` to cover at most ``max_span`` radians.
    Returns the ISCO lines whose name does not fit, for the legend.
    """
    from ._markers import marker_path

    ax.figure.canvas.draw()  # text widths need the final size
    arrows = style.isco_arrows
    unlabeled = []
    for line in ax.lines:
        sense = ISCO_SENSE.get(line.get_gid())
        if sense is None:
            continue
        radius = float(np.max(np.hypot(*line.get_xydata().T)))
        if radius == 0:
            continue
        top = sense > 0
        free_start, free_length = 0.0, 2 * np.pi
        if write_names:
            text = style.name(line.get_gid())
            size = style.target.small_size
            low, high = curved_text(ax, text, radius, top, None, style, size, draw=False)
            while high - low > max_span and size > min_fontsize:
                size = max(size - 0.5, min_fontsize)
                low, high = curved_text(ax, text, radius, top, None, style, size,
                                        draw=False)
            if high - low <= max_span:
                curved_text(ax, text, radius, top, line.get_color(), style, size)
                margin = 0.15
                free_start = high + margin
                free_length = 2 * np.pi - (high - low) - 2 * margin
            else:
                unlabeled.append(line)
        if arrows is None or free_length <= 0:
            continue
        count = arrows["count"]
        angles = free_start + free_length * (np.arange(count) + 0.5) / count
        chevrons = ax.scatter(radius * np.cos(angles), radius * np.sin(angles),
                              s=arrows["size"], facecolors="none",
                              edgecolors=line.get_color(), linewidths=arrows["linewidth"],
                              zorder=line.get_zorder() + 0.1, label="_nolegend_")
        chevrons.set_paths([marker_path("chevron", angle, tip=False)
                            for angle in np.degrees(angles) + 90.0 * sense])
    return unlabeled
