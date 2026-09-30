"""Posterior corner plot in the layout of GRAVITY Collaboration (2020), Fig. C.1.

Diagonal: histograms with the median (thick dashed) and the 16th and 84th
percentiles (thin dashed). Below: regions enclosing the ``Corner.levels``
fractions of the samples, filled from black inwards to light grey. Large
constants are subtracted where they shorten tick labels and are written in
the axis name; all names share one size, fitted to the panel side.
"""

from __future__ import annotations

import math
from typing import Any, Sequence

import numpy as np

from ._layout import locator, scaled_formatter
from .style import MM, Style, number, resolve_style


def display_offset(values: np.ndarray, ticks: int = 3) -> float:
    """Round constant to subtract so that tick labels are shortest.

    Candidates are 0 and the median truncated to steps of 10^k from the
    order of the central 99 % width upwards; ties keep the roundest, so an
    offset appears only when it shortens the labels.
    """
    low, median, high = np.percentile(values, (0.5, 50.0, 99.5))
    width = high - low
    if not width > 0:
        return 0.0
    ticker = locator(nbins=ticks)

    def label_length(offset: float) -> int:
        """Return the longest formatted tick label after subtracting ``offset``."""
        shown = ticker.tick_values(low - offset, high - offset)
        return max((len(number(value)) for value in shown), default=0)

    first = math.ceil(math.log10(width))
    candidates = [0.0] + [math.trunc(median / 10.0 ** k) * 10.0 ** k
                          for k in range(first + 4, first - 1, -1)]
    return min(candidates, key=label_length)


def credible_levels(density: np.ndarray, masses: Sequence[float]) -> list[float]:
    """Density thresholds enclosing the given fractions of the total."""
    ordered = np.sort(density.ravel())[::-1]
    cumulative = np.cumsum(ordered) / ordered.sum()
    return [float(ordered[min(np.searchsorted(cumulative, mass), ordered.size - 1)])
            for mass in masses]


def plot_corner(samples, labels: Sequence[str], *,
                style: Style | None = None) -> tuple[Any, np.ndarray]:
    """Lower-triangle corner plot of posterior samples.

    Parameters
    ----------
    samples : array-like
        Posterior samples with shape ``(n, k)``: ``n >= 2`` draws of
        ``k >= 2`` parameters, already in display units (e.g. degrees, mas).
    labels : sequence of str
        One axis name per parameter (``k`` names) including its unit, e.g.
        ``r"$a$ [mas]"``.
    style : Style or None, optional
        Appearance, mainly ``style.corner``.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure without an implicit show or save operation.
    axes : numpy.ndarray
        Object array of Matplotlib axes with shape ``(k, k)``. The lower
        triangle is filled; the other entries are ``None``.

    Raises
    ------
    ValueError
        If samples are not finite, have fewer than two samples or parameters,
        or ``labels`` does not have one name per parameter.
    TypeError
        If ``style`` is neither a :class:`~relatipy.plotting.Style` nor
        ``None``.
    ImportError
        If Matplotlib is unavailable.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> from relatipy.plotting import plot_corner
    >>> rng = np.random.default_rng(0)
    >>> samples = rng.normal(size=(500, 2))
    >>> fig, axes = plot_corner(samples, [r"$x$ [mas]", r"$y$ [mas]"])
    >>> axes.shape
    (2, 2)
    >>> axes[0, 1] is None
    True
    >>> plt.close(fig)
    """
    import matplotlib.pyplot as plt
    from scipy.ndimage import gaussian_filter

    style = resolve_style(style)
    settings = style.corner
    samples = np.asarray(samples, dtype=float)
    if samples.ndim != 2 or samples.shape[0] < 2 or samples.shape[1] < 2:
        raise ValueError("samples must have shape (n, k) with n, k >= 2")
    if not np.all(np.isfinite(samples)):
        raise ValueError("samples must be finite")
    count = samples.shape[1]
    labels = list(labels)
    if len(labels) != count:
        raise ValueError("labels must give one name per parameter")

    offsets = [display_offset(samples[:, i], settings.ticks) for i in range(count)]
    shown = samples - np.array(offsets)
    names = [label if offset == 0
             else f"{label} {'−' if offset > 0 else '+'} {number(abs(offset))}"
             for label, offset in zip(labels, offsets)]
    ranges = []
    for i in range(count):
        low, high = np.percentile(shown[:, i], (0.1, 99.9))
        extra = settings.margin * (high - low) or 0.5 * max(abs(low), 1.0)
        ranges.append((low - extra, high + extra))

    if settings.size_in is None:
        width, fixed_height = style.target.full_width_mm * MM, None
    else:
        width, fixed_height = settings.size_in
    fig = plt.figure(figsize=(width, fixed_height or width))
    fig.set_facecolor(style.background)
    grid = fig.add_gridspec(count, count, hspace=0.0, wspace=0.0)
    axes = np.full((count, count), None, dtype=object)
    frame = style.frame
    common = {"direction": frame.tick_direction, "width": frame.tick_width,
              "color": frame.color, "labelcolor": settings.text_color,
              "top": True, "right": True}
    for row in range(count):
        for column in range(row + 1):
            ax = fig.add_subplot(grid[row, column])
            axes[row, column] = ax
            ax.set_facecolor(style.background)
            for spine in ax.spines.values():
                spine.set_linewidth(settings.frame_width)
            x = shown[:, column]
            if row == column:
                ax.hist(x, bins=settings.hist_bins, range=ranges[column], density=True,
                        histtype="bar", rwidth=1.0, **settings.hist)
                low, median, high = np.percentile(x, (16.0, 50.0, 84.0))
                ax.axvline(median, **settings.median_line)
                for value in (low, high):
                    ax.axvline(value, **settings.quantile_line)
                ax.set_ylim(bottom=0.0)
            else:
                counts, x_edges, y_edges = np.histogram2d(
                    x, shown[:, row], bins=settings.bins,
                    range=(ranges[column], ranges[row]))
                density = gaussian_filter(counts.T, settings.smooth)
                thresholds = credible_levels(density, settings.levels)
                centers = (0.5 * (x_edges[1:] + x_edges[:-1]),
                           0.5 * (y_edges[1:] + y_edges[:-1]))
                bounds = [*thresholds[::-1], density.max() * 1.01]
                if np.all(np.diff(bounds) > 0):
                    ax.contourf(*centers, density, levels=bounds,
                                colors=list(settings.fills[:len(thresholds)][::-1]))
                    ax.contour(*centers, density, levels=thresholds[::-1][:-1],
                               colors=[settings.outline["color"]],
                               linewidths=[settings.outline["linewidth"]])
                ax.set_ylim(ranges[row])
            ax.set_xlim(ranges[column])
            ax.tick_params(which="major", length=0.7 * frame.tick_length,
                           labelsize=settings.tick_size, labelfontfamily=style.text_font,
                           **common)
            if settings.minor_ticks:
                ax.minorticks_on()
                ax.tick_params(which="minor", length=0.7 * frame.minor_tick_length,
                               **common)
            for axis in (ax.xaxis, ax.yaxis):
                axis.set_major_locator(locator(nbins=settings.ticks,
                                               inset=settings.tick_inset))
                axis.set_major_formatter(scaled_formatter())
            if row == column:  # a density axis has no useful scale
                from matplotlib.ticker import NullLocator

                ax.yaxis.set_major_locator(NullLocator())
                ax.yaxis.set_minor_locator(NullLocator())
            if row < count - 1:
                ax.tick_params(labelbottom=False)
            else:
                ax.set_xlabel(names[column], fontsize=settings.label_size,
                              fontfamily=style.text_font, color=settings.text_color)
                ax.xaxis.label.set_math_fontfamily(style.math_fontset)
            if column > 0 or row == 0:
                ax.tick_params(labelleft=False)
            else:
                ax.set_ylabel(names[row], fontsize=settings.label_size,
                              fontfamily=style.text_font, color=settings.text_color)
                ax.yaxis.label.set_math_fontfamily(style.math_fontset)
    fig.align_ylabels([axes[row, 0] for row in range(1, count)])
    fig.align_xlabels([axes[count - 1, column] for column in range(count)])

    # Outer margins measured in millimetres; the name size and tick thinning
    # depend on the panel size, so they are redone in each pass.
    pad = style.outer_margin_mm * MM
    left, bottom = 0.1, 0.1
    for _ in range(6):
        side = (width - left * width - pad) / count
        height = fixed_height or bottom * width + count * side + pad
        fig.set_size_inches(width, height)
        grid.update(left=left, right=1 - pad / width, bottom=bottom * width / height,
                    top=1 - pad / height)
        _fit_names(fig, axes, settings)
        _thin_ticks(fig, axes, settings)
        fig.canvas.draw()
        tight = fig.get_tightbbox(fig.canvas.get_renderer())
        change = (pad - tight.x0, pad - tight.y0)
        if max(abs(value) for value in change) < 0.002:
            break
        left += change[0] / width
        bottom += change[1] / width
    return fig, axes


def _fit_names(fig, axes, settings) -> None:
    """One size for every axis name: the largest at which the longest name
    fits its panel side (at most ``label_size``)."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    count = len(axes)
    side = axes[count - 1][0].get_window_extent(renderer).width
    texts = [axes[count - 1][column].xaxis.label for column in range(count)]
    texts += [axes[row][0].yaxis.label for row in range(1, count)]
    longest = max(max(text.get_window_extent(renderer).width,
                      text.get_window_extent(renderer).height) for text in texts)
    size = min(settings.label_size, texts[0].get_fontsize() * 0.96 * side / longest)
    for text in texts:
        text.set_fontsize(size)


def _thin_ticks(fig, axes, settings) -> None:
    """Horizontal bottom-row tick labels: fewer ticks where they would touch."""
    renderer = fig.canvas.get_renderer()
    gap = fig.dpi / 72
    for ax in axes[len(axes) - 1]:
        formatter = ax.xaxis.get_major_formatter()
        labels = ax.xaxis.get_ticklabels()
        properties = labels[0].get_fontproperties() if labels else None
        low, high = ax.get_xlim()
        for bins in range(settings.ticks, 0, -1):
            ticker = locator(nbins=bins, inset=settings.tick_inset)
            spans = []
            for value in ticker.tick_values(low, high):
                text_width = renderer.get_text_width_height_descent(
                    formatter(value, None), properties, ismath=False)[0]
                middle = ax.transData.transform((value, 0.0))[0]
                spans.append((middle - text_width / 2, middle + text_width / 2))
            if all(a[1] + gap <= b[0] for a, b in zip(spans, spans[1:])):
                break
        ax.xaxis.set_major_locator(ticker)
