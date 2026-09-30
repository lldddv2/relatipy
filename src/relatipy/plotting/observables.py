"""Astrometric and spectroscopic report of a stellar orbit.

Layout of GRAVITY Collaboration (2020), Fig. 2: the sky-projected orbit on
the left, and RA, Dec and line-of-sight velocity against time on the right,
stacked with one shared time axis. Data and model are supplied; nothing is
fitted or evaluated here.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from astropy import units as u

from ._layout import clear_joint_labels, locator
from .style import MM, Style, format_axes, legend_properties, resolve_style, set_axis_label


def _values(value, unit: u.UnitBase, name: str) -> np.ndarray:
    """One-dimensional finite array in ``unit`` (plain numbers are taken as such)."""
    if isinstance(value, u.Quantity):
        value = value.to_value(unit)
    array = np.asarray(value, dtype=float)
    if array.ndim != 1 or array.size == 0 or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be a nonempty finite one-dimensional array")
    return array


def plot_observables(
    *, astrometry_epochs, ra, ra_err, dec, dec_err,
    spectroscopy_epochs, velocity, velocity_err,
    model_epochs, model_ra, model_dec, model_velocity,
    center=(0.0, 0.0), residual_ra=None, residual_dec=None,
    names: dict[str, str] | None = None, style: Style | None = None,
) -> tuple[Any, tuple]:
    """Plot astrometry and line-of-sight velocity against a model.

    Parameters
    ----------
    astrometry_epochs, spectroscopy_epochs, model_epochs : array-like or Quantity
        Epochs in Julian years. Every array argument is one-dimensional,
        nonempty and finite; plain numbers are taken to be in the stated
        unit and quantities are converted to it.
    ra, ra_err, dec, dec_err : array-like or Quantity
        Sky offsets and their errors, in mas.
    velocity, velocity_err : array-like or Quantity
        Line-of-sight velocity and its error, in km/s.
    model_ra, model_dec, model_velocity : array-like or Quantity
        Model on ``model_epochs`` (a dense grid gives continuous curves).
    center : tuple, optional
        Two sky offsets ``(ra, dec)`` of the central mass, each a float in
        mas or an angle quantity. Default ``(0.0, 0.0)``.
    residual_ra, residual_dec : array-like or Quantity or None, optional
        Data minus model at ``astrometry_epochs``, in mas; drawn around a
        zero line in the RA and Dec panels.
    names : dict or None, optional
        Legend names for ``"data"``, ``"model"``, ``"residual"`` and
        ``"center"``; missing keys use ``style.observables.names``.
    style : Style or None, optional
        Appearance; the figure has the target's full width.

    Returns
    -------
    fig : matplotlib.figure.Figure
    axes : tuple of Axes
        ``(sky, ra_panel, dec_panel, velocity_panel)``.

    Raises
    ------
    ValueError
        If arrays are empty, not one-dimensional, not finite, or of
        mismatched lengths, or if only one of ``residual_ra`` and
        ``residual_dec`` is given.
    astropy.units.UnitConversionError
        If a quantity cannot be converted to the stated unit.
    TypeError
        If ``style`` is neither a :class:`~relatipy.plotting.Style` nor
        ``None``.
    ImportError
        If Matplotlib is unavailable.

    Notes
    -----
    RA grows to the east, so it increases to the left on the sky. Both sky
    axes have the same mas per millimetre; the Dec axis is extended upwards
    so the legend clears the highest point.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> from relatipy.plotting import plot_observables
    >>> epochs = np.array([2000.0, 2001.0, 2002.0])
    >>> model_epochs = np.linspace(2000.0, 2002.0, 50)
    >>> fig, (sky, ra_panel, dec_panel, v_panel) = plot_observables(
    ...     astrometry_epochs=epochs, ra=[1.0, 2.0, 1.5], ra_err=[0.1, 0.1, 0.1],
    ...     dec=[0.5, 1.0, 2.0], dec_err=[0.1, 0.1, 0.1],
    ...     spectroscopy_epochs=epochs, velocity=[100.0, -50.0, 20.0],
    ...     velocity_err=[5.0, 5.0, 5.0], model_epochs=model_epochs,
    ...     model_ra=np.interp(model_epochs, epochs, [1.0, 2.0, 1.5]),
    ...     model_dec=np.interp(model_epochs, epochs, [0.5, 1.0, 2.0]),
    ...     model_velocity=np.interp(model_epochs, epochs, [100.0, -50.0, 20.0]))
    >>> plt.close(fig)
    """
    import matplotlib.pyplot as plt

    style = resolve_style(style)
    settings = style.observables
    labels = {**settings.names, **(names or {})}
    yr, mas, kms = u.yr, u.mas, u.km / u.s
    t_pos = _values(astrometry_epochs, yr, "astrometry_epochs")
    t_rv = _values(spectroscopy_epochs, yr, "spectroscopy_epochs")
    t_model = _values(model_epochs, yr, "model_epochs")
    data = {"ra": _values(ra, mas, "ra"), "ra_err": _values(ra_err, mas, "ra_err"),
            "dec": _values(dec, mas, "dec"), "dec_err": _values(dec_err, mas, "dec_err"),
            "v": _values(velocity, kms, "velocity"),
            "v_err": _values(velocity_err, kms, "velocity_err")}
    model = {"ra": _values(model_ra, mas, "model_ra"),
             "dec": _values(model_dec, mas, "model_dec"),
             "v": _values(model_velocity, kms, "model_velocity")}
    for key in ("ra", "ra_err", "dec", "dec_err"):
        if data[key].shape != t_pos.shape:
            raise ValueError(f"{key} must match astrometry_epochs")
    if data["v"].shape != t_rv.shape or data["v_err"].shape != t_rv.shape:
        raise ValueError("velocity and velocity_err must match spectroscopy_epochs")
    if any(values.shape != t_model.shape for values in model.values()):
        raise ValueError("model arrays must match model_epochs")
    center = tuple(float(c.to_value(mas)) if isinstance(c, u.Quantity) else float(c)
                   for c in center)
    residuals = None
    if residual_ra is not None or residual_dec is not None:
        if residual_ra is None or residual_dec is None:
            raise ValueError("give both residual_ra and residual_dec, or neither")
        residuals = {"ra": _values(residual_ra, mas, "residual_ra"),
                     "dec": _values(residual_dec, mas, "residual_dec")}
        if any(values.shape != t_pos.shape for values in residuals.values()):
            raise ValueError("residuals must match astrometry_epochs")

    width = style.target.full_width_mm * MM
    height = settings.height_ratio * width
    fig = plt.figure(figsize=(width, height))
    fig.set_facecolor(style.background)
    grid = fig.add_gridspec(1, 2, width_ratios=(settings.sky_width_ratio, 1.0))
    sky = fig.add_subplot(grid[0, 0])
    column = grid[0, 1].subgridspec(3, 1, hspace=0.0)
    series = [fig.add_subplot(column[row, 0]) for row in range(3)]
    for ax in series[1:]:
        ax.sharex(series[0])

    points = {"fmt": "o", "color": settings.data_color, "markersize": settings.marker_size,
              "elinewidth": settings.error_width, "capsize": 0.0, "linestyle": "none",
              "zorder": 3}
    line = {"color": settings.model_color, "linewidth": settings.model_width, "zorder": 2}
    sky.plot(model["ra"], model["dec"], label=labels["model"], **line)
    sky.errorbar(data["ra"], data["dec"], xerr=data["ra_err"], yerr=data["dec_err"],
                 label=labels["data"], **points)
    mark = style.center_marker or {"marker": "x", "color": "black", "size": 18,
                                   "linewidth": 0.9}
    sky.scatter(*center, marker=mark["marker"], color=mark["color"], s=mark["size"],
                linewidths=mark["linewidth"], zorder=4, label=labels["center"])
    format_axes(sky, style)
    set_axis_label(sky, "x", "RA [mas]", style)
    set_axis_label(sky, "y", "Dec [mas]", style)

    rows = (("ra", t_pos, "ra_err", "RA [mas]"), ("dec", t_pos, "dec_err", "Dec [mas]"),
            ("v", t_rv, "v_err", r"$v_\mathrm{los}$ [km/s]"))
    for ax, (key, epochs, error, name) in zip(series, rows):
        ax.plot(t_model, model[key], **line)
        ax.errorbar(epochs, data[key], yerr=data[error], **points)
        ax.set_xlim(t_model[0], t_model[-1])
        format_axes(ax, style)
        ax.yaxis.set_major_locator(locator(nbins=5))
        set_axis_label(ax, "y", name, style)
    if settings.zero_line is not None:
        for ax in series[:2]:
            ax.axhline(0.0, zorder=1, **settings.zero_line)
    if residuals is not None:
        red = {**points, "color": settings.residual_color}
        for ax, key in ((series[0], "ra"), (series[1], "dec")):
            ax.errorbar(t_pos, residuals[key], yerr=data[f"{key}_err"], **red)
        sky.errorbar([], [], yerr=[], label=labels["residual"], **red)
    for ax in series[:-1]:
        ax.tick_params(labelbottom=False)
    set_axis_label(series[-1], "x", "Epoch [yr]", style)
    fig.align_ylabels(series)

    handles, texts = sky.get_legend_handles_labels()
    order = [texts.index(labels[key]) for key in ("data", "model", "residual", "center")
             if labels[key] in texts]
    legend = sky.legend([handles[i] for i in order], [texts[i] for i in order],
                        loc="upper left", handletextpad=0.4, borderaxespad=0.6,
                        **legend_properties(style))
    for text in legend.texts:
        text.set_math_fontfamily(style.math_fontset)

    _layout(fig, grid, sky, series, legend, data, model, style, width, height)
    return fig, (sky, *series)


def _layout(fig, grid, sky, series, legend, data, model, style, width, height) -> None:
    """Measured margins, equal sky scale, and Dec room for the legend."""
    settings = style.observables
    pad = style.outer_margin_mm * MM
    gap = settings.column_gap_mm * MM
    dec_top = max(model["dec"].max(), (data["dec"] + data["dec_err"]).max())
    dec_bottom = min(model["dec"].min(), (data["dec"] - data["dec_err"]).min())
    ra_high = max(model["ra"].max(), (data["ra"] + data["ra_err"]).max())
    ra_low = min(model["ra"].min(), (data["ra"] - data["ra_err"]).min())
    dec_pad = 0.04 * (dec_top - dec_bottom)
    ra_pad = 0.04 * (ra_high - ra_low)
    headroom = 0.0

    def sky_limits() -> None:
        """Equal mas per millimetre on both sky axes, filling the panel box."""
        box = sky.get_position().transformed(fig.transFigure)
        ratio = box.width / box.height
        low, high = dec_bottom - dec_pad, dec_top + dec_pad + headroom
        if (high - low) * ratio < ra_high - ra_low + 2 * ra_pad:
            low = high - (ra_high - ra_low + 2 * ra_pad) / ratio
        middle, half = 0.5 * (ra_high + ra_low), 0.5 * (high - low) * ratio
        sky.set_ylim(low, high)
        sky.set_xlim(middle + half, middle - half)  # RA grows to the left

    left, bottom, right, top, wspace = 0.08, 0.08, 0.98, 0.98, 0.25
    for _ in range(10):
        grid.update(left=left, bottom=bottom, right=right, top=top, wspace=wspace)
        sky_limits()
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        legend_bottom = sky.transData.inverted().transform(
            (0.0, legend.get_window_extent(renderer).y0))[1]
        per_mas = abs(sky.transData.transform((0, 1))[1] - sky.transData.transform((0, 0))[1])
        gap_mas = settings.legend_gap_mm * MM * fig.dpi / per_mas
        headroom = max(0.0, headroom + dec_top + gap_mas - legend_bottom)
        tight = fig.get_tightbbox(renderer)
        sky_right = sky.get_tightbbox(renderer).x1 / fig.dpi
        series_left = min(ax.get_tightbbox(renderer).x0 for ax in series) / fig.dpi
        change = (pad - tight.x0, pad - tight.y0, width - pad - tight.x1,
                  height - pad - tight.y1, gap - (series_left - sky_right))
        if (max(abs(value) for value in change) < 0.002
                and abs(dec_top + gap_mas - legend_bottom) < 0.01 * dec_pad):
            break
        left += change[0] / width
        bottom += change[1] / height
        right += change[2] / width
        top += change[3] / height
        wspace += change[4] / ((right - left) * width / (2 + wspace))
    sky_limits()
    for lower, upper in zip(series[1:], series[:-1]):
        clear_joint_labels(lower.yaxis, upper.yaxis)
