"""Time series of stored coordinates and osculating elements."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np
from astropy import units as u

from ._common import trajectory_status
from ._layout import locator, scaled_formatter
from .style import MM, Style, format_axes, publication_style, resolve_style, set_axis_label

if TYPE_CHECKING:
    from ..geodesic import Solution


_FAMILIES = {
    "cartesian": ("x", "y", "z"),
    "spherical": ("r", "theta", "phi"),
    "bl": ("R", "Theta", "Phi"),
    "elements": ("a", "e", "f"),
}
_COMPONENTS = tuple(name for family in _FAMILIES.values() for name in family)
_ANGLE_COMPONENTS = frozenset(("theta", "phi", "Theta", "Phi", "f"))
_LATEX_NAMES = {
    "theta": r"\theta",
    "Theta": r"\Theta",
    "phi": r"\phi",
    "Phi": r"\Phi",
    "tau": r"\tau",
}


def _tick_factor(values: np.ndarray) -> tuple[float, str]:
    """Select an axis-label factor while keeping plotted values unchanged."""
    finite = np.abs(values[np.isfinite(values)])
    maximum = float(np.max(finite)) if finite.size else 0.0
    if maximum == 0.0:
        return 1.0, ""
    exponent = int(np.floor(np.log10(maximum)))
    if -3 < exponent < 4:
        return 1.0, ""
    return 10.0 ** exponent, rf" $\times 10^{{{exponent}}}$"


def plot_evolution(
    solution: Solution, *, time: str = "t", coords: str = "cartesian",
    style: Style | None = None, **components: bool,
) -> tuple[Any, np.ndarray]:
    """Plot selected stored positions or elements against time.

    ``coords`` selects the three components enabled by default. Explicit
    boolean component keywords override those defaults, including ``False``.
    Other coordinate families can be added with explicit ``True`` keywords.
    The accepted component names are ``x``, ``y``, ``z``, ``r``, ``theta``,
    ``phi``, ``R``, ``Theta``, ``Phi``, ``a``, ``e`` and ``f``.

    Parameters
    ----------
    solution : relatipy.geodesic.Solution
        Source of stored samples; no interpolation or integration is done.
    time : {"t", "tau"}, optional
        Horizontal coordinate. Default is coordinate time ``"t"``.
    coords : {"cartesian", "cartesial", "spherical", "bl", "elements"}, optional
        Default family. ``"elements"`` selects stored osculating ``a``,
        ``e`` and ``f``. ``"cartesial"`` is an alias for
        ``"cartesian"``. Default is Cartesian.
    style : Style or None, optional
        Figure format; ``None`` uses the default RelatiPy style.
    **components : bool
        Named component switches. Omitted switches inherit the selected
        family's default. Explicit ``False`` excludes a component.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure without an implicit show or save operation.
    axes : numpy.ndarray
        One-dimensional axes array in component order, one panel per series.

    Raises
    ------
    ValueError
        If ``time`` or ``coords`` is unknown, a component name is unknown,
        all components are disabled, or the stored ``time`` samples are not
        finite.
    TypeError
        If a component switch is not a boolean or ``style`` is invalid.
    ImportError
        If Matplotlib is unavailable.

    Notes
    -----
    Lengths retain their stored units, eccentricity is dimensionless, and
    angles are displayed in degrees. Time retains its stored unit. Tick labels
    use a power-of-ten factor in the axis label when useful, without rescaling
    plotted values. Panels share borders without vertical gaps.
    Lines join stored samples only. Undefined element values leave gaps.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from astropy import units as u
    >>> from astropy.constants import c
    >>> from relatipy import Kerr
    >>> from relatipy.plotting import plot_evolution
    >>> bh = Kerr(mass=1 * u.Msun, spin=0.5)
    >>> solution = bh.orbit(x=12 * bh.r_g, vy=0.1 * c).solve(
    ...     tau_span=(0 * u.s, 2e-8 * u.s), method="dp45")
    >>> fig, axes = plot_evolution(solution, time="tau", coords="bl", Phi=False)
    >>> axes.shape
    (2,)
    >>> plt.close(fig)
    """
    if time not in ("t", "tau"):
        raise ValueError("time must be 't' or 'tau'")
    if coords == "cartesial":
        coords = "cartesian"
    if coords not in _FAMILIES:
        raise ValueError(
            "coords must be 'cartesian', 'cartesial', 'spherical', 'bl' or 'elements'"
        )
    unknown = components.keys() - _COMPONENTS
    if unknown:
        raise ValueError(f"unknown coordinate components: {', '.join(sorted(unknown))}")
    for name, enabled in components.items():
        if not isinstance(enabled, bool):
            raise TypeError(f"{name} must be a bool")
    selected = tuple(
        name for name in _COMPONENTS
        if components.get(name, name in _FAMILIES[coords])
    )
    if not selected:
        raise ValueError("at least one coordinate component must be enabled")
    style = resolve_style(style)
    horizontal = _horizontal(solution, time)

    elements = solution.orbital_elements() if any(
        name in ("a", "e", "f") for name in selected
    ) else None
    series = []
    for name in selected:
        quantity = getattr(elements, name) if name in ("a", "e", "f") else getattr(solution, name)
        unit = u.deg if name in _ANGLE_COMPONENTS else (u.one if name == "e" else quantity.unit)
        values = np.asarray(quantity if name == "e" else quantity.to_value(unit), dtype=float)
        if name in ("a", "e", "f"):
            values = np.where(np.isfinite(values), values, np.nan)
        elif not np.all(np.isfinite(values)):
            raise ValueError(f"stored {name} samples must be finite")
        unit_label = f" [{unit.to_string()}]" if unit != u.one else ""
        series.append((f"${_LATEX_NAMES.get(name, name)}${unit_label}", values))
    return _draw_panels(solution, time, horizontal, series, style)


def _horizontal(solution: Solution, time: str) -> u.Quantity:
    """Return finite stored ``t`` or ``tau`` samples for the horizontal axis."""
    if time not in ("t", "tau"):
        raise ValueError("time must be 't' or 'tau'")
    horizontal = getattr(solution, time)
    if not np.all(np.isfinite(horizontal.to_value(horizontal.unit))):
        raise ValueError(f"stored {time} samples must be finite")
    return horizontal


def _draw_panels(
    solution: Solution, time: str, horizontal: u.Quantity,
    series: list[tuple[str, np.ndarray]], style: Style, panel_mm: float = 30.0,
) -> tuple[Any, np.ndarray]:
    """Draw labelled series in stacked panels sharing the time axis."""
    import matplotlib.pyplot as plt

    horizontal_values = horizontal.to_value(horizontal.unit)

    width = style.target.width_mm * MM
    height = max(64.0, panel_mm * len(series)) * MM
    with publication_style(style):
        fig, grid = plt.subplots(
            len(series), 1, sharex=True, figsize=(width, height),
            gridspec_kw={"hspace": 0.0}, squeeze=False,
        )
    fig.set_facecolor(style.background)
    fig.subplots_adjust(left=0.26, right=0.97, bottom=0.13, top=0.96, hspace=0.0)
    axes = grid.ravel()
    for ax in axes:
        format_axes(ax, style)
    axes[-1].xaxis.set_major_locator(locator(nbins=5))
    time_factor, time_suffix = _tick_factor(np.asarray(horizontal_values))
    if time_suffix:
        axes[-1].xaxis.set_major_formatter(scaled_formatter(time_factor))
    line = dict(style.lines["trajectory"])
    if len(solution) == 1:
        line.setdefault("marker", "o")
    for ax, (label, values) in zip(axes, series):
        ax.plot(horizontal_values, values, **line)
        factor, suffix = _tick_factor(values)
        if suffix:
            ax.yaxis.set_major_formatter(scaled_formatter(factor))
        set_axis_label(ax, "y", f"{label}{suffix}", style)
    for ax in axes[:-1]:
        ax.tick_params(labelbottom=False)
    set_axis_label(
        axes[-1], "x",
        f"${_LATEX_NAMES.get(time, time)}$ [{horizontal.unit.to_string()}]{time_suffix}",
        style,
    )
    status = trajectory_status(solution)
    if status is not None:
        fig.suptitle(status, color=style.title_color)
    return fig, axes


_CONSTANTS = (
    ("E", "get_E", "E", "E_0"),
    ("Lz", "get_Lz", "L_z", "L_{z,0}"),
    ("Q", "get_Q", "Q", "Q_0"),
)
# Normalized magnitude below which drift is absolute, e.g. Q of an equatorial orbit.
_ABSOLUTE_DRIFT_BELOW = 1e-12
# Relative spread below which tick labels cannot resolve the values themselves.
_RESOLVED_SPREAD = 1e-4


def plot_constants(
    solution: Solution, *, time: str = "t", drift: bool = False,
    style: Style | None = None, **constants: bool,
) -> tuple[Any, np.ndarray]:
    """Plot the stored energy, axial angular momentum and Carter constant.

    Parameters
    ----------
    solution : relatipy.geodesic.Solution
        Solution produced by :meth:`relatipy.geodesic.Orbit.solve`; no
        interpolation or integration is done.
    time : {"t", "tau"}, optional
        Horizontal coordinate. Default is coordinate time ``"t"``.
    drift : bool, optional
        If ``True``, plot ``(X - X_0) / |X_0|`` for each constant ``X``
        relative to the first stored sample, or ``X - X_0`` when
        ``|X_0| < 1e-12`` (for example ``Q`` of an equatorial orbit). Default
        ``False`` plots the values themselves.
    style : Style or None, optional
        Figure format; ``None`` uses the default RelatiPy style.
    **constants : bool
        Switches named ``E``, ``Lz`` and ``Q``, all enabled by default.
        Explicit ``False`` removes a panel.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure without an implicit show or save operation.
    axes : numpy.ndarray
        One-dimensional axes array, one panel per constant in the order
        ``E``, ``Lz``, ``Q``.

    Raises
    ------
    ValueError
        If ``time`` or a constant name is unknown, all constants are
        disabled, stored times are not finite, or the solution carries no
        native Kerr samples.
    TypeError
        If ``drift`` or a switch is not a boolean, or ``style`` is invalid.
    ImportError
        If Matplotlib is unavailable.

    Notes
    -----
    Values use ``G = c = M = mu = 1``, as in :meth:`~relatipy.geodesic.Solution.get_E`,
    :meth:`~relatipy.geodesic.Solution.get_Lz` and :meth:`~relatipy.geodesic.Solution.get_Q`. Without ``drift``,
    a series whose spread is below ``1e-4`` of its magnitude is plotted as
    ``X - X_0``, with ``X_0`` (the first sample) printed inside the panel,
    because tick labels could not resolve the values. Time retains its
    stored unit. Lines join stored samples only, so the panels show the
    numerical drift of the integration.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from astropy import units as u
    >>> from astropy.constants import c
    >>> from relatipy import Kerr
    >>> from relatipy.plotting import plot_constants
    >>> bh = Kerr(mass=1 * u.Msun, spin=0.5)
    >>> solution = bh.orbit(x=12 * bh.r_g, vy=0.1 * c).solve(
    ...     tau_span=(0 * u.s, 2e-8 * u.s), method="dp45")
    >>> fig, axes = plot_constants(solution, time="tau", drift=True, Q=False)
    >>> axes.shape
    (2,)
    >>> plt.close(fig)
    """
    names = tuple(entry[0] for entry in _CONSTANTS)
    unknown = constants.keys() - set(names)
    if unknown:
        raise ValueError(f"unknown constants: {', '.join(sorted(unknown))}")
    for name, enabled in (("drift", drift), *constants.items()):
        if not isinstance(enabled, bool):
            raise TypeError(f"{name} must be a bool")
    selected = [entry for entry in _CONSTANTS if constants.get(entry[0], True)]
    if not selected:
        raise ValueError("at least one constant must be enabled")
    style = resolve_style(style)
    horizontal = _horizontal(solution, time)
    series = []
    notes = []
    for _, method, latex, initial in selected:
        values = np.asarray(getattr(solution, method)(), dtype=float)
        label = f"${latex}$"
        spread = float(np.ptp(values))
        if not drift and 0.0 < spread < _RESOLVED_SPREAD * float(np.max(np.abs(values))):
            notes.append(rf"${initial} \approx {values[0]:.12g}$")
            values = values - values[0]
            label = f"${latex} - {initial}$"
        else:
            notes.append(None)
        if drift:
            reference = values[0]
            if abs(reference) < _ABSOLUTE_DRIFT_BELOW:
                values = values - reference
                label = rf"$\Delta {latex}$"
            else:
                values = (values - reference) / abs(reference)
                label = rf"$\Delta {latex} / |{initial}|$"
        series.append((label, values))
    # Taller panels leave room for the power-of-ten factors in the labels.
    fig, axes = _draw_panels(solution, time, horizontal, series, style, panel_mm=40.0)
    for ax, note in zip(axes, notes):
        if note is not None:
            ax.text(0.03, 0.92, note, transform=ax.transAxes, ha="left", va="top",
                    fontsize=style.target.tick_size, family=style.text_font,
                    color=style.frame.color, zorder=5,
                    bbox={"facecolor": style.background, "edgecolor": "none",
                          "alpha": 0.85, "pad": 1.5})
    return fig, axes
