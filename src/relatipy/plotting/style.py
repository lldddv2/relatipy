"""Publication format shared by every RelatiPy figure.

A :class:`Style` holds everything that controls appearance: the document the
figure goes into (:class:`Target`), fonts, frame and ticks, line and marker
looks by role, legends and the settings of each figure family. Sizes are in
points or millimetres at the final printed size, so figures are designed at
the size they are printed and are never rescaled in LaTeX.

The same format is available for figures that RelatiPy does not draw::

    from relatipy.plotting import new_figure, save_figure

    fig, ax = new_figure()            # thesis column width, styled axes
    ax.plot([0, 1, 2], [0, 1, 4])
    ax.set_xlabel("$t$ [yr]")
    save_figure(fig, "figure")        # writes figure.pdf and figure.png

Matplotlib is imported only when a function that draws is called.
"""

from __future__ import annotations

import contextlib
import logging
import subprocess
from dataclasses import dataclass, field, replace
from functools import lru_cache
from pathlib import Path
from types import MappingProxyType
from typing import Any, Iterator, Mapping

MM = 1.0 / 25.4
"""Inches per millimetre."""

MINUS = "\N{MINUS SIGN}"


def _frozen(mapping: Mapping | None) -> Mapping | None:
    """Read-only copy of a mapping, with nested mappings also read-only."""
    if mapping is None:
        return None
    return MappingProxyType({key: _frozen(value) if isinstance(value, Mapping) else value
                             for key, value in mapping.items()})


@dataclass(frozen=True)
class Target:
    """Document a figure is printed in.

    Parameters
    ----------
    name : str
        Short identifier.
    width_mm, full_width_mm : float
        Width of a single figure and of a full-width figure, in millimetres.
    font_file : str
        TeX font file used when ``Style.font`` is ``None`` and TeX is found.
    font_fallback : tuple of str
        Families tried, in order, when the font file is not available.
    math_font : str
        Matplotlib math font set for this document.
    label_size, tick_size, small_size : float
        Axis names, tick labels, and legends and annotations, in points.
    """

    name: str
    width_mm: float
    full_width_mm: float
    font_file: str
    font_fallback: tuple[str, ...]
    math_font: str
    label_size: float
    tick_size: float
    small_size: float


TARGETS: Mapping[str, Target] = MappingProxyType({
    # Thesis: book class, 12 pt, A4 with 2.5 cm margins (text width 160 mm);
    # a single figure is 0.675 of the text width.
    "thesis": Target("thesis", 108.0, 160.0, "lmroman10-regular.otf",
                     ("Latin Modern Roman", "CMU Serif", "DejaVu Serif"), "cm",
                     10.0, 9.0, 8.0),
    # A&A: one column is 88 mm and two columns 180 mm. The lettering size
    # has not been confirmed from the official author guide.
    "aanda": Target("aanda", 88.0, 180.0, "texgyretermes-regular.otf",
                    ("TeX Gyre Termes", "Nimbus Roman", "Times New Roman",
                     "DejaVu Serif"), "stix", 8.0, 7.0, 6.0),
})
"""Known documents, by name."""


@dataclass(frozen=True)
class Frame:
    """Axes frame and ticks of two-dimensional panels.

    Parameters
    ----------
    color : str, optional
        Spine and tick colour.
    width : float, optional
        Spine width in points.
    tick_direction : {"in", "out", "inout"}, optional
        Direction passed to Matplotlib for major and minor ticks.
    tick_length, tick_width, minor_tick_length : float, optional
        Major-tick length, major-tick width, and minor-tick length in points.
    mirror_ticks, minor_ticks : bool, optional
        Whether to draw ticks on the opposite spines and enable minor ticks.
    """

    color: str = "black"
    width: float = 0.6
    tick_direction: str = "in"
    tick_length: float = 3.5
    tick_width: float = 0.6
    mirror_ticks: bool = True
    minor_ticks: bool = True
    minor_tick_length: float = 1.8


@dataclass(frozen=True)
class Views:
    """Layout settings for :func:`plot_views`.

    Parameters
    ----------
    width_factor : float, optional
        Figure width relative to the target's single-figure width.
    gap_mm : float, optional
        Gap between panels in millimetres.
    padding : float, optional
        Fraction of the largest data span added around each projection.
    """

    width_factor: float = 1.1        # width relative to a single figure
    gap_mm: float = 0.0              # between panels; 0 makes them touch
    padding: float = 0.06            # data margin, fraction of the largest span


@dataclass(frozen=True)
class Static3D:
    """Settings for static Matplotlib three-dimensional figures.

    Parameters
    ----------
    elevation, azimuth : float, optional
        Initial camera angles in degrees.
    pane_color, pane_edge_color : str, optional
        Face and edge colours of the coordinate panes.
    pane_edge_width : float, optional
        Pane-edge width in points.
    grid : mapping or None, optional
        Matplotlib grid keyword arguments; ``None`` suppresses the grid.
    tick_inward, tick_outward : float, optional
        Fractions of an axis used for inward and outward three-dimensional
        ticks.
    ticks : int, optional
        Maximum number of tick intervals per axis.
    label_pad : float, optional
        Distance between a label and its axis in points.
    depth_sorting : bool, optional
        Whether Matplotlib depth sorting may place references over a path.
    padding : float, optional
        Fraction of the data span added around the scene.
    height_factor : float, optional
        Initial figure height relative to its width.
    """

    elevation: float = 25.0
    azimuth: float = -60.0
    pane_color: str = "white"
    pane_edge_color: str = "#8C8C8C"
    pane_edge_width: float = 0.5
    grid: Mapping | None = field(default_factory=lambda: _frozen(
        {"color": "#E3E3E3", "linewidth": 0.4}))
    tick_inward: float = 0.25        # fraction of the axis used by Matplotlib
    tick_outward: float = 0.0
    ticks: int = 4                   # at most this many tick intervals
    label_pad: float = 5.0
    depth_sorting: bool = False      # False: trajectory above the references
    padding: float = 0.04
    height_factor: float = 0.9       # first layout only; the page is trimmed


@dataclass(frozen=True)
class Interactive:
    """Settings for interactive Plotly figures.

    Parameters
    ----------
    height : int, optional
        Figure height in CSS pixels.
    font_size, tick_size, small_size : float, optional
        Main, tick-label, and auxiliary font sizes in CSS pixels.
    trajectory_width, reference_width : float, optional
        Line widths in CSS pixels.
    marker_sizes : mapping, optional
        Marker sizes keyed by ``"initial"``, ``"end"``, and ``"current"``.
    equatorial_plane, horizon_surface, ergosurface : mapping or None, optional
        Trace styling dictionaries; ``None`` omits that reference.
    ergo_gap : float, optional
        Fraction of the horizon radius omitted near the ergosurface poles.
    surface_resolution : int, optional
        Number of angular samples used for reference surfaces.
    isco_cone, end_cone : float, optional
        Arrow lengths as fractions of the plotted span.
    bh_size : float, optional
        Black-hole marker size.
    grid_color, line_color : str, optional
        Grid and axis-line colours.
    ticks : int, optional
        Requested number of axis tick intervals.
    camera_distance : float, optional
        Initial Plotly camera distance from the origin.
    """

    height: int = 560
    font_size: float = 13
    tick_size: float = 12
    small_size: float = 12
    trajectory_width: float = 2.5
    reference_width: float = 2.0
    marker_sizes: Mapping = field(default_factory=lambda: _frozen(
        {"initial": 3, "end": 1, "current": 5}))
    equatorial_plane: Mapping | None = field(default_factory=lambda: _frozen(
        {"color": "#9E9E9E", "opacity": 0.12, "padding": 0.04}))
    horizon_surface: Mapping | None = field(default_factory=lambda: _frozen(
        {"color": "#222222", "opacity": 1.0}))
    ergosurface: Mapping | None = field(default_factory=lambda: _frozen(
        {"color": "#44AA99", "opacity": 0.35}))
    # The ergosurface meets the horizon at the poles; closer than this
    # fraction of the horizon radius it is not drawn (depth fighting).
    ergo_gap: float = 0.04
    surface_resolution: int = 72
    isco_cone: float = 0.025         # arrow length, fraction of the span
    end_cone: float = 0.045
    bh_size: float = 2
    grid_color: str = "#E3E3E3"
    line_color: str = "black"
    ticks: int = 5
    camera_distance: float = 2.2


@dataclass(frozen=True)
class Observables:
    """Appearance settings for :func:`plot_observables`.

    Parameters
    ----------
    data_color, model_color, residual_color : str, optional
        Colours of observed data, model samples, and residuals.
    marker_size, error_width, model_width : float, optional
        Data-marker size and error-bar/model line widths in points.
    zero_line : mapping or None, optional
        Matplotlib line keyword arguments for the residual zero line.
    sky_width_ratio, height_ratio : float, optional
        Relative widths and height of the panel arrangement.
    column_gap_mm, legend_gap_mm : float, optional
        Inter-column and legend gaps in millimetres.
    names : mapping, optional
        Labels keyed by data, model, residual, and centre roles.
    """

    data_color: str = "black"
    marker_size: float = 1.6
    error_width: float = 0.5
    model_color: str = "#7F7F7F"
    model_width: float = 0.8
    residual_color: str = "#CC3311"
    zero_line: Mapping | None = field(default_factory=lambda: _frozen(
        {"color": "#CC3311", "linewidth": 0.6}))
    sky_width_ratio: float = 0.9
    height_ratio: float = 0.8
    column_gap_mm: float = 3.0
    legend_gap_mm: float = 1.5
    names: Mapping = field(default_factory=lambda: _frozen(
        {"data": "Data", "model": "Model", "residual": "Data − model",
         "center": "Central mass"}))


@dataclass(frozen=True)
class Corner:
    """Appearance settings for :func:`plot_corner`.

    Parameters
    ----------
    size_in : tuple of float or None, optional
        Explicit ``(width, height)`` in inches; ``None`` uses a square,
        full-width figure.
    tick_size, label_size : float, optional
        Tick and largest axis-label sizes in points.
    text_color : str, optional
        Colour of annotations.
    frame_width : float, optional
        Width of subplot frames in points.
    minor_ticks : bool, optional
        Whether to show minor ticks.
    bins, hist_bins : int, optional
        Numbers of bins for two-dimensional contours and marginal histograms.
    smooth : float, optional
        Gaussian smoothing standard deviation in two-dimensional bins.
    levels : tuple of float, optional
        Enclosed posterior probabilities drawn as contours.
    fills : tuple of str, optional
        Fill colours ordered from inner to outer contour.
    outline, hist, median_line, quantile_line : mapping, optional
        Matplotlib keyword arguments for contour outlines, histograms, and
        summary lines.
    margin : float, optional
        Fractional margin around plotted samples.
    ticks : int, optional
        Target number of tick intervals on each panel.
    tick_inset : float, optional
        Fraction of a panel reserved between ticks and its edge.
    """

    size_in: tuple[float, float] | None = (8.0, 9.0)  # None: full width, square
    tick_size: float = 7
    label_size: float = 8            # largest; reduced to fit the panel side
    text_color: str = "#666666"
    frame_width: float = 1.0
    minor_ticks: bool = True
    bins: int = 40
    smooth: float = 1.0              # Gaussian sigma, in 2D bins
    hist_bins: int = 60
    levels: tuple[float, ...] = (0.683, 0.954, 0.997)
    fills: tuple[str, ...] = ("black", "#BFBFBF", "#F2F2F2")
    outline: Mapping = field(default_factory=lambda: _frozen(
        {"color": "#7F7F7F", "linewidth": 0.4}))
    margin: float = 0.25
    hist: Mapping = field(default_factory=lambda: _frozen(
        {"color": "#7F7F7F", "edgecolor": "#3F3F3F", "linewidth": 0.3}))
    median_line: Mapping = field(default_factory=lambda: _frozen(
        {"color": "black", "linewidth": 1.0, "linestyle": (0, (3, 1.5))}))
    quantile_line: Mapping = field(default_factory=lambda: _frozen(
        {"color": "black", "linewidth": 0.5, "linestyle": (0, (3, 1.5))}))
    ticks: int = 3
    tick_inset: float = 0.08


# Paul Tol "muted" green and purple for the prograde and retrograde ISCO:
# distinguishable with colour-vision deficiencies, not a red/blue pair, at
# least 3:1 against white and 11.7 L* apart in grayscale. Line styles and
# names repeat the cue.
_LINES = {
    "trajectory": {"color": "black", "linewidth": 0.8},
    "preview": {"color": "black", "linewidth": 0.8},
    "horizon": {"color": "black", "linewidth": 0.6},
    "isco_prograde": {"color": "#117733", "linewidth": 0.6, "linestyle": "--"},
    "isco_retrograde": {"color": "#882255", "linewidth": 0.6, "linestyle": "--"},
}
# "shape": "dot", "chevron" (open ">") or "triangle"; directional shapes
# point along the last step. "size" is the area in points². color None
# reuses the trajectory color.
_MARKERS = {
    "initial": {"shape": "dot", "color": None, "size": 1.5},
    "end": {"shape": "chevron", "color": None, "size": 22, "linewidth": 0.8},
    "current": {"shape": "dot", "color": "black", "size": 8},
}
# Colours of overlaid trajectories (plot_sols), in order: black, then Paul
# Tol "muted" hues other than the ISCO green and purple.
_TRAJECTORY_COLORS = ("black", "#332288", "#CC6677", "#44AA99", "#999933",
                      "#88CCEE", "#DDCC77", "#AA4499")
_NAMES = {
    "trajectory": "Trajectory", "preview": "Osculating Kepler\npreview",
    "horizon": "Outer horizon", "isco_prograde": "Prograde ISCO",
    "isco_retrograde": "Retrograde ISCO", "initial": "Initial position",
    "end": "End position", "current": "Current position", "ergosurface": "Ergosphere",
}


@dataclass(frozen=True)
class Style:
    """Complete appearance of RelatiPy figures.

    Build variants with :meth:`Style.replace`, e.g.
    ``DEFAULT_STYLE.replace(target="aanda")``. Mappings are read-only.
    Every field is optional; the defaults below give :data:`DEFAULT_STYLE`.

    Parameters
    ----------
    target : Target or str, optional
        Document, or the name of one in :data:`TARGETS`. Default
        ``"thesis"``.
    font : str or None, optional
        Text family, default ``"DejaVu Serif"``; ``None`` uses the
        target's font.
    math_font : str or None, optional
        Matplotlib math font set, default ``"dejavuserif"``; ``None`` uses
        the target's.
    background : str, optional
        Matplotlib figure and axes face colour, default ``"white"``.
    frame : Frame, optional
        Two-dimensional axes and tick settings, default ``Frame()``.
    outer_margin_mm : float, optional
        Minimum figure-edge margin in millimetres, default ``2.0``.
    length_unit : {None, "r_g"}, optional
        Orbit lengths in the stored data unit by default (``None``). Use
        ``"r_g"`` to display gravitational radii ``GM/c²`` explicitly.
    show_title : bool, optional
        Whether supported figures show their supplied title, default
        ``False``.
    title_color : str, optional
        Colour used for displayed titles, default ``"#555555"``.
    lines, markers : mapping, optional
        Looks by role (see the module constants for the keys). Defaults are
        the module's role tables.
    trajectory_colors : tuple of str, optional
        Colours of overlaid trajectories in :func:`plot_sols`, used in
        order and repeated when there are more trajectories. Default:
        black, then colour-blind-safe hues distinct from the ISCOs.
    names : mapping, optional
        Legend names by role. Default is the module's role-name table.
    isco_arrows, center_marker : mapping or None, optional
        Drawing settings for ISCO arrows and the central-mass marker;
        ``None`` omits the corresponding decoration. Both are set by
        default.
    isco_labels : bool, optional
        Whether two-dimensional views write ISCO names beside their circles,
        default ``True``.
    isco_label_gap : float, optional
        Separation of those ISCO labels from their circles in points,
        default ``1.5``.
    bh_name : str or None, optional
        Label for the central mass, default ``"BH"``; ``None`` omits it.
    bh_name_offset : tuple of float, optional
        Offset of the central-mass label in points, default ``(3.0, 3.0)``.
    spin_text : str or None, optional
        Format string for the spin note, default ``"spin = {spin}"``, or
        ``None`` to omit it.
    legend_frame : bool, optional
        Whether legends draw a frame, default ``False``.
    legend_hide : tuple of str, optional
        Roles never listed in a legend, default ``("initial", "end")``.
    legend_min_entries : int, optional
        A legend with fewer entries is dropped; the caption names the curve.
        Default ``2``.
    views : Views, optional
        Orthogonal-projection layout settings, default ``Views()``.
    static_3d : Static3D, optional
        Static three-dimensional rendering settings, default ``Static3D()``.
    interactive : Interactive, optional
        Plotly rendering settings, default ``Interactive()``.
    observables : Observables, optional
        Astrometric-observable plot settings, default ``Observables()``.
    corner : Corner, optional
        Posterior corner-plot settings, default ``Corner()``.

    Raises
    ------
    ValueError
        If ``target`` is a name not in :data:`TARGETS`, ``length_unit``
        is neither ``"r_g"`` nor ``None``, or ``trajectory_colors`` is not
        a nonempty sequence of colour strings.
    TypeError
        If ``target`` is neither a :class:`Target` nor a string.

    Examples
    --------
    >>> from relatipy.plotting import DEFAULT_STYLE, Style
    >>> style = DEFAULT_STYLE.replace(length_unit="r_g")
    >>> style.length_unit
    'r_g'
    >>> Style(target="thesis").target.name
    'thesis'
    """

    target: Target | str = "thesis"
    font: str | None = "DejaVu Serif"
    math_font: str | None = "dejavuserif"
    background: str = "white"
    frame: Frame = Frame()
    outer_margin_mm: float = 2.0
    length_unit: str | None = None
    show_title: bool = False
    title_color: str = "#555555"
    lines: Mapping = field(default_factory=lambda: _frozen(_LINES))
    markers: Mapping = field(default_factory=lambda: _frozen(_MARKERS))
    trajectory_colors: tuple[str, ...] = _TRAJECTORY_COLORS
    names: Mapping = field(default_factory=lambda: _frozen(_NAMES))
    isco_arrows: Mapping | None = field(default_factory=lambda: _frozen(
        {"count": 4, "size": 14, "linewidth": 0.6}))
    isco_labels: bool = True         # names along the circles in xy views
    isco_label_gap: float = 1.5      # points outside the circle
    center_marker: Mapping | None = field(default_factory=lambda: _frozen(
        {"marker": "x", "color": "black", "size": 18, "linewidth": 0.9}))
    bh_name: str | None = "BH"
    bh_name_offset: tuple[float, float] = (3.0, 3.0)
    spin_text: str | None = "spin = {spin}"
    legend_frame: bool = False
    legend_hide: tuple[str, ...] = ("initial", "end")
    legend_min_entries: int = 2
    views: Views = Views()
    static_3d: Static3D = Static3D()
    interactive: Interactive = Interactive()
    observables: Observables = Observables()
    corner: Corner = Corner()

    def __post_init__(self) -> None:
        """Resolve the target and freeze mutable mappings after validation."""
        target = self.target
        if isinstance(target, str):
            if target not in TARGETS:
                raise ValueError(f"unknown target {target!r}; known: {sorted(TARGETS)}")
            object.__setattr__(self, "target", TARGETS[target])
        elif not isinstance(target, Target):
            raise TypeError("target must be a Target or the name of one")
        if self.length_unit not in ("r_g", None):
            raise ValueError("length_unit must be 'r_g' or None")
        colors = self.trajectory_colors
        colors = () if isinstance(colors, str) else tuple(colors)
        if not colors or not all(isinstance(color, str) for color in colors):
            raise ValueError("trajectory_colors must be a nonempty sequence of colours")
        object.__setattr__(self, "trajectory_colors", colors)
        for name in ("lines", "markers", "names", "isco_arrows", "center_marker"):
            object.__setattr__(self, name, _frozen(getattr(self, name)))

    def replace(self, **changes: Any) -> Style:
        """Return a style copy with selected dataclass fields changed.

        Parameters
        ----------
        **changes
            Field names and replacement values accepted by :class:`Style`.

        Returns
        -------
        Style
            A newly validated immutable style.

        Raises
        ------
        TypeError
            If a keyword is not a :class:`Style` field.
        ValueError
            If a replacement violates :class:`Style` validation.
        """
        return replace(self, **changes)

    @property
    def text_font(self) -> str:
        """str: Resolved Matplotlib family used for text."""
        return self.font or resolve_font(self.target)

    @property
    def math_fontset(self) -> str:
        """str: Resolved Matplotlib math-text font set."""
        return self.math_font or self.target.math_font

    def name(self, role: str) -> str:
        """Return the configured legend label for a drawing role.

        Parameters
        ----------
        role : str
            Role key, such as ``"trajectory"`` or ``"horizon"``.

        Returns
        -------
        str
            Configured label, or ``role`` when it has no configured label.
        """
        return self.names.get(role, role)


DEFAULT_STYLE = Style()
"""Style used when a function receives ``style=None``.

Derive variants with :meth:`Style.replace`; the default itself is immutable.
"""


def resolve_style(style: Style | None) -> Style:
    """Resolve an optional style to an immutable :class:`Style`.

    Parameters
    ----------
    style : Style or None
        Style to use; ``None`` selects :data:`DEFAULT_STYLE`.

    Returns
    -------
    Style
        The supplied style or the default style.

    Raises
    ------
    TypeError
        If ``style`` is neither :class:`Style` nor ``None``.
    """
    if style is None:
        return DEFAULT_STYLE
    if not isinstance(style, Style):
        raise TypeError("style must be a relatipy.plotting.Style")
    return style


@lru_cache(maxsize=None)
def resolve_font(target: Target) -> str:
    """Resolve a usable Matplotlib text family for a publication target.

    The target's TeX font is registered when ``kpsewhich`` locates it. The
    configured fallback families are tested in order otherwise.

    Parameters
    ----------
    target : Target
        Document target that supplies the preferred font and fallbacks.

    Returns
    -------
    str
        A registered font family, a configured fallback, or ``"serif"``.

    Raises
    ------
    ImportError
        If Matplotlib is unavailable.
    """
    from matplotlib import font_manager

    try:
        path = subprocess.run(["kpsewhich", target.font_file], capture_output=True,
                              text=True, check=False, timeout=10).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        path = ""
    if path:
        # Some TeX font headers carry old timestamps; fontTools warns.
        logging.getLogger("fontTools").setLevel(logging.ERROR)
        font_manager.fontManager.addfont(path)
        return font_manager.FontProperties(fname=path).get_name()
    for name in target.font_fallback:
        try:
            font_manager.findfont(font_manager.FontProperties(family=name),
                                  fallback_to_default=False)
        except ValueError:
            continue
        return name
    return "serif"


def number(value: float) -> str:
    """Format a numeric tick value with a typographic minus sign.

    Parameters
    ----------
    value : float
        Value to format using Python's general numeric format.

    Returns
    -------
    str
        Compact decimal representation whose negative sign is ``U+2212``.
    """
    text = f"{value:g}"
    return text.replace("-", MINUS)


def rc_params(style: Style | None = None) -> dict[str, Any]:
    """Return Matplotlib rcParams matching a RelatiPy style.

    Parameters
    ----------
    style : Style or None, optional
        Source format. ``None`` uses :data:`DEFAULT_STYLE`.

    Returns
    -------
    dict of str to object
        Values suitable for :func:`matplotlib.rc_context` or ``rcParams``.

    Raises
    ------
    TypeError
        If ``style`` is neither :class:`Style` nor ``None``.
    """
    style = resolve_style(style)
    target, frame = style.target, style.frame
    return {
        "font.family": style.text_font,
        "mathtext.fontset": style.math_fontset,
        "font.size": target.small_size,
        "axes.labelsize": target.label_size,
        "axes.titlesize": target.small_size,
        "xtick.labelsize": target.tick_size,
        "ytick.labelsize": target.tick_size,
        "legend.fontsize": target.small_size,
        "legend.frameon": style.legend_frame,
        "axes.linewidth": frame.width,
        "axes.edgecolor": frame.color,
        "axes.unicode_minus": True,
        "axes.facecolor": style.background,
        "figure.facecolor": style.background,
        "xtick.direction": frame.tick_direction,
        "ytick.direction": frame.tick_direction,
        "xtick.top": frame.mirror_ticks,
        "ytick.right": frame.mirror_ticks,
        "xtick.major.size": frame.tick_length,
        "ytick.major.size": frame.tick_length,
        "xtick.major.width": frame.tick_width,
        "ytick.major.width": frame.tick_width,
        "xtick.minor.visible": frame.minor_ticks,
        "ytick.minor.visible": frame.minor_ticks,
        "xtick.minor.size": frame.minor_tick_length,
        "ytick.minor.size": frame.minor_tick_length,
        "xtick.minor.width": 0.8 * frame.tick_width,
        "ytick.minor.width": 0.8 * frame.tick_width,
        "lines.linewidth": style.lines["trajectory"]["linewidth"],
        "savefig.dpi": 300,
        # TrueType fonts keep PDF text editable and searchable.
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }


@contextlib.contextmanager
def publication_style(style: Style | None = None) -> Iterator[Style]:
    """Temporarily apply a RelatiPy Matplotlib publication format.

    Parameters
    ----------
    style : Style or None, optional
        Format to install for the context. ``None`` uses :data:`DEFAULT_STYLE`.

    Yields
    ------
    Style
        The resolved immutable style.

    Raises
    ------
    TypeError
        If ``style`` is neither :class:`Style` nor ``None``.
    ImportError
        If Matplotlib is unavailable.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from relatipy.plotting import MM, publication_style
    >>> with publication_style() as style:
    ...     fig, ax = plt.subplots(figsize=(style.target.width_mm * MM, 2))
    >>> plt.close(fig)
    """
    import matplotlib

    style = resolve_style(style)
    # The first pyplot figure of a session may select the backend inside this
    # context; inline backends then enable interactive mode, which rc_context
    # would undo on exit and later figures would not display. Keep that state.
    state = {}
    try:
        with matplotlib.rc_context(rc_params(style)):
            try:
                yield style
            finally:
                state["interactive"] = matplotlib.is_interactive()
    finally:
        if "interactive" in state:
            matplotlib.interactive(state["interactive"])


def format_axes(ax: Any, style: Style | None = None) -> Any:
    """Apply a RelatiPy frame and tick format to existing 2D axes.

    Tick labels use a typographic minus.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Existing two-dimensional axes to modify in place.
    style : Style or None, optional
        Format to apply. ``None`` uses :data:`DEFAULT_STYLE`.

    Returns
    -------
    matplotlib.axes.Axes
        The same ``ax`` after formatting.

    Raises
    ------
    TypeError
        If ``style`` is neither :class:`Style` nor ``None``.
    ImportError
        If Matplotlib is unavailable.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from relatipy.plotting import format_axes
    >>> fig, ax = plt.subplots()
    >>> format_axes(ax) is ax
    True
    >>> plt.close(fig)
    """
    from matplotlib.ticker import FuncFormatter

    style = resolve_style(style)
    target, frame = style.target, style.frame
    ax.set_facecolor(style.background)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color(frame.color)
        spine.set_linewidth(frame.width)
    common = {"direction": frame.tick_direction, "width": frame.tick_width,
              "colors": frame.color, "top": frame.mirror_ticks,
              "right": frame.mirror_ticks}
    ax.tick_params(which="major", length=frame.tick_length, labelsize=target.tick_size,
                   labelfontfamily=style.text_font, **common)
    if frame.minor_ticks:
        ax.minorticks_on()
        ax.tick_params(which="minor", length=frame.minor_tick_length,
                       **{**common, "width": 0.8 * frame.tick_width})
    else:
        ax.minorticks_off()
    for axis in (ax.xaxis, ax.yaxis):
        if axis.get_scale() == "linear":
            axis.set_major_formatter(FuncFormatter(lambda value, _: number(value)))
    return ax


def set_axis_label(ax: Any, which: str, text: str, style: Style | None = None) -> None:
    """Set one two-dimensional axis label using the style's fonts.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes whose label is set in place.
    which : {"x", "y"}
        Axis to label.
    text : str
        Label text. Matplotlib math markup is rendered in the style's math
        font set.
    style : Style or None, optional
        Format to apply. ``None`` uses :data:`DEFAULT_STYLE`.

    Raises
    ------
    ValueError
        If ``which`` is not ``"x"`` or ``"y"``.
    TypeError
        If ``style`` is neither :class:`Style` nor ``None``.
    """
    style = resolve_style(style)
    if which not in ("x", "y"):
        raise ValueError("which must be 'x' or 'y'")
    axis = ax.xaxis if which == "x" else ax.yaxis
    axis.set_label_text(text, fontsize=style.target.label_size,
                        fontfamily=style.text_font)
    axis.label.set_math_fontfamily(style.math_fontset)


def legend_properties(style: Style | None = None, size: float | None = None) -> dict:
    """Build Matplotlib ``Axes.legend`` keyword arguments for a style.

    Parameters
    ----------
    style : Style or None, optional
        Source format. ``None`` uses :data:`DEFAULT_STYLE`.
    size : float or None, optional
        Legend text size in points. ``None`` uses the target's auxiliary size.

    Returns
    -------
    dict
        ``frameon`` and font properties suitable for ``Axes.legend``.

    Raises
    ------
    TypeError
        If ``style`` is neither :class:`Style` nor ``None``.
    """
    style = resolve_style(style)
    return {"frameon": style.legend_frame,
            "prop": {"family": style.text_font, "size": size or style.target.small_size}}


def figure_width(width: str | float, style: Style | None = None) -> float:
    """Return a publication width in inches.

    Parameters
    ----------
    width : {"single", "full"} or float
        Named target width, or a positive width in millimetres.
    style : Style or None, optional
        Target dimensions. ``None`` uses :data:`DEFAULT_STYLE`.

    Returns
    -------
    float
        Width in inches.

    Raises
    ------
    ValueError
        If ``width`` is neither a supported name nor a positive real value.
    TypeError
        If ``style`` is neither :class:`Style` nor ``None``.
    """
    style = resolve_style(style)
    if width == "single":
        return style.target.width_mm * MM
    if width == "full":
        return style.target.full_width_mm * MM
    if isinstance(width, bool) or not isinstance(width, (int, float)) or width <= 0:
        raise ValueError("width must be 'single', 'full' or a positive width in mm")
    return float(width) * MM


def new_figure(width: str | float = "single", height: float | None = None,
               nrows: int = 1, ncols: int = 1, style: Style | None = None,
               **subplots: Any) -> tuple[Any, Any]:
    """Create a figure at printed size with axes in the publication format.

    Parameters
    ----------
    width : {"single", "full"} or float, optional
        Target figure width, or a width in millimetres.
    height : float or None, optional
        Height in millimetres; ``None`` uses 3/4 of the width.
    nrows, ncols : int, optional
        Grid of axes, as in ``plt.subplots``.
    style : Style or None, optional
        Format; ``None`` uses :data:`DEFAULT_STYLE`.
    **subplots
        Passed to ``plt.subplots`` (e.g. ``sharex=True``).

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure, with constrained layout and the style's outer margin.
    axes : matplotlib.axes.Axes or numpy.ndarray
        The axes returned by ``plt.subplots``, formatted with
        :func:`format_axes`: one axes object for a 1 x 1 grid, otherwise an
        array of axes.

    Raises
    ------
    ValueError
        If ``width`` is unsupported or non-positive.
    TypeError
        If ``style`` is invalid, or Matplotlib rejects ``**subplots``.
    ImportError
        If Matplotlib is unavailable.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from relatipy.plotting import new_figure
    >>> fig, axes = new_figure("full", height=40, ncols=2, sharey=True)
    >>> axes.shape
    (2,)
    >>> plt.close(fig)
    """
    import matplotlib.pyplot as plt
    import numpy as np

    style = resolve_style(style)
    inches = figure_width(width, style)
    tall = 0.75 * inches if height is None else float(height) * MM
    with publication_style(style):
        fig, axes = plt.subplots(nrows, ncols, figsize=(inches, tall),
                                 layout="constrained", **subplots)
    margin = style.outer_margin_mm * MM
    fig.get_layout_engine().set(w_pad=margin, h_pad=margin)
    fig.set_facecolor(style.background)
    for ax in np.ravel(axes):
        format_axes(ax, style)
    return fig, axes


def save_figure(fig: Any, path: str | Path, formats: tuple[str, ...] = ("pdf", "png"),
                dpi: int = 300, **savefig: Any) -> list[Path]:
    """Save a Matplotlib figure in one or more formats.

    PDF text is embedded as TrueType (Type 42) and stays editable. ``path``
    without suffix gets one file per format.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure whose ``savefig`` method is called.
    path : str or pathlib.Path
        Base path. A suffix overrides ``formats`` and selects one output.
    formats : tuple of str, optional
        Extensions used when ``path`` has no suffix, default
        ``("pdf", "png")``.
    dpi : int, optional
        Resolution passed to Matplotlib, default 300.
    **savefig
        Additional keyword arguments forwarded unchanged to ``fig.savefig``.

    Returns
    -------
    list of pathlib.Path
        Paths written in format order.

    Raises
    ------
    OSError
        If Matplotlib cannot write an output path.
    ImportError
        If Matplotlib is unavailable.

    Examples
    --------
    >>> import tempfile
    >>> from pathlib import Path
    >>> import matplotlib.pyplot as plt
    >>> from relatipy.plotting import new_figure, save_figure
    >>> fig, ax = new_figure()
    >>> with tempfile.TemporaryDirectory() as folder:
    ...     written = save_figure(fig, Path(folder) / "orbit")
    ...     [path.name for path in written]
    ['orbit.pdf', 'orbit.png']
    >>> plt.close(fig)
    """
    import matplotlib

    base = Path(path)
    if base.suffix:
        formats = (base.suffix.lstrip("."),)
        base = base.with_suffix("")
    written = []
    with matplotlib.rc_context({"pdf.fonttype": 42, "ps.fonttype": 42}):
        for extension in formats:
            target = base.with_suffix(f".{extension}")
            fig.savefig(target, dpi=dpi, **savefig)
            written.append(target)
    return written
