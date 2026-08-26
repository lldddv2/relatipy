"""
Matplotlib subplot styling for publication-quality scientific figures.

This module provides :class:`SciSubplot`, which applies a consistent
LaTeX-like serif look, inward ticks on all four sides, optional minor
ticks, and a thin frame—similar to figures in GRAVITY/ESO-style papers
(Gillessen et al.).

Examples
--------
>>> from relatipy.visualization._2D.sci import SciSubplot
>>> ws = SciSubplot()
>>> _ = ws.ax.plot([0, 1], [0, 1])
"""
import matplotlib.pyplot as plt


class SciSubplot:
    """
    Scientific-style subplot with GRAVITY/ESO-inspired aesthetics.

    The style uses LaTeX-like serif fonts (STIXGeneral / Computer Modern
    family stack), inward ticks on all four sides, optional automatic minor
    ticks, a thin frame, and no grid by default (common in astrophysics
    publications).

    Parameters
    ----------
    grid : bool, default False
        If True, draw a very subtle dashed grid on all axes.
    minor_ticks : bool, default True
        If True, enable automatic minor tick locators on every axis.
    subplot : tuple, list, dict, or None, default None
        If None, create a single figure and one axes via ``add_subplot(111)``.
        If a tuple or list, it is unpacked as positional arguments to
        :func:`matplotlib.pyplot.subplots` (e.g. ``(2, 2)`` for a 2×2 grid).
        If a dict, keys are passed as keyword arguments to ``subplots``.
    *args
        Extra positional arguments for :func:`matplotlib.pyplot.figure`
        (when ``subplot`` is None) or :func:`matplotlib.pyplot.subplots`.
    **kwargs
        Extra keyword arguments for ``figure`` or ``subplots``.

    Attributes
    ----------
    fig : matplotlib.figure.Figure
        The matplotlib figure instance.
    ax : matplotlib.axes.Axes
        Primary axes; the first flattened entry when ``axs`` is an array.
    axs : matplotlib.axes.Axes or numpy.ndarray
        Single axes or array of axes from ``subplots``, depending on layout.
    COLORS : dict
        Named colors aligned with GRAVITY/SINFONI-style palettes (keys:
        ``"black"``, ``"cyan"``, ``"red"``, ``"blue"``, ``"gray"``).

    Raises
    ------
    ValueError
        If ``subplot`` is neither None, a tuple/list, nor a dict.

    Examples
    --------
    Single axes (default):

    >>> ws = SciSubplot()
    >>> ws.ax.plot([0, 1], [0, 1])  # doctest: +ELLIPSIS
    [...]

    Grid of subplots:

    >>> ws = SciSubplot(subplot=(2, 2), figsize=(6, 6))
    >>> hasattr(ws.axs, "flat")
    True
    """

    # ── GRAVITY-style color palette ───────────────────────────────────────
    COLORS = {
        "black":   "black",
        "cyan":    "#4DB8FF",   # GRAVITY
        "red":     "#CC3300",   # flares / reference
        "blue":    "#1A3E8C",   # SINFONI
        "gray":    "#999999",   # fit curves
    }

    def __init__(self, grid=False, minor_ticks=True, subplot=None, *args, **kwargs):
        import matplotlib as mpl
        from matplotlib.ticker import AutoMinorLocator

        # ── LaTeX-like typography ─────────────────────────────────────────
        mpl.rcParams.update({
            "font.family":          "serif",
            "font.serif":           ["TeX Gyre Termes", "Times New Roman",
                                     "STIXGeneral", "DejaVu Serif"],
            "mathtext.fontset":     "stix",
            "font.size":            10,
            "axes.labelsize":       10,
            "xtick.labelsize":      9,
            "ytick.labelsize":      9,
            "legend.fontsize":      9,

            # ── Ticks ────────────────────────────────────────────────────
            "xtick.direction":      "in",
            "ytick.direction":      "in",
            "xtick.top":            True,
            "ytick.right":          True,
            "xtick.major.size":     3,
            "ytick.major.size":     3,
            "xtick.minor.size":     1,
            "ytick.minor.size":     1,
            "xtick.major.width":    0.6,
            "ytick.major.width":    0.6,
            "xtick.minor.width":    0.6,
            "ytick.minor.width":    0.6,

            # ── Lines / markers ───────────────────────────────────────────
            "lines.linewidth":      1.2,
            "lines.markersize":     4,
            "errorbar.capsize":     0,      # no caps on error bars

            # ── Frame and background ────────────────────────────────────
            "axes.linewidth":       0.6,
            "axes.edgecolor":       "black",
            "axes.facecolor":       "white",
            "figure.facecolor":     "white",

            # ── Legend ────────────────────────────────────────────────────
            "legend.frameon":       False,
            "legend.handlelength":  1.5,
            "legend.handletextpad": 0.5,    # space between handle and text
            "legend.labelspacing":  0.4,    # vertical space between entries

            # ── Resolution ─────────────────────────────────────────────
            "figure.dpi":           125,
            "savefig.dpi":          300,
            "savefig.bbox":         "tight",
        })

        # ── Figure and axes (subplot support) ───────────────────────────
        if subplot is None:
            # Original behavior: one figure and one axes
            self.fig = plt.figure(*args, **kwargs)
            self.ax = self.fig.add_subplot(111)
            self.axs = self.ax
        else:
            # Subplots via matplotlib API
            # Example:
            #   workspace = SciSubplot(subplot=(2, 2), figsize=(6, 6))
            #   fig, axs = workspace.fig, workspace.axs
            if isinstance(subplot, (tuple, list)):
                # subplot=(nrows, ncols, ...) + normal subplots args/kwargs
                self.fig, axs = plt.subplots(*subplot, *args, **kwargs)
            elif isinstance(subplot, dict):
                # subplot={"nrows": ..., "ncols": ..., ...} + args
                # Positional args must precede keyword-only dict unpacking.
                self.fig, axs = plt.subplots(*args, **subplot, **kwargs)
            else:
                raise ValueError(
                    "The 'subplot' parameter must be a tuple/list "
                    "(nrows, ncols, ...) or a dict of keyword arguments for "
                    "matplotlib.pyplot.subplots."
                )

            self.axs = axs
            # Primary axes for backward compatibility
            if hasattr(axs, "flat"):
                self.ax = axs.flat[0]
            else:
                self.ax = axs

        # ── Minor ticks ───────────────────────────────────────────────────
        if minor_ticks:
            axes_iter = (
                list(self.axs.flat)
                if hasattr(self.axs, "flat")
                else [self.axs]
            )
            for _ax in axes_iter:
                _ax.xaxis.set_minor_locator(AutoMinorLocator())
                _ax.yaxis.set_minor_locator(AutoMinorLocator())

        # ── Optional subtle grid ──────────────────────────────────────────
        if grid:
            axes_iter = (
                list(self.axs.flat)
                if hasattr(self.axs, "flat")
                else [self.axs]
            )
            for _ax in axes_iter:
                _ax.grid(True, linestyle="--", color="black",
                         alpha=0.12, lw=0.5, zorder=0)

    # ── Helpers ───────────────────────────────────────────────────────────

    def add_grid(self, linestyle="--", color="black", alpha=0.12):
        """
        Enable a subtle grid on the primary axes.

        Parameters
        ----------
        linestyle : str, default "--"
            Grid line style (matplotlib linestyle).
        color : str, default "black"
            Grid line color.
        alpha : float, default 0.12
            Opacity of the grid lines.

        Examples
        --------
        >>> ws = SciSubplot()
        >>> ws.add_grid()
        """
        self.ax.grid(True, linestyle=linestyle, color=color,
                     alpha=alpha, lw=0.5, zorder=0)

    def add_second_scale_horizontal(self, func, label=None):
        """
        Add a synchronized secondary *x* scale on the top spine.

        Tick positions match the primary axes; labels are computed by
        applying ``func`` to the primary tick values.

        Parameters
        ----------
        func : callable or tuple
            Forward mapping from primary axis values to displayed values.
            If a tuple ``(forward, inverse)``, only ``forward`` is used.
        label : str or None, default None
            Label for the secondary *x* axis (top).

        Returns
        -------
        matplotlib.axes.Axes
            The twin axes created with :meth:`~matplotlib.axes.Axes.twiny`.

        Examples
        --------
        >>> ws = SciSubplot()
        >>> ax2 = ws.add_second_scale_horizontal(lambda x: x * 2, label="2x")
        >>> isinstance(ax2.xaxis.get_major_formatter(), object)
        True
        """
        from matplotlib.ticker import FixedLocator, FuncFormatter
        forward = func[0] if isinstance(func, tuple) else func

        ax2 = self.ax.twiny()
        ax2.tick_params(which="major", direction="in", length=5, width=0.6)
        ax2.tick_params(which="minor", direction="in", length=2.5, width=0.6)
        for spine in ax2.spines.values():
            spine.set_linewidth(0.6)
        if label:
            ax2.set_xlabel(label)

        def _sync(event):
            ax2.set_xlim(self.ax.get_xlim())
            major_locs = self.ax.get_xticks()
            ax2.xaxis.set_major_locator(FixedLocator(major_locs))
            ax2.xaxis.set_major_formatter(
                FuncFormatter(lambda x, _: f"{forward(x):.4g}"))
            minor_locs = self.ax.xaxis.get_minorticklocs()
            if len(minor_locs):
                ax2.xaxis.set_minor_locator(FixedLocator(minor_locs))

        self.fig.canvas.mpl_connect("draw_event", _sync)
        return ax2

    def add_second_scale_vertical(self, func, label=None):
        """
        Add a synchronized secondary *y* scale on the right spine.

        Tick positions match the primary axes; labels are computed by
        applying ``func`` to the primary tick values.

        Parameters
        ----------
        func : callable or tuple
            Forward mapping from primary axis values to displayed values.
            If a tuple ``(forward, inverse)``, only ``forward`` is used.
        label : str or None, default None
            Label for the secondary *y* axis (right).

        Returns
        -------
        matplotlib.axes.Axes
            The twin axes created with :meth:`~matplotlib.axes.Axes.twinx`.

        Examples
        --------
        >>> ws = SciSubplot()
        >>> ax2 = ws.add_second_scale_vertical(lambda y: y ** 2, label="y²")
        >>> ax2.get_ylabel() == "y²"
        True
        """
        from matplotlib.ticker import FixedLocator, FuncFormatter
        forward = func[0] if isinstance(func, tuple) else func

        ax2 = self.ax.twinx()
        ax2.tick_params(which="major", direction="in", length=5, width=0.6)
        ax2.tick_params(which="minor", direction="in", length=2.5, width=0.6)
        for spine in ax2.spines.values():
            spine.set_linewidth(0.6)
        if label:
            ax2.set_ylabel(label)

        def _sync(event):
            ax2.set_ylim(self.ax.get_ylim())
            major_locs = self.ax.get_yticks()
            ax2.yaxis.set_major_locator(FixedLocator(major_locs))
            ax2.yaxis.set_major_formatter(
                FuncFormatter(lambda y, _: f"{forward(y):.4g}"))
            minor_locs = self.ax.yaxis.get_minorticklocs()
            if len(minor_locs):
                ax2.yaxis.set_minor_locator(FixedLocator(minor_locs))

        self.fig.canvas.mpl_connect("draw_event", _sync)
        return ax2
