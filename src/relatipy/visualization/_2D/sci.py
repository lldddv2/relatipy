import matplotlib.pyplot as plt


class SciSubplot:
    """
    Subplot con aspecto científico estilo GRAVITY/ESO (Gillessen et al.)
    - Fuente serif tipo LaTeX (STIXGeneral / Computer Modern)
    - Ticks hacia adentro en los 4 lados
    - Minor ticks automáticos
    - Marco delgado
    - Sin grid por defecto (papers de astrofísica raramente lo usan)
    """

    # ── Paleta de colores estilo GRAVITY ──────────────────────────────────
    COLORS = {
        "black":   "black",
        "cyan":    "#4DB8FF",   # GRAVITY
        "red":     "#CC3300",   # flares / referencia
        "blue":    "#1A3E8C",   # SINFONI
        "gray":    "#999999",   # curvas de ajuste
    }

    def __init__(self, grid=False, minor_ticks=True, *args, **kwargs):
        import matplotlib as mpl
        from matplotlib.ticker import AutoMinorLocator

        # ── Tipografía estilo LaTeX ────────────────────────────────────────
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

            # ── Líneas / marcadores ───────────────────────────────────────
            "lines.linewidth":      1.2,
            "lines.markersize":     4,
            "errorbar.capsize":     0,      # sin caps en barras de error

            # ── Marco y fondo ─────────────────────────────────────────────
            "axes.linewidth":       0.6,
            "axes.edgecolor":       "black",
            "axes.facecolor":       "white",
            "figure.facecolor":     "white",

            # ── Leyenda ───────────────────────────────────────────────────
            "legend.frameon":       False,
            "legend.handlelength":  1.5,
            "legend.handletextpad": 0.5,    # espacio entre handle y texto
            "legend.labelspacing":  0.4,    # espacio vertical entre entradas

            # ── Resolución ────────────────────────────────────────────────
            "figure.dpi":           125,
            "savefig.dpi":          300,
            "savefig.bbox":         "tight",
        })

        self.fig = plt.figure(*args, **kwargs)
        self.ax  = self.fig.add_subplot(111)

        # ── Minor ticks ───────────────────────────────────────────────────
        if minor_ticks:
            self.ax.xaxis.set_minor_locator(AutoMinorLocator())
            self.ax.yaxis.set_minor_locator(AutoMinorLocator())

        # ── Grid (opcional, muy sutil) ────────────────────────────────────
        if grid:
            self.add_grid()

    # ── Helpers ───────────────────────────────────────────────────────────

    def add_grid(self, linestyle="--", color="black", alpha=0.12):
        self.ax.grid(True, linestyle=linestyle, color=color,
                     alpha=alpha, lw=0.5, zorder=0)

    def add_second_scale_horizontal(self, func, label=None):
        """
        Segunda escala en el eje superior.
        func: callable o tupla (forward, inverse).
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
        Segunda escala en el eje derecho.
        func: callable o tupla (forward, inverse).
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