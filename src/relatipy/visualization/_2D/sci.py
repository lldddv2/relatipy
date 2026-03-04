import matplotlib.pyplot as plt

class SciSubplot:
    """
    Subplot con aspecto científico.
    - Ticks inside
    - Ticks ambos lados
    - Letra monospace
    """
    def __init__(self, grid=True, minor_ticks=True, *args, **kwargs):
        from matplotlib import rcParams
        rcParams['font.family'] = 'monospace'
        rcParams['font.monospace'] = ['JetBrains Mono', 'DejaVu Sans Mono']

        self.fig = plt.figure(*args, **kwargs)
        self.ax = self.fig.add_subplot(111)

        # Marco más delgado
        for spine in self.ax.spines.values():
            spine.set_linewidth(0.3)

        # Ticks inside, más largos
        self.ax.tick_params(which='major', direction='in', top=True, right=True,
                            length=6, width=0.4)

        # Ticks intermedias (minor)
        if minor_ticks:
            from matplotlib.ticker import AutoMinorLocator
            self.ax.xaxis.set_minor_locator(AutoMinorLocator())
            self.ax.yaxis.set_minor_locator(AutoMinorLocator())
            self.ax.tick_params(which='minor', direction='in', top=True, right=True,
                                length=3, width=0.4)

        if grid:
            self.add_grid()

    def add_grid(self, linestyle='--', color='black', alpha=0.15):
        self.ax.grid(True, linestyle=linestyle, color=color, alpha=alpha, lw=0.7)

    def add_second_scale_horizontal(self, func, label=None):
        """
        Agrega una segunda escala en el eje superior.
        func: función o tupla (forward, inverse) que convierte las unidades del eje x.
        """
        from matplotlib.ticker import FixedLocator, FuncFormatter
        forward = func[0] if isinstance(func, tuple) else func

        ax2 = self.ax.twiny()
        ax2.tick_params(which='major', direction='in', length=6, width=0.4)
        ax2.tick_params(which='minor', direction='in', length=3, width=0.4)
        for spine in ax2.spines.values():
            spine.set_linewidth(0.4)
        if label:
            ax2.set_xlabel(label)

        def _sync(event):
            ax2.set_xlim(self.ax.get_xlim())
            major_locs = self.ax.get_xticks()
            ax2.xaxis.set_major_locator(FixedLocator(major_locs))
            ax2.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{forward(x):.4g}'))
            minor_locs = self.ax.xaxis.get_minorticklocs()
            if len(minor_locs):
                ax2.xaxis.set_minor_locator(FixedLocator(minor_locs))

        self.fig.canvas.mpl_connect('draw_event', _sync)
        return ax2

    def add_second_scale_vertical(self, func, label=None):
        """
        Agrega una segunda escala en el eje derecho.
        func: función o tupla (forward, inverse) que convierte las unidades del eje y.
        """
        from matplotlib.ticker import FixedLocator, FuncFormatter
        forward = func[0] if isinstance(func, tuple) else func

        ax2 = self.ax.twinx()
        ax2.tick_params(which='major', direction='in', length=6, width=0.4)
        ax2.tick_params(which='minor', direction='in', length=3, width=0.4)
        for spine in ax2.spines.values():
            spine.set_linewidth(0.4)
        if label:
            ax2.set_ylabel(label)

        def _sync(event):
            ax2.set_ylim(self.ax.get_ylim())
            major_locs = self.ax.get_yticks()
            ax2.yaxis.set_major_locator(FixedLocator(major_locs))
            ax2.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f'{forward(y):.4g}'))
            minor_locs = self.ax.yaxis.get_minorticklocs()
            if len(minor_locs):
                ax2.yaxis.set_minor_locator(FixedLocator(minor_locs))

        self.fig.canvas.mpl_connect('draw_event', _sync)
        return ax2