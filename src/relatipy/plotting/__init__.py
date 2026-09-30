"""Publication-ready figures of orbits, observables and posteriors.

Orbit figures (:func:`plot_solution`, :func:`preview_orbit`,
:func:`plot_views`) show stored samples or the native osculating preview as
figures: ``"views"`` is the static Matplotlib default. Explicit 3D views
default to interactive Plotly; ``interactive=`` chooses the rendering mode.
:func:`plot_evolution` shows coordinate time series. :func:`plot_observables`
and :func:`plot_corner` report an orbit fit. Every figure follows a :class:`Style`, whose 2D format is also
available for custom figures through :func:`new_figure`,
:func:`format_axes`, :func:`publication_style` and :func:`save_figure`.

Matplotlib and Plotly are imported only when a figure is drawn. Figures are
returned without being shown or saved.
"""

from .corner import plot_corner
from .evolution import plot_evolution
from .interactive import plot_solution_interactive, preview_orbit_interactive
from .observables import plot_observables
from .orbits import plot_solution, preview_orbit
from .style import (
    DEFAULT_STYLE, MM, TARGETS, Corner, Frame, Interactive, Observables, Static3D,
    Style, Target, Views, figure_width, format_axes, new_figure, publication_style,
    rc_params, save_figure, set_axis_label,
)
from .views import plot_views

__all__ = [
    "DEFAULT_STYLE", "MM", "TARGETS", "Corner", "Frame", "Interactive", "Observables",
    "Static3D", "Style", "Target", "Views", "figure_width", "format_axes",
    "new_figure", "plot_corner", "plot_evolution", "plot_observables", "plot_solution",
    "plot_solution_interactive", "plot_views", "preview_orbit",
    "preview_orbit_interactive", "publication_style", "rc_params", "save_figure",
    "set_axis_label",
]
