.. _plotting-reference:

Plotting
========

.. py:module:: relatipy.plotting

:mod:`relatipy.plotting` creates figures from stored public data. Plotting
does not integrate a trajectory or interpolate a solution. Matplotlib is
imported only when a static figure is drawn; Plotly is optional for interactive
three-dimensional figures. Functions return figures without showing or saving
them.

Orbit and time-series figures
-----------------------------

``plot_solution`` draws stored Cartesian samples of a
:class:`relatipy.geodesic.Solution`. ``preview_orbit`` draws the native
osculating preview of an :class:`relatipy.geodesic.Orbit`; it is not an
integrated Kerr trajectory. ``plot_views`` renders three orthogonal Cartesian
projections at a shared physical scale. ``plot_evolution`` draws stored
coordinate or osculating-element time series. Cartesian coordinates are
expressed in the selected displayed length unit, while axes retain their
physical units.

``plot_solution_interactive`` and ``preview_orbit_interactive`` create Plotly
three-dimensional figures. Static plotting functions raise :class:`ValueError`
for invalid projections, source data, units, or style choices. Missing optional
plotting dependencies raise :class:`ImportError` when the affected function is
called.

Fit and publication figures
---------------------------

``plot_observables`` plots observed and predicted astrometry or velocity
series. ``plot_corner`` expects finite posterior samples with shape ``(n, k)``,
where ``n >= 2`` and ``k >= 2``; labels contain one display name, including any
unit, per parameter. ``Style`` and its immutable nested settings control
figure sizes, fonts, axes, roles, and plot-specific appearance. ``Target``
contains dimensions in millimetres and text sizes in points. ``MM`` converts
millimetres to inches. ``new_figure``, ``format_axes``,
``publication_style``, and ``save_figure`` support consistently formatted
custom Matplotlib figures.

API
---

Style configuration
~~~~~~~~~~~~~~~~~~~

.. autodata:: relatipy.plotting.MM

.. autodata:: relatipy.plotting.TARGETS

.. autodata:: relatipy.plotting.DEFAULT_STYLE

.. autoclass:: relatipy.plotting.Target
   :show-inheritance:

.. autoclass:: relatipy.plotting.Frame
   :show-inheritance:

.. autoclass:: relatipy.plotting.Views
   :show-inheritance:

.. autoclass:: relatipy.plotting.Static3D
   :show-inheritance:

.. autoclass:: relatipy.plotting.Interactive
   :show-inheritance:

.. autoclass:: relatipy.plotting.Observables
   :show-inheritance:

.. autoclass:: relatipy.plotting.Corner
   :show-inheritance:

.. autoclass:: relatipy.plotting.Style
   :show-inheritance:

.. autofunction:: relatipy.plotting.figure_width

.. autofunction:: relatipy.plotting.format_axes

.. autofunction:: relatipy.plotting.new_figure

.. autofunction:: relatipy.plotting.publication_style

.. autofunction:: relatipy.plotting.rc_params

.. autofunction:: relatipy.plotting.save_figure

.. autofunction:: relatipy.plotting.set_axis_label

Figure functions
~~~~~~~~~~~~~~~~

.. autofunction:: relatipy.plotting.plot_solution

.. autofunction:: relatipy.plotting.preview_orbit

.. autofunction:: relatipy.plotting.plot_views

.. autofunction:: relatipy.plotting.plot_solution_interactive

.. autofunction:: relatipy.plotting.preview_orbit_interactive

.. autofunction:: relatipy.plotting.plot_evolution

.. autofunction:: relatipy.plotting.plot_observables

.. autofunction:: relatipy.plotting.plot_corner
