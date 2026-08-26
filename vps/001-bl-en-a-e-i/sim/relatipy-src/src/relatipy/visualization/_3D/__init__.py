"""
Three-dimensional visualization for relativistic orbits and black-hole spacetimes.

This subpackage re-exports helpers for building Plotly 3D figures: a convenience
routine that combines an orbit with a compact object and reference geometry,
plus higher-level plot types for Schwarzschild and Kerr spacetimes and
composable scene elements (paths and planes).

Exports
-------
``construct_basic_path_plot``
    Build a :class:`plotly.graph_objects.Figure` for a Cartesian trajectory,
    optional equatorial reference plane, and sphere of radius ``R_s``.
``PlotSchwarzschild``, ``PlotKerr``
    Plot constructors for static black-hole visualizations.
``OrbitPath``, ``EquatorialPlane``, ``SquarePlane``
    Reusable building blocks for paths and planes in 3D scenes.

Notes
-----
Rendering uses Plotly (:mod:`plotly.graph_objects`). Implementation details live
in :mod:`~relatipy.visualization._3D.orbits` and
:mod:`~relatipy.visualization._3D.plot_black_hole`.

Examples
--------
>>> from relatipy.visualization._3D import __all__
>>> sorted(__all__)
['EquatorialPlane', 'OrbitPath', 'PlotKerr', 'PlotSchwarzschild', 'SquarePlane', 'construct_basic_path_plot']
"""

from .orbits import construct_basic_path_plot
from .plot_black_hole import (
    PlotSchwarzschild,
    PlotKerr,
    OrbitPath,
    EquatorialPlane,
    SquarePlane,
)

__all__ = [
    "construct_basic_path_plot",
    "PlotSchwarzschild",
    "PlotKerr",
    "OrbitPath",
    "EquatorialPlane",
    "SquarePlane",
]
