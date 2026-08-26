"""
Three-dimensional black-hole visualization with Plotly.

This subpackage provides plot constructors for Schwarzschild and Kerr spacetimes
and small helper types for orbit paths and reference planes. Plotters build on a
shared base that attaches metric instances, orbit paths, and optional equatorial
geometry, then produce :class:`plotly.graph_objects.Figure` scenes.

Exports
-------
``PlotSchwarzschild``, ``PlotKerr``
    High-level plotters for the static Schwarzschild horizon and the rotating
    Kerr horizon (and related scene elements).
``OrbitPath``, ``EquatorialPlane``, ``SquarePlane``
    Lightweight containers for styling and geometry passed into the base
    plotter API.

Notes
-----
Figures are built with Plotly. ``PlotSchwarzschild`` expects a
:class:`~relatipy.numeric.metrics.schwarzschild_metric.Schwarzschild` metric;
``PlotKerr`` expects a compatible Kerr metric instance. The base class
:class:`~relatipy.visualization._3D.plot_black_hole.base.BasePlotBlackHole` is
not re-exported here but is used internally.

Examples
--------
>>> from relatipy.visualization._3D.plot_black_hole import (
...     PlotSchwarzschild,
...     PlotKerr,
...     OrbitPath,
... )
>>> PlotKerr.__name__
'PlotKerr'
>>> PlotSchwarzschild.__name__
'PlotSchwarzschild'
>>> OrbitPath.__name__
'OrbitPath'
"""

from .base import OrbitPath, EquatorialPlane, SquarePlane, BasePlotBlackHole
from .kerr import PlotKerr
from .schwarzschild import PlotSchwarzschild

__all__ = [
    "PlotSchwarzschild",
    "PlotKerr",
    "OrbitPath",
    "EquatorialPlane",
    "SquarePlane",
]
