"""
Three-dimensional visualization of a Schwarzschild black hole.

This module provides a Plotly-based plotter that draws the spherically symmetric
event horizon at :math:`r = R_s` (Schwarzschild radius). The horizon mesh is built
with the same oblate-spheroidal parametrization used for Kerr (see
``_build_kerr_surface`` in ``kerr.py``), specialized to zero spin (:math:`a = 0`)
and constant radial coordinate.

Notes
-----
The outer horizon of the Schwarzschild geometry is a coordinate sphere of radius
:math:`R_s`. The implementation reuses ``_build_kerr_surface`` with
:math:`r(\\theta) = R_s` everywhere and :math:`a = 0`, which reduces to Cartesian
coordinates on a sphere.

Examples
--------
>>> from relatipy.numeric.metrics.schwarzschild_metric import Schwarzschild
>>> from relatipy.visualization._3D.plot_black_hole.schwarzschild import PlotSchwarzschild
>>> metric = Schwarzschild(5.972e24)
>>> plotter = PlotSchwarzschild(metric)
>>> isinstance(plotter.metric, Schwarzschild)
True
"""
import numpy as np

from .base import BasePlotBlackHole, EquatorialPlane
from ..base_elements import construct_surface
from .kerr import _build_kerr_surface


class PlotSchwarzschild(BasePlotBlackHole):
    """
    Plotly 3D plotter for a Schwarzschild black hole.

    The event horizon is rendered as a sphere of radius ``metric.R_s``. An
    equatorial reference plane is added by default (see :class:`EquatorialPlane`).

    Parameters
    ----------
    metric : Schwarzschild
        Schwarzschild metric instance providing ``R_s`` for the horizon radius.

    Attributes
    ----------
    metric : Schwarzschild
        The Schwarzschild metric (same as ``metric`` argument).
    """

    def __init__(self, metric):
        """
        Construct a Schwarzschild black-hole plotter.

        Parameters
        ----------
        metric : Schwarzschild
            Metric instance used for the horizon radius and related quantities.

        Raises
        ------
        TypeError
            If ``metric`` is not a :class:`~relatipy.numeric.metrics.schwarzschild_metric.Schwarzschild`
            instance.

        Examples
        --------
        >>> from relatipy.numeric.metrics.schwarzschild_metric import Schwarzschild
        >>> from relatipy.visualization._3D.plot_black_hole.schwarzschild import PlotSchwarzschild
        >>> plotter = PlotSchwarzschild(Schwarzschild(5.972e24))
        >>> plotter.metric.R_s > 0
        True
        """
        from relatipy.numeric.metrics.schwarzschild_metric import Schwarzschild

        if not isinstance(metric, Schwarzschild):
            raise TypeError("metric must be a Schwarzschild instance")
        super().__init__(metric)
        self.add_plane(EquatorialPlane())

    def _build_black_hole_elements(self) -> list:
        """
        Build Plotly surfaces for the event horizon.

        The horizon uses oblate-spheroidal coordinates with spin :math:`a = 0` and
        constant :math:`r = R_s`, producing a spherical surface in Cartesian space.

        Returns
        -------
        list
            A one-element list containing a ``plotly.graph_objects.Surface`` for the
            event horizon (black, semi-transparent).

        Notes
        -----
        Grid resolution is fixed at ``n_theta = n_phi = 60`` for the Kerr surface
        builder.

        Examples
        --------
        >>> from relatipy.numeric.metrics.schwarzschild_metric import Schwarzschild
        >>> from relatipy.visualization._3D.plot_black_hole.schwarzschild import PlotSchwarzschild
        >>> plotter = PlotSchwarzschild(Schwarzschild(5.972e24))
        >>> elements = plotter._build_black_hole_elements()
        >>> len(elements) == 1
        True
        """
        # Oblate-spheroidal with a=0, r = R_s constant
        x, y, z = _build_kerr_surface(
            lambda theta: np.full_like(theta, self.metric.R_s),
            0,
            n_theta=60,
            n_phi=60,
        )
        horizon = construct_surface(
            x, y, z,
            color='black',
            opacity=0.6,
            name="Event horizon",
        )
        return [horizon]
