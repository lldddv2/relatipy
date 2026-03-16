"""
Schwarzschild black hole 3D plotter (event horizon via oblate-spheroidal parametrization with a=0).
"""
import numpy as np

from .base import BasePlotBlackHole, EquatorialPlane
from ..base_elements import construct_surface
from .kerr import _build_kerr_surface


class PlotSchwarzschild(BasePlotBlackHole):
    """Plotter for Schwarzschild black hole (event horizon as sphere r = R_s)."""

    def __init__(self, metric):
        from relatipy.numeric.metrics.schwarzschild_metric import Schwarzschild

        if not isinstance(metric, Schwarzschild):
            raise TypeError("metric must be a Schwarzschild instance")
        super().__init__(metric)
        self.add_plane(EquatorialPlane())

    def _build_black_hole_elements(self) -> list:
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
