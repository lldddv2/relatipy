"""
Kerr black hole 3D plotter (event horizon + optional ergosphere).
"""
import numpy as np

from .base import BasePlotBlackHole, EquatorialPlane
from ..base_elements import construct_surface, construct_arrow_scatter, construct_path
import plotly.graph_objects as go
from ..orbits import construct_isco


def _build_kerr_surface(r_func, a, n_theta=60, n_phi=60) -> tuple:
    """
    Build (x, y, z) mesh using oblate-spheroidal parametrization:
      x = √(r² + a²) sin θ cos φ
      y = √(r² + a²) sin θ sin φ
      z = r cos θ
    with r = r_func(theta_mesh), θ ∈ [0, π], φ ∈ [0, 2π].

    Parameters
    ----------
    r_func : callable
        r_func(theta_array) -> array (same shape as theta).
    a : float
        Kerr spin parameter in geometric units.
    n_theta, n_phi : int
        Grid resolution.

    Returns
    -------
    (x_mesh, y_mesh, z_mesh) : tuple of 2D np.ndarray
    """
    theta = np.linspace(0, np.pi, n_theta)
    phi = np.linspace(0, 2 * np.pi, n_phi)
    theta_mesh, phi_mesh = np.meshgrid(theta, phi, indexing='ij')
    r_mesh = r_func(theta_mesh)
    r2_a2 = r_mesh**2 + a**2
    sqrt_r2_a2 = np.sqrt(np.maximum(r2_a2, 0))
    x_mesh = sqrt_r2_a2 * np.sin(theta_mesh) * np.cos(phi_mesh)
    y_mesh = sqrt_r2_a2 * np.sin(theta_mesh) * np.sin(phi_mesh)
    z_mesh = r_mesh * np.cos(theta_mesh)
    return (x_mesh, y_mesh, z_mesh)


class PlotKerr(BasePlotBlackHole):
    """Plotter for Kerr black hole (event horizon + optional ergosphere)."""

    def __init__(self, metric, show_ergosphere=True, show_isco_prograde=True, show_isco_retrograde=True):
        from relatipy.numeric.metrics.kerr_metric import Kerr

        if not isinstance(metric, Kerr):
            raise TypeError("metric must be a Kerr instance")
        super().__init__(metric)
        self.show_ergosphere = show_ergosphere
        self.show_isco_prograde = show_isco_prograde
        self.show_isco_retrograde = show_isco_retrograde
        self.add_plane(EquatorialPlane())

    def _build_black_hole_elements(self) -> list:
        M = self.metric.R_s / 2
        a = self.metric.a
        r_plus = M + np.sqrt(M**2 - a**2)

        x_h, y_h, z_h = _build_kerr_surface(
            lambda theta: np.full_like(theta, r_plus),
            a,
        )
        horizon = construct_surface(
            x_h, y_h, z_h,
            color='black',
            opacity=0.6,
            name="Event horizon",
        )
        elements = [horizon]

        if self.show_ergosphere:
            def r_ergo(theta):
                return M + np.sqrt(np.maximum(M**2 - a**2 * np.cos(theta)**2, 0))

            x_e, y_e, z_e = _build_kerr_surface(r_ergo, a)
            ergo = construct_surface(
                x_e, y_e, z_e,
                color='royalblue',
                opacity=0.15,
                name="Ergosphere",
            )
            elements.append(ergo)

        # Spin arrow: only if a != 0. Scatter-style so size scales with zoom (not fixed like Cone).
        R_s = self.metric.R_s
        if abs(a) > 1e-12:
            tail_z = np.sign(a) * R_s
            spin_traces = construct_arrow_scatter(
                0, 0, tail_z,
                u=0, v=0, w=10 * a,
                color='blue',
                line_width=4,
                opacity=0.9,
                label=f"Spin (a: {a:.1f})",
            )
            elements.extend(spin_traces)
        
        # ISCOs
        if self.show_isco_prograde:
            isco_prograde = construct_isco(self.metric.isco_prograde, True, 'blue', 8, 4, 0.6)
            elements.extend(isco_prograde)
        if self.show_isco_retrograde:
            isco_retrograde = construct_isco(self.metric.isco_retrograde, False, 'green', 8, 4, 0.6)
            elements.extend(isco_retrograde)

        return elements
