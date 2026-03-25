"""
Three-dimensional Kerr black-hole visualization.

This module provides :class:`PlotKerr`, a Plotly-based helper that draws the
outer event horizon, an optional ergosphere, spin-direction cues, and optional
prograde/retrograde innermost stable circular orbit (ISCO) rings in the
equatorial plane. Surfaces use an oblate-spheroidal mesh consistent with
Boyer--Lindquist-style radius :math:`r` and Kerr spin parameter :math:`a`.

Notes
-----
:class:`PlotKerr` requires a :class:`~relatipy.numeric.metrics.kerr_metric.Kerr`
metric instance; passing any other metric type raises :exc:`TypeError`.

References
----------
Bardeen, J. M., Press, W. H., & Teukolsky, S. A. (1972). Rotating black holes:
locally nonrotating frames, energy extraction, and scalar synchrotron radiation.
*The Astrophysical Journal*, 178, 347--370.

Examples
--------
>>> from relatipy.numeric.metrics import Kerr
>>> from relatipy.visualization._3D.plot_black_hole import PlotKerr
>>> plotter = PlotKerr(Kerr(1.0, 0.5), show_ergosphere=True)
>>> plotter.show_isco_prograde
True
"""
from __future__ import annotations

from collections.abc import Callable

import numpy as np

from .base import BasePlotBlackHole, EquatorialPlane
from ..base_elements import construct_arrow_scatter, construct_surface
from ..orbits import construct_isco


def _build_kerr_surface(
    r_func: Callable[[np.ndarray], np.ndarray],
    a: float,
    n_theta: int = 60,
    n_phi: int = 60,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Build Cartesian coordinate meshes for a Kerr-related surface.

    The mapping uses an oblate-spheroidal parametrization with polar angle
    :math:`\\theta \\in [0, \\pi]` and azimuth :math:`\\phi \\in [0, 2\\pi)`:

    .. math::

        x &= \\sqrt{r^2 + a^2}\\, \\sin\\theta \\cos\\phi,\\\\
        y &= \\sqrt{r^2 + a^2}\\, \\sin\\theta \\sin\\phi,\\\\
        z &= r \\cos\\theta,

    where :math:`r = r(\\theta)` is supplied by ``r_func`` on the internal
    :math:`\\theta` grid. Non-positive values under the square root are clamped
    before taking :func:`numpy.sqrt` for numerical safety.

    Parameters
    ----------
    r_func : callable
        Maps a :class:`numpy.ndarray` of :math:`\\theta` values to an array of
        radii :math:`r` with the same shape (typically Boyer--Lindquist :math:`r`).
    a : float
        Kerr spin parameter in length units consistent with the metric (same
        convention as :attr:`~relatipy.numeric.metrics.kerr_metric.Kerr.a`).
    n_theta : int, optional
        Number of samples along :math:`\\theta`. Default is ``60``.
    n_phi : int, optional
        Number of samples along :math:`\\phi`. Default is ``60``.

    Returns
    -------
    x_mesh, y_mesh, z_mesh : tuple of numpy.ndarray
        Two-dimensional meshgrid arrays of shape ``(n_theta, n_phi)`` with
        Cartesian coordinates.

    Examples
    --------
    Constant-radius sphere slice (fixed :math:`r`) with :math:`a = 0`:

    >>> import numpy as np
    >>> from relatipy.visualization._3D.plot_black_hole.kerr import _build_kerr_surface
    >>> x, y, z = _build_kerr_surface(lambda th: np.ones_like(th), 0.0, n_theta=8, n_phi=8)
    >>> x.shape
    (8, 8)
    """
    theta = np.linspace(0, np.pi, n_theta)
    phi = np.linspace(0, 2 * np.pi, n_phi)
    theta_mesh, phi_mesh = np.meshgrid(theta, phi, indexing="ij")
    r_mesh = r_func(theta_mesh)
    r2_a2 = r_mesh**2 + a**2
    sqrt_r2_a2 = np.sqrt(np.maximum(r2_a2, 0))
    x_mesh = sqrt_r2_a2 * np.sin(theta_mesh) * np.cos(phi_mesh)
    y_mesh = sqrt_r2_a2 * np.sin(theta_mesh) * np.sin(phi_mesh)
    z_mesh = r_mesh * np.cos(theta_mesh)
    return (x_mesh, y_mesh, z_mesh)


class PlotKerr(BasePlotBlackHole):
    """
    Plotly-based 3D figure builder for a Kerr black hole.

    Renders the event horizon as a closed surface, optionally the ergosphere,
    a scatter-style spin arrow when :math:`a \\neq 0`, and optional ISCO rings
    in the equatorial plane. An :class:`EquatorialPlane` reference disk is added
    in the constructor.

    Parameters
    ----------
    metric : relatipy.numeric.metrics.kerr_metric.Kerr
        Kerr metric instance supplying :math:`R_s`, spin :math:`a`, and ISCO
        radii used for geometry and orbit overlays.
    show_ergosphere : bool, optional
        If ``True``, draw the ergosphere surface. Default is ``True``.
    show_isco_prograde : bool, optional
        If ``True``, draw the prograde equatorial ISCO. Default is ``True``.
    show_isco_retrograde : bool, optional
        If ``True``, draw the retrograde equatorial ISCO. Default is ``True``.

    Attributes
    ----------
    metric : relatipy.numeric.metrics.kerr_metric.Kerr
        Kerr metric bound to this plotter (see :class:`BasePlotBlackHole`).
    show_ergosphere : bool
        Whether the ergosphere is included in the figure.
    show_isco_prograde : bool
        Whether the prograde ISCO ring is drawn.
    show_isco_retrograde : bool
        Whether the retrograde ISCO ring is drawn.

    Examples
    --------
    >>> from relatipy.numeric.metrics import Kerr
    >>> from relatipy.visualization._3D.plot_black_hole import PlotKerr
    >>> fig = PlotKerr(Kerr(1.0, 0.5)).plot()
    >>> fig is not None
    True
    """

    def __init__(
        self,
        metric,
        show_ergosphere: bool = True,
        show_isco_prograde: bool = True,
        show_isco_retrograde: bool = True,
    ) -> None:
        """
        Construct a Kerr black-hole plotter bound to a ``Kerr`` metric.

        Parameters
        ----------
        metric : relatipy.numeric.metrics.kerr_metric.Kerr
            Kerr metric instance.
        show_ergosphere : bool, optional
            Include the ergosphere surface. Default is ``True``.
        show_isco_prograde : bool, optional
            Include the prograde ISCO. Default is ``True``.
        show_isco_retrograde : bool, optional
            Include the retrograde ISCO. Default is ``True``.

        Raises
        ------
        TypeError
            If ``metric`` is not an instance of
            :class:`~relatipy.numeric.metrics.kerr_metric.Kerr`.

        Examples
        --------
        >>> from relatipy.numeric.metrics import Kerr
        >>> from relatipy.visualization._3D.plot_black_hole import PlotKerr
        >>> PlotKerr(Kerr(1.0, 0.9)).metric.a > 0
        True
        """
        from relatipy.numeric.metrics.kerr_metric import Kerr

        if not isinstance(metric, Kerr):
            raise TypeError("metric must be a Kerr instance")
        super().__init__(metric)
        self.show_ergosphere = show_ergosphere
        self.show_isco_prograde = show_isco_prograde
        self.show_isco_retrograde = show_isco_retrograde
        self.add_plane(EquatorialPlane())

    def _build_black_hole_elements(self) -> list:
        """
        Assemble Plotly traces for the Kerr horizon and optional overlays.

        Computes the outer horizon radius :math:`r_+ = M + \\sqrt{M^2 - a^2}`
        with :math:`M = R_s/2`, builds the horizon mesh, optionally adds the
        ergosphere :math:`r_{\\mathrm{ergo}}(\\theta)`, a spin arrow along
        :math:`\\pm \\hat{z}` when :math:`a \\neq 0`, and ISCO paths when
        enabled.

        Returns
        -------
        list
            Plotly graph objects (e.g. surface and scatter traces) to append to
            the figure.

        Examples
        --------
        >>> from relatipy.numeric.metrics import Kerr
        >>> from relatipy.visualization._3D.plot_black_hole import PlotKerr
        >>> elements = PlotKerr(Kerr(1.0, 0.5))._build_black_hole_elements()
        >>> len(elements) >= 1
        True
        """
        M = self.metric.R_s / 2
        a = self.metric.a
        r_plus = M + np.sqrt(M**2 - a**2)

        x_h, y_h, z_h = _build_kerr_surface(
            lambda theta: np.full_like(theta, r_plus),
            a,
        )
        horizon = construct_surface(
            x_h,
            y_h,
            z_h,
            color="black",
            opacity=0.6,
            name="Event horizon",
        )
        elements = [horizon]

        if self.show_ergosphere:

            def r_ergo(theta):
                return M + np.sqrt(np.maximum(M**2 - a**2 * np.cos(theta) ** 2, 0))

            x_e, y_e, z_e = _build_kerr_surface(r_ergo, a)
            ergo = construct_surface(
                x_e,
                y_e,
                z_e,
                color="royalblue",
                opacity=0.15,
                name="Ergosphere",
            )
            elements.append(ergo)

        # Spin arrow: only if a != 0. Scatter-style so size scales with zoom (not fixed like Cone).
        R_s = self.metric.R_s
        if abs(a) > 1e-12:
            tail_z = np.sign(a) * R_s
            spin_traces = construct_arrow_scatter(
                0,
                0,
                tail_z,
                u=0,
                v=0,
                w=10 * a,
                color="blue",
                line_width=4,
                opacity=0.9,
                label=f"Spin (a: {a:.1f})",
            )
            elements.extend(spin_traces)

        # ISCOs
        if self.show_isco_prograde:
            isco_prograde = construct_isco(
                self.metric.isco_prograde, True, "blue", 8, 4, 0.6
            )
            elements.extend(isco_prograde)
        if self.show_isco_retrograde:
            isco_retrograde = construct_isco(
                self.metric.isco_retrograde, False, "green", 8, 4, 0.6
            )
            elements.extend(isco_retrograde)

        return elements
