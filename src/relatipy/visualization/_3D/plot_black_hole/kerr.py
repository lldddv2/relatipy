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
import plotly.graph_objects as go

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


def _arrow_with_label(
    p0: np.ndarray,
    p1: np.ndarray,
    perp: np.ndarray,
    color: str,
    label: str,
    line_width: int = 3,
    opacity: float = 0.85,
    legendgroup: str = "visual_line",
) -> list:
    """
    Shaft + V-head + text label; V wings constrained to the plane via *perp*.

    Parameters
    ----------
    p0, p1 : ndarray, shape (3,)
        Start and tip of the arrow.
    perp : ndarray, shape (3,)
        Unit vector perpendicular to ``p1 - p0`` for the V-head wings.
    color : str
        Plotly color string.
    label : str
        Text placed at the tip and used as trace name.
    line_width, opacity, legendgroup :
        Visual style; all traces share the same legendgroup.

    Returns
    -------
    list of plotly.graph_objects.Scatter3d
        Four traces: shaft, two wing segments, text marker.
    """
    dx = p1 - p0
    length = float(np.linalg.norm(dx))
    if length < 1e-12:
        return []
    d_hat = dx / length
    back = p1 - 0.2 * length * d_hat
    left = back + 0.12 * length * perp
    right = back - 0.12 * length * perp

    def _seg(a, b):
        return go.Scatter3d(
            x=[float(a[0]), float(b[0])],
            y=[float(a[1]), float(b[1])],
            z=[float(a[2]), float(b[2])],
            mode="lines",
            line=dict(color=color, width=line_width),
            opacity=opacity,
            legendgroup=legendgroup,
            showlegend=False,
            name=label,
        )

    return [
        _seg(p0, p1),
        _seg(p1, left),
        _seg(p1, right),
        go.Scatter3d(
            x=[float(p1[0])], y=[float(p1[1])], z=[float(p1[2])],
            mode="text",
            text=[label],
            textfont=dict(color=color, size=11),
            legendgroup=legendgroup,
            showlegend=False,
            name=label,
        ),
    ]


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
        self._poe_orientation = None
        self.add_plane(EquatorialPlane())

    def add_path(self, path, color='red', opacity=0.6, size=5, show_start=True, show_end=True, label="orbit") -> None:
        """
        Append an orbit path; if ``path`` is a ``ProperOrbitalElements`` instance,
        store its orientation angles to activate the visual-line indicator.

        Parameters
        ----------
        path : object
            Trajectory with ``convert_to("Cartesian")`` support.
        color, opacity, size, show_start, show_end, label :
            Forwarded to :class:`~relatipy.visualization._3D.plot_black_hole.base.OrbitPath`.

        Examples
        --------
        >>> from relatipy.numeric.metrics import Kerr
        >>> from relatipy.visualization._3D.plot_black_hole import PlotKerr
        >>> pk = PlotKerr(Kerr(1.0, 0.5))
        >>> pk.add_path(None)
        >>> pk._poe_orientation is None
        True
        """
        super().add_path(path, color=color, opacity=opacity, size=size,
                         show_start=show_start, show_end=show_end, label=label)
        from relatipy.numeric.coordinates.proper_orbital_elements import ProperOrbitalElements
        if isinstance(path, ProperOrbitalElements):
            self._poe_orientation = (path.zeta, path.eta)

    def _build_visual_indicator(
        self,
        zeta_deg: float,
        eta_deg: float,
        half_length: float,
        plane_size: float,
    ) -> list:
        """
        Build visual-line, sky-plane, and Dec=0 arrow traces.

        The visual line runs along the observer's line-of-sight direction
        :math:`\\hat{d} = (-\\sin\\zeta, 0, \\cos\\zeta)` in the BH frame,
        from :math:`-2u\\hat{d}` to :math:`+2u\\hat{d}` (``half_length`` =
        :math:`2u`).  A cone arrowhead marks the :math:`+2u` end (toward the
        observer).  A translucent square patch of half-side ``plane_size``
        spans the sky plane (perpendicular to :math:`\\hat{d}` through the
        origin).  An arrow in the sky plane indicates the Dec = 0 direction
        (:math:`\\hat{x}_{\\mathrm{sky}}` mapped to the BH frame via the
        rotation :math:`R = R_y(-\\zeta)\\,R_z(-\\eta)`).

        Parameters
        ----------
        zeta_deg : float
            Polar angle of the BH spin from the line of sight, in degrees.
        eta_deg : float
            Position angle of the spin projection on the sky plane, in degrees.
        half_length : float
            Half-extent of the dashed line in coordinate units.
        plane_size : float
            Half-side of the sky-plane square patch (matches equatorial plane).

        Returns
        -------
        list of plotly.graph_objects traces

        Examples
        --------
        >>> from relatipy.numeric.metrics import Kerr
        >>> from relatipy.visualization._3D.plot_black_hole import PlotKerr
        >>> pk = PlotKerr(Kerr(1.0, 0.5))
        >>> traces = pk._build_visual_indicator(45.0, 90.0, 10.0, 20.0)
        >>> len(traces) >= 3
        True
        """
        zeta_rad = np.radians(zeta_deg)
        eta_rad = np.radians(eta_deg)
        cz, sz = np.cos(zeta_rad), np.sin(zeta_rad)
        cn, sn = np.cos(eta_rad), np.sin(eta_rad)

        # LOS direction in BH frame: R @ ẑ_sky = (-sin ζ, 0, cos ζ)
        d = np.array([-sz, 0.0, cz])
        p0 = -half_length * d
        p1 = half_length * d

        traces = []

        # Dashed line through BH along line of sight
        traces.append(go.Scatter3d(
            x=[p0[0], p1[0]],
            y=[p0[1], p1[1]],
            z=[p0[2], p1[2]],
            mode="lines",
            line=dict(color="green", width=4, dash="dash"),
            opacity=0.9,
            name="Visual line",
            legendgroup="visual_line",
            showlegend=True,
        ))

        # Cone arrowhead at p1 (toward observer)
        cone_len = 0.15 * half_length
        traces.append(go.Cone(
            x=[p1[0] - cone_len * d[0]],
            y=[p1[1] - cone_len * d[1]],
            z=[p1[2] - cone_len * d[2]],
            u=[cone_len * d[0]],
            v=[cone_len * d[1]],
            w=[cone_len * d[2]],
            anchor="tail",
            showscale=False,
            sizemode="absolute",
            colorscale=[[0, "green"], [1, "green"]],
            opacity=0.9,
            sizeref=cone_len / 3,
            legendgroup="visual_line",
            showlegend=False,
            name="Visual line",
        ))

        # Sky plane perpendicular to LOS through origin (size = equatorial plane)
        ref = np.array([0.0, 1.0, 0.0]) if abs(d[1]) < 0.9 else np.array([1.0, 0.0, 0.0])
        t1 = np.cross(d, ref)
        t1 /= np.linalg.norm(t1)
        t2 = np.cross(d, t1)
        t2 /= np.linalg.norm(t2)

        vals = np.linspace(-plane_size, plane_size, 20)
        S, Q = np.meshgrid(vals, vals)
        X = S * t1[0] + Q * t2[0]
        Y = S * t1[1] + Q * t2[1]
        Z = S * t1[2] + Q * t2[2]

        traces.append(go.Surface(
            x=X, y=Y, z=Z,
            colorscale=[[0, "green"], [1, "green"]],
            opacity=0.1,
            showscale=False,
            name="Sky plane",
            showlegend=True,
        ))

        # Dec arrow: R @ x̂_sky = (cos ζ cos η, −sin η, sin ζ cos η)
        dec_dir = np.array([cz * cn, -sn, sz * cn])
        # RA arrow: R @ ŷ_sky = (cos ζ sin η, cos η, sin ζ sin η)
        ra_dir = np.array([cz * sn, cn, sz * sn])
        origin = np.zeros(3)
        arrow_len = 0.35 * plane_size
        # Wings constrained to sky plane: Dec uses RA as perp, RA uses Dec as perp
        traces.extend(_arrow_with_label(origin, arrow_len * dec_dir, ra_dir,  "green", "DEC 0"))
        traces.extend(_arrow_with_label(origin, arrow_len * ra_dir,  dec_dir, "green", "AR 0"))

        return traces

    def plot(self, show_center=True, show_visual_indicator=True, show_equatorial_axes=True):
        """
        Build the figure; adds visual-line indicator when a
        ``ProperOrbitalElements`` path has been registered via :meth:`add_path`.

        Also sets the initial camera to the observer's point of view (looking
        from the line-of-sight direction toward the BH) with Dec north up.

        Parameters
        ----------
        show_center : bool, optional
            Draw a marker at the coordinate origin. Default is ``True``.
        show_visual_indicator : bool, optional
            Draw the visual line, sky plane, and Dec = 0 arrow when a
            ``ProperOrbitalElements`` path is registered. Default is ``True``.
        show_equatorial_axes : bool, optional
            Draw x and y axis lines on the equatorial plane (opacity 0.2).
            Default is ``True``.

        Returns
        -------
        plotly.graph_objects.Figure

        Examples
        --------
        >>> from relatipy.numeric.metrics import Kerr
        >>> from relatipy.visualization._3D.plot_black_hole import PlotKerr
        >>> fig = PlotKerr(Kerr(1.0, 0.5)).plot()
        >>> fig is not None
        True
        """
        fig = super().plot(show_center)

        # Scene half-range (mirrors base-class computation, reused below)
        max_range = self.metric.R_s * 10
        for op in self._paths:
            cart = op.path.convert_to("Cartesian")
            sv = cart.state_vector
            spatial = sv[1:4] if sv.ndim == 2 else np.asarray(sv[1:4]).reshape(3, -1)
            max_range = max(max_range, float(np.max(np.abs(spatial))))
        plane_size = max_range * 1.2

        if show_equatorial_axes:
            _ax = np.array([plane_size, 0.0, 0.0])
            _ay = np.array([0.0, plane_size, 0.0])
            for p1, perp, label in [
                (_ax, np.array([0., 1., 0.]), "x"),
                (_ay, np.array([1., 0., 0.]), "y"),
            ]:
                for t in _arrow_with_label(
                    -p1, p1, perp, "gray", label,
                    line_width=2, opacity=0.2, legendgroup="equatorial_axes",
                ):
                    fig.add_trace(t)

        if self._poe_orientation is not None and show_visual_indicator:
            zeta_deg, eta_deg = self._poe_orientation
            u = self.metric.R_s + 10.0 * abs(self.metric.a)
            half_length = 2.0 * u

            for trace in self._build_visual_indicator(zeta_deg, eta_deg, half_length, plane_size):
                fig.add_trace(trace)

            # Camera: observer sits at +d (LOS direction), Dec north up
            zeta_rad = np.radians(zeta_deg)
            eta_rad = np.radians(eta_deg)
            cz, sz = np.cos(zeta_rad), np.sin(zeta_rad)
            cn, sn = np.cos(eta_rad), np.sin(eta_rad)
            d = np.array([-sz, 0.0, cz])
            dec_zero = np.array([cz * cn, -sn, sz * cn])
            eye_scale = 2.0
            fig.update_layout(scene=dict(camera=dict(
                eye=dict(x=float(d[0]) * eye_scale,
                         y=float(d[1]) * eye_scale,
                         z=float(d[2]) * eye_scale),
                up=dict(x=float(dec_zero[0]),
                        y=float(dec_zero[1]),
                        z=float(dec_zero[2])),
            )))

        return fig

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
