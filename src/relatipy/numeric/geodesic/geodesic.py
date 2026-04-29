"""
Geodesic integration for test particles in a given spacetime metric.

This module defines :class:`Geodesic`, which integrates the geodesic equation
using SciPy ODE solvers, a specialized Kerr ``Radau2`` path, or a 6th-order
symplectic Yoshida integrator. Paths can be returned as coordinate objects in
the metric's native chart or converted back to the caller's original system.

Notes
-----
The geodesic acceleration in coordinate form is

.. math::

    \\frac{\\mathrm{d} u^\\sigma}{\\mathrm{d}\\tau}
    = -\\Gamma^\\sigma_{\\mu\\nu} u^\\mu u^\\nu,

where :math:`u^\\mu = \\mathrm{d} x^\\mu / \\mathrm{d}\\tau` and
:math:`\\Gamma^\\sigma_{\\mu\\nu}` are the Christoffel symbols of the metric.

See Also
--------
relatipy.numeric.metrics :
    Metric implementations used with :class:`Geodesic`.

Examples
--------
The class is constructed with a metric instance; integration is done via
:meth:`Geodesic.get_path` or :meth:`Geodesic.get_path_periodic`:

>>> from relatipy.numeric.geodesic import Geodesic  # doctest: +SKIP
>>> # g = SchwarzschildMetric(...); geo = Geodesic(g)  # doctest: +SKIP
"""

import numpy
from itertools import product
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

from ..coordinates import coordinate_systems


class Geodesic:
    """
    Integrate timelike geodesics for a fixed spacetime metric.

    The right-hand side uses Christoffel symbols from ``metric`` to advance
    :math:`(x^\\mu, u^\\mu)` in proper time. Optional Kerr-specific projection
    and integrators assume Boyer-Lindquist coordinates when selected.

    Parameters
    ----------
    metric : object
        A metric object providing ``valid_coordinate``, ``kwargs``,
        :meth:`get_christoffel_symbols`, :meth:`metric`,
        :meth:`get_4state_vector`, and :meth:`get_dxs_dt_from_4velocity`, and
        (for Kerr) attributes such as ``a`` and ``R_s`` as required by the
        chosen integrator.

    Attributes
    ----------
    metric : object
        The spacetime metric used for connection coefficients and state
        assembly.
    valid_coordinate : str
        Name of the coordinate system accepted by this metric (same as
        ``metric.valid_coordinate``).

    Examples
    --------
    >>> from relatipy.numeric.geodesic import Geodesic  # doctest: +SKIP
    >>> # geo = Geodesic(my_metric)  # doctest: +SKIP
    """

    def __init__(self, metric):
        self.metric = metric
        self.valid_coordinate = self.metric.valid_coordinate

    def model_geodesic(self, tau, ys0):
        """
        Right-hand side of the geodesic ODE in first-order form.

        State vector ``ys0`` is ``[q^0, q^1, q^2, q^3, u^0, u^1, u^2, u^3]``.
        The acceleration components are

        .. math::

            a^\\sigma = -\\Gamma^\\sigma_{\\mu\\nu} u^\\mu u^\\nu.

        Parameters
        ----------
        tau : float
            Affine parameter (proper time) at the current step. Not used in the
            acceleration for stationary metrics because Christoffel symbols
            depend only on position.
        ys0 : ndarray, shape (8,)
            Concatenated position :math:`q^\\mu` and four-velocity :math:`u^\\mu`.

        Returns
        -------
        ndarray, shape (8,)
            Time derivative of the state,
            :math:`[\\mathrm{d}q^\\mu/\\mathrm{d}\\tau, \\mathrm{d}u^\\mu/\\mathrm{d}\\tau]`.

        Examples
        --------
        >>> import numpy as np
        >>> # rhs = Geodesic(metric).model_geodesic(0.0, y8)  # doctest: +SKIP
        """
        xs0 = ys0[:4]
        us0 = ys0[4:]

        as_ = numpy.zeros(4)

        chris = self.metric.get_christoffel_symbols(xs0)

        for sigma, mu, nu in product(range(4), repeat=3):
            as_[sigma] -= chris[sigma, mu, nu] * us0[mu] * us0[nu]

        return numpy.concatenate([us0, as_])

    def _project_constraints(self, y, E0, Lz0, Q0, tol=1e-12, max_iter=20):
        """
        Newton iteration to project four-velocity onto Kerr constraints.

        For Boyer-Lindquist Kerr geodesics, the four-velocity is adjusted so
        that (approximately) the normalization, conserved energy :math:`E`,
        azimuthal angular momentum :math:`L_z`, and Carter constant :math:`Q`
        match the targets ``E0``, ``Lz0``, and ``Q0``.

        The residuals are

        .. math::

            C_1 &= g_{\\mu\\nu} u^\\mu u^\\nu - 1,\\\\
            C_2 &= -g_{0\\mu} u^\\mu - E_0,\\\\
            C_3 &= g_{3\\mu} u^\\mu - L_{z,0},\\\\
            C_4 &= p_\\theta^2 + \\cos^2\\theta\\left[
                a^2(1-E^2) + \\frac{L_z^2}{\\sin^2\\theta}\\right] - Q_0,

        with :math:`p_\\mu = g_{\\mu\\nu} u^\\nu`.

        Parameters
        ----------
        y : ndarray, shape (8,)
            State ``[q^0, q^1, q^2, q^3, u^0, u^1, u^2, u^3]`` in Boyer-Lindquist
            coordinates.
        E0 : float
            Target conserved energy parameter :math:`E_0`.
        Lz0 : float
            Target conserved :math:`z`-component of angular momentum :math:`L_{z,0}`.
        Q0 : float
            Target Carter constant :math:`Q_0`.
        tol : float, optional
            Convergence tolerance on :math:`\\max_i |C_i|`. Default is ``1e-12``.
        max_iter : int, optional
            Maximum Newton iterations. Default is ``20``.

        Returns
        -------
        ndarray, shape (8,)
            State with the same positions as ``y`` and an updated four-velocity
            ``u^\\mu``.

        Notes
        -----
        If the iteration does not reach ``tol`` within ``max_iter`` steps, the
        last iterate is returned.

        Examples
        --------
        >>> import numpy as np
        >>> # y_new = Geodesic(kerr_metric)._project_constraints(y, E0, Lz0, Q0)  # doctest: +SKIP
        """
        q = y[:4]
        u = y[4:].copy()
        g = self.metric.metric(q)
        a_spin = self.metric.a
        theta = q[2]
        s = numpy.sin(theta)
        c = numpy.cos(theta)
        s2 = s * s
        c2 = c * c
        g00, g03, g30, g33, g22 = g[0, 0], g[0, 3], g[3, 0], g[3, 3], g[2, 2]

        for _ in range(max_iter):
            p = g @ u
            E    = -p[0]
            Lz   =  p[3]
            p_th =  p[2]
            C = numpy.array([
                numpy.dot(u, p) - 1.0,
                E  - E0,
                Lz - Lz0,
                p_th**2 + c2 * (a_spin**2 * (1.0 - E**2) + Lz**2 / s2) - Q0,
            ])
            if numpy.max(numpy.abs(C)) < tol:
                break
            dQ0 = c2 * (2.0 * a_spin**2 * E * g00 + 2.0 * Lz * g30 / s2)
            dQ3 = c2 * (2.0 * a_spin**2 * E * g03 + 2.0 * Lz * g33 / s2)
            J = numpy.array([
                2.0 * p,
                [-g00, 0.0, 0.0, -g03],
                [ g30, 0.0, 0.0,  g33],
                [ dQ0, 0.0, 2.0 * g22**2 * u[2], dQ3],
            ])
            u += numpy.linalg.solve(J, -C)

        return numpy.concatenate([q, u])

    def get_path(self, initial_conditions, taus, integrator="Radau",
                 adaptative=True, steps_per_period=100):
        """
        Compute a geodesic and return it as a coordinate object.

        Integrates the geodesic in the metric's valid coordinate system,
        optionally using SciPy ``solve_ivp``, a Kerr ``Radau2`` specialized
        integrator, or a Yoshida symplectic scheme. The result is wrapped in
        the appropriate registered coordinate class and may be converted back
        to the caller's original system.

        Parameters
        ----------
        initial_conditions : CoordinateSystem
            Initial conditions in an arbitrary registered coordinate system;
            converted internally to ``self.valid_coordinate`` when needed.
        taus : array_like
            Strictly increasing proper-time samples spanning the integration
            interval. At least two values are required (see
            :meth:`_get_path_from_4state_vector`).
        integrator : str, optional
            ``"Radau"``, ``"DOP853"``, another ``solve_ivp`` method name,
            ``"Yoshida6"``, or ``"Radau2"`` (Kerr Boyer-Lindquist only). Default
            is ``"Radau"``.
        adaptative : bool, optional
            If ``True``, the integrator chooses its own output times (dense
            trajectory for SciPy; full internal steps for Yoshida6/Radau2). If
            ``False``, the solution is interpolated to ``taus`` (cubic
            interpolation for Yoshida6/Radau2; ``t_eval`` for SciPy).
        steps_per_period : int, optional
            Used with ``"Yoshida6"`` and ``"Radau2"`` when an orbital period is
            available: steps per orbit, scaling the total step count. Also used
            as a fallback scale for Radau2 when no period is known. Default is
            ``100``.

        Returns
        -------
        CoordinateSystem
            Trajectory in native coordinates, or converted to the original
            coordinate name when that differs from the metric's chart (except
            when the original system is ``"OrbitalElements"``, which is not
            converted back automatically).

        Raises
        ------
        ValueError
            If ``taus`` has fewer than two points, or if ``integrator`` is
            ``"Radau2"`` while ``self.valid_coordinate`` is not
            ``"BoyerLindquist"``.

        Examples
        --------
        >>> import numpy as np
        >>> # taus = np.linspace(0.0, 1.0, 51)
        >>> # path = Geodesic(metric).get_path(initial_conditions, taus)  # doctest: +SKIP
        """
        original_coordinate = initial_conditions.name_metric
        original_kwargs = initial_conditions.kwargs

        # Grab period before any coordinate conversion
        period = None
        if integrator.lower() in ("yoshida6", "radau2", "mino"):
            if hasattr(initial_conditions, "_get_period"):
                period = initial_conditions._get_period()

        if initial_conditions.name_metric != self.valid_coordinate:
            initial_conditions = initial_conditions.convert_to(
                self.valid_coordinate, **self.metric.kwargs
            )

        ys0 = self.metric.get_4state_vector(initial_conditions)

        sol = self._get_path_from_4state_vector(
            ys0, taus, integrator=integrator, adaptative=adaptative,
            steps_per_period=steps_per_period, period=period,
        )

        dxs_dt = self.metric.get_dxs_dt_from_4velocity(sol[4:])
        ys = coordinate_systems[initial_conditions.name_metric](
            sol[:4], vels=dxs_dt[1:], from_dxs_dt=True, **initial_conditions.kwargs
        )

        if (original_coordinate != initial_conditions.name_metric
                and original_coordinate != "OrbitalElements"):
            ys = ys.convert_to(original_coordinate, **original_kwargs)

        return ys

    def _get_path_from_4state_vector(
        self, ys0, taus, integrator="Radau", adaptative=True,
        steps_per_period=100, period=None,
    ):
        """
        Integrate the 8-dimensional geodesic state from a four-state vector.

        Dispatches to Yoshida6, Kerr Radau2, or generic ``solve_ivp`` based on
        ``integrator``. When ``adaptative`` is ``False``, Yoshida6 and Radau2
        results are interpolated to the requested ``taus``.

        Parameters
        ----------
        ys0 : ndarray, shape (8,)
            Initial state ``[q^\\mu, u^\\mu]`` in the metric's valid coordinates.
        taus : array_like
            Monotonic proper-time grid; must contain at least two values.
        integrator : str, optional
            Same meaning as in :meth:`get_path`. Default is ``"Radau"``.
        adaptative : bool, optional
            If ``True``, return the integrator's native time sampling. If
            ``False``, interpolate (Yoshida6/Radau2) or use ``t_eval`` (SciPy).
        steps_per_period : int, optional
            Steps per orbital period for Yoshida6/Radau2 when ``period`` is
            positive; otherwise a fallback step count is used.
        period : float or None, optional
            Orbital period from ``initial_conditions._get_period()`` when
            available; used to set the total number of steps for Yoshida6 and
            Radau2.

        Returns
        -------
        ndarray, shape (8, N)
            Integrated trajectory. Column :math:`k` is the state at the
            :math:`k`-th output time (or interpolated time when
            ``adaptative`` is ``False`` for Yoshida6/Radau2).

        Raises
        ------
        ValueError
            If ``len(taus) < 2`` or if ``integrator`` is ``"radau2"`` but the
            metric is not in Boyer-Lindquist coordinates.

        Notes
        -----
        If SciPy integration fails (``sol.status == -1``), a warning is printed
        and the partial solution in ``sol.y`` is still returned.

        Examples
        --------
        >>> import numpy as np
        >>> # y_traj = Geodesic(metric)._get_path_from_4state_vector(y0, np.linspace(0,1,11))  # doctest: +SKIP
        """
        if len(taus) < 2:
            raise ValueError("taus must contain at least 2 values.")

        taus_np = numpy.asarray(taus, dtype=float)
        t_span = (taus_np[0], taus_np[-1])
        span = t_span[1] - t_span[0]

        if integrator.lower() == "yoshida6":
            from .integrators.kerr.yoshida6 import Yoshida6Integrator

            if period is not None and period > 0:
                # steps proportional to the number of orbital periods
                n_periods_float = span / period
                n_steps = max(int(numpy.ceil(n_periods_float * steps_per_period)), 100)
            else:
                # fallback: 50 steps per output point
                n_steps = max(len(taus_np) * 50, 1000)

            result = Yoshida6Integrator(self.metric).integrate(
                ys0, t_span, n_steps, stride=1
            )

            if adaptative:
                return result.y  # (8, N_steps+1)

            # Interpolate to the requested taus
            out = numpy.zeros((8, len(taus_np)))
            for i in range(8):
                f = interp1d(result.t, result.y[i], kind="cubic")
                out[i] = f(taus_np)
            return out

        if integrator.lower() == "radau2":
            if self.valid_coordinate != "BoyerLindquist":
                raise ValueError(
                    "Radau2 requires a Kerr metric (Boyer-Lindquist coordinates)."
                )
            from .integrators.kerr.radau import _integrate_kerr_radau2

            g0 = self.metric.metric(ys0[:4])
            p0 = g0 @ ys0[4:]
            E0 = -p0[0]
            Lz0 = p0[3]
            p_th0 = p0[2]
            theta0 = ys0[2]
            c2 = numpy.cos(theta0) ** 2
            s2 = numpy.sin(theta0) ** 2
            Q0 = p_th0 ** 2 + c2 * (
                self.metric.a ** 2 * (1.0 - E0 ** 2) + Lz0 ** 2 / s2
            )

            if period is not None and period > 0:
                n_periods_float = span / period
                n_steps = max(int(numpy.ceil(n_periods_float * steps_per_period)), 100)
            else:
                n_steps = max(len(taus_np) * steps_per_period, 1000)

            result = _integrate_kerr_radau2(
                self.metric.R_s, self.metric.a,
                ys0, t_span, n_steps,
                E0=E0, Lz0=Lz0, Q0=Q0,
            )

            if adaptative:
                return result.y

            # Interpolate to requested taus
            y_out = numpy.zeros((8, len(taus_np)))
            for i in range(8):
                f = interp1d(result.t, result.y[i], kind="cubic")
                y_out[i] = f(taus_np)
            return y_out

        if integrator.lower() == "mino":
            if self.valid_coordinate != "BoyerLindquist":
                raise ValueError(
                    "Mino requires a Kerr metric (Boyer-Lindquist coordinates)."
                )
            from .integrators.kerr.mino import _integrate_kerr_mino

            g0    = self.metric.metric(ys0[:4])
            p0    = g0 @ ys0[4:]
            E0    = -p0[0]
            Lz0   = p0[3]
            p_th0 = p0[2]
            theta0 = ys0[2]
            c2 = numpy.cos(theta0) ** 2
            s2 = numpy.sin(theta0) ** 2
            Q0 = p_th0 ** 2 + c2 * (
                self.metric.a ** 2 * (1.0 - E0 ** 2) + Lz0 ** 2 / s2
            )

            if period is not None and period > 0:
                n_periods_float = span / period
                n_steps = max(int(numpy.ceil(n_periods_float * steps_per_period)), 100)
            else:
                n_steps = max(len(taus_np) * steps_per_period, 1000)

            result = _integrate_kerr_mino(
                self.metric.R_s, self.metric.a,
                ys0, t_span, n_steps,
                E0=E0, Lz0=Lz0, Q0=Q0,
            )

            if adaptative:
                return result.y

            y_out = numpy.zeros((8, len(taus_np)))
            for i in range(8):
                f = interp1d(result.t, result.y[i], kind="cubic")
                y_out[i] = f(taus_np)
            return y_out

        # scipy integrators
        if adaptative:
            sol = solve_ivp(self.model_geodesic, t_span, ys0, method=integrator)
        else:
            sol = solve_ivp(
                self.model_geodesic, t_span, ys0, t_eval=taus_np, method=integrator
            )

        if sol.status == -1:
            print("WARNING: Integration failed.")

        return sol.y

    def get_path_periodic(
        self,
        initial_conditions,
        n_periods=500,
        steps_per_period=100,
        output_per_period=10,
    ):
        """
        Long-horizon integration with the Yoshida 6th-order symplectic scheme.

        Integrates for ``n_periods`` orbital periods using a step count derived
        from ``steps_per_period`` and sub-samples the trajectory with ``stride``
        so that roughly ``output_per_period`` points are kept per period.

        Parameters
        ----------
        initial_conditions : OrbitalElements
            Initial conditions; must provide :meth:`_get_period` for the orbital
            period used to set the time span and step scaling.
        n_periods : int, optional
            Number of orbital periods to integrate. Default is ``500``.
        steps_per_period : int, optional
            Yoshida steps per period (accuracy control). Relative energy error per
            orbit scales roughly as :math:`(2\\pi / N)^6` with :math:`N` equal to
            ``steps_per_period``. Default is ``100``.
        output_per_period : int, optional
            Target number of stored output points per period (via ``stride``).
            Default is ``10``.

        Returns
        -------
        CoordinateSystem
            Path expressed in the metric's native coordinate system, optionally
            converted back from ``OrbitalElements`` to the original chart.

        Examples
        --------
        >>> # path = Geodesic(metric).get_path_periodic(orbital_ic, n_periods=100)  # doctest: +SKIP
        """
        from .integrators.kerr.yoshida6 import Yoshida6Integrator

        period = initial_conditions._get_period()

        original_coordinate = initial_conditions.name_metric
        original_kwargs = initial_conditions.kwargs

        if initial_conditions.name_metric != self.valid_coordinate:
            initial_conditions = initial_conditions.convert_to(
                self.valid_coordinate, **self.metric.kwargs
            )

        state0 = self.metric.get_4state_vector(initial_conditions)

        tau_end = float(n_periods) * period
        n_steps = n_periods * steps_per_period
        stride = max(1, n_steps // (n_periods * output_per_period))

        result = Yoshida6Integrator(self.metric).integrate(
            state0, (0.0, tau_end), n_steps, stride=stride
        )

        sol = result.y  # (8, N_out)

        dxs_dt = self.metric.get_dxs_dt_from_4velocity(sol[4:])
        ys = coordinate_systems[initial_conditions.name_metric](
            sol[:4], vels=dxs_dt[1:], from_dxs_dt=True, **initial_conditions.kwargs
        )

        if (original_coordinate != initial_conditions.name_metric
                and original_coordinate != "OrbitalElements"):
            ys = ys.convert_to(original_coordinate, **original_kwargs)

        return ys
