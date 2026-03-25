import numpy
from itertools import product
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

from ..coordinates import coordinate_systems


class Geodesic:
    def __init__(self, metric):
        self.metric = metric
        self.valid_coordinate = self.metric.valid_coordinate

    def model_geodesic(self, tau, ys0):
        xs0 = ys0[:4]
        us0 = ys0[4:]

        as_ = numpy.zeros(4)

        chris = self.metric.get_christoffel_symbols(xs0)

        for sigma, mu, nu in product(range(4), repeat=3):
            as_[sigma] -= chris[sigma, mu, nu] * us0[mu] * us0[nu]

        return numpy.concatenate([us0, as_])

    def _project_constraints(self, y, E0, Lz0, Q0, tol=1e-12, max_iter=20):
        """
        Newton-project u^μ onto {C1=0, C2=0, C3=0, C4=0} for Kerr geodesics.

        Constraints:
          C1 = g_{μν} u^μ u^ν − 1  (normalization)
          C2 = −g_{0μ} u^μ − E0     (energy)
          C3 =  g_{3μ} u^μ − Lz0   (angular momentum)
          C4 = p_θ² + cos²θ[a²(1−E²) + Lz²/sin²θ] − Q0  (Carter)

        Parameters
        ----------
        y           : ndarray(8,) — [q^0..q^3, u^0..u^3] in Boyer-Lindquist
        E0, Lz0, Q0 : conserved quantities from initial conditions
        tol         : convergence tolerance on max|C_i|

        Returns
        -------
        ndarray(8,) with corrected u^μ
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
        Returns the geodesic path for a test particle in the given metric.

        Parameters
        ----------
        initial_conditions : CoordinateSystem
            Initial conditions in a specified coordinate system.
        taus : list
            List of proper time values that define the integration span.
        integrator : str
            "Radau", "DOP853" (or any scipy method), or "Yoshida6".
        adaptative : bool
            If True, the integrator chooses its own output times.
            If False, the solution is evaluated/interpolated at `taus`.
        steps_per_period : int
            Only used with Yoshida6. Number of Yoshida steps per orbital
            period. Requires initial_conditions with a _get_period() method
            (e.g. OrbitalElements). Default 100.
        """
        original_coordinate = initial_conditions.name_metric
        original_kwargs = initial_conditions.kwargs

        # Grab period before any coordinate conversion
        period = None
        if integrator.lower() in ("yoshida6", "radau2"):
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
        if len(taus) < 2:
            raise ValueError("taus must contain at least 2 values.")

        taus_np = numpy.asarray(taus, dtype=float)
        t_span = (taus_np[0], taus_np[-1])
        span = t_span[1] - t_span[0]

        if integrator.lower() == "yoshida6":
            from .integrators.yoshida6 import Yoshida6Integrator

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
                    "Radau2 requires a Kerr metric (BoyerLindquist coordinates)."
                )
            from .integrators.radau import _integrate_kerr_radau2

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
        Long-term integration using Yoshida 6th-order symplectic integrator.

        Parameters
        ----------
        initial_conditions : OrbitalElements
            Initial conditions. Period is auto-computed via _get_period().
        n_periods : int
            Number of orbital periods to integrate.
        steps_per_period : int
            Yoshida steps per period (controls accuracy, default 100).
            Relative energy error per orbit ≈ (2π/steps_per_period)^6.
        output_per_period : int
            Output points stored per period (default 10).

        Returns
        -------
        CoordinateSystem
            Path in the metric's native coordinate system.
        """
        from .integrators.yoshida6 import Yoshida6Integrator

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
