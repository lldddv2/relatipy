import numpy
from itertools import product
from scipy.integrate import solve_ivp

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

    def get_path(self, initial_conditions, taus):
        """
        Returns the geodesic equations for a test particle in the given metric.

        Parameters
        ----------
        initial_conditions : CoordinateSystem
            Initial conditions in a specified coordinate system.
        taus : list
            List of proper time values where the solution is evaluated.
        """
        original_coordinate = initial_conditions.name_metric
        original_kwargs = initial_conditions.kwargs
        if initial_conditions.name_metric != self.valid_coordinate:
            initial_conditions = initial_conditions.convert_to(self.valid_coordinate, **self.metric.kwargs)

        ys0 = self.metric.get_4state_vector(initial_conditions)
        sol = self._get_path_from_4state_vector(ys0, taus)

        dxs_dt = self.metric.get_dxs_dt_from_4velocity(sol[4:])
        ys = coordinate_systems[initial_conditions.name_metric](
            sol[:4], vels=dxs_dt[1:], from_dxs_dt=True, **initial_conditions.kwargs
        )

        if original_coordinate != initial_conditions.name_metric:
            ys = ys.convert_to(original_coordinate, **original_kwargs)

        return ys

    def _get_path_from_4state_vector(self, ys0, taus):
        """
        Returns the geodesic equations for a test particle in the given metric.

        Parameters
        ----------
        ys0 : list
            Initial conditions [x0, x1, x2, x3, u0, u1, u2, u3].
        taus : list
            List of proper time values where the solution is evaluated.
        """
        t_span = (taus[0], taus[-1])

        sol = solve_ivp(self.model_geodesic, t_span, ys0, t_eval=taus, method="Radau")
        if sol.status == -1:
            print("WARNING: Integration failed.")

        return sol.y
