from curses import echo
from numpy import array, zeros_like, concatenate
from itertools import product

from ..constants import _c
from ..utils.dimensions import validator

class CoordinateBase:
    def __init__(
        self, xs, vels=None, from_dxs_dt=False, system_name="CoordinateBase", **kwargs
    ):
        self.name_metric = system_name

        xs = validator.validate_vector(xs)
        vels = validator.validate_vector(vels)

        self.kwargs = kwargs
        for key, value in kwargs.items():
            setattr(self, key, value)

        if len(xs) != 4:
            raise ValueError(
                "xs must be a list or array of length 4 representing (t, x1, x2, x3), not {}".format(
                    len(xs)
                )
            )

        self.xs = array(xs)

        if vels is None:
            self.vels = zeros_like(xs[1:])
            self.dxs_dt = zeros_like(xs[1:])
        else:
            if not from_dxs_dt:
                self.vs = vels
                self.dxs_dt = self._get_dxs_dt_from_vs()
            else:
                self.dxs_dt = vels
                self.vs = self._get_vs_from_dxs_dt()

        self.state_vector = concatenate((self.xs, self.vs))

    def convert_to_cartesian(self):
        from .cartesian import Cartesian

        xs_p, vs_p = self._convert_to_cartesian(self.xs, self.vs, **self.kwargs)
        coordinate = Cartesian(xs_p, vels=vs_p, from_dxs_dt=False)
        return coordinate

    def convert_to(self, target_system, **kwargs):
        from . import coordinate_systems

        if target_system not in coordinate_systems:
            raise ValueError(
                f"Unsupported target coordinate system: {target_system}. Supported systems are: {list(coordinate_systems.keys())}"
            )

        cartesian = self.convert_to_cartesian()
        if target_system == "Cartesian":
            return cartesian

        target_class = coordinate_systems[target_system]
        return target_class._convert_from_cartesian(cartesian.xs, cartesian.vs, **kwargs)

    # Numerical methods
    def _get_vs_from_dxs_dt(self):
        raise NotImplementedError(
            "Subclasses must implement get_vs_from_dxs_dt method."
        )

    def _get_dxs_dt_from_vs(self):
        raise NotImplementedError(
            "Subclasses must implement get_dxs_dt_from_vs method."
        )

    # Statics & numerical methods
    @staticmethod
    def _convert_to_cartesian(xs, vs, **kwargs):
        raise NotImplementedError(
            "Subclasses must implement _convert_to_cartesian static method."
        )

    @staticmethod
    def _convert_from_cartesian(xs_p, vs_p, **kwargs):
        raise NotImplementedError(
            "Subclasses must implement _convert_from_cartesian static method."
        )

    # Magic methods
    def __getitem__(self, index):
        # Permite indexación simple y múltiple (tupla)
        if isinstance(index, tuple):
            return self.state_vector[index]
        if index > 6 or index < 0:
            raise IndexError(
                "Index out of range. Must be between 0 and 6. [x0,x1,x2,x3,v1,v2,v3]"
            )
        return self.state_vector[index]

    def _get_ds_dtau(self, metric):
        """
        Calcula g_{mu nu} u^mu u^nu. Para una geodésica tipo tiempo debe ser -c^2.
        
        Parameters
        ----------
        metric : callable
            Función metric(xs, **kwargs) que devuelve g_{mu nu} (4x4).
        
        Returns
        -------
        float
            g_{mu nu} u^mu u^nu (debería ser -c^2)
        """
        from numpy import einsum, ones

        g = metric.metric(self.xs)  # (N, 4, 4)
        
        u = ones((4, len(self.dxs_dt[0])))
        u[1:, :] = self.dxs_dt  # (4, N)

        return einsum('ijn,in,jn->n', g, u, u)   # (N,)