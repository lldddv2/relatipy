"""Private orbit integrated through the forced Kerr layer.

:class:`ForcedOrbit` is an :class:`~relatipy.geodesic.Orbit` whose native
integration calls ``relatipy._force_core`` with one bound external force
instead of ``relatipy._core``. Initial states, sampled solutions,
interpolation and coordinate views are the homogeneous ones; only the
equation of motion differs. The orbit's metric is the background
:class:`~relatipy.metrics.Kerr`, which solutions keep for interpolation and
plotting.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from astropy import units as u

from .exceptions import IntegrationError
from .orbit import Orbit

if TYPE_CHECKING:
    from ..forces._bound import BoundForce
    from ..metrics.kerr import Kerr
    from ._options import IntegrationOptions
    from .solution import Solution
    from .state import InitialState

_PROJECTION_MESSAGE = (
    "projection_radau is not available with an external force: its projector "
    "restores Carter's constant, which a force does not conserve; use radau, "
    "dop853 or dp45"
)


def _reject_projection(method: object) -> None:
    """Raise ``ValueError`` for the projected Radau method."""
    if isinstance(method, str) and method == "projection_radau":
        raise ValueError(_PROJECTION_MESSAGE)


class ForcedOrbit(Orbit):
    """Orbit whose integrations include one bound external force.

    Created by :meth:`relatipy.metrics.PerturbedKerr.orbit`; not a public
    entry point.
    """

    __slots__ = ("_bound",)

    def __init__(
        self,
        metric: Kerr,
        initial_state: InitialState,
        initial_y: np.ndarray,
        length_unit: u.UnitBase,
        time_unit: u.UnitBase,
        input_family: str,
        bound: BoundForce | None = None,
    ) -> None:
        """Store the orbit data and the bound force it integrates with."""
        super().__init__(
            metric, initial_state, initial_y, length_unit, time_unit, input_family
        )
        self._bound = bound

    def copy(self) -> ForcedOrbit:
        """Return an independently evolving copy with the same force.

        Returns
        -------
        Orbit
            Copy retaining the saved initial state, current state and force.
        """
        clone = ForcedOrbit(
            self._metric,
            self._initial,
            self._initial_y,
            self._length_unit,
            self._time_unit,
            self._input_family,
            self._bound,
        )
        clone._current = self._current
        clone._current_y = np.frombuffer(self._current_y.tobytes(), dtype=np.float64)
        return clone

    def _call_native(
        self,
        y0: np.ndarray,
        tau0: float,
        tau1: float,
        options: IntegrationOptions,
        *,
        store_steps: bool,
    ) -> tuple[np.ndarray, np.ndarray, dict, int]:
        """Integrate one normalized span through the forced native solve."""
        from .. import _force_core

        _reject_projection(options.method)
        bound = self._bound
        if bound is None:
            raise IntegrationError("forced orbit has no bound force")
        scale = self._metric._time_scale
        result = _force_core.integrate_kerr(
            self._metric.spin,
            bound.callback,
            bound.model,
            np.ascontiguousarray(y0),
            tau0,
            tau1,
            options.method,
            options.rtol,
            options.atol,
            None if options.first_step is None
            else (options.first_step / scale).to_value(u.one),
            None if options.max_step is None
            else (options.max_step / scale).to_value(u.one),
            store_steps,
        )
        taus, states, stats, status = result
        taus = np.asarray(taus, dtype=np.float64)
        states = np.asarray(states, dtype=np.float64)
        if (
            taus.ndim != 1
            or states.shape != (len(taus), 8)
            or len(taus) == 0
            or not np.all(np.isfinite(taus))
            or not np.all(np.isfinite(states))
        ):
            raise IntegrationError("native integration returned invalid states")
        if status not in (-1, 0, 1):
            raise IntegrationError("native integration returned invalid status")
        return taus, states, stats, status

    def integrate(self, tau: u.Quantity, *, method: str = "radau", **options) -> None:
        """Advance under the force; see :meth:`relatipy.geodesic.Orbit.integrate`.

        Raises
        ------
        ValueError
            Also for ``method="projection_radau"``.
        """
        _reject_projection(method)
        return super().integrate(tau, method=method, **options)

    def solve(self, *, method: str = "radau", **options) -> Solution:
        """Solve under the force; see :meth:`relatipy.geodesic.Orbit.solve`.

        Raises
        ------
        ValueError
            Also for ``method="projection_radau"``.
        """
        _reject_projection(method)
        return super().solve(method=method, **options)
