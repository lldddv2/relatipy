"""Read-only coordinate-time solutions for Kerr null geodesics."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from astropy import units as u
from scipy.interpolate import CubicSpline

from . import _native as native
from ._null_state import NullState, NullTermination, delegate_null_properties
from .integration import IntegrationInfo
from .interpolation import _exact_indices, _time_query
from .solution import _validate_index


def _restore_null_solution(
    state: NullState, integration: IntegrationInfo, status: int,
    message: str, termination: NullTermination | None,
) -> NullSolution:
    """Restore a solution through validation, including diagnostic arrays."""
    info = IntegrationInfo(
        method=integration.method, rtol=integration.rtol, atol=integration.atol,
        max_step=integration.max_step, first_step=integration.first_step,
        n_steps=integration.n_steps, nfev=integration.nfev,
    )
    return NullSolution(
        state=state, integration=info, status=status, message=message,
        termination=termination,
    )


@dataclass(frozen=True, slots=True, init=False, eq=False)
class NullSolution:
    """Read-only sampled result from a future null-geodesic integration.

    Parameters
    ----------
    state : NullState
        Nonempty series of stored states with strictly increasing ``t``.
    integration : IntegrationInfo
        Effective method, tolerances, and integration counters. Step controls
        are ``None`` because the affine parameter is internal.
    status : int
        ``0`` for reaching the requested time, ``1`` for a terminal event,
        or ``-1`` for a numerical failure.
    message : str
        Human-readable outcome description.
    termination : NullTermination or None, optional
        Last valid terminal state, required exactly when ``status == 1``.

    Raises
    ------
    TypeError
        If an argument has the wrong type.
    ValueError
        If samples, status, or termination information are inconsistent.

    Notes
    -----
    Indexing selects samples. :meth:`at` interpolates Cartesian position and
    coordinate velocity separately in coordinate time; C converts the result
    to the stored BL chart. Interpolation does not enforce the null constraint.
    Observer frames, redshift, past integration, and the direction root choice
    in the ergoregion remain pending. Proper time and affine time are absent.

    Examples
    --------
    >>> from astropy import units as u
    >>> from astropy.constants import c
    >>> from relatipy import Kerr
    >>> metric = Kerr(mass=1 * u.M_sun, spin=0)
    >>> ray = metric.null(R=10 * metric.r_g, Theta=90 * u.deg,
    ...                   Phi=0 * u.rad, vR=1 * u.km / u.s)
    >>> solution = ray.solve(t_eval=[0, 0.01, 0.02] * (metric.r_g / c),
    ...                      method="dp45")
    >>> len(solution), solution.success
    (3, True)
    >>> solution[0].t
    <Quantity 0. s>
    >>> bool(solution.at(solution.t[1]).R == solution.R[1])
    True
    """

    _state: NullState
    integration: IntegrationInfo
    status: int
    message: str
    termination: NullTermination | None

    def __init__(
        self, *, state: NullState, integration: IntegrationInfo, status: int,
        message: str, termination: NullTermination | None = None,
    ) -> None:
        if not isinstance(state, NullState):
            raise TypeError("state must be a NullState")
        if state._scalar or state.t.ndim != 1:
            raise ValueError("NullSolution state must be a one-dimensional state series")
        if len(state.t) == 0:
            raise ValueError("NullSolution must contain at least one stored sample")
        if not np.all(np.isfinite(state.t.value)) or np.any(np.diff(state.t.value) <= 0):
            raise ValueError("NullSolution coordinate times must be finite and strictly increasing")
        if not isinstance(integration, IntegrationInfo):
            raise TypeError("integration must be an IntegrationInfo")
        if isinstance(status, (bool, np.bool_)) or status not in (-1, 0, 1):
            raise ValueError("status must be one of -1, 0 or 1")
        if not isinstance(message, str):
            raise TypeError("message must be a string")
        if status == 1:
            if not isinstance(termination, NullTermination):
                raise ValueError("termination is required when status is 1")
            if termination.reason not in ("horizon", "escape"):
                raise ValueError("null termination reason must be horizon or escape")
        elif termination is not None:
            raise ValueError("termination must be None unless status is 1")
        for name, value in (
            ("_state", state), ("integration", integration), ("status", int(status)),
            ("message", message), ("termination", termination),
        ):
            object.__setattr__(self, name, value)

    def __reduce_ex__(self, protocol: int):
        """Restore validated read-only arrays across copy and pickle."""
        return _restore_null_solution, (
            self._state, self.integration, self.status, self.message,
            self.termination,
        )

    def __len__(self) -> int:
        """Return the number of stored samples.

        Returns
        -------
        int
            Number of coordinate-time samples.
        """
        return len(self._state.t)

    def __getitem__(self, index: Any) -> NullState:
        """Select stored samples without interpolation.

        Parameters
        ----------
        index : int, slice, array-like of int, or array-like of bool
            Sample-axis selection; a boolean mask must match solution length.

        Returns
        -------
        NullState
            Immutable scalar state for integer indexing, or a state series.

        Raises
        ------
        IndexError
            If the index shape, dtype, mask length, or bounds are invalid.
        """
        return self._state._select(_validate_index(index, len(self)))

    @property
    def success(self) -> bool:
        """Report whether integration ended without numerical failure.

        Returns
        -------
        bool
            ``True`` for status ``0`` or ``1``.
        """
        return self.status in (0, 1)

    def at(self, t: u.Quantity) -> NullState:
        """Select or interpolate states at coordinate times.

        Parameters
        ----------
        t : astropy.units.Quantity
            Finite scalar or one-dimensional coordinate-time query, within
            the stored domain. Query order need not be increasing.

        Returns
        -------
        NullState
            Immutable scalar or series matching the query shape.

        Raises
        ------
        TypeError, astropy.units.UnitConversionError
            If the query lacks time units or has incompatible units.
        ValueError
            If the query shape or values are invalid, it is outside the
            domain, or C cannot reconstruct an interpolated state.

        Notes
        -----
        Exact samples copy the canonical stored row. Other rows use separate
        cubic splines for Cartesian position and coordinate velocity, with
        extrapolation disabled. C converts each row to BL; its affine scale
        is set to ``k^t = 1``. Stored azimuth selects the continuous branch.
        One stored sample supports only exact selection.

        Examples
        --------
        >>> from astropy import units as u
        >>> from astropy.constants import c
        >>> from relatipy import Kerr
        >>> metric = Kerr(mass=1 * u.M_sun, spin=0)
        >>> ray = metric.null(R=10 * metric.r_g, Theta=90 * u.deg,
        ...                   Phi=0 * u.rad, vR=1 * u.km / u.s)
        >>> solution = ray.solve(t_eval=[0, 0.01, 0.02] * (metric.r_g / c),
        ...                      method="dp45")
        >>> state = solution.at(solution.t[1] / 2)
        >>> state.t.shape, state.xyz.shape
        ((), (3,))
        """
        source = self._state
        query, scalar = _time_query(t, source._time_unit, "t")
        samples = np.asarray(source.t.to_value(source._time_unit))
        if np.any((query < samples[0]) | (query > samples[-1])):
            raise ValueError("t lies outside the stored solution domain")
        exact = _exact_indices(samples, query)
        if np.all(exact >= 0):
            return source._select(int(exact[0]) if scalar else exact)
        if len(samples) < 2:
            raise ValueError("a single stored sample only supports exact selection")

        from .. import _core

        matched = exact >= 0
        canonical = np.empty((len(query), 8), dtype=np.float64)
        canonical[matched] = source._canonical[exact[matched]]
        unknown = query[~matched]
        xyz = CubicSpline(
            samples, source.xyz.to_value(source._length_unit), axis=0,
            extrapolate=False,
        )(unknown)
        speed_unit = source._length_unit / source._time_unit
        vxyz = CubicSpline(
            samples, source.vxyz.to_value(speed_unit), axis=0,
            extrapolate=False,
        )(unknown)
        rows = np.empty((len(unknown), 7), dtype=np.float64)
        rows[:, 0] = (unknown * source._time_unit / source._time_scale).to_value(u.one)
        rows[:, 1:4] = (xyz * source._length_unit / source._length_scale).to_value(u.one)
        rows[:, 4:7] = (
            vxyz * speed_unit * source._time_scale / source._length_scale
        ).to_value(u.one)
        for index, row in zip(np.flatnonzero(~matched), rows):
            bl, status = _core.null_family_to_bl(source._spin, "cartesian", row)
            if status != 0:
                raise ValueError(
                    f"interpolated null state at query row {index} cannot be reconstructed "
                    f"({native.describe_reconstruction_status(status)})"
                )
            canonical[index, :4] = bl[:4]
            canonical[index, 4] = 1.0
            canonical[index, 5:] = bl[4:]
        reference = np.interp(
            unknown, samples, source.Phi.to_value(u.rad),
        )
        canonical[~matched, 3] += (2 * np.pi) * np.rint(
            (reference - canonical[~matched, 3]) / (2 * np.pi)
        )
        result = NullState(
            canonical=canonical, spin=source._spin,
            length_scale=source._length_scale, time_scale=source._time_scale,
            length_unit=source._length_unit, time_unit=source._time_unit,
            t_physical=query * source._time_unit,
        )
        return result._select(0) if scalar else result


delegate_null_properties(NullSolution, "_state")
