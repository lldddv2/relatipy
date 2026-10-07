"""Scalar Kerr null geodesics with coordinate time as the public clock."""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import numpy as np
from astropy import units as u

from .._validation import immutable_array, readonly_quantity
from . import _native as native
from ._null_state import (
    NullState, NullTermination, _NullIntegrationTerminated,
    delegate_null_properties,
)
from ._options import integration_options
from .exceptions import IntegrationError, IntegrationWarning
from .integration import IntegrationInfo
from .null_solution import NullSolution
from .orbit import (
    _AUTO_RTOL, _MAX_FIT_TOL, _MIN_FIT_RTOL, _classify_inputs,
    _initial_components_from_scales, _scalar_quantity, _state_scales,
)

if TYPE_CHECKING:
    from ..metrics.kerr import Kerr

_POSITION_KEYS = (("x", "y", "z"), ("r", "theta", "phi"), ("R", "Theta", "Phi"))
_VELOCITY_KEYS = ("vx", "vy", "vz", "vr", "vtheta", "vphi", "vR", "vTheta", "vPhi")
_CONSTANT_KEYS = ("b", "eta", "radial_sign", "polar_sign")


def _check_status(status: int, *, direction: bool = False) -> None:
    """Translate native initial-state rejection without selecting a root."""
    if status == 0:
        return
    if direction and status == 3:
        raise ValueError("ergoregion has two future null roots; root choice is pending")
    if status == 8:
        raise ValueError("no future null tangent")
    raise ValueError(native.describe_reconstruction_status(status))


def _sign(value: object, name: str) -> int:
    """Validate an integer branch sign, excluding booleans."""
    if isinstance(value, (bool, np.bool_, str)):
        raise TypeError(f"{name} must be an integer sign")
    if not isinstance(value, (int, np.integer)) or value not in (-1, 1):
        raise ValueError(f"{name} sign must be an integer -1 or +1")
    return int(value)


def _resolve_options(y0: np.ndarray, *, method: object, rtol: object, atol: object):
    """Resolve the same automatic tolerance scales as the timelike frontend."""
    if not isinstance(method, str) or method not in ("radau", "dop853", "dp45"):
        raise ValueError("null method must be one of radau, dop853, dp45")
    options = integration_options(
        method=method, rtol=_AUTO_RTOL if rtol is None else rtol,
        atol=0.0 if atol is None else atol, first_step=None, max_step=None,
    )
    if isinstance(options.atol, np.ndarray) and options.atol.shape != (8,):
        raise ValueError("atol array must have shape (8,)")
    scales = _state_scales(y0)
    if atol is None:
        options = options._replace(atol=immutable_array(options.rtol * scales))
    problems = []
    if rtol is not None and options.rtol > _MAX_FIT_TOL:
        problems.append(f"rtol={options.rtol:g} is looser than {_MAX_FIT_TOL:g}")
    if rtol is not None and 0 < options.rtol < _MIN_FIT_RTOL:
        problems.append(f"rtol={options.rtol:g} is below {_MIN_FIT_RTOL:.1g}, unattainable in double precision")
    if atol is not None:
        ratio = np.broadcast_to(options.atol, (8,))[:5] / scales[:5]
        names = ("t", "R", "Theta", "Phi", "k^t")
        loose = [names[i] for i in np.flatnonzero(ratio > _MAX_FIT_TOL)]
        if loose:
            problems.append(f"atol exceeds {_MAX_FIT_TOL:g} times the characteristic scale of {', '.join(loose)}")
    if problems:
        warnings.warn(
            "unfit integration tolerances for this null geodesic: "
            + "; ".join(problems)
            + ". The trajectory may be inaccurate or the integration may fail; "
            "omit rtol and atol to choose them automatically",
            IntegrationWarning, stacklevel=3,
        )
    return options


def _failure_message(stats: dict) -> str:
    """Describe a numerical failure using the shared native decoder."""
    return "Kerr null integration failed before the target time (" + native.describe_integrator_status(stats.get("native_status")) + ")"


def _restore_null(
    cls: type, metric: Kerr, initial: NullState, initial_y: np.ndarray,
    current: NullState, current_y: np.ndarray, length_unit: u.UnitBase,
    time_unit: u.UnitBase, b: u.Quantity, eta: u.Quantity,
) -> Null:
    """Rebuild photon data through validation without sharing mutable caches."""
    def restore_state(state: NullState) -> NullState:
        return NullState(
            canonical=state._canonical, spin=metric.spin,
            length_scale=metric._length_scale, time_scale=metric._time_scale,
            length_unit=length_unit, time_unit=time_unit,
            t_physical=state.t, scalar=True, xyz_physical=state._xyz_physical,
        )

    photon = cls(metric, restore_state(initial), initial_y, length_unit, time_unit)
    values = np.asarray(current_y, dtype=np.float64)
    if values.shape != (8,) or not np.all(np.isfinite(values)):
        raise ValueError("current native null state must be finite with shape (8,)")
    photon._current = restore_state(current)
    photon._current_y = immutable_array(values)
    photon._b = readonly_quantity(_scalar_quantity(b, u.m, "b"), u.m, "b", ndim=(0,))
    photon._eta = readonly_quantity(
        _scalar_quantity(eta, u.m**2, "eta"), u.m**2, "eta", ndim=(0,),
    )
    return photon


class Null:
    """One mutable photon state with immutable saved initial conditions.

    Parameters
    ----------
    metric, initial_state, initial_y, length_unit, time_unit
        Internal validated construction data. Use :meth:`Kerr.null
        <relatipy.metrics.Kerr.null>` to create a photon.

    Notes
    -----
    Public velocities are coordinate derivatives with respect to ``t``.
    The affine parameter remains internal. Past integration, observer frames,
    redshift, and selection between two ergoregion roots are pending.

    Examples
    --------
    >>> from astropy import units as u
    >>> from astropy.constants import c
    >>> from relatipy import Kerr
    >>> bh = Kerr(mass=1 * u.Msun, spin=0)
    >>> photon = bh.null(x=10 * bh.r_g, vx=c)
    >>> photon.integrate(1e-8 * u.s, method="dp45")
    >>> bool(photon.t == 1e-8 * u.s)
    True
    """

    __slots__ = ("_metric", "_initial", "_current", "_initial_y", "_current_y",
                 "_length_unit", "_time_unit", "_b", "_eta")

    def __init__(self, metric: Kerr, initial_state: NullState, initial_y: np.ndarray,
                 length_unit: u.UnitBase, time_unit: u.UnitBase) -> None:
        """Store validated scalar initial data and native invariants."""
        from .. import _core

        self._metric = metric
        self._initial = self._current = initial_state
        self._initial_y = immutable_array(np.asarray(initial_y, dtype=np.float64))
        self._current_y = self._initial_y
        self._length_unit, self._time_unit = length_unit, time_unit
        invariants, status = _core.null_invariants(metric.spin, self._initial_y)
        _check_status(status)
        self._b = readonly_quantity(
            (invariants["impact_parameter"] * metric._length_scale).to(length_unit), u.m, "b", ndim=(0,))
        self._eta = readonly_quantity(
            (invariants["eta"] * metric._length_scale**2).to(length_unit**2), u.m**2, "eta", ndim=(0,))

    def __reduce_ex__(self, protocol: int):
        """Restore immutable data and exact presentation times after copying."""
        return _restore_null, (
            type(self), self._metric, self._initial, self._initial_y,
            self._current, self._current_y, self._length_unit, self._time_unit,
            self._b, self._eta,
        )

    @classmethod
    def _from_kwargs(cls, metric: Kerr, **values: object) -> Null:
        """Validate one public family and delegate null construction to C."""
        from .. import _core

        if not any(values.get(key) is not None for keys in _POSITION_KEYS for key in keys):
            raise ValueError("one position family is required")
        has_constants = any(values.get(key) is not None for key in _CONSTANT_KEYS)
        has_velocity = any(values.get(key) is not None for key in _VELOCITY_KEYS)
        if has_constants and has_velocity:
            raise ValueError("provide either a direction or b, eta and signs")
        if has_constants and any(values.get(key) is None for key in _CONSTANT_KEYS):
            raise ValueError("constants input requires b, eta, radial_sign and polar_sign")
        # Reuse Orbit's family validation without introducing its element API.
        family_values = {key: values.get(key) for keys in _POSITION_KEYS for key in keys}
        family_values.update({key: values.get(key) for key in _VELOCITY_KEYS})
        family_values.update({key: None for key in ("a", "e", "inc", "Omega", "omega", "f")})
        family = _classify_inputs(family_values)
        t = _scalar_quantity(values["t"], u.s, "t")
        row, length_unit = _initial_components_from_scales(
            metric._length_scale, metric._time_scale, family, family_values, t)
        if has_constants:
            b = _scalar_quantity(values["b"], u.m, "b")
            eta = _scalar_quantity(values["eta"], u.m**2, "eta")
            radial_sign = _sign(values["radial_sign"], "radial_sign")
            polar_sign = _sign(values["polar_sign"], "polar_sign")
        bl, status = _core.null_family_to_bl(metric.spin, family, row)
        if status == 3:
            raise ValueError("position is outside the exterior chart or at the Kerr horizon")
        _check_status(status)
        if bl[1] <= _core.null_horizon_threshold(metric.spin):
            raise ValueError("initial radius must exceed the null horizon threshold")
        if has_constants:
            canonical, status = _core.null_state_from_constants(
                metric.spin, np.ascontiguousarray(bl[:4]),
                float((b / metric._length_scale).to_value(u.one)),
                float((eta / metric._length_scale**2).to_value(u.one)), radial_sign, polar_sign)
            _check_status(status)
        else:
            if not np.any(row[4:] != 0):
                raise ValueError("a nonzero spatial direction is required")
            canonical, status = _core.null_state_from_direction(metric.spin, family, row)
            _check_status(status, direction=True)
        canonical = np.asarray(canonical, dtype=np.float64)
        if canonical.shape != (8,) or not np.all(np.isfinite(canonical)):
            raise ValueError("native null initial state is invalid")
        # Preserve supplied Cartesian positions, including exact default zeros,
        # at the presentation boundary; C still determines all velocities.
        xyz_physical = None
        if family == "cartesian":
            xyz_physical = np.array([
                0.0 if values[key] is None else values[key].to_value(length_unit)
                for key in ("x", "y", "z")
            ]) * length_unit
        initial = NullState(
            canonical=canonical.reshape(1, 8), spin=metric.spin,
            length_scale=metric._length_scale, time_scale=metric._time_scale,
            length_unit=length_unit, time_unit=t.unit, t_physical=t, scalar=True,
            xyz_physical=xyz_physical)
        return cls(metric, initial, canonical, length_unit, t.unit)

    def _state(self, canonical: np.ndarray, *, scalar: bool, t_physical=None) -> NullState:
        """Wrap native rows with presentation units and optional exact times."""
        return NullState(
            canonical=np.asarray(canonical).reshape(-1, 8), spin=self._metric.spin,
            length_scale=self._metric._length_scale, time_scale=self._metric._time_scale,
            length_unit=self._length_unit, time_unit=self._time_unit,
            t_physical=t_physical, scalar=scalar)

    def _escape_radius(self, value: object, y0: np.ndarray) -> float:
        """Validate and normalize the optional escape threshold."""
        if value is None:
            return 0.0
        radius = _scalar_quantity(value, u.m, "r_escape")
        normalized = float((radius / self._metric._length_scale).to_value(u.one))
        if normalized <= y0[1]:
            raise ValueError("r_escape must exceed the starting BL radius")
        return normalized

    def _call_native(self, y0, end, evaluation, options, escape, *, store_steps):
        """Call the coarse native operation and validate its result shapes."""
        from .. import _core

        scale = self._metric._time_scale
        final_time = float((end / scale).to_value(u.one))
        sample_times = None
        if evaluation is not None:
            sample_times = np.array((evaluation / scale).to_value(u.one), dtype=np.float64)
            # Python already validated t0 <= t_eval <= t_final in physical
            # units; absorb unit-conversion roundoff at both ends for C.
            sample_times[0] = max(sample_times[0], float(y0[0]))
            sample_times[-1] = min(sample_times[-1], final_time)
        result = _core.integrate_kerr_null(
            self._metric.spin, y0, final_time, sample_times,
            escape, options.method, options.rtol, options.atol, store_steps)
        states, final, stats, status, reason = result
        states, final = np.asarray(states), np.asarray(final)
        if (states.ndim != 2 or states.shape[1] != 8 or final.shape != (8,)
                or not np.all(np.isfinite(states)) or not np.all(np.isfinite(final))
                or status not in (-1, 0, 1) or (status == 1 and reason not in ("horizon", "escape"))):
            raise IntegrationError("native null integration returned invalid states or status")
        return states, final, stats, status, reason

    @property
    def initial(self) -> NullState:
        """Return the immutable saved initial state.

        Returns
        -------
        NullState
            Scalar coordinate and coordinate-velocity views.
        """
        return self._initial

    @property
    def b(self) -> u.Quantity:
        """Return the initial impact parameter.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar length in the presentation unit.
        """
        return self._b

    @property
    def eta(self) -> u.Quantity:
        """Return the initial Carter parameter.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar squared length in the presentation unit squared.
        """
        return self._eta

    def reset(self) -> None:
        """Restore the exact initial state.

        Returns
        -------
        None
            This photon is reset without integration.
        """
        self._current, self._current_y = self._initial, self._initial_y

    def copy(self) -> Null:
        """Return an independently evolving copy.

        Returns
        -------
        Null
            Copy with independent state storage, invariant quantities and
            lazy view caches. The immutable source metric is shared.
        """
        restore, arguments = self.__reduce_ex__(4)
        return restore(*arguments)

    def integrate(self, t: u.Quantity, *, method="radau", rtol=None, atol=None,
                  r_escape=None) -> None:
        """Advance the current point to a future coordinate time.

        Parameters
        ----------
        t : astropy.units.Quantity
            Finite scalar absolute time greater than the current time.
        method : {"radau", "dop853", "dp45"}, optional
            Native method, default ``"radau"``.
        rtol : float or None, optional
            Relative tolerance; ``None`` selects ``1e-10``.
        atol : float, array_like or None, optional
            Scalar or shape ``(8,)`` normalized absolute tolerance. ``None``
            selects ``rtol * max(abs(y0), floors)`` using Orbit's floors.
            State order is ``(t/T0, R/r_g, Theta, Phi, k^t, k^R, k^Theta, k^Phi)``.
        r_escape : astropy.units.Quantity or None, optional
            Escape radius greater than the current BL radius; default disabled.

        Returns
        -------
        None
            The current point is updated to the last valid native state.

        Raises
        ------
        IntegrationTerminated
            Horizon or escape terminates integration; carries ``reason``,
            ``t``, ``state`` and ``termination``.
        IntegrationError
            A numerical failure occurs after retaining the last valid state.
        TypeError, astropy.units.UnitConversionError, ValueError
            Invalid time, units, method, tolerances, or escape radius.

        Warns
        -----
        IntegrationWarning
            Explicit tolerances meet the unfit criteria of architecture §3.3.1.
            The warning leaves the supplied tolerances unchanged.

        Notes
        -----
        Past integration and null projection are pending. No affine step
        controls are exposed. Tolerances control local error, not global error.
        """
        target = _scalar_quantity(t, u.s, "t")
        if target <= self._current.t:
            raise ValueError("integration requires a future coordinate time")
        options = _resolve_options(self._current_y, method=method, rtol=rtol, atol=atol)
        escape = self._escape_radius(r_escape, self._current_y)
        _, final, stats, status, reason = self._call_native(
            self._current_y, target, None, options, escape, store_steps=False)
        last = self._state(final, scalar=True, t_physical=target if status == 0 else None)
        self._current, self._current_y = last, immutable_array(final)
        if status == 1:
            raise _NullIntegrationTerminated(reason, last.t, last,
                                            f"Kerr null integration reached {reason}")
        if status == -1:
            raise IntegrationError(_failure_message(stats))

    def solve(self, *, t_span=None, t_eval=None, method="radau", rtol=None,
              atol=None, r_escape=None) -> NullSolution:
        """Solve from saved initial conditions without mutating this photon.

        Parameters
        ----------
        t_span : tuple or list of astropy.units.Quantity, optional
            Two scalar absolute times; start equals the saved initial time,
            end is strictly later. Either ``t_span`` or ``t_eval`` is required.
        t_eval : astropy.units.Quantity, optional
            Nonempty finite one-dimensional strictly increasing sample times
            within the span. Without a span, the last sample is the target.
        method : {"radau", "dop853", "dp45"}, optional
            Native method, default ``"radau"``.
        rtol : float or None, optional
            Relative local tolerance; default ``None`` selects ``1e-10``.
        atol : float, array_like or None, optional
            Scalar or shape ``(8,)`` normalized tolerance, in order
            ``(t/T0, R/r_g, Theta, Phi, k^t, k^R, k^Theta, k^Phi)``.
            ``None`` selects the automatic initial-state scales.
        r_escape : astropy.units.Quantity or None, optional
            Escape radius greater than the initial BL radius; default disabled.

        Returns
        -------
        NullSolution
            Read-only samples reached. Status is 0 at the target, 1 at a
            horizon or escape event, and -1 on numerical failure.

        Raises
        ------
        IntegrationTerminated, IntegrationError
            With ``t_eval``, an event or numerical failure occurs before the
            first requested sample. Without ``t_eval``, numerical failure
            returns a partial solution with at least the initial state and
            ``status=-1``.
        TypeError, astropy.units.UnitConversionError, ValueError
            Invalid time span, samples, units, or numerical options.

        Warns
        -----
        IntegrationWarning
            Explicit tolerances meet the criteria described in ``integrate``.

        Notes
        -----
        Without ``t_eval``, stores initial state and accepted steps. With
        ``t_eval``, C samples cubic Hermite interpolation in its internal
        affine parameter. Its interpolation error can dominate local error.
        Past integration, observer frames, redshift and the ergoregion root
        choice remain pending.

        Examples
        --------
        >>> from astropy import units as u
        >>> from astropy.constants import c
        >>> from relatipy import Kerr
        >>> photon = Kerr(mass=1 * u.Msun, spin=0).null(x=20 * u.km, vx=c)
        >>> solution = photon.solve(t_eval=[0, 1e-8, 2e-8] * u.s, method="dp45")
        >>> len(solution), solution.status, solution.success
        (3, 0, True)
        """
        if t_span is None and t_eval is None:
            raise ValueError("t_span or t_eval is required")
        start, end = self._initial.t, None
        if t_span is not None:
            if not isinstance(t_span, (tuple, list)) or len(t_span) != 2:
                raise ValueError("t_span requires exactly two time quantities")
            start = _scalar_quantity(t_span[0], u.s, "t_span[0]")
            end = _scalar_quantity(t_span[1], u.s, "t_span[1]")
            if start != self._initial.t:
                raise ValueError("t_span[0] must equal the initial coordinate time")
            if end <= start:
                raise ValueError("t_span requires a future end time")
        evaluation = None
        if t_eval is not None:
            evaluation = readonly_quantity(t_eval, u.s, "t_eval", ndim=(1,))
            numbers = evaluation.to_value(self._time_unit)
            if evaluation.size == 0 or not np.all(np.isfinite(numbers)) or np.any(np.diff(numbers) <= 0):
                raise ValueError("t_eval must be nonempty, finite and strictly increasing")
            if evaluation[0] < start:
                raise ValueError("t_eval precedes the initial time")
            if end is None:
                end = evaluation[-1]
            elif evaluation[-1] > end:
                raise ValueError("t_eval exceeds t_span")
        if end <= start:
            raise ValueError("solve requires a future end time")
        options = _resolve_options(self._initial_y, method=method, rtol=rtol, atol=atol)
        escape = self._escape_radius(r_escape, self._initial_y)
        states, final, stats, status, reason = self._call_native(
            self._initial_y, end, evaluation, options, escape, store_steps=True)
        message = (_failure_message(stats) if status == -1 else
                   f"Kerr null integration reached {reason}" if status == 1 else
                   "Kerr null integration reached the target time")
        terminal = None
        if status == 1:
            last = self._state(final, scalar=True)
            terminal = NullTermination(reason, last.t, last)
        if not len(states):
            if status == 1:
                raise _NullIntegrationTerminated(reason, terminal.t, terminal.state,
                                                message + " before the first t_eval sample")
            raise IntegrationError(message + "; no t_eval sample was reached")
        if evaluation is not None:
            physical_times = evaluation[:len(states)].to(self._time_unit)
        else:
            numbers = (states[:, 0] * self._metric._time_scale).to_value(self._time_unit)
            numbers[0] = self._initial.t.to_value(self._time_unit)
            if status == 0:
                numbers[-1] = end.to_value(self._time_unit)
            physical_times = numbers * self._time_unit
        state = self._state(states, scalar=False, t_physical=physical_times)
        info = IntegrationInfo(options.method, options.rtol, options.atol, None, None,
                               int(stats["n_steps"]), int(stats["nfev"]))
        return NullSolution(state=state, integration=info, status=status,
                            message=message, termination=terminal)


delegate_null_properties(Null, "_current")
