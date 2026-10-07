"""Scalar Kerr orbit frontend over the native initial-state and solver kernels."""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import numpy as np
from astropy import units as u
from astropy.constants import G, c

from .._validation import dimensionless_scalar, immutable_array, readonly_quantity
from ..coordinates import KerrOrbitalElements, OrbitalElements
from . import _native as native
from ._deferred import DeferredInitialState, DeferredState
from ._delegation import delegate_state_properties
from ._options import IntegrationOptions, integration_options
from .exceptions import IntegrationError, IntegrationTerminated, IntegrationWarning
from .integration import IntegrationInfo, Termination
from .solution import Solution, constants_of_motion, select_state
from .state import InitialState, State

if TYPE_CHECKING:
    from ..metrics.kerr import Kerr


_FAMILIES = {
    "elements": ("a", "e", "inc", "Omega", "omega", "f"),
    "cartesian": ("x", "y", "z", "vx", "vy", "vz"),
    "spherical": ("r", "theta", "phi", "vr", "vtheta", "vphi"),
    "bl": ("R", "Theta", "Phi", "vR", "vTheta", "vPhi"),
}
_METHODS = ("radau", "dop853", "dp45", "projection_radau")
_BOUND_KEYS = ("p", "e", "x", "q_r0", "q_theta0", "q_phi0")
# Automatic relative tolerance used when ``rtol`` is omitted.
_AUTO_RTOL = 1e-10
# Loosest relative tolerance, and loosest atol per characteristic scale of the
# t, R, Theta, Phi and u^t components, before warning about unfit tolerances.
_MAX_FIT_TOL = 1e-6
# Tightest relative tolerance attainable in double precision.
_MIN_FIT_RTOL = 100 * np.finfo(float).eps
_STATE_NAMES = ("t", "R", "Theta", "Phi", "u^t", "u^R", "u^Theta", "u^Phi")
# Largest tau_eval spacing, as a fraction of the Kepler period, before warning.
_SAMPLES_PER_PERIOD = 10


def _scalar_quantity(value: object, unit: u.UnitBase, name: str) -> u.Quantity:
    """Validate one finite scalar quantity compatible with ``unit``."""
    result = readonly_quantity(value, unit, name, ndim=(0,))
    if not np.isfinite(result.to_value(unit)):
        raise ValueError(f"{name} must be finite")
    return result


def _eccentricity(value: object) -> float:
    """Validate and return one finite dimensionless eccentricity value."""
    number = dimensionless_scalar(value, "e")
    if not np.isfinite(number):
        raise ValueError("e must be finite")
    return number


def _classify_inputs(values: dict[str, object]) -> str:
    """Identify the sole supported initial-condition family in ``values``.

    Raises ``TypeError`` for a non-``KerrOrbitalElements`` bound container and
    ``ValueError`` for missing or incompatible families and required fields.
    """
    # Discriminate bound input before looking at the shared e/x names.
    # Without p/phases, x retains its original Cartesian meaning.
    bound = values.get("elements")
    raw_bound = any(
        values.get(key) is not None
        for key in ("p", "q_r0", "q_theta0", "q_phi0")
    )
    if bound is not None:
        if not isinstance(bound, KerrOrbitalElements):
            raise TypeError("elements must be a KerrOrbitalElements instance")
        if any(
            value is not None for key, value in values.items()
            if key not in ("elements", "tau", "t")
        ):
            raise ValueError("incompatible initial-condition families with elements")
        return "bound"
    if raw_bound:
        incompatible = [
            key for key, value in values.items()
            if value is not None and key not in (*_BOUND_KEYS, "tau", "t")
        ]
        if incompatible:
            raise ValueError(
                "incompatible initial-condition families: bound and "
                + ", ".join(incompatible)
            )
        missing = [key for key in ("p", "e", "x") if values.get(key) is None]
        if missing:
            raise ValueError("bound input requires " + ", ".join(missing))
        return "bound"
    supplied = {
        family: {key for key in keys if values[key] is not None}
        for family, keys in _FAMILIES.items()
    }
    active = [family for family, names in supplied.items() if names]
    if len(active) != 1:
        if not active:
            raise ValueError("one initial-condition family is required")
        raise ValueError(
            "incompatible initial-condition families: " + ", ".join(active)
        )
    family = active[0]
    required = {
        "elements": ("a",),
        "cartesian": ("x", "y", "z"),
        "spherical": ("r", "theta", "phi"),
        "bl": ("R", "Theta", "Phi"),
    }[family]
    if family == "cartesian":
        if not any(values[key] is not None for key in required):
            raise ValueError("Cartesian input requires x, y, or z")
    else:
        missing = [key for key in required if values[key] is None]
        if missing:
            raise ValueError(
                f"{family} input requires " + ", ".join(missing)
            )
    return family


def _initial_components(
    metric: Kerr, family: str, values: dict[str, object], t: u.Quantity
) -> tuple[np.ndarray, u.UnitBase]:
    """Normalize one initial-condition family using a metric's scales."""
    return _initial_components_from_scales(
        metric._length_scale, metric._time_scale, family, values, t
    )


def _initial_components_from_scales(
    length_scale: u.Quantity,
    time_scale: u.Quantity,
    family: str,
    values: dict[str, object],
    t: u.Quantity,
) -> tuple[np.ndarray, u.UnitBase]:
    """Convert one orbit input family without constructing a Kerr object."""
    speed_scale = c.to(u.m / u.s)
    if family == "bound":
        elements = values.get("elements")
        if elements is None:
            elements = KerrOrbitalElements(**{
                key: values[key] for key in _BOUND_KEYS
                if values.get(key) is not None
            })
        components = [
            (t / time_scale).to_value(u.one),
            (elements.p / length_scale).to_value(u.one),
            elements.e, elements.x,
            elements.q_r0.to_value(u.rad),
            elements.q_theta0.to_value(u.rad),
            elements.q_phi0.to_value(u.rad),
        ]
        return np.asarray(components, dtype=np.float64), elements.p.unit
    if family == "elements":
        length = _scalar_quantity(values["a"], u.m, "a")
        components = [
            (t / time_scale).to_value(u.one),
            (length / length_scale).to_value(u.one),
            _eccentricity(0 if values["e"] is None else values["e"]),
        ]
        for key in ("inc", "Omega", "omega", "f"):
            angle = _scalar_quantity(
                0 * u.rad if values[key] is None else values[key], u.rad, key
            )
            components.append(angle.to_value(u.rad))
        return np.asarray(components, dtype=np.float64), length.unit

    first_name = next(
        key
        for key in ("x", "y", "z", "r", "R")
        if values.get(key) is not None
    )
    length_unit = _scalar_quantity(values[first_name], u.m, first_name).unit
    components = [(t / time_scale).to_value(u.one)]
    keys = _FAMILIES[family]
    for key in keys:
        raw = values[key]
        if key in ("x", "y", "z"):
            val = _scalar_quantity(
                0 * length_unit if raw is None else raw, u.m, key
            )
            components.append((val / length_scale).to_value(u.one))
        elif key in ("r", "R"):
            val = _scalar_quantity(raw, u.m, key)
            components.append((val / length_scale).to_value(u.one))
        elif key in ("theta", "phi", "Theta", "Phi"):
            val = _scalar_quantity(raw, u.rad, key)
            components.append(val.to_value(u.rad))
        elif key in ("vx", "vy", "vz", "vr", "vR"):
            val = _scalar_quantity(
                0 * u.m / u.s if raw is None else raw, u.m / u.s, key
            )
            components.append((val / speed_scale).to_value(u.one))
        else:
            val = _scalar_quantity(
                0 * u.rad / u.s if raw is None else raw,
                u.rad / u.s,
                key,
            )
            components.append((val * time_scale).to_value(u.rad))
    return np.asarray(components, dtype=np.float64), length_unit


def _state_from_native(
    metric: Kerr,
    canonical: np.ndarray,
    tau: np.ndarray,
    *,
    length_unit: u.UnitBase,
    time_unit: u.UnitBase,
    initial: bool = False,
    scalar: bool = False,
    tau_physical: u.Quantity | None = None,
) -> State:
    """Wrap one native batch in read-only physical state containers."""
    from .. import _core

    canonical = np.asarray(canonical, dtype=np.float64)
    tau = np.asarray(tau, dtype=np.float64)
    if canonical.ndim != 2 or canonical.shape[1] != 8:
        raise RuntimeError("native canonical states must have shape (n, 8)")
    if tau.shape != (canonical.shape[0],):
        raise RuntimeError("native proper times do not match canonical states")
    rows, cartesian, row_status = _core.reconstruct_canonical_batch(
        metric.spin, np.ascontiguousarray(canonical)
    )
    rows = np.asarray(rows, dtype=np.float64)
    cartesian = np.asarray(cartesian, dtype=np.float64)
    row_status = np.asarray(row_status)
    if (
        rows.shape != (len(tau), native.ROW_WIDTH)
        or cartesian.shape != (len(tau), 7)
        or row_status.shape != (len(tau),)
    ):
        raise RuntimeError("native state reconstruction returned invalid shape")
    failed = np.flatnonzero(row_status != 0)
    if failed.size:
        row = int(failed[0])
        raise ValueError(
            f"state row {row} cannot be reconstructed "
            f"({native.describe_reconstruction_status(row_status[row])})"
        )
    native.validate_reconstructed_rows(rows)
    if not np.all(np.isfinite(cartesian)):
        raise ValueError("native state reconstruction returned non-finite values")

    length_scale = metric._length_scale
    time_scale = metric._time_scale
    speed_scale = c.to(u.m / u.s)
    velocity_unit = length_unit / time_unit
    angular_rate_scale = u.rad / time_scale
    angular_velocity_unit = u.rad / time_unit

    def length(values: np.ndarray) -> u.Quantity:
        """Convert native normalized lengths to the requested length unit."""
        return (values * length_scale).to(length_unit)

    def time(values: np.ndarray) -> u.Quantity:
        """Convert native normalized times to the requested time unit."""
        return (values * time_scale).to(time_unit)

    def speed(values: np.ndarray) -> u.Quantity:
        """Convert native normalized speeds to the requested velocity unit."""
        return (values * speed_scale).to(velocity_unit)

    converters = {
        "length": length,
        "angle": lambda values: values * u.rad,
        "speed": speed,
        "rate": lambda values: (values * angular_rate_scale).to(
            angular_velocity_unit
        ),
        "one": lambda values: values * u.one,
    }

    if initial or scalar:
        if len(tau) != 1:
            raise RuntimeError("scalar state reconstruction requires one row")
        rows = rows[0]
        cartesian = cartesian[0]
        tau = tau[0]
    proper_time = time(tau) if tau_physical is None else tau_physical.to(time_unit)
    if proper_time.shape != np.shape(tau):
        raise RuntimeError("physical proper time does not match native state shape")
    return native.assemble_state(
        lambda column: converters[native.COLUMN_KIND[column]](rows[..., column]),
        tau=proper_time,
        t=time(cartesian[..., 0]),
        x=length(cartesian[..., 1]),
        y=length(cartesian[..., 2]),
        z=length(cartesian[..., 3]),
        vxyz=speed(cartesian[..., 4:7]),
        state_type=InitialState if initial else State,
    )


def _integration_options(
    *,
    method: object,
    rtol: object,
    atol: object,
    first_step: object,
    max_step: object,
) -> IntegrationOptions:
    """Validate settings and apply the Kerr solver's method and shape limits."""
    options = integration_options(
        method=method, rtol=rtol, atol=atol, first_step=first_step, max_step=max_step
    )
    if options.method not in _METHODS:
        raise ValueError("method must be one of radau, dop853, dp45, projection_radau")
    if isinstance(options.atol, np.ndarray) and options.atol.shape != (8,):
        raise ValueError("atol array must have shape (8,)")
    return options


def _failure_message(stats: dict) -> str:
    """Build the public numerical-failure message from native diagnostics."""
    reason = native.describe_integrator_status(stats.get("native_status"))
    return f"Kerr integration failed before the target time ({reason})"


def _state_scales(y0: np.ndarray) -> np.ndarray:
    """Characteristic magnitude of each normalized native state component.

    Uses the initial radius ``R0`` (in ``r_g``) and Newtonian orders of
    magnitude: dynamical time ``R0**1.5`` for ``t``, ``R0`` for ``R``, one
    radian for the angles, one for ``u^t``, circular speed ``R0**-0.5`` for
    ``u^R`` and angular rate ``R0**-1.5`` for ``u^Theta`` and ``u^Phi``. Each
    scale is at least the magnitude of the initial component.
    """
    radius = float(y0[1])
    floor = np.array([
        radius**1.5, radius, 1.0, 1.0, 1.0,
        radius**-0.5, radius**-1.5, radius**-1.5,
    ])
    return np.maximum(np.abs(np.asarray(y0, dtype=float)), floor)


def _resolve_options(
    y0: np.ndarray,
    *,
    method: object,
    rtol: object,
    atol: object,
    first_step: object,
    max_step: object,
) -> IntegrationOptions:
    """Validate settings and fill omitted tolerances automatically.

    An omitted ``rtol`` becomes ``1e-10``. An omitted ``atol`` becomes the
    shape ``(8,)`` array ``rtol * _state_scales(y0)``, so each component is
    controlled relative to its own characteristic magnitude. Explicit values
    that are unfit for the state emit :class:`IntegrationWarning`.
    """
    options = _integration_options(
        method=method,
        rtol=_AUTO_RTOL if rtol is None else rtol,
        atol=0.0 if atol is None else atol,
        first_step=first_step,
        max_step=max_step,
    )
    scales = _state_scales(y0)
    if atol is None:
        options = options._replace(atol=immutable_array(options.rtol * scales))
    problems = []
    if rtol is not None and options.rtol > _MAX_FIT_TOL:
        problems.append(f"rtol={options.rtol:g} is looser than {_MAX_FIT_TOL:g}")
    if rtol is not None and 0.0 < options.rtol < _MIN_FIT_RTOL:
        problems.append(
            f"rtol={options.rtol:g} is below {_MIN_FIT_RTOL:.1g}, "
            "unattainable in double precision"
        )
    if atol is not None:
        ratio = np.broadcast_to(options.atol, (8,))[:5] / scales[:5]
        loose = [_STATE_NAMES[i] for i in np.flatnonzero(ratio > _MAX_FIT_TOL)]
        if loose:
            problems.append(
                f"atol exceeds {_MAX_FIT_TOL:g} times the characteristic "
                f"scale of {', '.join(loose)}"
            )
    if problems:
        warnings.warn(
            "unfit integration tolerances for this orbit: "
            + "; ".join(problems)
            + ". The trajectory may be inaccurate or the integration may "
            "fail; omit rtol and atol to choose them automatically",
            IntegrationWarning,
            stacklevel=3,
        )
    return options


def _keplerian_period(elements: OrbitalElements, mass: u.Quantity) -> u.Quantity | None:
    """Newtonian period ``2*pi*sqrt(a**3 / (G*M))`` of an elliptic osculating conic.

    Returns ``None`` unless ``0 <= e < 1`` and ``a`` is finite and positive.
    """
    eccentricity = float(elements.e)
    semi_major = elements.a
    if not (0.0 <= eccentricity < 1.0) or not np.isfinite(semi_major.value):
        return None
    if not semi_major.value > 0.0:
        return None
    return 2 * np.pi * np.sqrt(semi_major**3 / (G * mass))


def _warn_sparse_sampling(
    initial: InitialState, mass: u.Quantity, evaluation: u.Quantity
) -> None:
    """Warn when tau_eval is too sparse for the initial Kepler period."""
    if evaluation.size < 2:
        return
    period = _keplerian_period(initial.orbital_elements(), mass)
    if period is None:
        return
    period = period.to(evaluation.unit)
    spacing = np.max(np.diff(evaluation))
    if spacing <= period / _SAMPLES_PER_PERIOD:
        return
    warnings.warn(
        f"tau_eval spacing up to {spacing:.4g} exceeds 1/{_SAMPLES_PER_PERIOD} "
        f"of the estimated Kepler period {period:.4g} of the initial state; "
        "the sampled trajectory may alias the orbit. Use denser tau_eval "
        "samples",
        IntegrationWarning,
        stacklevel=3,
    )


class Orbit:
    """One mutable point state with immutable saved Kerr initial conditions.

    Create orbits with :meth:`relatipy.metrics.Kerr.orbit`.  The constructor
    receives already validated native state and is not a public entry point.
    """

    __slots__ = (
        "_metric",
        "_initial",
        "_current",
        "_initial_y",
        "_current_y",
        "_length_unit",
        "_time_unit",
        "_input_family",
    )

    def __init__(
        self,
        metric: Kerr,
        initial_state: InitialState,
        initial_y: np.ndarray,
        length_unit: u.UnitBase,
        time_unit: u.UnitBase,
        input_family: str,
    ) -> None:
        """Store already validated metric, state, and native orbit data.

        This internal constructor receives a scalar initial state and its
        canonical eight-component representation. Public callers construct an
        orbit through :meth:`relatipy.metrics.Kerr.orbit`.
        """
        self._metric = metric
        self._initial = initial_state
        self._current = initial_state
        self._initial_y = np.frombuffer(
            np.asarray(initial_y, dtype=np.float64).tobytes(), dtype=np.float64
        )
        self._current_y = self._initial_y
        self._length_unit = length_unit
        self._time_unit = time_unit
        self._input_family = input_family

    @classmethod
    def _from_kwargs(cls, metric: Kerr, **values: object) -> Orbit:
        """Build an orbit from the keyword interface of ``Kerr.orbit``.

        The method classifies the requested input family, converts it to the
        native normalized convention, and rejects a non-exterior initial
        Boyer--Lindquist radius.
        """
        from .. import _core

        family = _classify_inputs(values)
        tau = _scalar_quantity(values["tau"], u.s, "tau")
        t = _scalar_quantity(values["t"], u.s, "t")
        components, length_unit = _initial_components(metric, family, values, t)
        canonical = np.asarray(
            _core.initial_kerr(metric.spin, family, np.ascontiguousarray(components)),
            dtype=np.float64,
        )
        if canonical.shape != (8,) or not np.all(np.isfinite(canonical)):
            raise ValueError("native Kerr initial state is invalid")
        if canonical[1] <= metric._radii[0]:
            raise ValueError("initial state must have R > the outer Kerr horizon")
        normalized_tau = float((tau / metric._time_scale).to_value(u.one))
        initial = DeferredInitialState(
            canonical=canonical.reshape(1, 8),
            tau=np.array([normalized_tau], dtype=np.float64),
            spin=metric.spin,
            length_scale=metric._length_scale,
            time_scale=metric._time_scale,
            length_unit=length_unit,
            time_unit=tau.unit,
            input_family=family,
            tau_physical=tau,
            scalar=True,
        )
        return cls(metric, initial, canonical, length_unit, tau.unit, family)

    @property
    def initial(self) -> InitialState:
        """Return the immutable state saved when the orbit was constructed.

        Returns
        -------
        InitialState
            Scalar read-only initial state. It is unchanged by
            :meth:`integrate` and is the start point for :meth:`solve`.
        """
        return self._initial

    def reset(self) -> None:
        """Restore the exact saved initial state without integrating.

        Returns
        -------
        None
        """
        self._current = self._initial
        self._current_y = self._initial_y

    def copy(self) -> Orbit:
        """Return an independently evolving copy of this orbit.

        Returns
        -------
        Orbit
            Copy retaining the same saved initial state and current state.
        """
        clone = Orbit(
            self._metric,
            self._initial,
            self._initial_y,
            self._length_unit,
            self._time_unit,
            self._input_family,
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
        """Integrate one normalized span and validate native result shapes.

        The returned proper-time vector has shape ``(n,)`` and the canonical
        state array has shape ``(n, 8)``. Raises ``IntegrationError`` for
        invalid native output or an unsupported status code.
        """
        from .. import _core

        scale = self._metric._time_scale
        result = _core.integrate_kerr(
            self._metric.spin,
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

    def integrate(
        self,
        tau: u.Quantity,
        *,
        method: str = "radau",
        rtol: float | None = None,
        atol: float | np.ndarray | None = None,
        first_step: u.Quantity | None = None,
        max_step: u.Quantity | None = None,
    ) -> None:
        """Advance to an absolute proper time, mutating this orbit.

        Parameters
        ----------
        tau : astropy.units.Quantity
            Scalar absolute proper time at or after the saved initial time.
        method : {"radau", "dop853", "dp45", "projection_radau"}, optional
            Native integration method, default ``"radau"``.
        rtol : float or None, optional
            Relative local-error tolerance for the normalized native state.
            The default ``None`` selects ``1e-10``.
        atol : float, array_like or None, optional
            Non-negative absolute tolerance: a scalar or a shape ``(8,)``
            array in native state order
            ``(t/T0, R/r_g, Theta, Phi, u^t, u^R, u^Theta, u^Phi)``. The
            default ``None`` selects ``rtol`` times the characteristic scale
            of each component at the starting state.
        first_step, max_step : astropy.units.Quantity or None, optional
            Positive proper-time initial and maximum step sizes. The default
            ``None`` lets the native integrator choose.

        Returns
        -------
        None

        Raises
        ------
        IntegrationTerminated
            An internal terminal event occurred. The last valid state remains.
        IntegrationError
            Numerical failure occurred. The last valid state remains.
        TypeError, astropy.units.UnitConversionError
            The target or an option has the wrong type or unit.
        ValueError
            The target or numerical options are invalid.

        Warns
        -----
        IntegrationWarning
            An explicit ``rtol`` is looser than ``1e-6`` or below
            ``100 * eps``, or an explicit ``atol`` exceeds ``1e-6`` times the
            characteristic scale of ``t``, ``R``, ``Theta``, ``Phi`` or
            ``u^t``.

        Notes
        -----
        A target before the current state restarts from saved initial
        conditions and integrates forward; it does not reverse the solver.
        """
        target = _scalar_quantity(tau, u.s, "tau")
        if target < self._initial.tau:
            raise ValueError("tau precedes the saved initial proper time")
        restart = target < self._current.tau
        source = self._initial if restart else self._current
        y0 = self._initial_y if restart else self._current_y
        options = _resolve_options(
            y0,
            method=method,
            rtol=rtol,
            atol=atol,
            first_step=first_step,
            max_step=max_step,
        )
        if target == self._current.tau:
            return None
        scale = self._metric._time_scale
        tau0 = float((source.tau / scale).to_value(u.one))
        tau1 = float((target / scale).to_value(u.one))
        taus, states, stats, status = self._call_native(
            y0, tau0, tau1, options, store_steps=False
        )
        last = DeferredState(
            canonical=states[-1:],
            tau=taus[-1:],
            spin=self._metric.spin,
            length_scale=self._metric._length_scale,
            time_scale=self._metric._time_scale,
            length_unit=self._length_unit,
            time_unit=self._time_unit,
            input_family=self._input_family,
            scalar=True,
            tau_physical=target if status == 0 else None,
        )
        self._current = last
        self._current_y = np.frombuffer(states[-1].tobytes(), dtype=np.float64)
        if status == 1:
            horizon = bool(stats.get("crossed_outer_horizon", False))
            raise IntegrationTerminated(
                "outer_horizon" if horizon else "internal_terminal_event",
                last.tau,
                last,
                "integration reached the outer Kerr horizon"
                if horizon else "integration reached an internal terminal event",
            )
        if status == -1:
            raise IntegrationError(_failure_message(stats))
        return None

    def solve(
        self,
        *,
        tau_span: tuple[u.Quantity, u.Quantity] | None = None,
        tau_eval: u.Quantity | None = None,
        method: str = "radau",
        rtol: float | None = None,
        atol: float | np.ndarray | None = None,
        first_step: u.Quantity | None = None,
        max_step: u.Quantity | None = None,
    ) -> Solution:
        """Integrate from saved initial conditions without mutating this orbit.

        Parameters
        ----------
        tau_span : tuple of astropy.units.Quantity, optional
            Initial and final absolute proper time. The first equals the
            saved initial proper time.
        tau_eval : astropy.units.Quantity, optional
            Nonempty, finite, strictly increasing one-dimensional sample
            times within the span. Either this or ``tau_span`` is required.
        method : {"radau", "dop853", "dp45", "projection_radau"}, optional
            Native integration method, default ``"radau"``.
        rtol : float or None, optional
            Relative local-error tolerance for the normalized native state.
            The default ``None`` selects ``1e-10``.
        atol : float, array_like or None, optional
            Non-negative absolute tolerance: a scalar or a shape ``(8,)``
            array in native state order
            ``(t/T0, R/r_g, Theta, Phi, u^t, u^R, u^Theta, u^Phi)``. The
            default ``None`` selects ``rtol`` times the characteristic scale
            of each component at the starting state.
        first_step, max_step : astropy.units.Quantity or None, optional
            Positive proper-time initial and maximum step sizes. The default
            ``None`` lets the native integrator choose.

        Returns
        -------
        Solution
            Immutable sampled trajectory, including any partial result.  With
            ``tau_eval``, it holds the requested samples reached before a
            numerical failure or terminal event.

        Raises
        ------
        IntegrationError
            Numerical failure occurred before the first ``tau_eval`` sample.
        IntegrationTerminated
            A terminal event occurred before the first ``tau_eval`` sample.
        TypeError, astropy.units.UnitConversionError
            A time or option has the wrong type or unit.
        ValueError
            The time span, samples, or numerical options are invalid.

        Warns
        -----
        IntegrationWarning
            The largest ``tau_eval`` spacing exceeds one tenth of the Kepler
            period estimated from the initial osculating elements (bound
            orbits only), or an explicit ``rtol`` or ``atol`` is unfit for
            the initial state (see :meth:`integrate`).

        Notes
        -----
        The native integrator works in the mass-normalized state; proper
        times are divided by ``T0 = G M / c^3`` before integration and
        restored afterwards. ``Solution.status`` is ``0`` when the final
        time was reached (its last stored proper time is then exactly the
        requested end), ``1`` for an internal terminal event such as a
        confirmed outer-horizon crossing, and ``-1`` for a numerical failure.
        A failure or terminal event does not raise when at least one sample
        can be returned: the partial :class:`~relatipy.geodesic.Solution`
        has ``success == False`` and keeps the last valid state.

        Without ``tau_eval`` the solution stores the initial state and every
        accepted native step. With ``tau_eval`` the requested samples reached
        before the end of integration are obtained from those steps with
        :meth:`Solution.at`, so they are interpolated rather than integrated.
        Tolerances control local error only and do not bound the global
        trajectory error.

        Examples
        --------
        >>> from astropy import units as u
        >>> from astropy.constants import c
        >>> from relatipy import Kerr
        >>> bh = Kerr(mass=1 * u.Msun, spin=0.5)
        >>> orbit = bh.orbit(x=12 * bh.r_g, vy=0.1 * c)
        >>> solution = orbit.solve(tau_span=(0 * u.s, 2e-8 * u.s), method="dp45")
        >>> solution.status, solution.success
        (0, True)
        >>> sampled = orbit.solve(tau_eval=[0, 1e-8, 2e-8] * u.s, method="dp45")
        >>> len(sampled)
        3
        """
        options = _resolve_options(
            self._initial_y,
            method=method,
            rtol=rtol,
            atol=atol,
            first_step=first_step,
            max_step=max_step,
        )
        if tau_span is None and tau_eval is None:
            raise ValueError("tau_span or tau_eval is required")
        if tau_span is not None:
            if not isinstance(tau_span, (tuple, list)) or len(tau_span) != 2:
                raise ValueError("tau_span must contain exactly two time quantities")
            start = _scalar_quantity(tau_span[0], u.s, "tau_span[0]")
            end = _scalar_quantity(tau_span[1], u.s, "tau_span[1]")
            if start != self._initial.tau:
                raise ValueError("tau_span[0] must equal the initial proper time")
            if end < start:
                raise ValueError("tau_span end must not precede its start")
        else:
            start = self._initial.tau
            end = None
        if tau_eval is not None:
            evaluation = readonly_quantity(tau_eval, u.s, "tau_eval", ndim=(1,))
            if evaluation.size == 0:
                raise ValueError("tau_eval must be a nonempty one-dimensional array")
            eval_values = np.asarray(evaluation.to_value(self._time_unit), dtype=float)
            if (
                not np.all(np.isfinite(eval_values))
                or np.any(np.diff(eval_values) <= 0)
            ):
                raise ValueError("tau_eval must be finite and strictly increasing")
            if evaluation[0] < start:
                raise ValueError("tau_eval precedes tau_span")
            if end is None:
                end = evaluation[-1]
            elif evaluation[-1] > end:
                raise ValueError("tau_eval exceeds tau_span")
            _warn_sparse_sampling(self._initial, self._metric.mass, evaluation)
        else:
            evaluation = None
        assert end is not None
        scale = self._metric._time_scale
        tau0 = float((start / scale).to_value(u.one))
        tau1 = float((end / scale).to_value(u.one))
        taus, states, stats, status = self._call_native(
            self._initial_y, tau0, tau1, options, store_steps=True
        )
        physical_taus = (taus * scale).to(self._time_unit)
        if status == 0:
            exact_taus = np.array(physical_taus.value, copy=True)
            exact_taus[-1] = end.to_value(self._time_unit)
            physical_taus = exact_taus * self._time_unit
        full_state = DeferredState(
            canonical=states,
            tau=taus,
            spin=self._metric.spin,
            length_scale=self._metric._length_scale,
            time_scale=self._metric._time_scale,
            length_unit=self._length_unit,
            time_unit=self._time_unit,
            input_family=self._input_family,
            tau_physical=physical_taus,
        )
        info = IntegrationInfo(
            method=options.method,
            rtol=options.rtol,
            atol=options.atol,
            max_step=options.max_step,
            first_step=options.first_step,
            n_steps=int(stats["n_steps"]),
            nfev=int(stats["nfev"]),
        )
        horizon = bool(stats.get("crossed_outer_horizon", False))
        message = {
            -1: _failure_message(stats),
            0: "Kerr integration reached the target time",
            1: ("Kerr integration reached the outer horizon" if horizon
                else "Kerr integration reached an internal terminal event"),
        }[status]
        terminal = None
        reason = "outer_horizon" if horizon else "internal_terminal_event"
        if status == 1:
            final_state = select_state(full_state, -1)
            terminal = Termination(reason, final_state.tau, final_state)
        full = Solution(
            state=full_state,
            integration=info,
            status=status,
            message=message,
            termination=terminal,
            _metric=self._metric,
        )
        if evaluation is None:
            return full
        available = evaluation[evaluation <= full.tau[-1]]
        if not available.size:
            if status == 1:
                raise IntegrationTerminated(
                    reason, terminal.tau, terminal.state,
                    f"{message} before the first tau_eval sample",
                )
            raise IntegrationError(f"{message}; no tau_eval sample was reached")
        return Solution(
            state=full.at(tau=available),
            integration=info,
            status=status,
            message=message,
            termination=terminal,
            _metric=self._metric,
        )

    def preview(
        self,
        *,
        projection: str = "views",
        interactive: bool | None = None,
        show_horizon: bool = True,
        show_isco: bool = True,
        style=None,
        fig=None,
        ax=None,
    ):
        """Plot the current osculating Kepler path and equatorial references.

        Parameters
        ----------
        projection : {"3d", "xy", "xz", "yz", "views"}, optional
            Spatial view, default ``"views"`` (three orthogonal planes).
        interactive : bool or None, optional
            ``True``: Plotly figure; ``False``: Matplotlib ``(fig, ax)``.
            ``None`` (default) is ``True`` for ``"3d"`` and ``False`` otherwise.
        show_horizon, show_isco : bool, optional
            Include equatorial reference circles, default ``True``.
        style : relatipy.plotting.Style or None, optional
            Appearance; ``None`` uses the default publication style and
            displays lengths in the orbit's stored position unit.
        fig : matplotlib.figure.Figure or plotly.graph_objects.Figure, optional
            Existing figure to draw into. Leave unset for ``"views"``.
        ax : matplotlib.axes.Axes, optional
            Existing static axes belonging to ``fig``. Leave unset for
            interactive figures and for ``"views"``.

        Returns
        -------
        plotly.graph_objects.Figure or tuple
            Interactive Plotly figure, static ``(fig, ax)``, or
            ``(fig, (main, top, right))`` for ``"views"``. See
            :func:`relatipy.plotting.preview_orbit`.

        Raises
        ------
        ImportError
            If Plotly is required for an interactive figure or Matplotlib is
            required for a static figure but is unavailable.
        TypeError
            If a plotting argument has an invalid type.
        ValueError
            If the current osculating elements have undefined angles or a
            non-finite semi-major axis, or plotting arguments are invalid.
        RuntimeError
            If the native preview returns invalid positions or reference radii.

        Notes
        -----
        The instantaneous Kepler conic is not an integrated Kerr trajectory.
        The horizon and ISCO circles are equatorial references; their plotted
        Cartesian radii are supplied by the native Kerr kernel.
        """
        from ..plotting import preview_orbit

        return preview_orbit(
            self, projection=projection, interactive=interactive, show_horizon=show_horizon,
            show_isco=show_isco, style=style, fig=fig, ax=ax,
        )

    def _osculating_preview(
        self, samples: int = 360
    ) -> tuple[u.Quantity, u.Quantity]:
        """Return the native osculating conic and equatorial reference radii.

        Returns
        -------
        path : astropy.units.Quantity
            Cartesian conic points of shape ``(samples, 3)``.
        radii : astropy.units.Quantity
            Cartesian equatorial radii of the outer horizon and the prograde
            and retrograde ISCO.  Both results use the orbit length unit.

        Raises
        ------
        ValueError
            If elements are undefined or the semi-major axis is infinite.
        RuntimeError
            If the native preview returns invalid positions or radii.
        """
        from .. import _core

        elements = self.orbital_elements()
        if not np.all(elements.defined):
            raise ValueError(
                "osculating preview requires defined orbital elements; "
                "near-radial states have undefined angles"
            )
        if not np.isfinite(elements.a.value):
            raise ValueError(
                "osculating preview requires a finite semi-major axis; "
                "exactly parabolic states have a = +inf"
            )
        path, radii = _core.osculating_preview(
            self._metric.spin, np.ascontiguousarray(self._current_y), samples,
        )
        path = np.asarray(path, dtype=np.float64)
        radii = np.asarray(radii, dtype=np.float64)
        if path.ndim != 2 or path.shape[1] != 3 or path.shape[0] < 2:
            raise RuntimeError("native osculating preview returned invalid path")
        if radii.shape != (3,) or not np.all(np.isfinite(radii)):
            raise RuntimeError("native osculating preview returned invalid radii")
        if not np.all(np.isfinite(path)):
            raise RuntimeError("native osculating preview returned non-finite points")
        scale = self._metric._length_scale
        return (path * scale).to(self._length_unit), (radii * scale).to(
            self._length_unit
        )

    def orbital_elements(self) -> OrbitalElements:
        """Return classical osculating elements of the current point state.

        Returns
        -------
        OrbitalElements
            Immutable scalar osculating elements. The ``defined`` mask
            identifies available elements; positive infinite ``a`` represents
            an exactly parabolic conic.
        """
        return self._current.orbital_elements()

    def get_keplerian_period(self) -> u.Quantity:
        """Return the Keplerian period of the current osculating conic.

        The period is ``2*pi*sqrt(a**3 / (G*M))``, with ``a`` the osculating
        semi-major axis of the current state and ``M`` the black-hole mass.

        Returns
        -------
        astropy.units.Quantity
            Scalar period in the orbit's time unit.

        Raises
        ------
        ValueError
            If the osculating conic is not elliptic (``e`` outside
            ``[0, 1)`` or ``a`` not finite and positive).

        Notes
        -----
        This is a Newtonian estimate from instantaneous elements, not the
        radial, azimuthal or polar period of the Kerr geodesic. It changes as
        the orbit advances.

        Examples
        --------
        >>> from astropy import units as u
        >>> from relatipy import Kerr
        >>> bh = Kerr(mass=1 * u.Msun, spin=0.0)
        >>> orbit = bh.orbit(a=100 * bh.r_g, e=0.0, inc=0 * u.rad,
        ...                  Omega=0 * u.rad, omega=0 * u.rad, f=0 * u.rad)
        >>> orbit.get_keplerian_period().unit == orbit.tau.unit
        True
        """
        period = _keplerian_period(self.orbital_elements(), self._metric.mass)
        if period is None:
            raise ValueError(
                "the Keplerian period needs an elliptic osculating conic "
                "(0 <= e < 1 and finite positive a)"
            )
        return period.to(self._time_unit)


    def get_E(self) -> float:
        """Return the specific energy of the current point state.

        Returns
        -------
        float
            ``E = -u_t`` normalized with ``G = c = M = mu = 1``, i.e.
            ``E / (mu c^2)``. Bound orbits have ``E < 1``.

        Raises
        ------
        ValueError
            If the current state lies on a Boyer--Lindquist coordinate
            singularity.

        Notes
        -----
        The value follows the current state, so after :meth:`integrate` it
        includes the numerical drift of the integration.
        """
        return float(constants_of_motion(self._current)[0, 0])

    def get_Lz(self) -> float:
        """Return the specific axial angular momentum of the current state.

        Returns
        -------
        float
            ``Lz = u_phi`` in units of ``mu G M / c`` (normalized with
            ``G = c = M = mu = 1``).

        Raises
        ------
        ValueError
            If the current state lies on a Boyer--Lindquist coordinate
            singularity.
        """
        return float(constants_of_motion(self._current)[0, 1])

    def get_Q(self) -> float:
        """Return the Carter constant of the current point state.

        Returns
        -------
        float
            ``Q = u_theta**2 + cos(theta)**2 * (a**2 (1 - E**2) + Lz**2 / sin(theta)**2)``
            in units of ``(mu G M / c)**2``, with fixed rest mass ``mu = 1``.
            It is not ``K = Q + (Lz - a E)**2``; equatorial orbits have
            ``Q = 0``.

        Raises
        ------
        ValueError
            If the current state lies on a Boyer--Lindquist coordinate
            singularity.
        """
        return float(constants_of_motion(self._current)[0, 2])

delegate_state_properties(Orbit, "_current")
