"""Define immutable integration diagnostic and termination records."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from astropy import units as u

from .._validation import non_negative_integer, readonly_quantity, require_instance
from ._options import integration_options
from .state import State


@dataclass(frozen=True, slots=True, eq=False)
class IntegrationInfo:
    """Immutable effective integration configuration and diagnostics.

    Parameters
    ----------
    method : str
        Non-empty integration-method identifier, stored exactly as supplied.
    rtol : float
        Effective scalar relative tolerance.
    atol : float or array_like
        Effective unitless absolute tolerance: a non-negative scalar or a
        nonempty, one-dimensional, read-only array of non-negative values.
    max_step, first_step : astropy.units.Quantity or None
        Optional finite, positive proper-time step limits.
    n_steps, nfev : int
        Non-negative step and function-evaluation counts.

    Attributes
    ----------
    method, rtol, atol, max_step, first_step, n_steps, nfev
        Immutable normalized values described in ``Parameters``.

    Raises
    ------
    TypeError
        If ``method``, a tolerance, or a counter has the wrong type.
    astropy.units.UnitConversionError
        If a step quantity is not compatible with time.
    ValueError
        If ``method`` is empty, a tolerance is negative, non-finite, or has an
        invalid shape, a step is non-positive or non-finite, or a counter is
        negative.

    Notes
    -----
    Orbit integrates the normalized eight-component Boyer--Lindquist state
    ``(t/T0, R/r_g, Theta, Phi, u^t, u^R, u^Theta, u^Phi)``. A vector
    ``atol`` therefore has eight entries. This container stores the
    effective tolerances used, including the automatic ones chosen by
    :meth:`Orbit.solve <relatipy.geodesic.Orbit.solve>` when ``rtol`` or
    ``atol`` is omitted; it sets no defaults itself. The local error scale
    is ``atol + rtol * abs(y)`` for each component of ``y``.

    Examples
    --------
    >>> info = IntegrationInfo("EXAMPLE_METHOD", 1e-9, 1e-12, None, None, 42, 301)
    >>> info.n_steps
    42
    """

    method: str
    rtol: float
    atol: float | np.ndarray
    max_step: u.Quantity | None
    first_step: u.Quantity | None
    n_steps: int
    nfev: int

    def __post_init__(self) -> None:
        """Validate and normalize integration settings and diagnostics."""
        options = integration_options(
            method=self.method,
            rtol=self.rtol,
            atol=self.atol,
            first_step=self.first_step,
            max_step=self.max_step,
        )
        object.__setattr__(self, "rtol", options.rtol)
        object.__setattr__(self, "atol", options.atol)
        object.__setattr__(self, "max_step", options.max_step)
        object.__setattr__(self, "first_step", options.first_step)
        object.__setattr__(
            self, "n_steps", non_negative_integer(self.n_steps, "n_steps")
        )
        object.__setattr__(self, "nfev", non_negative_integer(self.nfev, "nfev"))


@dataclass(frozen=True, slots=True, eq=False)
class Termination:
    """Store structured information about an internal terminal event.

    Parameters
    ----------
    reason : str
        Non-empty stable identifier for the terminal condition.
    tau : astropy.units.Quantity
        Scalar proper time of the event, compatible with time and equal to
        ``state.tau`` after unit conversion.
    state : State
        Last valid scalar state at the event.

    Attributes
    ----------
    reason : str
        Stored terminal-condition identifier.
    tau : astropy.units.Quantity
        Immutable scalar proper time.
    state : State
        Immutable scalar state supplied at construction.

    Raises
    ------
    TypeError
        If ``reason`` is not a string, ``tau`` is not a quantity, or ``state``
        is not a :class:`State`.
    astropy.units.UnitConversionError
        If ``tau`` is not compatible with time.
    ValueError
        If ``reason`` is empty, ``tau`` is not scalar, ``state`` is
        vectorized, or ``tau`` does not match ``state.tau``.
    """

    reason: str
    tau: u.Quantity
    state: State

    def __post_init__(self) -> None:
        """Validate and freeze the terminal-event record."""
        if not isinstance(self.reason, str):
            raise TypeError("reason must be a string")
        if not self.reason:
            raise ValueError("reason must not be empty")
        tau = readonly_quantity(self.tau, u.s, "tau", ndim=(0,))
        require_instance(self.state, State, "state")
        if self.state.tau.shape != ():
            raise ValueError("termination state must be scalar")
        if not np.array_equal(tau.to_value(self.state.tau.unit), self.state.tau.value):
            raise ValueError("termination tau must match state.tau")
        object.__setattr__(self, "tau", tau)
