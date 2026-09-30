"""Internal validation of integration options.

:class:`~relatipy.geodesic.IntegrationInfo` and the orbit integration methods
share these rules, so the same input is accepted or rejected with the same
message in both places.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np
from astropy import units as u

from .._validation import (
    immutable_array,
    non_negative_real,
    readonly_quantity,
)

class IntegrationOptions(NamedTuple):
    """Validated integration settings in their normalized form."""

    method: str
    rtol: float
    atol: float | np.ndarray
    first_step: u.Quantity | None
    max_step: u.Quantity | None


def absolute_tolerance(value: object) -> float | np.ndarray:
    """Validate a scalar or one-dimensional absolute tolerance.

    Parameters
    ----------
    value : object
        Unitless real scalar or nonempty one-dimensional sequence of real
        numbers. A zero-dimensional NumPy array is a scalar.

    Returns
    -------
    float or numpy.ndarray
        A scalar float or an independent, immutable float array.

    Raises
    ------
    TypeError
        If values have units or are not real numbers.
    ValueError
        If the array has an invalid shape or any value is negative or
        non-finite.
    """
    if isinstance(value, u.Quantity):
        raise TypeError("atol must contain unitless real numbers")
    if isinstance(value, (list, tuple)):
        if any(isinstance(item, u.Quantity) for item in value):
            raise TypeError("atol must contain unitless real numbers")
        if any(isinstance(item, (bool, np.bool_)) for item in value):
            raise TypeError("atol must contain real numbers, excluding booleans")
    try:
        array = np.asarray(value)
    except (TypeError, ValueError) as exc:
        raise TypeError("atol must contain real numbers") from exc
    if array.dtype.kind not in "iuf":
        raise TypeError("atol must contain real numbers")
    if array.ndim == 0:
        return non_negative_real(array.item(), "atol")
    if array.ndim != 1 or array.size == 0:
        raise ValueError("atol must be scalar or a nonempty one-dimensional array")
    with np.errstate(over="ignore", invalid="ignore"):
        result = np.asarray(array, dtype=float)
    if not np.all(np.isfinite(result)) or np.any(result < 0):
        raise ValueError("atol must be finite and non-negative")
    return immutable_array(result)


def positive_step(value: object, name: str) -> u.Quantity:
    """Validate a finite, positive, scalar proper-time step.

    Parameters
    ----------
    value : astropy.units.Quantity
        Scalar step compatible with time.
    name : str
        Field name used in validation errors.

    Returns
    -------
    astropy.units.Quantity
        Read-only copy in the input unit.

    Raises
    ------
    TypeError
        If ``value`` is not an Astropy quantity.
    astropy.units.UnitConversionError
        If ``value`` is not compatible with time.
    ValueError
        If ``value`` is not scalar, finite, and positive.

    Examples
    --------
    >>> from astropy import units as u
    >>> positive_step(2 * u.s, "max_step")
    <Quantity 2. s>
    """
    result = readonly_quantity(value, u.s, name, ndim=(0,))
    if not np.isfinite(result.value) or result.value <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return result


def integration_options(
    *,
    method: object,
    rtol: object,
    atol: object,
    first_step: object,
    max_step: object,
) -> IntegrationOptions:
    """Validate and normalize integration settings.

    Parameters
    ----------
    method : str
        Non-empty method identifier.
    rtol : float
        Finite, non-negative relative tolerance.
    atol : float or array_like
        Finite, non-negative scalar or one-dimensional absolute tolerance.
    first_step, max_step : astropy.units.Quantity or None
        Optional finite, positive proper-time steps.

    Returns
    -------
    IntegrationOptions
        Normalized settings.

    Raises
    ------
    TypeError, astropy.units.UnitConversionError, ValueError
        If a setting has the wrong type, unit, shape, or value.
    """
    if not isinstance(method, str):
        raise TypeError("method must be a string")
    if not method:
        raise ValueError("method must not be empty")
    return IntegrationOptions(
        method=method,
        rtol=non_negative_real(rtol, "rtol"),
        atol=absolute_tolerance(atol),
        first_step=(
            None if first_step is None else positive_step(first_step, "first_step")
        ),
        max_step=None if max_step is None else positive_step(max_step, "max_step"),
    )
