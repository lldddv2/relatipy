"""Internal array, quantity, number, and type validation shared across subpackages."""

from __future__ import annotations

import math

import numpy as np
from astropy import units as u

_NDIM_DESCRIPTIONS = {
    (0,): "a scalar",
    (1,): "a one-dimensional series",
    (0, 1): "a scalar or a one-dimensional series",
}
REAL_TYPES = (int, float, np.integer, np.floating)


def immutable_array(value: object) -> np.ndarray:
    """Copy numeric data into immutable, contiguous storage.

    Parameters
    ----------
    value : object
        Scalar or array-like numeric data to copy.

    Returns
    -------
    numpy.ndarray
        A C-contiguous array backed by immutable bytes.  Its shape and dtype
        match the copied input.

    Raises
    ------
    TypeError
        If the array does not have a numeric dtype.

    Examples
    --------
    >>> array = immutable_array([1.0, 2.0])
    >>> array.flags.writeable
    False
    """
    array = np.asarray(value)
    if not np.issubdtype(array.dtype, np.number):
        raise TypeError("quantity values must be numeric")
    # ``tobytes`` always yields one C-ordered copy, so no other copy is needed.
    result = np.frombuffer(array.tobytes(), dtype=array.dtype).reshape(array.shape)
    result.setflags(write=False)
    return result


def readonly_quantity(
    value: u.Quantity,
    unit: u.UnitBase,
    name: str,
    *,
    ndim: tuple[int, ...] = (0, 1),
) -> u.Quantity:
    """Validate a quantity and copy it into immutable storage.

    Parameters
    ----------
    value : astropy.units.Quantity
        Quantity to validate and copy.
    unit : astropy.units.UnitBase
        Unit representing the required physical dimension.
    name : str
        Field name used in validation errors.
    ndim : tuple of int, optional
        Allowed numbers of array dimensions.  The default accepts a scalar or
        a one-dimensional series.

    Returns
    -------
    astropy.units.Quantity
        A read-only copy in the input unit.

    Raises
    ------
    TypeError
        If ``value`` is not an Astropy quantity or its data are not numeric.
    astropy.units.UnitConversionError
        If ``value`` is not dimensionally compatible with ``unit``.
    ValueError
        If the number of dimensions is not listed in ``ndim``.

    Examples
    --------
    >>> from astropy import units as u
    >>> quantity = readonly_quantity([1, 2] * u.km, u.m, "distance")
    >>> quantity.unit
    Unit("km")
    """
    if not isinstance(value, u.Quantity):
        raise TypeError(f"{name} must be an astropy.units.Quantity")
    if not value.unit.is_equivalent(unit):
        raise u.UnitConversionError(
            f"{name} must have units compatible with {unit}; got {value.unit}"
        )

    array = immutable_array(value.value)
    if array.ndim not in ndim:
        expected = _NDIM_DESCRIPTIONS.get(ndim, f"ndim in {ndim}")
        raise ValueError(f"{name} must be {expected}; got shape {array.shape}")
    return u.Quantity(array, value.unit, copy=False)


def require_instance(value: object, expected: type, name: str) -> None:
    """Require a field to have its declared runtime type.

    Parameters
    ----------
    value : object
        Value to inspect.
    expected : type
        Required runtime type.
    name : str
        Field name used in the error message.

    Raises
    ------
    TypeError
        If ``value`` is not an instance of ``expected``.

    Examples
    --------
    >>> require_instance(3, int, "value")
    """
    if not isinstance(value, expected):
        raise TypeError(f"{name} must be an instance of {expected.__name__}")


def dimensionless_scalar(value: object, name: str) -> float:
    """Convert a real number or dimensionless scalar quantity to ``float``.

    Parameters
    ----------
    value : object
        Python or NumPy real number, excluding booleans, or a scalar
        dimensionless Astropy quantity.
    name : str
        Field name used in validation errors.

    Returns
    -------
    float
        Converted value.  Finiteness and range are left to the caller.

    Raises
    ------
    TypeError
        If ``value`` is a boolean, not a real number, or a quantity with
        physical dimensions.
    ValueError
        If ``value`` is a non-scalar quantity.

    Examples
    --------
    >>> from astropy import units as u
    >>> dimensionless_scalar(50 * u.percent, "spin")
    0.5
    """
    if isinstance(value, u.Quantity):
        if not value.unit.is_equivalent(u.one):
            raise TypeError(f"{name} must be dimensionless; got unit {value.unit}")
        if value.ndim != 0:
            raise ValueError(f"{name} must be a scalar; got shape {value.shape}")
        return float(value.to_value(u.one))
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, REAL_TYPES):
        raise TypeError(f"{name} must be a real dimensionless scalar")
    return float(value)


def non_negative_real(value: object, name: str) -> float:
    """Convert a real scalar to a finite, non-negative float.

    Parameters
    ----------
    value : object
        Python or NumPy real number, excluding booleans.
    name : str
        Field name used in validation errors.

    Returns
    -------
    float
        Validated finite value.

    Raises
    ------
    TypeError
        If ``value`` is a boolean or not a real number.
    ValueError
        If ``value`` is non-finite or negative.

    Examples
    --------
    >>> non_negative_real(1.5, "rtol")
    1.5
    """
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, REAL_TYPES):
        raise TypeError(f"{name} must be a real number")
    result = float(value)
    if not math.isfinite(result) or result < 0.0:
        raise ValueError(f"{name} must be finite and non-negative")
    return result


def non_negative_integer(value: object, name: str) -> int:
    """Validate and normalize a non-negative integer diagnostic.

    Parameters
    ----------
    value : object
        Python or NumPy integer, excluding booleans.
    name : str
        Field name used in validation errors.

    Returns
    -------
    int
        The value normalized to a Python integer.

    Raises
    ------
    TypeError
        If ``value`` is not an integer or is a boolean.
    ValueError
        If ``value`` is negative.

    Examples
    --------
    >>> non_negative_integer(3, "n_steps")
    3
    """
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise TypeError(f"{name} must be an integer")
    result = int(value)
    if result < 0:
        raise ValueError(f"{name} must be non-negative")
    return result
