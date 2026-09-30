"""Shared immutable quantity validation for coordinate value objects."""

from __future__ import annotations

from typing import Callable

import numpy as np
from astropy import units as u

from .._validation import immutable_array, readonly_quantity


def _require_same_shape(**values: u.Quantity) -> tuple[int, ...]:
    """Require all named quantities to have exactly the same shape.

    Parameters
    ----------
    **values : astropy.units.Quantity
        Named quantities to compare.  At least one value is required.

    Returns
    -------
    tuple of int
        The common shape.

    Raises
    ------
    ValueError
        If any quantity has a different shape.
    StopIteration
        If no named quantity is supplied.  Internal callers always provide at
        least one value.

    Examples
    --------
    >>> from astropy import units as u
    >>> _require_same_shape(x=[1, 2] * u.m, y=[3, 4] * u.s)
    (2,)
    """
    iterator = iter(values.items())
    first_name, first_value = next(iterator)
    expected = first_value.shape
    for name, value in iterator:
        if value.shape != expected:
            raise ValueError(
                f"{name} has shape {value.shape}, but {first_name} has shape "
                f"{expected}"
            )
    return expected


def _stack_quantities(
    values: tuple[u.Quantity, u.Quantity, u.Quantity],
    *,
    name: str,
) -> u.Quantity:
    """Stack three compatible quantities along a new final axis.

    Parameters
    ----------
    values : tuple of astropy.units.Quantity
        Exactly three quantities with mutually compatible units and shapes.
    name : str
        Group name used in a unit-conversion error.

    Returns
    -------
    astropy.units.Quantity
        Read-only stacked data expressed in the first component's unit, with
        shape ``(3,)`` or ``(n, 3)``.

    Raises
    ------
    astropy.units.UnitConversionError
        If a component cannot be converted to the first component's unit.
    ValueError
        If NumPy cannot stack the component shapes.

    Examples
    --------
    >>> from astropy import units as u
    >>> _stack_quantities((1 * u.km, 2000 * u.m, 3 * u.km), name="xyz")
    <Quantity [1., 2., 3.] km>
    """
    unit = values[0].unit
    try:
        array = immutable_array(
            np.stack([value.to_value(unit) for value in values], axis=-1)
        )
    except u.UnitConversionError as exc:
        raise u.UnitConversionError(
            f"all components of {name} must have compatible units"
        ) from exc
    return u.Quantity(array, unit, copy=False)


def frozen_quantities(
    *,
    hook: Callable[[object], None] | None = None,
    **specs: u.UnitBase | tuple[u.UnitBase, tuple[int, ...]],
) -> Callable[[type], type]:
    """Class decorator validating, converting, and freezing quantity fields.

    Wraps ``__init__`` so that, after the underlying dataclass constructor
    runs, every field named in ``specs`` is passed through
    :func:`~relatipy._validation.readonly_quantity` against its expected
    physical unit, then all of them are checked with
    :func:`_require_same_shape`.  Must be applied
    *after* ``@dataclass(frozen=True, eq=False)`` (i.e. written above it in
    source, so it decorates last), since it relies on the dataclass-generated
    ``__init__`` already existing on the class.  Use ``eq=False``: generated
    field-by-field equality cannot compare array quantities, so instances
    compare and hash by identity.

    Parameters
    ----------
    hook : callable, optional
        Called with the instance once all named fields are validated,
        converted, and frozen.  Use it to derive and freeze extra fields
        (e.g. a stacked vector) or to enforce checks beyond a shared shape.
    **specs : astropy.units.UnitBase or tuple
        Mapping of field name to its required physical unit, or to a
        ``(unit, ndim)`` pair overriding the default allowed number of
        dimensions accepted by :func:`~relatipy._validation.readonly_quantity`.

    Returns
    -------
    Callable[[type], type]
        A decorator that mutates and returns the given class.

    Raises
    ------
    TypeError
        At instantiation, if a named field is not an Astropy quantity.
    astropy.units.UnitConversionError
        At instantiation, if a named field is dimensionally incompatible
        with its declared unit.
    ValueError
        At instantiation, if the named fields do not all share the same
        shape, or if ``hook`` raises its own validation error.

    Examples
    --------
    >>> from dataclasses import dataclass
    >>> from astropy import units as u
    >>> @frozen_quantities(r=u.m, phi=u.rad)
    ... @dataclass(frozen=True, slots=True, eq=False)
    ... class Point:
    ...     r: u.Quantity
    ...     phi: u.Quantity
    >>> Point(4 * u.km, 1 * u.rad).r
    <Quantity 4. km>
    """

    def decorator(cls: type) -> type:
        """Install quantity validation on one frozen dataclass type."""
        original_init = cls.__init__

        def __init__(self, *args: object, **kwargs: object) -> None:
            """Initialize fields, then validate, convert, and freeze them."""
            original_init(self, *args, **kwargs)
            values = {}
            for name, spec in specs.items():
                unit, ndim = spec if isinstance(spec, tuple) else (spec, (0, 1))
                values[name] = readonly_quantity(
                    getattr(self, name), unit, name, ndim=ndim
                )
            _require_same_shape(**values)
            for name, value in values.items():
                object.__setattr__(self, name, value)
            if hook is not None:
                hook(self)

        cls.__init__ = __init__
        return cls

    return decorator
