"""Store immutable classical orbital elements.

This module validates and stores orbital elements that were computed
elsewhere.  It does not derive elements from a state or select a dynamical
convention.

Examples
--------
>>> from astropy import units as u
>>> from relatipy.coordinates import OrbitalElements
>>> elements = OrbitalElements(
...     10 * u.km, 0.2, 1 * u.rad, 2 * u.rad, 3 * u.rad, 4 * u.rad
... )
>>> elements.e
0.2
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from astropy import units as u

from .._validation import immutable_array
from .core import frozen_quantities


def _readonly_eccentricity(value: object) -> float | np.ndarray:
    """Validate and copy dimensionless eccentricity values.

    Parameters
    ----------
    value : object
        A numeric scalar or one-dimensional series.  An Astropy quantity must
        be dimensionless.

    Returns
    -------
    float or numpy.ndarray
        A ``float`` for scalar input or a read-only array for a series.

    Raises
    ------
    astropy.units.UnitConversionError
        If an Astropy quantity is not dimensionless.
    TypeError
        If a quantity holds non-numeric data.
    ValueError
        If the input cannot be converted to floating-point data or has more
        than one dimension.

    Examples
    --------
    >>> _readonly_eccentricity([0.1, 0.2]).flags.writeable
    False
    """
    if isinstance(value, u.Quantity):
        if not value.unit.is_equivalent(u.one):
            raise u.UnitConversionError(
                f"e must be dimensionless; got units of {value.unit}"
            )
        array = immutable_array(value.to_value(u.one))
    else:
        array = immutable_array(np.asarray(value, dtype=float))

    if array.ndim not in (0, 1):
        raise ValueError(
            f"e must be a scalar or one-dimensional series; got shape {array.shape}"
        )
    if array.ndim == 0:
        return float(array)
    return array


def _set_eccentricity(self: OrbitalElements) -> None:
    """Validate ``e`` against the already-frozen angular shape, then freeze it."""
    eccentricity = _readonly_eccentricity(self.e)
    shape = self.a.shape
    if np.shape(eccentricity) != shape:
        raise ValueError(
            f"e has shape {np.shape(eccentricity)}, but a has shape {shape}"
        )
    object.__setattr__(self, "e", eccentricity)


@frozen_quantities(
    a=u.m,
    inc=u.rad,
    Omega=u.rad,
    omega=u.rad,
    f=u.rad,
    hook=_set_eccentricity,
)
@dataclass(frozen=True, slots=True, eq=False)
class OrbitalElements:
    """Classical orbital elements for one state or a state series.

    Parameters
    ----------
    a : astropy.units.Quantity
        Semi-major axis as a length quantity, scalar or shape ``(n,)``.
        Positive infinity represents an exactly parabolic osculating conic.
    e : float, numpy.ndarray, or astropy.units.Quantity
        Dimensionless eccentricity.  A scalar remains a ``float`` and a
        series is stored as a read-only NumPy array.  A quantity must have
        dimensionless units.
    inc, Omega, omega, f : astropy.units.Quantity
        Inclination, longitude of ascending node, argument of periapsis and
        true anomaly as angular quantities with the same shape as ``a``.
        Native reconstruction uses ``NaN`` for all four angles when angular
        momentum is too small to define an orbital plane.

    Attributes
    ----------
    a, e, inc, Omega, omega, f
        Immutable stored values.  Compatible input units are preserved except
        that a dimensionless quantity supplied for ``e`` is stored as a
        ``float`` or NumPy array.
    defined : numpy.ndarray
        Read-only boolean mask with shape ``(6,)`` or ``(n, 6)`` in the order
        ``(a, e, inc, Omega, omega, f)``. Finite values and positive infinity
        for ``a`` are defined; ``NaN`` and other infinities are not.

    Raises
    ------
    TypeError
        If a dimensional field is not an Astropy quantity.
    astropy.units.UnitConversionError
        If a field has incompatible units.
    ValueError
        If fields do not all have one common scalar or series shape.

    Notes
    -----
    This object stores elements computed elsewhere.  It does not derive them
    from a state vector and does not impose dynamical conventions beyond
    units, shape and immutability.

    Examples
    --------
    >>> from astropy import units as u
    >>> elements = OrbitalElements(
    ...     10 * u.km, 0.2, 1 * u.deg, 2 * u.deg, 3 * u.deg, 4 * u.deg
    ... )
    >>> elements.a
    <Quantity 10. km>
    >>> elements.e
    0.2
    """

    a: u.Quantity
    e: float | np.ndarray | u.Quantity
    inc: u.Quantity
    Omega: u.Quantity
    omega: u.Quantity
    f: u.Quantity

    @property
    def defined(self) -> np.ndarray:
        """Return which stored elements are defined, in element order.

        Returns
        -------
        numpy.ndarray
            Immutable boolean mask of shape ``(6,)`` for scalar elements or
            ``(n, 6)`` for a series. Columns are ``a``, ``e``, ``inc``,
            ``Omega``, ``omega``, and ``f``. Positive infinite ``a`` is
            defined for an exactly parabolic conic.

        Notes
        -----
        The mask reports stored values; it does not validate a state or
        impose orbital conventions. Its immutable bytes backing prevents
        enabling writes with ``setflags(write=True)``.
        """
        axis = self.a.value
        mask = np.stack(
            (
                np.isfinite(axis) | np.isposinf(axis),
                np.isfinite(self.e),
                np.isfinite(self.inc.value),
                np.isfinite(self.Omega.value),
                np.isfinite(self.omega.value),
                np.isfinite(self.f.value),
            ),
            axis=-1,
        )
        return np.frombuffer(mask.tobytes(), dtype=np.bool_).reshape(mask.shape)