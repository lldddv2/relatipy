"""Immutable parameters for a stable bound timelike Kerr geodesic."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from astropy import units as u

from .._validation import dimensionless_scalar
from .core import frozen_quantities


_PHASES = ("q_r0", "q_theta0", "q_phi0")


def _bounded_scalar(value: object, name: str, lower: float, upper: float) -> float:
    """Return a finite dimensionless scalar in the closed range ``[lower, upper]``."""
    number = dimensionless_scalar(value, name)
    if not (np.isfinite(number) and lower <= number <= upper):
        raise ValueError(f"{name} must be finite and in [{lower:g}, {upper:g}]")
    return number


def _validate_bound_elements(elements: KerrOrbitalElements) -> None:
    """Check parameter domains and store ``e`` and ``x`` as floats.

    Runs after the dimensional fields were converted to read-only scalar
    quantities. Kerr dynamics, stability and chart limits are checked in C.
    """
    for name in ("p", *_PHASES):
        if np.iscomplexobj(getattr(elements, name).value):
            raise TypeError(f"{name} must be real")

    p = elements.p.to_value(u.m)
    if not (np.isfinite(p) and p > 0):
        raise ValueError("p must be finite and positive")

    e = _bounded_scalar(elements.e, "e", 0.0, 1.0)
    if e == 1.0:
        raise ValueError("bound eccentricity e must satisfy 0 <= e < 1")
    x = _bounded_scalar(elements.x, "x", -1.0, 1.0)

    for name in _PHASES:
        if not np.isfinite(getattr(elements, name).to_value(u.rad)):
            raise ValueError(f"{name} must be finite")

    object.__setattr__(elements, "e", e)
    object.__setattr__(elements, "x", x)


@frozen_quantities(
    p=(u.m, (0,)),
    q_r0=(u.rad, (0,)),
    q_theta0=(u.rad, (0,)),
    q_phi0=(u.rad, (0,)),
    hook=_validate_bound_elements,
)
@dataclass(frozen=True, slots=True, eq=False)
class KerrOrbitalElements:
    """Scalar bound-geodesic parameters in the spin-aligned Kerr frame.

    Parameters
    ----------
    p : astropy.units.Quantity
        Positive semilatus rectum as a physical length. The native value is
        ``p / bh.r_g``; the radial turning points are ``p / (1 +/- e)``.
    e : float or astropy.units.Quantity
        Dimensionless eccentricity in ``[0, 1)``.
    x : float or astropy.units.Quantity
        Dimensionless inclination parameter in ``[-1, 1]``, using the
        KerrGeoPy convention ``sign(Lz) * sin(theta_min)`` away from polar
        orbits. Positive/negative values select prograde/retrograde motion;
        ``abs(x) == 1`` is equatorial and ``x == 0`` is polar.
    q_r0, q_theta0, q_phi0 : astropy.units.Quantity, optional
        Initial radial, polar and azimuthal Mino phases as finite angles,
        default zero. Radial/polar zero denotes the inner/northern turning
        point. These are not Keplerian true anomalies. The radial phase is
        redundant when ``e == 0``; the polar phase is redundant for an
        equatorial orbit. ``q_phi0`` is the initial Boyer--Lindquist azimuth.

    Attributes
    ----------
    p, q_r0, q_theta0, q_phi0 : astropy.units.Quantity
        Immutable scalar copies, preserving compatible input units.
    e, x : float
        Validated dimensionless values stored as plain floats.

    Raises
    ------
    TypeError
        If a dimensional field is not a real Astropy quantity, or ``e`` or
        ``x`` is neither a real number nor a dimensionless quantity.
    astropy.units.UnitConversionError
        If a field has incompatible units.
    ValueError
        If a field is non-scalar, non-finite or outside its parameter range.

    Notes
    -----
    This immutable value object checks units, shapes and parameter ranges.
    Stability and the Boyer--Lindquist chart domain are checked by the native
    constructor when passed to ``bh.orbit(elements=...)``. It does not store
    the mass, spin, initial time or osculating Kepler elements.
    Polar orbits can start away from the axis, but integration through the
    axis is not supported by the current Boyer--Lindquist chart.

    Examples
    --------
    >>> from astropy import units as u
    >>> elements = KerrOrbitalElements(10 * u.km, 0.3, 0.8)
    >>> elements.e, elements.q_r0
    (0.3, <Quantity 0. rad>)
    """

    p: u.Quantity
    e: float | u.Quantity
    x: float | u.Quantity
    q_r0: u.Quantity = field(default_factory=lambda: 0 * u.rad)
    q_theta0: u.Quantity = field(default_factory=lambda: 0 * u.rad)
    q_phi0: u.Quantity = field(default_factory=lambda: 0 * u.rad)
