"""Immutable Boyer--Lindquist coordinate and velocity containers."""

from __future__ import annotations

from dataclasses import dataclass

from astropy import units as u

from .core import frozen_quantities


@frozen_quantities(t=u.s, R=u.m, Theta=u.rad, Phi=u.rad)
@dataclass(frozen=True, slots=True, eq=False)
class BoyerLindquistCoordinates:
    """Store coordinate time and Boyer--Lindquist position.

    Parameters
    ----------
    t : astropy.units.Quantity
        Coordinate time, scalar or shape ``(n,)``.
    R : astropy.units.Quantity
        Boyer--Lindquist radial length with the same shape as ``t``.
    Theta, Phi : astropy.units.Quantity
        Boyer--Lindquist polar and azimuthal angles with the same shape as
        ``t``.

    Raises
    ------
    TypeError
        If a component is not an Astropy quantity.
    astropy.units.UnitConversionError
        If a component has an incompatible physical dimension.
    ValueError
        If component shapes differ or are not scalar/one-dimensional.

    Notes
    -----
    The container stores supplied values without transforming them. Views
    produced by :class:`~relatipy.geodesic.Orbit` and
    :class:`~relatipy.geodesic.Solution` hold the Boyer--Lindquist chart in
    which the native Kerr geodesic is integrated; ``R`` is a coordinate
    radius, not a Euclidean distance. For nonzero spin it differs from the
    spherical ``r`` of :class:`SphericalCoordinates`.

    Examples
    --------
    >>> from astropy import units as u
    >>> point = BoyerLindquistCoordinates(
    ...     0 * u.s, 4 * u.km, 1 * u.rad, 2 * u.rad
    ... )
    >>> point.R
    <Quantity 4. km>
    """

    t: u.Quantity
    R: u.Quantity
    Theta: u.Quantity
    Phi: u.Quantity


@frozen_quantities(vR=u.m / u.s, vTheta=u.rad / u.s, vPhi=u.rad / u.s)
@dataclass(frozen=True, slots=True, eq=False)
class BoyerLindquistVelocity:
    """Store Boyer--Lindquist coordinate-velocity components.

    Parameters
    ----------
    vR : astropy.units.Quantity
        Radial velocity, scalar or shape ``(n,)``, with units compatible with
        length per time.
    vTheta, vPhi : astropy.units.Quantity
        Polar and azimuthal angular velocities with the same shape as ``vR``
        and units compatible with angle per time.

    Raises
    ------
    TypeError
        If a component is not an Astropy quantity.
    astropy.units.UnitConversionError
        If a component has an incompatible physical dimension.
    ValueError
        If component shapes differ or are not scalar/one-dimensional.

    Examples
    --------
    >>> from astropy import units as u
    >>> velocity = BoyerLindquistVelocity(
    ...     3 * u.km / u.s, 0.1 * u.rad / u.s, 0.2 * u.rad / u.s
    ... )
    >>> velocity.vR
    <Quantity 3. km / s>
    """

    vR: u.Quantity
    vTheta: u.Quantity
    vPhi: u.Quantity


@frozen_quantities(uR=u.m / u.s, uTheta=u.rad / u.s, uPhi=u.rad / u.s)
@dataclass(frozen=True, slots=True, eq=False)
class BoyerLindquistFourVelocity:
    """Store spatial Boyer--Lindquist four-velocity components.

    Parameters
    ----------
    uR : astropy.units.Quantity
        Proper-time derivative of radial length, scalar or shape ``(n,)``.
    uTheta, uPhi : astropy.units.Quantity
        Proper-time derivatives of the polar and azimuthal angles with the
        same shape as ``uR``.

    Raises
    ------
    TypeError
        If a component is not an Astropy quantity.
    astropy.units.UnitConversionError
        If a component has an incompatible physical dimension.
    ValueError
        If component shapes differ or are not scalar/one-dimensional.

    Examples
    --------
    >>> from astropy import units as u
    >>> velocity = BoyerLindquistFourVelocity(
    ...     3 * u.km / u.s, 0.1 * u.rad / u.s, 0.2 * u.rad / u.s
    ... )
    >>> velocity.uR
    <Quantity 3. km / s>
    """

    uR: u.Quantity
    uTheta: u.Quantity
    uPhi: u.Quantity