"""Immutable spherical coordinate and velocity containers."""

from __future__ import annotations

from dataclasses import dataclass

from astropy import units as u

from .core import frozen_quantities


@frozen_quantities(t=u.s, r=u.m, theta=u.rad, phi=u.rad)
@dataclass(frozen=True, slots=True, eq=False)
class SphericalCoordinates:
    """Store coordinate time and spherical position.

    Parameters
    ----------
    t : astropy.units.Quantity
        Coordinate time, scalar or shape ``(n,)``.
    r : astropy.units.Quantity
        Radial length with the same shape as ``t``.
    theta, phi : astropy.units.Quantity
        Polar and azimuthal angles with the same shape as ``t``.

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
    :class:`~relatipy.geodesic.Solution` are Euclidean spherical coordinates
    of the oblate Cartesian position (see
    :class:`CartesianCoordinates`): ``r = sqrt(x**2 + y**2 + z**2)``,
    ``theta = atan2(sqrt(x**2 + y**2), z)`` and ``phi`` equal to the
    Boyer--Lindquist azimuth. For nonzero spin, ``r`` and ``theta`` differ
    from the Boyer--Lindquist ``R`` and ``Theta`` of
    :class:`BoyerLindquistCoordinates`; the two are not interchangeable.

    Examples
    --------
    >>> from astropy import units as u
    >>> point = SphericalCoordinates(0 * u.s, 4 * u.km, 1 * u.rad, 2 * u.rad)
    >>> point.r
    <Quantity 4. km>
    """

    t: u.Quantity
    r: u.Quantity
    theta: u.Quantity
    phi: u.Quantity


@frozen_quantities(vr=u.m / u.s, vtheta=u.rad / u.s, vphi=u.rad / u.s)
@dataclass(frozen=True, slots=True, eq=False)
class SphericalVelocity:
    """Store spherical coordinate-velocity components.

    Parameters
    ----------
    vr : astropy.units.Quantity
        Radial velocity, scalar or shape ``(n,)``.
    vtheta, vphi : astropy.units.Quantity
        Polar and azimuthal angular velocities with the same shape as ``vr``.

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
    >>> velocity = SphericalVelocity(3 * u.km / u.s, 0.1 * u.rad / u.s,
    ...                              0.2 * u.rad / u.s)
    >>> velocity.vr
    <Quantity 3. km / s>
    """

    vr: u.Quantity
    vtheta: u.Quantity
    vphi: u.Quantity


@frozen_quantities(ur=u.m / u.s, utheta=u.rad / u.s, uphi=u.rad / u.s)
@dataclass(frozen=True, slots=True, eq=False)
class SphericalFourVelocity:
    """Store spatial spherical four-velocity components.

    Parameters
    ----------
    ur : astropy.units.Quantity
        Proper-time derivative of radial length, scalar or shape ``(n,)``.
    utheta, uphi : astropy.units.Quantity
        Proper-time derivatives of the polar and azimuthal angles.

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
    >>> velocity = SphericalFourVelocity(
    ...     3 * u.km / u.s, 0.1 * u.rad / u.s, 0.2 * u.rad / u.s
    ... )
    >>> velocity.ur
    <Quantity 3. km / s>
    """

    ur: u.Quantity
    utheta: u.Quantity
    uphi: u.Quantity