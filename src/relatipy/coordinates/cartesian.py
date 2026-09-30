"""Immutable Cartesian coordinate and velocity containers."""

from __future__ import annotations

from dataclasses import dataclass, field

from astropy import units as u

from .core import _stack_quantities, frozen_quantities


def _set_xyz(self: CartesianCoordinates) -> None:
    """Stack and freeze the derived Cartesian position."""
    object.__setattr__(
        self, "_xyz", _stack_quantities((self.x, self.y, self.z), name="xyz")
    )


@frozen_quantities(t=u.s, x=u.m, y=u.m, z=u.m, hook=_set_xyz)
@dataclass(frozen=True, slots=True, eq=False)
class CartesianCoordinates:
    """Coordinate time and Cartesian position.

    Parameters
    ----------
    t : astropy.units.Quantity
        Coordinate time, scalar or shape ``(n,)``, with units compatible with
        time.
    x, y, z : astropy.units.Quantity
        Cartesian lengths with the same shape as ``t`` and units compatible
        with length. Mutually compatible units may differ and are preserved.

    Attributes
    ----------
    t, x, y, z : astropy.units.Quantity
        Immutable copies of the supplied coordinate components, preserving
        compatible input units.
    xyz : astropy.units.Quantity
        Stacked Cartesian position with shape ``(3,)`` or ``(n, 3)``. Values
        are expressed in the unit of ``x``.

    Raises
    ------
    TypeError
        If a component is not an Astropy quantity.
    astropy.units.UnitConversionError
        If a component has incompatible units.
    ValueError
        If component shapes differ or are not scalar/one-dimensional.

    Notes
    -----
    The container stores supplied values without transforming them.
    Positions produced by :class:`~relatipy.geodesic.Orbit` and
    :class:`~relatipy.geodesic.Solution` use right-handed, spin-aligned
    *oblate* Cartesian axes centered on the black hole. They are related to
    Boyer--Lindquist coordinates ``(R, Theta, Phi)`` by
    ``x = sqrt(R**2 + a**2) sin(Theta) cos(Phi)``,
    ``y = sqrt(R**2 + a**2) sin(Theta) sin(Phi)`` and ``z = R cos(Theta)``,
    where ``a = spin * r_g`` (see :class:`~relatipy.metrics.Kerr`). They
    reduce to ordinary spherical-to-Cartesian axes only for ``spin = 0``.

    Examples
    --------
    >>> from astropy import units as u
    >>> point = CartesianCoordinates(0 * u.s, 1 * u.km, 2 * u.km, 3 * u.km)
    >>> point.xyz
    <Quantity [1., 2., 3.] km>
    """

    t: u.Quantity
    x: u.Quantity
    y: u.Quantity
    z: u.Quantity
    _xyz: u.Quantity = field(init=False, repr=False)

    @property
    def xyz(self) -> u.Quantity:
        """Return the stacked Cartesian position.

        Returns
        -------
        astropy.units.Quantity
            Read-only position with shape ``(3,)`` or ``(n, 3)`` in the unit
            of :attr:`x`.
        """
        return self._xyz


def _check_three_components(self: CartesianStateVector) -> None:
    """Require ``xyz`` (and therefore ``vxyz``) to end in three components."""
    if self.xyz.shape[-1] != 3:
        raise ValueError(
            f"xyz and vxyz must end in three components; got {self.xyz.shape}"
        )


@frozen_quantities(
    xyz=(u.m, (1, 2)),
    vxyz=(u.m / u.s, (1, 2)),
    hook=_check_three_components,
)
@dataclass(frozen=True, slots=True, eq=False)
class CartesianStateVector:
    """Cartesian phase-space state split into homogeneous blocks.

    Parameters
    ----------
    xyz : astropy.units.Quantity
        Position with shape ``(3,)`` or ``(n, 3)`` and units compatible with
        length.
    vxyz : astropy.units.Quantity
        Coordinate velocity ``dx^i/dt`` (derivative with respect to
        coordinate time) with the same shape as ``xyz`` and units compatible
        with length per time. States produced by RelatiPy use the oblate
        Cartesian axes described in :class:`CartesianCoordinates`.

    Attributes
    ----------
    xyz, vxyz : astropy.units.Quantity
        Immutable copies of the supplied homogeneous vectors, preserving
        compatible input units.

    Raises
    ------
    TypeError
        If either component is not an Astropy quantity.
    astropy.units.UnitConversionError
        If either component has an incompatible physical dimension.
    ValueError
        If shapes differ, are not one- or two-dimensional, or do not end in
        exactly three Cartesian components.

    Examples
    --------
    >>> from astropy import units as u
    >>> vector = CartesianStateVector(
    ...     [1, 2, 3] * u.km, [4, 5, 6] * u.km / u.s
    ... )
    >>> vector.vxyz.shape
    (3,)
    """

    xyz: u.Quantity
    vxyz: u.Quantity