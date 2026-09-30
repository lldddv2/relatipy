"""Define immutable geodesic state records.

The objects in this module assemble already computed coordinate,
coordinate-velocity, four-velocity, and orbital-element representations into
one validated state.  They check units, sample cardinality, and consistency,
but do not evaluate physical transformations or equations of motion.

Notes
-----
One state uses scalar component quantities and Cartesian vectors of shape
``(3,)``.  A series of ``n`` states uses component shape ``(n,)`` and vector
shape ``(n, 3)``.  Inputs are copied into storage that is read-only.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from astropy import units as u

from .._validation import readonly_quantity, require_instance
from ..coordinates import (
    BoyerLindquistCoordinates,
    BoyerLindquistFourVelocity,
    BoyerLindquistVelocity,
    CartesianCoordinates,
    CartesianStateVector,
    OrbitalElements,
    SphericalCoordinates,
    SphericalFourVelocity,
    SphericalVelocity,
)


def _component_shape(value: u.Quantity) -> tuple[int, ...]:
    """Return the sample shape of a scalar component quantity.

    Parameters
    ----------
    value : astropy.units.Quantity
        Component whose shape is required.

    Returns
    -------
    tuple of int
        The quantity shape.

    Examples
    --------
    >>> from astropy import units as u
    >>> _component_shape([1, 2] * u.s)
    (2,)
    """
    return value.shape


def _require_state_shape(
    expected: tuple[int, ...],
    *,
    name: str,
    actual: tuple[int, ...],
) -> None:
    """Require one representation to match the state's sample shape.

    Parameters
    ----------
    expected : tuple of int
        Required sample shape.
    name : str
        Representation name used in the error message.
    actual : tuple of int
        Observed sample shape.

    Raises
    ------
    ValueError
        If ``actual`` and ``expected`` differ.

    Examples
    --------
    >>> _require_state_shape((2,), name="position", actual=(2,))
    """
    if actual != expected:
        raise ValueError(f"{name} has state shape {actual}, expected {expected}")


def _require_same_time(reference: u.Quantity, other: u.Quantity, name: str) -> None:
    """Require two coordinate-time quantities to represent identical values.

    Parameters
    ----------
    reference, other : astropy.units.Quantity
        Scalar or one-dimensional coordinate times.  Values are compared
        after expressing ``reference`` in the unit of ``other``.
    name : str
        Name of ``other`` used in the error message.

    Raises
    ------
    astropy.units.UnitConversionError
        If the coordinate-time units are incompatible.
    ValueError
        If the converted values differ.

    Examples
    --------
    >>> from astropy import units as u
    >>> _require_same_time(1 * u.s, 1000 * u.ms, "other")
    """
    if not np.array_equal(reference.to_value(other.unit), other.value):
        raise ValueError(f"{name} must describe the same coordinate time as txyz.t")


def _readonly_spatial_vector(
    value: u.Quantity,
    unit: u.UnitBase,
    name: str,
) -> u.Quantity:
    """Validate an immutable Cartesian vector or vector series.

    Parameters
    ----------
    value : astropy.units.Quantity
        Vector with shape ``(3,)`` or ``(n, 3)``.
    unit : astropy.units.UnitBase
        Unit representing the required physical dimension.
    name : str
        Field name used in validation errors.

    Returns
    -------
    astropy.units.Quantity
        A read-only copy in the input unit.

    Raises
    ------
    TypeError
        If ``value`` is not an Astropy quantity.
    astropy.units.UnitConversionError
        If ``value`` is dimensionally incompatible with ``unit``.
    ValueError
        If ``value`` is not one- or two-dimensional or its final dimension is
        not three.

    Examples
    --------
    >>> from astropy import units as u
    >>> _readonly_spatial_vector([1, 2, 3] * u.km, u.m, "xyz").shape
    (3,)
    """
    result = readonly_quantity(value, unit, name, ndim=(1, 2))
    if result.shape[-1] != 3:
        raise ValueError(f"{name} must end in three components; got {result.shape}")
    return result


def _vector_state_shape(value: u.Quantity) -> tuple[int, ...]:
    """Return the sample shape represented by a Cartesian vector array.

    Parameters
    ----------
    value : astropy.units.Quantity
        Validated vector with shape ``(3,)`` or ``(n, 3)``.

    Returns
    -------
    tuple of int
        ``()`` for one vector or ``(n,)`` for a vector series.

    Examples
    --------
    >>> from astropy import units as u
    >>> _vector_state_shape([[1, 2, 3], [4, 5, 6]] * u.km)
    (2,)
    """
    return () if value.ndim == 1 else value.shape[:-1]


@dataclass(frozen=True, slots=True, init=False, eq=False)
class State:
    """One immutable orbital state or a vectorized series of states.

    Parameters
    ----------
    tau : astropy.units.Quantity
        Proper time, scalar or shape ``(n,)``, with units compatible with
        time.
    txyz : relatipy.coordinates.CartesianCoordinates
    trqp : relatipy.coordinates.SphericalCoordinates
    tRQP : relatipy.coordinates.BoyerLindquistCoordinates
        Cartesian, spherical and Boyer--Lindquist coordinate views containing
        values that were calculated elsewhere.
    vxyz : astropy.units.Quantity
        Cartesian coordinate velocity ``dx^i/dt`` (derivative with respect
        to coordinate time ``t``) with shape ``(3,)`` or ``(n, 3)`` and
        units compatible with length per time.
    vrqp : relatipy.coordinates.SphericalVelocity
    vRQP : relatipy.coordinates.BoyerLindquistVelocity
        Spherical and Boyer--Lindquist coordinate-velocity views, also
        derivatives with respect to coordinate time ``t``.
    ut : astropy.units.Quantity
        Dimensionless temporal four-velocity component, scalar or shape
        ``(n,)``.
    uxyz : astropy.units.Quantity
        Cartesian spatial four-velocity ``dx^i/dtau`` (derivative with
        respect to proper time ``tau``) with shape ``(3,)`` or ``(n, 3)``
        and units compatible with length per time.
    urqp : relatipy.coordinates.SphericalFourVelocity
    uRQP : relatipy.coordinates.BoyerLindquistFourVelocity
        Spherical and Boyer--Lindquist spatial four-velocity views, also
        derivatives with respect to proper time ``tau``.
    orbital_elements : relatipy.coordinates.OrbitalElements
        Classical osculating orbital elements calculated elsewhere.

    Attributes
    ----------
    tau, txyz, trqp, tRQP, vxyz, vrqp, vRQP, ut, uxyz, urqp, uRQP
        Immutable validated views described in ``Parameters``.
    t, x, y, z, xyz, r, theta, phi, R, Theta, Phi
        Read-only coordinate delegates.
    vx, vy, vz, vr, vtheta, vphi, vR, vTheta, vPhi
        Read-only coordinate-velocity delegates (derivatives with respect to
        coordinate time ``t``).
    ux, uy, uz, ur, utheta, uphi, uR, uTheta, uPhi
        Read-only spatial four-velocity delegates (derivatives with respect
        to proper time ``tau``).
    state_vector : relatipy.coordinates.CartesianStateVector
        Read-only Cartesian position and coordinate velocity.

    Raises
    ------
    TypeError
        If a typed container has the wrong type or a dimensional value is not
        an Astropy quantity.
    astropy.units.UnitConversionError
        If ``tau``, ``ut`` or a vector has incompatible units.
    ValueError
        If representations have inconsistent scalar/vectorized shapes or
        coordinate times.

    Notes
    -----
    The constructor is a boundary for already computed representations.  It
    verifies units, cardinality and consistency but performs no physical
    coordinate or orbital-element transformations.

    Examples
    --------
    >>> from astropy import units as u
    >>> from relatipy.coordinates import (BoyerLindquistCoordinates,
    ...     BoyerLindquistFourVelocity, BoyerLindquistVelocity,
    ...     CartesianCoordinates, SphericalCoordinates,
    ...     SphericalFourVelocity, SphericalVelocity)
    >>> from relatipy.coordinates import OrbitalElements
    >>> t = 0 * u.s
    >>> state = State(
    ...     tau=t,
    ...     txyz=CartesianCoordinates(t, 1 * u.km, 2 * u.km, 3 * u.km),
    ...     trqp=SphericalCoordinates(t, 4 * u.km, 1 * u.rad, 2 * u.rad),
    ...     tRQP=BoyerLindquistCoordinates(t, 4 * u.km, 1 * u.rad, 2 * u.rad),
    ...     vxyz=[1, 2, 3] * u.km / u.s,
    ...     vrqp=SphericalVelocity(1 * u.km / u.s, 2 * u.rad / u.s,
    ...                               3 * u.rad / u.s),
    ...     vRQP=BoyerLindquistVelocity(1 * u.km / u.s, 2 * u.rad / u.s,
    ...                                  3 * u.rad / u.s),
    ...     ut=1 * u.one, uxyz=[1, 2, 3] * u.km / u.s,
    ...     urqp=SphericalFourVelocity(1 * u.km / u.s, 2 * u.rad / u.s,
    ...                                  3 * u.rad / u.s),
    ...     uRQP=BoyerLindquistFourVelocity(1 * u.km / u.s,
    ...                                     2 * u.rad / u.s,
    ...                                     3 * u.rad / u.s),
    ...     orbital_elements=OrbitalElements(10 * u.km, 0.2, 1 * u.rad,
    ...                                       2 * u.rad, 3 * u.rad, 4 * u.rad),
    ... )
    >>> state.xyz.shape
    (3,)
    """

    tau: u.Quantity
    txyz: CartesianCoordinates
    trqp: SphericalCoordinates
    tRQP: BoyerLindquistCoordinates
    vxyz: u.Quantity
    vrqp: SphericalVelocity
    vRQP: BoyerLindquistVelocity
    ut: u.Quantity
    uxyz: u.Quantity
    urqp: SphericalFourVelocity
    uRQP: BoyerLindquistFourVelocity
    _orbital_elements: OrbitalElements = field(repr=False)
    _state_vector: CartesianStateVector = field(repr=False)

    def __init__(
        self,
        *,
        tau: u.Quantity,
        txyz: CartesianCoordinates,
        trqp: SphericalCoordinates,
        tRQP: BoyerLindquistCoordinates,
        vxyz: u.Quantity,
        vrqp: SphericalVelocity,
        vRQP: BoyerLindquistVelocity,
        ut: u.Quantity,
        uxyz: u.Quantity,
        urqp: SphericalFourVelocity,
        uRQP: BoyerLindquistFourVelocity,
        orbital_elements: OrbitalElements,
    ) -> None:
        """Initialize a state from precomputed, mutually consistent views."""
        require_instance(txyz, CartesianCoordinates, "txyz")
        require_instance(trqp, SphericalCoordinates, "trqp")
        require_instance(tRQP, BoyerLindquistCoordinates, "tRQP")
        require_instance(vrqp, SphericalVelocity, "vrqp")
        require_instance(vRQP, BoyerLindquistVelocity, "vRQP")
        require_instance(urqp, SphericalFourVelocity, "urqp")
        require_instance(uRQP, BoyerLindquistFourVelocity, "uRQP")
        require_instance(orbital_elements, OrbitalElements, "orbital_elements")

        tau_value = readonly_quantity(tau, u.s, "tau")
        ut_value = readonly_quantity(ut, u.one, "ut")
        vxyz_value = _readonly_spatial_vector(vxyz, u.m / u.s, "vxyz")
        uxyz_value = _readonly_spatial_vector(uxyz, u.m / u.s, "uxyz")
        state_shape = tau_value.shape

        shapes = {
            "txyz": _component_shape(txyz.t),
            "trqp": _component_shape(trqp.t),
            "tRQP": _component_shape(tRQP.t),
            "vxyz": _vector_state_shape(vxyz_value),
            "vrqp": _component_shape(vrqp.vr),
            "vRQP": _component_shape(vRQP.vR),
            "ut": _component_shape(ut_value),
            "uxyz": _vector_state_shape(uxyz_value),
            "urqp": _component_shape(urqp.ur),
            "uRQP": _component_shape(uRQP.uR),
            "orbital_elements": _component_shape(orbital_elements.a),
        }
        for name, shape in shapes.items():
            _require_state_shape(state_shape, name=name, actual=shape)

        _require_same_time(txyz.t, trqp.t, "trqp.t")
        _require_same_time(txyz.t, tRQP.t, "tRQP.t")

        object.__setattr__(self, "tau", tau_value)
        object.__setattr__(self, "txyz", txyz)
        object.__setattr__(self, "trqp", trqp)
        object.__setattr__(self, "tRQP", tRQP)
        object.__setattr__(self, "vxyz", vxyz_value)
        object.__setattr__(self, "vrqp", vrqp)
        object.__setattr__(self, "vRQP", vRQP)
        object.__setattr__(self, "ut", ut_value)
        object.__setattr__(self, "uxyz", uxyz_value)
        object.__setattr__(self, "urqp", urqp)
        object.__setattr__(self, "uRQP", uRQP)
        object.__setattr__(self, "_orbital_elements", orbital_elements)
        object.__setattr__(
            self,
            "_state_vector",
            CartesianStateVector(xyz=txyz.xyz, vxyz=vxyz_value),
        )

    @property
    def t(self) -> u.Quantity:
        """Return coordinate time.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar or shape ``(n,)`` coordinate time.

        Examples
        --------
        >>> from types import SimpleNamespace
        >>> from astropy import units as u
        >>> proxy = SimpleNamespace(txyz=SimpleNamespace(t=1 * u.s))
        >>> bool(State.t.fget(proxy) == 1 * u.s)
        True
        """
        return self.txyz.t

    @property
    def x(self) -> u.Quantity:
        """Return the Cartesian x coordinate.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar or shape ``(n,)`` length.

        Examples
        --------
        >>> from types import SimpleNamespace
        >>> from astropy import units as u
        >>> proxy = SimpleNamespace(txyz=SimpleNamespace(x=1 * u.km))
        >>> bool(State.x.fget(proxy) == 1 * u.km)
        True
        """
        return self.txyz.x

    @property
    def y(self) -> u.Quantity:
        """Return the Cartesian y coordinate.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar or shape ``(n,)`` length.

        Examples
        --------
        >>> from types import SimpleNamespace
        >>> from astropy import units as u
        >>> proxy = SimpleNamespace(txyz=SimpleNamespace(y=2 * u.km))
        >>> bool(State.y.fget(proxy) == 2 * u.km)
        True
        """
        return self.txyz.y

    @property
    def z(self) -> u.Quantity:
        """Return the Cartesian z coordinate.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar or shape ``(n,)`` length.

        Examples
        --------
        >>> from types import SimpleNamespace
        >>> from astropy import units as u
        >>> proxy = SimpleNamespace(txyz=SimpleNamespace(z=3 * u.km))
        >>> bool(State.z.fget(proxy) == 3 * u.km)
        True
        """
        return self.txyz.z

    @property
    def xyz(self) -> u.Quantity:
        """Return the Cartesian position vector.

        Returns
        -------
        astropy.units.Quantity
            Read-only position with shape ``(3,)`` or ``(n, 3)``.

        Examples
        --------
        >>> from types import SimpleNamespace
        >>> from astropy import units as u
        >>> proxy = SimpleNamespace(txyz=SimpleNamespace(xyz=[1, 2, 3] * u.km))
        >>> State.xyz.fget(proxy).shape
        (3,)
        """
        return self.txyz.xyz

    @property
    def r(self) -> u.Quantity:
        """Return the spherical radial coordinate.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar or shape ``(n,)`` length.

        Examples
        --------
        >>> from types import SimpleNamespace
        >>> from astropy import units as u
        >>> proxy = SimpleNamespace(trqp=SimpleNamespace(r=4 * u.km))
        >>> bool(State.r.fget(proxy) == 4 * u.km)
        True
        """
        return self.trqp.r

    @property
    def theta(self) -> u.Quantity:
        """Return the spherical polar angle.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar or shape ``(n,)`` angle.

        Examples
        --------
        >>> from types import SimpleNamespace
        >>> from astropy import units as u
        >>> proxy = SimpleNamespace(trqp=SimpleNamespace(theta=1 * u.rad))
        >>> bool(State.theta.fget(proxy) == 1 * u.rad)
        True
        """
        return self.trqp.theta

    @property
    def phi(self) -> u.Quantity:
        """Return the spherical azimuthal angle.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar or shape ``(n,)`` angle.

        Examples
        --------
        >>> from types import SimpleNamespace
        >>> from astropy import units as u
        >>> proxy = SimpleNamespace(trqp=SimpleNamespace(phi=2 * u.rad))
        >>> bool(State.phi.fget(proxy) == 2 * u.rad)
        True
        """
        return self.trqp.phi

    @property
    def R(self) -> u.Quantity:
        """Return the Boyer--Lindquist radial coordinate.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar or shape ``(n,)`` length.

        Examples
        --------
        >>> from types import SimpleNamespace
        >>> from astropy import units as u
        >>> proxy = SimpleNamespace(tRQP=SimpleNamespace(R=5 * u.km))
        >>> bool(State.R.fget(proxy) == 5 * u.km)
        True
        """
        return self.tRQP.R

    @property
    def Theta(self) -> u.Quantity:
        """Return the Boyer--Lindquist polar angle.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar or shape ``(n,)`` angle.

        Examples
        --------
        >>> from types import SimpleNamespace
        >>> from astropy import units as u
        >>> proxy = SimpleNamespace(tRQP=SimpleNamespace(Theta=1 * u.rad))
        >>> bool(State.Theta.fget(proxy) == 1 * u.rad)
        True
        """
        return self.tRQP.Theta

    @property
    def Phi(self) -> u.Quantity:
        """Return the Boyer--Lindquist azimuthal angle.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar or shape ``(n,)`` angle.

        Examples
        --------
        >>> from types import SimpleNamespace
        >>> from astropy import units as u
        >>> proxy = SimpleNamespace(tRQP=SimpleNamespace(Phi=2 * u.rad))
        >>> bool(State.Phi.fget(proxy) == 2 * u.rad)
        True
        """
        return self.tRQP.Phi

    @property
    def vx(self) -> u.Quantity:
        """Return the Cartesian x coordinate velocity.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar or shape ``(n,)`` derivative with respect to
            coordinate time ``t``, in length per time.

        Examples
        --------
        >>> from types import SimpleNamespace
        >>> from astropy import units as u
        >>> proxy = SimpleNamespace(vxyz=[1, 2, 3] * u.km / u.s)
        >>> bool(State.vx.fget(proxy) == 1 * u.km / u.s)
        True
        """
        return self.vxyz[..., 0]

    @property
    def vy(self) -> u.Quantity:
        """Return the Cartesian y coordinate velocity.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar or shape ``(n,)`` derivative with respect to
            coordinate time ``t``, in length per time.

        Examples
        --------
        >>> from types import SimpleNamespace
        >>> from astropy import units as u
        >>> proxy = SimpleNamespace(vxyz=[1, 2, 3] * u.km / u.s)
        >>> bool(State.vy.fget(proxy) == 2 * u.km / u.s)
        True
        """
        return self.vxyz[..., 1]

    @property
    def vz(self) -> u.Quantity:
        """Return the Cartesian z coordinate velocity.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar or shape ``(n,)`` derivative with respect to
            coordinate time ``t``, in length per time.

        Examples
        --------
        >>> from types import SimpleNamespace
        >>> from astropy import units as u
        >>> proxy = SimpleNamespace(vxyz=[1, 2, 3] * u.km / u.s)
        >>> bool(State.vz.fget(proxy) == 3 * u.km / u.s)
        True
        """
        return self.vxyz[..., 2]

    @property
    def vr(self) -> u.Quantity:
        """Return the spherical radial coordinate velocity.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar or shape ``(n,)`` derivative with respect to
            coordinate time ``t``, in length per time.

        Examples
        --------
        >>> from types import SimpleNamespace
        >>> from astropy import units as u
        >>> proxy = SimpleNamespace(vrqp=SimpleNamespace(vr=1 * u.km / u.s))
        >>> bool(State.vr.fget(proxy) == 1 * u.km / u.s)
        True
        """
        return self.vrqp.vr

    @property
    def vtheta(self) -> u.Quantity:
        """Return the spherical polar coordinate velocity.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar or shape ``(n,)`` derivative with respect to
            coordinate time ``t``, in angle per time.

        Examples
        --------
        >>> from types import SimpleNamespace
        >>> from astropy import units as u
        >>> proxy = SimpleNamespace(vrqp=SimpleNamespace(vtheta=1 * u.rad / u.s))
        >>> bool(State.vtheta.fget(proxy) == 1 * u.rad / u.s)
        True
        """
        return self.vrqp.vtheta

    @property
    def vphi(self) -> u.Quantity:
        """Return the spherical azimuthal coordinate velocity.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar or shape ``(n,)`` derivative with respect to
            coordinate time ``t``, in angle per time.

        Examples
        --------
        >>> from types import SimpleNamespace
        >>> from astropy import units as u
        >>> proxy = SimpleNamespace(vrqp=SimpleNamespace(vphi=2 * u.rad / u.s))
        >>> bool(State.vphi.fget(proxy) == 2 * u.rad / u.s)
        True
        """
        return self.vrqp.vphi

    @property
    def vR(self) -> u.Quantity:
        """Return the Boyer--Lindquist radial coordinate velocity.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar or shape ``(n,)`` derivative with respect to
            coordinate time ``t``, in length per time.

        Examples
        --------
        >>> from types import SimpleNamespace
        >>> from astropy import units as u
        >>> proxy = SimpleNamespace(vRQP=SimpleNamespace(vR=1 * u.km / u.s))
        >>> bool(State.vR.fget(proxy) == 1 * u.km / u.s)
        True
        """
        return self.vRQP.vR

    @property
    def vTheta(self) -> u.Quantity:
        """Return the Boyer--Lindquist polar coordinate velocity.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar or shape ``(n,)`` derivative with respect to
            coordinate time ``t``, in angle per time.

        Examples
        --------
        >>> from types import SimpleNamespace
        >>> from astropy import units as u
        >>> proxy = SimpleNamespace(vRQP=SimpleNamespace(vTheta=1 * u.rad / u.s))
        >>> bool(State.vTheta.fget(proxy) == 1 * u.rad / u.s)
        True
        """
        return self.vRQP.vTheta

    @property
    def vPhi(self) -> u.Quantity:
        """Return the Boyer--Lindquist azimuthal coordinate velocity.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar or shape ``(n,)`` derivative with respect to
            coordinate time ``t``, in angle per time.

        Examples
        --------
        >>> from types import SimpleNamespace
        >>> from astropy import units as u
        >>> proxy = SimpleNamespace(vRQP=SimpleNamespace(vPhi=2 * u.rad / u.s))
        >>> bool(State.vPhi.fget(proxy) == 2 * u.rad / u.s)
        True
        """
        return self.vRQP.vPhi

    @property
    def ux(self) -> u.Quantity:
        """Return the Cartesian x component of four-velocity.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar or shape ``(n,)`` proper-time derivative of
            Cartesian length.

        Examples
        --------
        >>> from types import SimpleNamespace
        >>> from astropy import units as u
        >>> proxy = SimpleNamespace(uxyz=[1, 2, 3] * u.km / u.s)
        >>> bool(State.ux.fget(proxy) == 1 * u.km / u.s)
        True
        """
        return self.uxyz[..., 0]

    @property
    def uy(self) -> u.Quantity:
        """Return the Cartesian y component of four-velocity.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar or shape ``(n,)`` proper-time derivative of
            Cartesian length.

        Examples
        --------
        >>> from types import SimpleNamespace
        >>> from astropy import units as u
        >>> proxy = SimpleNamespace(uxyz=[1, 2, 3] * u.km / u.s)
        >>> bool(State.uy.fget(proxy) == 2 * u.km / u.s)
        True
        """
        return self.uxyz[..., 1]

    @property
    def uz(self) -> u.Quantity:
        """Return the Cartesian z component of four-velocity.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar or shape ``(n,)`` proper-time derivative of
            Cartesian length.

        Examples
        --------
        >>> from types import SimpleNamespace
        >>> from astropy import units as u
        >>> proxy = SimpleNamespace(uxyz=[1, 2, 3] * u.km / u.s)
        >>> bool(State.uz.fget(proxy) == 3 * u.km / u.s)
        True
        """
        return self.uxyz[..., 2]

    @property
    def ur(self) -> u.Quantity:
        """Return the spherical radial component of four-velocity.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar or shape ``(n,)`` length per proper time.

        Examples
        --------
        >>> from types import SimpleNamespace
        >>> from astropy import units as u
        >>> proxy = SimpleNamespace(urqp=SimpleNamespace(ur=1 * u.km / u.s))
        >>> bool(State.ur.fget(proxy) == 1 * u.km / u.s)
        True
        """
        return self.urqp.ur

    @property
    def utheta(self) -> u.Quantity:
        """Return the spherical polar component of four-velocity.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar or shape ``(n,)`` angle per proper time.

        Examples
        --------
        >>> from types import SimpleNamespace
        >>> from astropy import units as u
        >>> proxy = SimpleNamespace(urqp=SimpleNamespace(utheta=1 * u.rad / u.s))
        >>> bool(State.utheta.fget(proxy) == 1 * u.rad / u.s)
        True
        """
        return self.urqp.utheta

    @property
    def uphi(self) -> u.Quantity:
        """Return the spherical azimuthal component of four-velocity.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar or shape ``(n,)`` angle per proper time.

        Examples
        --------
        >>> from types import SimpleNamespace
        >>> from astropy import units as u
        >>> proxy = SimpleNamespace(urqp=SimpleNamespace(uphi=2 * u.rad / u.s))
        >>> bool(State.uphi.fget(proxy) == 2 * u.rad / u.s)
        True
        """
        return self.urqp.uphi

    @property
    def uR(self) -> u.Quantity:
        """Return the Boyer--Lindquist radial four-velocity component.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar or shape ``(n,)`` length per proper time.

        Examples
        --------
        >>> from types import SimpleNamespace
        >>> from astropy import units as u
        >>> proxy = SimpleNamespace(uRQP=SimpleNamespace(uR=1 * u.km / u.s))
        >>> bool(State.uR.fget(proxy) == 1 * u.km / u.s)
        True
        """
        return self.uRQP.uR

    @property
    def uTheta(self) -> u.Quantity:
        """Return the Boyer--Lindquist polar four-velocity component.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar or shape ``(n,)`` angle per proper time.

        Examples
        --------
        >>> from types import SimpleNamespace
        >>> from astropy import units as u
        >>> proxy = SimpleNamespace(uRQP=SimpleNamespace(uTheta=1 * u.rad / u.s))
        >>> bool(State.uTheta.fget(proxy) == 1 * u.rad / u.s)
        True
        """
        return self.uRQP.uTheta

    @property
    def uPhi(self) -> u.Quantity:
        """Return the Boyer--Lindquist azimuthal four-velocity component.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar or shape ``(n,)`` angle per proper time.

        Examples
        --------
        >>> from types import SimpleNamespace
        >>> from astropy import units as u
        >>> proxy = SimpleNamespace(uRQP=SimpleNamespace(uPhi=2 * u.rad / u.s))
        >>> bool(State.uPhi.fget(proxy) == 2 * u.rad / u.s)
        True
        """
        return self.uRQP.uPhi

    @property
    def state_vector(self) -> CartesianStateVector:
        """Return Cartesian position and coordinate velocity.

        Returns
        -------
        CartesianStateVector
            Immutable Cartesian phase-space state.

        Examples
        --------
        >>> from types import SimpleNamespace
        >>> marker = object()
        >>> State.state_vector.fget(SimpleNamespace(_state_vector=marker)) is marker
        True
        """
        return self._state_vector

    def orbital_elements(self) -> OrbitalElements:
        """Return the stored classical osculating orbital elements.

        Returns
        -------
        OrbitalElements
            Immutable scalar or vectorized elements aligned with this state.
        """
        return self._orbital_elements


class InitialState(State):
    """Represent immutable saved initial conditions.

    ``InitialState`` accepts the same keyword-only parameters and performs the
    same validation as :class:`State`.  Its distinct type is intended to
    identify the conditions preserved by an Orbit instance.

    Parameters
    ----------
    **fields
        The keyword-only fields documented in :class:`State`, with the same
        types, shapes, units and validation.

    Attributes
    ----------
    tau, txyz, trqp, tRQP, vxyz, vrqp, vRQP, ut, uxyz, urqp, uRQP
        Immutable validated state views inherited from :class:`State`.

    See Also
    --------
    State : Parameters, exceptions and conversion views.
    relatipy.geodesic.Orbit.initial : Saved initial conditions of an orbit.
    """
