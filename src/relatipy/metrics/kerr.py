"""Public Kerr parameters and mass-scaled characteristic radii."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from astropy import units as u
from astropy.constants import G, c

from .._validation import dimensionless_scalar, readonly_quantity
from ..coordinates import KerrOrbitalElements
from ..geodesic.orbit import Orbit

# Private codes matching ``rp_kerr_surface`` in the native core.
_SURFACE_CODES = {"outer_horizon": 0, "ergosurface": 1}


@dataclass(frozen=True, slots=True)
class Horizons:
    """Boyer--Lindquist radii of the Kerr horizons.

    Parameters
    ----------
    event, cauchy : astropy.units.Quantity
        Scalar outer and inner horizon coordinate radii, with units compatible
        with length.

    Attributes
    ----------
    event, cauchy : astropy.units.Quantity
        Read-only scalar radii in the input-compatible length units.

    Raises
    ------
    TypeError
        If either radius is not an Astropy quantity.
    astropy.units.UnitConversionError
        If either radius is not compatible with length.
    ValueError
        If either radius is not scalar.
    """

    event: u.Quantity
    cauchy: u.Quantity

    def __post_init__(self) -> None:
        """Validate and freeze the two scalar horizon radii."""
        object.__setattr__(
            self, "event", readonly_quantity(self.event, u.m, "event", ndim=(0,))
        )
        object.__setattr__(
            self, "cauchy", readonly_quantity(self.cauchy, u.m, "cauchy", ndim=(0,))
        )


class Kerr:
    """A Kerr black hole with immutable physical mass and dimensionless spin.

    Parameters
    ----------
    mass : astropy.units.Quantity
        Finite, positive scalar mass.
    spin : float or astropy.units.Quantity
        Finite dimensionless scalar in the closed interval ``[0, 1]``.

    Raises
    ------
    TypeError, astropy.units.UnitConversionError
        If ``mass`` or ``spin`` has an invalid physical type or unit.
    ValueError
        If mass or spin is non-scalar, non-finite, or outside its accepted
        range, or if mass cannot produce finite geometric scales.
    RuntimeError
        If the private native calculation returns invalid characteristic radii.

    Notes
    -----
    Instances are immutable values: two metrics are equal, and hash alike,
    when their masses in kilograms and their spins are equal.

    ``spin`` is the dimensionless ratio ``a / M`` of the Kerr parameter to
    the mass in geometric units; negative values are not accepted, so
    retrograde motion is selected by the orbit's initial conditions. The
    native core evaluates the metric in Boyer--Lindquist coordinates
    ``(t, r, theta, phi)`` with signature ``(-, +, +, +)`` and geometric
    units ``G = c = 1``, with lengths scaled by
    :attr:`r_g` ``= G M / c^2`` and times by ``G M / c^3``. Every public
    radius and time is converted back to an Astropy quantity; the returned
    radii are Boyer--Lindquist coordinate radii, not Cartesian distances.

    Examples
    --------
    >>> from astropy import units as u
    >>> Kerr(mass=1 * u.M_sun, spin=0.5)
    Kerr(mass=<Quantity 1. solMass>, spin=0.5)
    """

    __slots__ = ("_mass", "_spin", "_length_scale", "_time_scale", "_radii")

    def __init__(self, *, mass: u.Quantity, spin: float | u.Quantity) -> None:
        """Initialize one validated immutable Kerr value.

        Parameters
        ----------
        mass : astropy.units.Quantity
            Finite, positive scalar mass.
        spin : float or astropy.units.Quantity
            Finite dimensionless scalar in ``[0, 1]``.

        Raises
        ------
        TypeError, astropy.units.UnitConversionError
            If an input has an invalid physical type or unit.
        ValueError
            If an input is non-scalar, non-finite, or outside its accepted
            range, or if mass cannot produce finite geometric scales.
        RuntimeError
            If private native characteristic radii are invalid.
        """
        from .. import _core

        value = readonly_quantity(mass, u.kg, "mass", ndim=(0,))
        mass_kg = float(value.to_value(u.kg))
        if not np.isfinite(mass_kg) or mass_kg <= 0:
            raise ValueError("mass must be finite and positive")
        spin_value = dimensionless_scalar(spin, "spin")
        if not np.isfinite(spin_value) or not 0 <= spin_value <= 1:
            raise ValueError("spin must be finite and in [0, 1]")

        length_scale = (G * value / c**2).to(u.m)
        time_scale = (G * value / c**3).to(u.s)
        if (not np.isfinite(length_scale.to_value(u.m))
                or not np.isfinite(time_scale.to_value(u.s))
                or length_scale <= 0 * u.m or time_scale <= 0 * u.s):
            raise ValueError("mass produces nonrepresentable geometric scales")
        radii = np.asarray(_core.kerr_properties(spin_value), dtype=np.float64)
        if radii.shape != (6,) or not np.all(np.isfinite(radii)):
            raise RuntimeError("native Kerr characteristic radii are invalid")
        for name, attribute in (
            ("_mass", value),
            ("_spin", spin_value),
            ("_length_scale", readonly_quantity(length_scale, u.m, "r_g", ndim=(0,))),
            ("_time_scale", readonly_quantity(time_scale, u.s, "t_g", ndim=(0,))),
            ("_radii", tuple(float(radius) for radius in radii)),
        ):
            object.__setattr__(self, name, attribute)

    def __setattr__(self, name: str, value: object) -> None:
        """Reject attribute assignment after construction.

        Raises
        ------
        AttributeError
            Always; Kerr metric values are immutable.
        """
        raise AttributeError("Kerr metrics are immutable")

    def __delattr__(self, name: str) -> None:
        """Reject attribute deletion after construction.

        Raises
        ------
        AttributeError
            Always; Kerr metric values are immutable.
        """
        raise AttributeError("Kerr metrics are immutable")

    def __reduce__(self) -> tuple:
        """Return constructor arguments used by copying and pickling."""
        return _rebuild_kerr, (self._mass, self._spin)

    def __repr__(self) -> str:
        """Return an evaluable-style representation of mass and spin."""
        return f"Kerr(mass={self._mass!r}, spin={self._spin!r})"

    def __eq__(self, other: object) -> bool:
        """Compare metrics by mass in kilograms and dimensionless spin."""
        if not isinstance(other, Kerr):
            return NotImplemented
        return self._key() == other._key()

    def __hash__(self) -> int:
        """Hash the same normalized mass-and-spin key used for equality."""
        return hash(self._key())

    def _key(self) -> tuple[float, float]:
        """Return the normalized immutable value key for equality and hashing."""
        return float(self._mass.to_value(u.kg)), self._spin

    def _radius(self, index: int, name: str) -> u.Quantity:
        """Scale one native dimensionless radius and return it read-only."""
        return readonly_quantity(
            self._radii[index] * self._length_scale, u.m, name, ndim=(0,)
        )

    @property
    def mass(self) -> u.Quantity:
        """Return the immutable physical mass.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar mass in the unit supplied at construction.
        """
        return self._mass

    @property
    def spin(self) -> float:
        """Return the dimensionless spin parameter.

        Returns
        -------
        float
            Finite spin in the closed interval ``[0, 1]``.
        """
        return self._spin

    @property
    def r_g(self) -> u.Quantity:
        """Return the gravitational length scale ``G M / c²``.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar length in metres.
        """
        return self._length_scale

    @property
    def r_s(self) -> u.Quantity:
        """Return the Schwarzschild length scale ``2 G M / c²``.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar length in metres.
        """
        return readonly_quantity(2 * self._length_scale, u.m, "r_s", ndim=(0,))

    @property
    def horizons(self) -> Horizons:
        """Return the outer and inner horizon coordinate radii.

        Returns
        -------
        Horizons
            Immutable event and Cauchy radii in metres.
        """
        return Horizons(self._radius(0, "event"), self._radius(1, "cauchy"))

    @property
    def r_isco_prograde(self) -> u.Quantity:
        """Return the prograde equatorial ISCO radius.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar Boyer--Lindquist coordinate radius in metres.
        """
        return self._radius(2, "r_isco_prograde")

    @property
    def r_isco_retrograde(self) -> u.Quantity:
        """Return the retrograde equatorial ISCO radius.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar Boyer--Lindquist coordinate radius in metres.
        """
        return self._radius(3, "r_isco_retrograde")

    @property
    def r_photon_prograde(self) -> u.Quantity:
        """Return the prograde equatorial photon orbit radius.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar Boyer--Lindquist coordinate radius in metres.
        """
        return self._radius(4, "r_photon_prograde")

    @property
    def r_photon_retrograde(self) -> u.Quantity:
        """Return the retrograde equatorial photon orbit radius.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar Boyer--Lindquist coordinate radius in metres.
        """
        return self._radius(5, "r_photon_retrograde")

    def r_ergosurface(self, *, theta: u.Quantity) -> u.Quantity:
        """Evaluate the outer ergosurface radius at polar angle ``theta``.

        Parameters
        ----------
        theta : astropy.units.Quantity
            Scalar or one-dimensional polar angle in ``[0, pi]``.

        Returns
        -------
        astropy.units.Quantity
            Read-only outer stationary-limit radius in metres, with the same
            scalar or one-dimensional shape as ``theta``.

        Raises
        ------
        TypeError, astropy.units.UnitConversionError
            If ``theta`` is not an angle quantity.
        ValueError
            If angles are non-finite or outside the polar domain.
        RuntimeError
            If the private native calculation returns an invalid result.
        """
        from .. import _core

        angle = readonly_quantity(theta, u.rad, "theta")
        radians = np.asarray(angle.to_value(u.rad), dtype=np.float64)
        if (
            not np.all(np.isfinite(radians))
            or np.any(radians < 0)
            or np.any(radians > np.pi)
        ):
            raise ValueError("theta must be finite and in [0, pi]")
        flattened = np.ascontiguousarray(radians.reshape(-1))
        radii = np.asarray(
            _core.kerr_ergosurface(self._spin, flattened), dtype=np.float64
        )
        if radii.shape != flattened.shape or not np.all(np.isfinite(radii)):
            raise RuntimeError("native Kerr ergosurface result is invalid")
        result = radii.reshape(radians.shape) * self._length_scale
        return readonly_quantity(result, u.m, "r_ergosurface")

    def _surface_profile(
        self, surface: str, theta: u.Quantity
    ) -> tuple[u.Quantity, u.Quantity]:
        """Evaluate the Cartesian meridional profile of a reference surface.

        Private helper for plotting. The native core maps the Boyer--Lindquist
        radius of the surface to ``rho = hypot(r, a) sin(theta)`` and
        ``z = r cos(theta)``.

        Parameters
        ----------
        surface : {"outer_horizon", "ergosurface"}
            Outer event horizon or outer ergosurface (stationary limit).
        theta : astropy.units.Quantity
            Scalar or one-dimensional polar angle in ``[0, pi]``.

        Returns
        -------
        rho, z : astropy.units.Quantity
            Read-only cylindrical radius and height in metres, each with the
            shape of ``theta``.

        Raises
        ------
        TypeError, astropy.units.UnitConversionError
            If ``theta`` is not an angle quantity.
        ValueError
            If ``surface`` is unknown or angles are non-finite or outside the
            polar domain.
        """
        from .. import _core

        if not isinstance(surface, str) or surface not in _SURFACE_CODES:
            raise ValueError(
                "surface must be 'outer_horizon' or 'ergosurface'"
            )
        angle = readonly_quantity(theta, u.rad, "theta")
        radians = np.asarray(angle.to_value(u.rad), dtype=np.float64)
        if (
            not np.all(np.isfinite(radians))
            or np.any(radians < 0)
            or np.any(radians > np.pi)
        ):
            raise ValueError("theta must be finite and in [0, pi]")
        flattened = np.ascontiguousarray(radians.reshape(-1))
        rho, z = _core.kerr_surface_profile(
            self._spin, _SURFACE_CODES[surface], flattened
        )
        rho = np.asarray(rho, dtype=np.float64)
        z = np.asarray(z, dtype=np.float64)
        if (
            rho.shape != flattened.shape
            or z.shape != flattened.shape
            or not np.all(np.isfinite(rho))
            or not np.all(np.isfinite(z))
        ):
            raise RuntimeError("native Kerr surface profile result is invalid")
        return (
            readonly_quantity(
                rho.reshape(radians.shape) * self._length_scale, u.m, "rho"
            ),
            readonly_quantity(
                z.reshape(radians.shape) * self._length_scale, u.m, "z"
            ),
        )

    def orbit(
        self,
        *,
        elements: KerrOrbitalElements | None = None,
        p: u.Quantity | None = None,
        q_r0: u.Quantity | None = None,
        q_theta0: u.Quantity | None = None,
        q_phi0: u.Quantity | None = None,
        a: u.Quantity | None = None,
        e: float | u.Quantity | None = None,
        inc: u.Quantity | None = None,
        Omega: u.Quantity | None = None,
        omega: u.Quantity | None = None,
        f: u.Quantity | None = None,
        x: float | u.Quantity | None = None,
        y: u.Quantity | None = None,
        z: u.Quantity | None = None,
        vx: u.Quantity | None = None,
        vy: u.Quantity | None = None,
        vz: u.Quantity | None = None,
        r: u.Quantity | None = None,
        theta: u.Quantity | None = None,
        phi: u.Quantity | None = None,
        vr: u.Quantity | None = None,
        vtheta: u.Quantity | None = None,
        vphi: u.Quantity | None = None,
        R: u.Quantity | None = None,
        Theta: u.Quantity | None = None,
        Phi: u.Quantity | None = None,
        vR: u.Quantity | None = None,
        vTheta: u.Quantity | None = None,
        vPhi: u.Quantity | None = None,
        tau: u.Quantity = 0 * u.s,
        t: u.Quantity = 0 * u.s,
    ) -> Orbit:
        """Construct one scalar timelike orbit from one input family.

        Parameters
        ----------
        elements : relatipy.coordinates.KerrOrbitalElements, optional
            Immutable bound-geodesic parameters. Cannot be combined with
            individual coordinate or orbital-element arguments.
        p : astropy.units.Quantity, optional
            Bound Kerr semilatus rectum as a physical length. Selects the
            bound family, which requires dimensionless ``e`` and ``x``.
            Stability is checked in C; marginal and unstable orbits fail.
        q_r0, q_theta0, q_phi0 : astropy.units.Quantity, optional
            Bound Kerr initial Mino phases, default zero radians. Radial
            zero is periapsis; polar zero is the northern turning point.
            These phases are not Keplerian true anomalies.
        a : astropy.units.Quantity, optional
            Scalar classical semi-major axis as a physical length. It selects
            the classical-element family.
        e : float or astropy.units.Quantity, optional
            Scalar dimensionless classical eccentricity. It defaults to zero
            with the classical-element family.
        inc, Omega, omega, f : astropy.units.Quantity, optional
            Scalar inclination, longitude of ascending node, argument of
            periapsis, and true anomaly. Each is an angle and defaults to zero
            with the classical-element family.
        x, y, z, vx, vy, vz : optional
            Scalar Cartesian position lengths and coordinate-time velocities.
            At least one position component selects the Cartesian family;
            omitted Cartesian position and velocity components default to
            zero. With ``p``, ``x`` instead means the dimensionless Kerr
            inclination parameter: positive prograde, negative retrograde,
            absolute value one equatorial. See ``KerrOrbitalElements``.
        r, theta, phi, vr, vtheta, vphi : astropy.units.Quantity, optional
            Scalar spherical radial length, polar and azimuthal angles, and
            coordinate-time velocities. ``r``, ``theta``, and ``phi`` are
            required when this family is selected; omitted velocities default
            to zero with radial units of length per time and angular units of
            angle per time.
        R, Theta, Phi, vR, vTheta, vPhi : astropy.units.Quantity, optional
            Scalar Boyer--Lindquist radial length, polar and azimuthal angles,
            and coordinate-time velocities. ``R``, ``Theta``, and ``Phi`` are
            required when this family is selected; omitted velocities default
            to zero with radial units of length per time and angular units of
            angle per time.
        tau, t : astropy.units.Quantity, optional
            Scalar initial proper and coordinate times, respectively, each
            defaulting to zero seconds.

        Returns
        -------
        Orbit
            Independently evolving scalar orbit.

        Raises
        ------
        TypeError, ValueError, astropy.units.UnitConversionError
            If no single family is selected, a required value is missing, or a
            scalar shape, value, unit, or initial state is invalid. The initial
            Boyer--Lindquist radius must exceed the outer horizon.
        """
        arguments = locals().copy()
        arguments.pop("self")
        return Orbit._from_kwargs(self, **arguments)


def _rebuild_kerr(mass: u.Quantity, spin: float) -> Kerr:
    """Recreate a copied or pickled metric through its validating constructor.

    Parameters
    ----------
    mass : astropy.units.Quantity
        Stored scalar physical mass.
    spin : float
        Stored dimensionless spin.

    Returns
    -------
    Kerr
        A newly validated immutable metric.
    """
    return Kerr(mass=mass, spin=spin)
