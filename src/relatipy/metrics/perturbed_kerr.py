"""A Kerr black hole whose timelike geodesics feel one external force."""

from __future__ import annotations

from typing import TYPE_CHECKING

from astropy import units as u

from ..coordinates import KerrOrbitalElements
from ..forces import ExternalForce
from ..geodesic._forced_orbit import ForcedOrbit
from .kerr import Horizons, Kerr

if TYPE_CHECKING:
    from ..forces._bound import BoundForce
    from ..geodesic.orbit import Orbit


class PerturbedKerr:
    """A Kerr black hole plus one external Newtonian force.

    Parameters
    ----------
    mass : astropy.units.Quantity
        Finite, positive scalar black-hole mass.
    spin : float or astropy.units.Quantity
        Finite dimensionless spin in ``[0, 1]``.
    force : relatipy.forces.ExternalForce
        Force felt by every orbit, such as a :class:`~relatipy.forces.ThinDisk`
        or a :class:`~relatipy.forces.Potential`.

    Raises
    ------
    TypeError
        If ``force`` is not an external force, or as :class:`Kerr`.
    ValueError, astropy.units.UnitConversionError
        As :class:`Kerr`.

    Notes
    -----
    The background geometry is the Kerr metric exposed by :attr:`kerr`; the
    force is a test-particle perturbation that does not change the metric.
    Orbits solve

    .. math::

       \\frac{du^\\mu}{d\\tau} + \\Gamma^\\mu{}_{\\alpha\\beta} u^\\alpha u^\\beta
       = \\left(g^{\\mu\\nu} + u^\\mu u^\\nu\\right) F_\\nu ,

    where ``F_nu`` is the covector of the force's Newtonian acceleration in
    the black-hole Cartesian frame (``z`` along the spin), so
    ``u_mu u^mu = -1`` is preserved. Directions of a force, such as a disk
    normal, are given in that frame.

    Instances are immutable values: two are equal, and hash alike, when
    their :attr:`kerr` metrics and forces are equal. The force is bound to
    the mass, and compiled if needed, on the first call of :meth:`orbit`.
    Kerr quantities that a force changes, such as ISCO radii and constants
    of motion, are available only through :attr:`kerr` and refer to the
    unperturbed metric.

    Examples
    --------
    >>> from astropy import units as u
    >>> from relatipy.forces import ThinDisk
    >>> bh = PerturbedKerr(mass=4e6 * u.M_sun, spin=0.5,
    ...                    force=ThinDisk(mass=1e3 * u.M_sun, r_out=1362.5 * u.au))
    >>> bh.kerr
    Kerr(mass=<Quantity 4000000. solMass>, spin=0.5)
    """

    __slots__ = ("_kerr", "_force", "_bound")

    def __init__(
        self,
        *,
        mass: u.Quantity,
        spin: float | u.Quantity,
        force: ExternalForce,
    ) -> None:
        """Validate the Kerr parameters and the force."""
        if not isinstance(force, ExternalForce):
            raise TypeError("force must be a relatipy.forces.ExternalForce")
        object.__setattr__(self, "_kerr", Kerr(mass=mass, spin=spin))
        object.__setattr__(self, "_force", force)
        object.__setattr__(self, "_bound", None)

    def __setattr__(self, name: str, value: object) -> None:
        """Reject attribute assignment after construction.

        Raises
        ------
        AttributeError
            Always; perturbed metrics are immutable.
        """
        raise AttributeError("PerturbedKerr metrics are immutable")

    def __delattr__(self, name: str) -> None:
        """Reject attribute deletion after construction.

        Raises
        ------
        AttributeError
            Always; perturbed metrics are immutable.
        """
        raise AttributeError("PerturbedKerr metrics are immutable")

    def __reduce__(self) -> tuple:
        """Return constructor arguments used by copying and pickling."""
        return _rebuild_perturbed_kerr, (self._kerr.mass, self._kerr.spin, self._force)

    def __repr__(self) -> str:
        """Return an evaluable-style representation."""
        return (
            f"PerturbedKerr(mass={self._kerr.mass!r}, spin={self._kerr.spin!r}, "
            f"force={self._force!r})"
        )

    def __eq__(self, other: object) -> bool:
        """Compare the background metrics and the forces."""
        if not isinstance(other, PerturbedKerr):
            return NotImplemented
        return self._kerr == other._kerr and self._force == other._force

    def __hash__(self) -> int:
        """Hash the metric and force used for equality."""
        return hash((self._kerr, self._force))

    @property
    def kerr(self) -> Kerr:
        """Return the unperturbed background metric.

        Returns
        -------
        Kerr
            Immutable Kerr metric with the same mass and spin.
        """
        return self._kerr

    @property
    def force(self) -> ExternalForce:
        """Return the external force.

        Returns
        -------
        relatipy.forces.ExternalForce
            The immutable force given at construction.
        """
        return self._force

    @property
    def mass(self) -> u.Quantity:
        """Return the black-hole mass, as :attr:`Kerr.mass`.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar mass in the unit supplied at construction.
        """
        return self._kerr.mass

    @property
    def spin(self) -> float:
        """Return the dimensionless spin, as :attr:`Kerr.spin`.

        Returns
        -------
        float
            Finite spin in ``[0, 1]``.
        """
        return self._kerr.spin

    @property
    def r_g(self) -> u.Quantity:
        """Return the gravitational length scale ``G M / c²``.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar length in metres.
        """
        return self._kerr.r_g

    @property
    def horizons(self) -> Horizons:
        """Return the horizons of the background metric.

        Returns
        -------
        Horizons
            Immutable event and Cauchy radii in metres. Orbits stop at the
            outer horizon as for :class:`Kerr`.
        """
        return self._kerr.horizons

    def _binding(self) -> BoundForce:
        """Return the force bound to this mass, binding it on first use."""
        bound = self._bound
        if bound is None:
            bound = self._force._bind(self._kerr.mass)
            object.__setattr__(self, "_bound", bound)
        return bound

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
        """Construct one scalar timelike orbit that feels the force.

        Parameters
        ----------
        elements, p, q_r0, q_theta0, q_phi0, a, e, inc, Omega, omega, f, \
x, y, z, vx, vy, vz, r, theta, phi, vr, vtheta, vphi, R, Theta, Phi, \
vR, vTheta, vPhi, tau, t
            One initial-condition family, exactly as
            :meth:`relatipy.metrics.Kerr.orbit`. Initial states are those of
            the background metric; bound Kerr parameters describe the
            unperturbed geodesic through the initial point.

        Returns
        -------
        Orbit
            Independently evolving scalar orbit. :meth:`~relatipy.geodesic.Orbit.integrate`
            and :meth:`~relatipy.geodesic.Orbit.solve` behave as for
            :class:`Kerr` and return the same
            :class:`~relatipy.geodesic.Solution` type, integrated with the
            force; ``method`` is ``"radau"``, ``"dop853"`` or ``"dp45"``.

        Raises
        ------
        TypeError, ValueError, astropy.units.UnitConversionError
            As :meth:`relatipy.metrics.Kerr.orbit`.
        ImportError, RuntimeError
            If a run-time force needs SymPy or a C compiler that is missing.

        Notes
        -----
        ``projection_radau`` raises :class:`ValueError`: its projector
        restores Carter's constant, which a force does not conserve. Radau
        uses a finite-difference Jacobian. A failing force callback is a
        numerical failure (``Solution.status == -1`` or
        :class:`~relatipy.geodesic.IntegrationError`). ``get_E``, ``get_Lz``
        and ``get_Q`` report the Kerr quantities of the current state, which
        a force does not conserve.
        """
        arguments = locals().copy()
        arguments.pop("self")
        bound = self._binding()
        orbit = ForcedOrbit._from_kwargs(self._kerr, **arguments)
        orbit._bound = bound
        return orbit


def _rebuild_perturbed_kerr(
    mass: u.Quantity, spin: float, force: ExternalForce
) -> PerturbedKerr:
    """Recreate a copied or pickled perturbed metric through its constructor."""
    return PerturbedKerr(mass=mass, spin=spin, force=force)
