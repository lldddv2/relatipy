"""Built-in Plummer sphere centred on the black hole."""

from __future__ import annotations

import numpy as np
from astropy import units as u

from .core import ExternalForce


class Plummer(ExternalForce):
    """Plummer sphere centred on the black hole.

    Parameters
    ----------
    mass : astropy.units.Quantity
        Finite scalar sphere mass ``m``, zero or positive. Zero gives an
        identically zero acceleration.
    radius : astropy.units.Quantity
        Finite, positive scalar scale radius ``b``.

    Raises
    ------
    TypeError
        If a parameter is not a scalar quantity.
    astropy.units.UnitConversionError
        If ``mass`` is not a mass or ``radius`` not a length.
    ValueError
        If a parameter is not finite, ``mass`` is negative or ``radius`` is
        not positive.

    Notes
    -----
    The potential and acceleration in the black-hole Cartesian frame are

    .. math::

       \\Phi = -\\frac{G m}{\\sqrt{|x|^2 + b^2}}, \\qquad
       a = -\\frac{G m\\, x}{(|x|^2 + b^2)^{3/2}} ,

    evaluated by the native callback with ``m`` in black-hole masses and
    ``b`` in ``G M / c^2``. The potential is static and axisymmetric about
    the spin axis, so ``E' = -u_t exp(Phi)`` and ``L' = u_phi exp(Phi)`` are
    conserved. The parameters ``mass`` and ``radius`` form the force block
    of the perturbed MCMC models, in the units given here.

    Examples
    --------
    >>> from astropy import units as u
    >>> from relatipy.forces import Plummer
    >>> Plummer(mass=1e4 * u.M_sun, radius=0.01 * u.pc)
    Plummer(mass=<Quantity 10000. solMass>, radius=<Quantity 0.01 pc>)
    """

    __slots__ = ()

    def __init__(self, *, mass: u.Quantity, radius: u.Quantity) -> None:
        """Validate the sphere mass and scale radius."""
        for name, value, unit in (("mass", mass, u.kg), ("radius", radius, u.m)):
            if not isinstance(value, u.Quantity):
                raise TypeError(f"{name} must be an astropy quantity")
            value.to(unit)
        super().__init__({"mass": mass, "radius": radius})

    def _check_vector(self, values: np.ndarray) -> np.ndarray:
        values = super()._check_vector(values)
        if values[0] < 0:
            raise ValueError("Plummer mass must not be negative")
        if not values[1] > 0:
            raise ValueError("Plummer radius must be positive")
        return values

    def __repr__(self) -> str:
        """Return an evaluable-style representation."""
        return (
            f"Plummer(mass={self._parameters['mass']!r}, "
            f"radius={self._parameters['radius']!r})"
        )

    def __reduce__(self) -> tuple:
        """Return constructor arguments used by copying and pickling."""
        return _rebuild_plummer, (self._parameters["mass"], self._parameters["radius"])

    def _key(self) -> tuple:
        return self._parameter_key()

    def _replace(self, parameters: dict[str, object]) -> Plummer:
        return Plummer(mass=parameters["mass"], radius=parameters["radius"])

    def _native_callback(self) -> tuple[int, object]:
        from .. import _force_core

        return int(_force_core.plummer_address()), _force_core


def _rebuild_plummer(mass: u.Quantity, radius: u.Quantity) -> Plummer:
    """Recreate a pickled Plummer sphere through its validating constructor."""
    return Plummer(mass=mass, radius=radius)
