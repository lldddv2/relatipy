"""Private native force handles bound to a black-hole mass.

A bound force owns the model memory passed to one native
``rp_external_acceleration_fn`` and keeps alive whatever the callback
depends on (a compiled library, a table). Its parameters are written in
place, so a sampler can move between proposals without rebuilding tables or
recompiling. Raw addresses are never pickled: a pickled handle stores its
force and last update and binds again when it is loaded.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np
from astropy import constants as const
from astropy import units as u

if TYPE_CHECKING:
    from .core import ExternalForce

_BASES = (u.m, u.s, u.kg, u.rad)
_G = float(const.G.si.value)
_C = float(const.c.si.value)
# Length and time units per kilogram of black-hole mass: r_g = G M / c^2 and
# t_g = G M / c^3, so lengths divide by r_g and times by t_g.
_PER_LENGTH = _C**2 / _G
_PER_TIME = _C**3 / _G


def geometric_conversion(unit: u.UnitBase, name: str) -> tuple[float, float]:
    """Return ``(factor, power)`` mapping a value in ``unit`` to geometric units.

    A value ``v`` in ``unit`` is ``v * factor * M**(-power)`` in the units of
    a black hole of mass ``M`` kilograms, with ``G = c = M = 1``. ``unit``
    must decompose into metres, seconds, kilograms and radians; radians are
    dimensionless.

    Raises
    ------
    astropy.units.UnitConversionError
        If ``unit`` involves another base unit, such as a temperature.
    """
    try:
        decomposed = u.Unit(unit).decompose(bases=_BASES)
    except u.UnitConversionError:
        raise u.UnitConversionError(
            f"{name} must have a unit built from length, time, mass and angle"
        ) from None
    exponents = {base: float(power) for base, power in
                 zip(decomposed.bases, decomposed.powers)}
    length = exponents.get(u.m, 0.0)
    time = exponents.get(u.s, 0.0)
    mass = exponents.get(u.kg, 0.0)
    factor = float(decomposed.scale) * _PER_LENGTH**length * _PER_TIME**time
    if not math.isfinite(factor) or factor == 0.0:
        raise ValueError(f"{name} unit cannot be represented in geometric units")
    return factor, length + time + mass


class BoundForce:
    """One native force handle for one black-hole mass.

    Parameters
    ----------
    force : ExternalForce
        Force whose parameters this handle stores.
    callback : int
        Address of the native callback.
    size : int
        Number of doubles in the model data.
    keepalive : object
        Objects the callback depends on, kept referenced by this handle.

    Notes
    -----
    :attr:`callback` and :attr:`model` are valid while this object is
    referenced. Calls that pass them to native code must hold a reference
    and must not run :meth:`update` concurrently.
    """

    __slots__ = (
        "_force", "_callback", "_buffer", "_keepalive", "_mass", "_values",
        "_rotation",
    )

    def __init__(
        self, force: ExternalForce, callback: int, size: int, keepalive: object
    ) -> None:
        self._force = force
        self._callback = int(callback)
        self._buffer = np.zeros(max(int(size), 1), dtype=np.float64)
        self._keepalive = keepalive
        self._mass = math.nan
        self._values: np.ndarray | None = None
        self._rotation: np.ndarray | None = None

    def __reduce__(self) -> tuple:
        """Pickle the force and the last update, not the native addresses."""
        return _restore_bound, (self._force, self._mass, self._values, self._rotation)

    @property
    def force(self) -> ExternalForce:
        """Force whose callback and parameter layout this handle uses."""
        return self._force

    @property
    def callback(self) -> int:
        """Address of the native ``rp_external_acceleration_fn``."""
        return self._callback

    @property
    def model(self) -> int:
        """Address of the model data passed to the callback."""
        return int(self._buffer.ctypes.data)

    @property
    def mass(self) -> float:
        """Black-hole mass in kilograms of the last update."""
        return self._mass

    @property
    def model_data(self) -> np.ndarray:
        """Read-only copy of the geometric model data."""
        data = self._buffer.copy()
        data.flags.writeable = False
        return data

    def update(
        self,
        mass: float,
        values: np.ndarray,
        rotation: np.ndarray | None = None,
    ) -> None:
        """Write new parameters in place.

        Parameters
        ----------
        mass : float
            Black-hole mass in kilograms.
        values : numpy.ndarray
            Parameters in the force's :attr:`ExternalForce.parameter_units`,
            already validated by ``ExternalForce._check_vector``.
        rotation : numpy.ndarray or None
            Body-to-observer rotation whose columns are the black-hole axes
            in observer axes ``(Dec, RA, away)``. ``None`` means directions
            are given in the black-hole frame.
        """
        self._force._write(self._buffer, float(mass), values, rotation)
        self._mass = float(mass)
        self._values = np.array(values, dtype=np.float64)
        self._rotation = None if rotation is None else np.array(rotation, dtype=np.float64)


def _restore_bound(
    force: ExternalForce,
    mass: float,
    values: np.ndarray | None,
    rotation: np.ndarray | None,
) -> BoundForce:
    """Bind a pickled handle again and replay its last update."""
    bound = force._handle()
    if values is not None:
        bound.update(mass, values, rotation)
    return bound
