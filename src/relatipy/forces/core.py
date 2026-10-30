"""Base class of external forces acting on a perturbed Kerr geodesic."""

from __future__ import annotations

from collections.abc import Mapping
from types import MappingProxyType

import numpy as np
from astropy import units as u

from .._validation import REAL_TYPES
from ._bound import BoundForce, geometric_conversion


def _parameter_quantity(value: object, name: str) -> u.Quantity:
    """Validate one finite scalar parameter, a number or an Astropy quantity."""
    if isinstance(value, u.Quantity):
        if value.ndim != 0:
            raise ValueError(f"{name} must be a scalar")
        quantity = u.Quantity(value, copy=True)
    elif isinstance(value, (bool, np.bool_)) or not isinstance(value, REAL_TYPES):
        raise TypeError(f"{name} must be a real number or a scalar quantity")
    else:
        quantity = u.Quantity(float(value), u.one)
    if not np.isfinite(quantity.value):
        raise ValueError(f"{name} must be finite")
    quantity.flags.writeable = False
    return quantity


class ExternalForce:
    """Immutable external force with named run-time parameters.

    This base class is not instantiated directly; use the built-in
    :class:`~relatipy.forces.ThinDisk`.

    Notes
    -----
    Every force evaluates a Newtonian acceleration in the black-hole
    Cartesian frame (``z`` along the spin), in geometric units of the black
    hole with ``G = c = M = 1``. Each parameter keeps the unit it was given
    in; when a force is bound to a black hole the parameter is converted to
    geometric units through its dimensions in length ``L``, time ``T`` and
    mass ``M``: it is divided by ``r_g**L * t_g**T * M**M`` with
    ``r_g = G M / c^2`` and ``t_g = G M / c^3``. A plain number is
    dimensionless and is used as given.

    Forces are immutable values: two forces are equal, and hash alike, when
    they have the same type, definition and parameters in the same units.
    Native resources (a table) are created on first use and shared by every
    binding of the same object.
    """

    __slots__ = ("_parameters", "_conversions", "_native")

    def __init__(self, parameters: Mapping[str, object]) -> None:
        """Validate and store named parameters (called by subclasses)."""
        stored = {}
        conversions = []
        for name, value in parameters.items():
            if not isinstance(name, str) or not name.isidentifier():
                raise ValueError(f"parameter name {name!r} must be an identifier")
            quantity = _parameter_quantity(value, name)
            stored[name] = quantity
            conversions.append(geometric_conversion(quantity.unit, name))
        object.__setattr__(self, "_parameters", MappingProxyType(stored))
        factors = np.array([item[0] for item in conversions], dtype=np.float64)
        powers = np.array([item[1] for item in conversions], dtype=np.float64)
        factors.flags.writeable = False
        powers.flags.writeable = False
        object.__setattr__(self, "_conversions", (factors, powers))
        object.__setattr__(self, "_native", None)
        self._check_vector(self._parameter_vector())

    def __setattr__(self, name: str, value: object) -> None:
        """Reject attribute assignment after construction.

        Raises
        ------
        AttributeError
            Always; forces are immutable.
        """
        raise AttributeError("forces are immutable")

    def __delattr__(self, name: str) -> None:
        """Reject attribute deletion after construction.

        Raises
        ------
        AttributeError
            Always; forces are immutable.
        """
        raise AttributeError("forces are immutable")

    @property
    def parameters(self) -> Mapping[str, u.Quantity]:
        """Read-only mapping from parameter name to its scalar quantity.

        Returns
        -------
        mapping
            Parameters in declaration order; plain numbers appear as
            dimensionless quantities.
        """
        return self._parameters

    @property
    def parameter_names(self) -> tuple[str, ...]:
        """Parameter names in declaration order.

        Returns
        -------
        tuple of str
            The order of the force block of the perturbed MCMC models.
        """
        return tuple(self._parameters)

    @property
    def parameter_units(self) -> Mapping[str, u.UnitBase]:
        """Unit of each parameter, keyed by name.

        Returns
        -------
        mapping
            The unit each parameter was given in; perturbed MCMC models read
            the force block in these units.
        """
        return MappingProxyType(
            {name: value.unit for name, value in self._parameters.items()}
        )

    def with_parameters(self, **values: object) -> ExternalForce:
        """Return a copy with some parameters replaced.

        Parameters
        ----------
        **values
            New values of existing parameters, as numbers or scalar
            quantities.

        Returns
        -------
        ExternalForce
            A new force of the same type and definition. Native resources
            that do not depend on parameters are shared with this one.

        Raises
        ------
        ValueError, TypeError, astropy.units.UnitConversionError
            If a name is unknown or a value is invalid for the force.
        """
        unknown = sorted(set(values) - set(self._parameters))
        if unknown:
            raise ValueError(f"unknown parameters: {', '.join(unknown)}")
        merged = dict(self._parameters)
        merged.update(values)
        clone = self._replace(merged)
        object.__setattr__(clone, "_native", self._native)
        return clone

    def __eq__(self, other: object) -> bool:
        """Compare type, definition and parameters with their units."""
        if not isinstance(other, ExternalForce):
            return NotImplemented
        return type(self) is type(other) and self._key() == other._key()

    def __hash__(self) -> int:
        """Hash the key used for equality."""
        return hash((type(self).__name__, self._key()))

    def _parameter_key(self) -> tuple:
        return tuple(
            (name, float(value.value), value.unit.to_string())
            for name, value in self._parameters.items()
        )

    def _parameter_vector(self) -> np.ndarray:
        """Parameter values in their own units, in declaration order."""
        return np.array(
            [float(value.value) for value in self._parameters.values()],
            dtype=np.float64,
        )

    def _mcmc_parameters(self, prefix: str = "force") -> tuple[tuple[str, u.UnitBase], ...]:
        """``(name, unit)`` pairs of the force block of an MCMC vector."""
        return tuple(
            (f"{prefix}.{name}", value.unit) for name, value in self._parameters.items()
        )

    def _bind(
        self, mass: u.Quantity | float, rotation: np.ndarray | None = None
    ) -> BoundForce:
        """Bind this force to a black-hole mass (kilograms or a quantity)."""
        mass_kg = float(u.Quantity(mass, u.kg).value)
        bound = self._handle()
        bound.update(mass_kg, self._parameter_vector(), rotation)
        return bound

    def _handle(self) -> BoundForce:
        """Return a new handle whose parameters :meth:`BoundForce.update` writes."""
        callback, keepalive = self._native_callback()
        return BoundForce(self, callback, self._new_storage(), keepalive)

    def _check_vector(self, values: np.ndarray) -> np.ndarray:
        """Validate parameter values in their units; raise ``ValueError``."""
        if not np.all(np.isfinite(values)):
            raise ValueError("force parameters must be finite")
        return values

    # Subclass hooks.
    def _key(self) -> tuple:
        raise NotImplementedError

    def _replace(self, parameters: dict[str, object]) -> ExternalForce:
        raise NotImplementedError

    def _native_callback(self) -> tuple[int, object]:
        """Callback address and the native object it depends on."""
        raise NotImplementedError

    def _new_storage(self) -> object:
        """Model data of a new handle; default: one double per parameter."""
        return np.zeros(max(len(self._parameters), 1), dtype=np.float64)

    def _storage_address(self, storage: object) -> int:
        """Address of the model data passed to the native callback."""
        return int(storage.ctypes.data)

    def _storage_data(self, storage: object) -> np.ndarray:
        """Geometric model data, for inspection and tests."""
        return storage

    def _write(
        self,
        storage: object,
        mass: float,
        values: np.ndarray,
        rotation: np.ndarray | None,
    ) -> None:
        """Write geometric model data; default: one double per parameter."""
        factors, powers = self._conversions
        storage[: factors.size] = values * factors * mass ** (-powers)
