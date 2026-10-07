"""Private immutable states and termination records for null geodesics.

Coordinate views are reconstructed by C when first requested. No proper time,
affine parameter, four-velocity, or orbital elements are exposed.
"""

from __future__ import annotations

from dataclasses import dataclass
from operator import attrgetter

import numpy as np
from astropy import units as u

from .._validation import immutable_array, readonly_quantity
from ..coordinates import (
    BoyerLindquistCoordinates, BoyerLindquistVelocity, CartesianCoordinates,
    CartesianStateVector, SphericalCoordinates, SphericalVelocity,
)
from . import _native as native
from .exceptions import IntegrationTerminated


NULL_PROPERTIES = (
    "t", "txyz", "trqp", "tRQP", "x", "y", "z", "xyz", "r", "theta",
    "phi", "R", "Theta", "Phi", "vxyz", "vrqp", "vRQP", "vx", "vy", "vz",
    "vr", "vtheta", "vphi", "vR", "vTheta", "vPhi", "state_vector",
)


def _restore_null_state(
    canonical: np.ndarray, spin: float, length_scale: u.Quantity,
    time_scale: u.Quantity, length_unit: u.UnitBase, time_unit: u.UnitBase,
    t_physical: u.Quantity, scalar: bool, xyz_physical: u.Quantity | None,
) -> NullState:
    """Rebuild a copied or unpickled null state without eager views."""
    return NullState(
        canonical=canonical, spin=spin, length_scale=length_scale,
        time_scale=time_scale, length_unit=length_unit, time_unit=time_unit,
        t_physical=t_physical, scalar=scalar, xyz_physical=xyz_physical,
    )


class NullState:
    """Immutable scalar or series of native null states with lazy views.

    Parameters
    ----------
    canonical : numpy.ndarray
        Finite native states with shape ``(n, 8)``.
    spin : float
        Dimensionless spin of the source metric.
    length_scale, time_scale : astropy.units.Quantity
        Geometric length and time scales of the source metric.
    length_unit, time_unit : astropy.units.UnitBase
        Presentation units for lengths and coordinate time.
    t_physical : astropy.units.Quantity or None, optional
        Exact physical coordinate times, overriding unit roundoff.
    xyz_physical : astropy.units.Quantity or None, optional
        Exact initial Cartesian positions, with shape ``(3,)`` for a scalar
        or ``(n, 3)`` for a series. Native coordinate velocities remain lazy.
    scalar : bool, optional
        Whether the single stored row represents a scalar state.

    Raises
    ------
    ValueError
        If the native array or physical time has invalid shape or values.

    Notes
    -----
    Observer frames, redshift, past integration and the direction root choice
    in the ergoregion remain pending. Coordinate velocities use ``t``.
    """

    __slots__ = (
        "_canonical", "_scalar", "_spin", "_length_scale", "_time_scale",
        "_length_unit", "_time_unit", "_t", "_xyz_physical", "_views",
    )

    def __init__(
        self, *, canonical: np.ndarray, spin: float,
        length_scale: u.Quantity, time_scale: u.Quantity,
        length_unit: u.UnitBase, time_unit: u.UnitBase,
        t_physical: u.Quantity | None = None, scalar: bool = False,
        xyz_physical: u.Quantity | None = None,
    ) -> None:
        values = np.asarray(canonical, dtype=np.float64)
        if values.ndim != 2 or values.shape[1] != 8:
            raise ValueError("canonical null states must have shape (n, 8)")
        if scalar and len(values) != 1:
            raise ValueError("a scalar null state requires one canonical row")
        if not np.all(np.isfinite(values)):
            raise ValueError("canonical null states must be finite")
        physical = (
            (values[:, 0] * time_scale).to(time_unit)
            if t_physical is None else t_physical.to(time_unit)
        )
        if scalar and t_physical is None:
            physical = physical[0]
        expected = () if scalar else (len(values),)
        if physical.shape != expected or not np.all(np.isfinite(physical.value)):
            raise ValueError("physical coordinate time must match null state shape")
        xyz = None
        if xyz_physical is not None:
            xyz = readonly_quantity(
                xyz_physical, u.m, "xyz_physical", ndim=(1, 2),
            ).to(length_unit)
            xyz = readonly_quantity(xyz, u.m, "xyz_physical", ndim=(1, 2))
            if xyz.shape != expected + (3,) or not np.all(np.isfinite(xyz.value)):
                raise ValueError("physical Cartesian position must match null state shape")
        for name, value in (
            ("_canonical", immutable_array(values)), ("_scalar", bool(scalar)),
            ("_spin", float(spin)), ("_length_scale", length_scale),
            ("_time_scale", time_scale), ("_length_unit", length_unit),
            ("_time_unit", time_unit),
            ("_t", readonly_quantity(physical, u.s, "t")),
            ("_xyz_physical", xyz), ("_views", {}),
        ):
            object.__setattr__(self, name, value)

    def __reduce_ex__(self, protocol: int):
        """Preserve stored rows and exact presentation data across copy and pickle."""
        return _restore_null_state, (
            np.array(self._canonical, copy=True), self._spin, self._length_scale,
            self._time_scale, self._length_unit, self._time_unit, self._t,
            self._scalar, self._xyz_physical,
        )

    def __setattr__(self, name: str, value: object) -> None:
        raise AttributeError("null states are immutable")

    def __delattr__(self, name: str) -> None:
        raise AttributeError("null states are immutable")

    def _length(self, values: np.ndarray) -> u.Quantity:
        return (values * self._length_scale).to(self._length_unit)

    def _speed(self, values: np.ndarray) -> u.Quantity:
        return (values * self._length_scale / self._time_scale).to(
            self._length_unit / self._time_unit
        )

    def _rate(self, values: np.ndarray) -> u.Quantity:
        return (values * u.rad / self._time_scale).to(u.rad / self._time_unit)

    def _ensure(self, family: str) -> None:
        if family in self._views:
            return
        from .. import _core

        rows, statuses = _core.null_views_batch(
            self._spin, family, np.ascontiguousarray(self._canonical)
        )
        rows, statuses = np.asarray(rows), np.asarray(statuses)
        if rows.shape != (len(self._canonical), 7) or statuses.shape != (len(rows),):
            raise RuntimeError("native null view reconstruction returned invalid shape")
        failed = np.flatnonzero(statuses != 0)
        if failed.size:
            row = int(failed[0])
            raise ValueError(
                f"null state row {row} cannot be reconstructed "
                f"({native.describe_reconstruction_status(statuses[row])})"
            )
        if not np.all(np.isfinite(rows)):
            raise ValueError("native null reconstruction returned non-finite values")
        if self._scalar:
            rows = rows[0]
        if family == "cartesian":
            xyz = (
                self._length(rows[..., 1:4])
                if self._xyz_physical is None else self._xyz_physical
            )
            coordinates = CartesianCoordinates(
                self.t, xyz[..., 0], xyz[..., 1], xyz[..., 2],
            )
            velocity = readonly_quantity(
                self._speed(rows[..., 4:7]), u.m / u.s, "vxyz", ndim=(1, 2)
            )
            self._views[family] = (
                coordinates, velocity,
                CartesianStateVector(xyz=coordinates.xyz, vxyz=velocity),
            )
        elif family in ("spherical", "bl"):
            coordinates_type = (
                SphericalCoordinates if family == "spherical" else BoyerLindquistCoordinates
            )
            velocity_type = (
                SphericalVelocity if family == "spherical" else BoyerLindquistVelocity
            )
            self._views[family] = (
                coordinates_type(
                    self.t, self._length(rows[..., 1]),
                    rows[..., 2] * u.rad, rows[..., 3] * u.rad,
                ),
                velocity_type(
                    self._speed(rows[..., 4]), self._rate(rows[..., 5]),
                    self._rate(rows[..., 6]),
                ),
            )
        else:
            raise ValueError(f"unknown null view family {family!r}")

    def _view(self, family: str, index: int):
        self._ensure(family)
        return self._views[family][index]

    def _select(self, index: object) -> NullState:
        return NullState(
            canonical=np.asarray(self._canonical[index]).reshape(-1, 8),
            spin=self._spin, length_scale=self._length_scale,
            time_scale=self._time_scale, length_unit=self._length_unit,
            time_unit=self._time_unit, t_physical=self.t[index],
            scalar=isinstance(index, (int, np.integer)),
            xyz_physical=(
                None if self._xyz_physical is None else self._xyz_physical[index]
            ),
        )

    @property
    def t(self) -> u.Quantity:
        """Return coordinate time.

        Returns
        -------
        astropy.units.Quantity
            Read-only scalar or shape ``(n,)`` time in the presentation unit.
        """
        return self._t


def _view_property(family: str, index: int, label: str, result: str) -> property:
    return property(
        lambda self: self._view(family, index),
        doc=f"""Return {label}, reconstructed and cached on first access.

        Returns
        -------
        {result}
            Read-only scalar or series view in the presentation units.
        """,
    )


for _name, _family, _index, _label, _result in (
    ("txyz", "cartesian", 0, "Cartesian coordinates", "CartesianCoordinates"),
    ("trqp", "spherical", 0, "spherical coordinates", "SphericalCoordinates"),
    ("tRQP", "bl", 0, "Boyer--Lindquist coordinates", "BoyerLindquistCoordinates"),
    ("vxyz", "cartesian", 1, "Cartesian coordinate velocity", "astropy.units.Quantity"),
    ("vrqp", "spherical", 1, "spherical coordinate velocity", "SphericalVelocity"),
    ("vRQP", "bl", 1, "Boyer--Lindquist coordinate velocity", "BoyerLindquistVelocity"),
    ("state_vector", "cartesian", 2, "Cartesian position and velocity", "CartesianStateVector"),
):
    setattr(NullState, _name, _view_property(_family, _index, _label, _result))

for _container, _components in (
    ("txyz", ("x", "y", "z", "xyz")),
    ("trqp", ("r", "theta", "phi")), ("tRQP", ("R", "Theta", "Phi")),
    ("vrqp", ("vr", "vtheta", "vphi")), ("vRQP", ("vR", "vTheta", "vPhi")),
):
    for _name in _components:
        setattr(NullState, _name, property(
            attrgetter(f"{_container}.{_name}"),
            doc=f"""Return coordinate component ``{_name}``.

            Returns
            -------
            astropy.units.Quantity
                Read-only component in the presentation units.
            """,
        ))

for _axis, _name in enumerate(("vx", "vy", "vz")):
    setattr(NullState, _name, property(
        lambda self, axis=_axis: self.vxyz[..., axis],
        doc=f"""Return Cartesian coordinate velocity ``{_name}``.

        Returns
        -------
        astropy.units.Quantity
            Read-only component in presentation length per time units.
        """,
    ))


def delegate_null_properties(owner: type, attribute: str) -> None:
    """Install the approved coordinate-time views on a wrapper class."""
    for name in NULL_PROPERTIES:
        setattr(owner, name, property(
            attrgetter(f"{attribute}.{name}"), doc=getattr(NullState, name).__doc__,
        ))


@dataclass(frozen=True, slots=True, eq=False)
class NullTermination:
    """Immutable terminal-event information for a null geodesic.

    Parameters
    ----------
    reason : str
        Non-empty internal event identifier, normally ``horizon`` or ``escape``.
    t : astropy.units.Quantity
        Scalar coordinate time equal to ``state.t``.
    state : NullState
        Last valid scalar null state; it may be outside the sampled solution.

    Raises
    ------
    TypeError, astropy.units.UnitConversionError, ValueError
        If the reason, time, or state is invalid or inconsistent.

    Notes
    -----
    Copying and pickling rebuild the record through validation and preserve
    the read-only coordinate time.
    """

    reason: str
    t: u.Quantity
    state: NullState

    def __post_init__(self) -> None:
        if not isinstance(self.reason, str):
            raise TypeError("reason must be a string")
        if not self.reason:
            raise ValueError("reason must not be empty")
        if not isinstance(self.state, NullState):
            raise TypeError("state must be a NullState")
        time = readonly_quantity(self.t, u.s, "t", ndim=(0,))
        if not self.state._scalar:
            raise ValueError("termination state must be scalar")
        if not np.array_equal(time.to_value(self.state.t.unit), self.state.t.value):
            raise ValueError("termination t must match state.t")
        object.__setattr__(self, "t", time)

    def __reduce_ex__(self, protocol: int):
        """Rebuild through validation after copying or unpickling."""
        return type(self), (self.reason, self.t, self.state)


class _NullIntegrationTerminated(IntegrationTerminated):
    """Coordinate-time termination compatible with IntegrationTerminated."""

    __slots__ = ()

    def __init__(
        self, reason: str, t: u.Quantity, state: NullState, message: str | None = None,
    ) -> None:
        self._termination = NullTermination(reason, t, state)
        RuntimeError.__init__(self, message if message is not None else reason)

    def __reduce_ex__(self, protocol: int):
        """Preserve the validated terminal payload and message when restored."""
        return type(self), (self.reason, self.t, self.state, str(self))

    @property
    def t(self) -> u.Quantity:
        """Return the last valid state's coordinate time."""
        return self._termination.t

    @property
    def tau(self):
        """Reject proper-time access for null geodesics."""
        raise AttributeError("null geodesics do not expose proper time")

    @property
    def termination(self) -> NullTermination:
        """Return immutable coordinate-time termination information."""
        return self._termination
