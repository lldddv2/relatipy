"""Lazy, read-only views of stored canonical Kerr states."""

from __future__ import annotations

import numpy as np
from astropy import units as u
from astropy.constants import c

from .._validation import readonly_quantity
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
from . import _native as native
from .state import InitialState, State


def _restore_deferred_state(
    cls: type[DeferredState], canonical: np.ndarray, tau: np.ndarray,
    spin: float, length_scale: u.Quantity, time_scale: u.Quantity,
    length_unit: u.UnitBase, time_unit: u.UnitBase, input_family: str,
    tau_physical: u.Quantity, scalar: bool,
) -> DeferredState:
    """Rebuild a copied or unpickled lazy state without eager extra views."""
    return cls(
        canonical=canonical, tau=tau, spin=spin,
        length_scale=length_scale, time_scale=time_scale,
        length_unit=length_unit, time_unit=time_unit,
        input_family=input_family, tau_physical=tau_physical, scalar=scalar,
    )


class DeferredState(State):
    """Expose canonical Boyer--Lindquist samples and cache other views on demand.

    This is an internal construction path. The public ``State(...)`` constructor
    retains its explicit, fully validated input contract.
    """

    __slots__ = ("_canonical", "_scalar", "_spin", "_length_scale", "_time_scale",
                 "_length_unit", "_time_unit", "_input_family", "_views")

    def __init__(
        self,
        *,
        canonical: np.ndarray,
        tau: np.ndarray,
        spin: float,
        length_scale: u.Quantity,
        time_scale: u.Quantity,
        length_unit: u.UnitBase,
        time_unit: u.UnitBase,
        input_family: str,
        tau_physical: u.Quantity | None = None,
        scalar: bool = False,
    ) -> None:
        """Store canonical samples and initialize the required BL views.

        Parameters
        ----------
        canonical : numpy.ndarray
            Finite native state array of shape ``(n, 8)``.
        tau : numpy.ndarray
            Finite normalized proper times of shape ``(n,)``.
        spin : float
            Dimensionless Kerr spin used by native reconstruction.
        length_scale, time_scale : astropy.units.Quantity
            Physical geometric scales for the source metric.
        length_unit, time_unit : astropy.units.UnitBase
            Presentation units for reconstructed lengths and times.
        input_family : str
            Initial-condition family whose view is populated eagerly.
        tau_physical : astropy.units.Quantity or None, optional
            Already reconstructed proper time. If omitted, it is reconstructed
            from ``tau`` and ``time_scale``.
        scalar : bool, optional
            Whether the stored row represents one scalar state.

        Raises
        ------
        RuntimeError
            If state or proper-time shapes, or reconstructed time shape, are
            inconsistent.
        ValueError
            If native states or normalized proper times are non-finite.
        """
        values = np.asarray(canonical, dtype=np.float64)
        proper = np.asarray(tau, dtype=np.float64)
        if values.ndim != 2 or values.shape[1] != 8:
            raise RuntimeError("canonical states must have shape (n, 8)")
        if proper.shape != (len(values),):
            raise RuntimeError("proper times do not match canonical states")
        if not np.all(np.isfinite(values)) or not np.all(np.isfinite(proper)):
            raise ValueError("canonical states and proper times must be finite")
        frozen = np.frombuffer(np.ascontiguousarray(values).tobytes(), dtype=np.float64)
        frozen = frozen.reshape(values.shape)
        object.__setattr__(self, "_canonical", frozen)
        object.__setattr__(self, "_scalar", bool(scalar))
        object.__setattr__(self, "_spin", float(spin))
        object.__setattr__(self, "_length_scale", length_scale)
        object.__setattr__(self, "_time_scale", time_scale)
        object.__setattr__(self, "_length_unit", length_unit)
        object.__setattr__(self, "_time_unit", time_unit)
        object.__setattr__(self, "_input_family", input_family)
        object.__setattr__(self, "_views", {})

        physical = (proper * time_scale).to(time_unit) if tau_physical is None else tau_physical.to(time_unit)
        if scalar and tau_physical is None:
            physical = physical[0]
        expected_shape = () if scalar else proper.shape
        if physical.shape != expected_shape:
            raise RuntimeError("physical proper time does not match canonical state shape")
        object.__setattr__(self, "tau", readonly_quantity(physical, u.s, "tau"))
        self._build_bl()
        if input_family in ("cartesian", "spherical", "elements"):
            self._ensure(input_family)

    def _length(self, values: np.ndarray) -> u.Quantity:
        """Convert normalized lengths to the stored presentation unit."""
        return (values * self._length_scale).to(self._length_unit)

    def _time(self, values: np.ndarray) -> u.Quantity:
        """Convert normalized times to the stored presentation unit."""
        return (values * self._time_scale).to(self._time_unit)

    def _speed(self, values: np.ndarray) -> u.Quantity:
        """Convert normalized speeds to the stored length-per-time unit."""
        return (values * c).to(self._length_unit / self._time_unit)

    def _rate(self, values: np.ndarray) -> u.Quantity:
        """Convert normalized angular rates to radians per stored time unit."""
        return (values * u.rad / self._time_scale).to(u.rad / self._time_unit)

    def _build_bl(self) -> None:
        """Populate cached Boyer--Lindquist coordinate and velocity views."""
        row = self._canonical[0] if self._scalar else self._canonical
        ut = row[..., 4] * u.one
        t = self._time(row[..., 0])
        self._views["bl"] = (
            BoyerLindquistCoordinates(t, self._length(row[..., 1]), row[..., 2] * u.rad, row[..., 3] * u.rad),
            BoyerLindquistVelocity(
                self._speed(row[..., 5] / row[..., 4]),
                self._rate(row[..., 6] / row[..., 4]),
                self._rate(row[..., 7] / row[..., 4]),
            ),
            BoyerLindquistFourVelocity(
                self._speed(row[..., 5]), self._rate(row[..., 6]), self._rate(row[..., 7])
            ),
            readonly_quantity(ut, u.one, "ut"),
        )

    def _ensure(self, family: str) -> None:
        """Reconstruct and cache one non-Boyer--Lindquist state family.

        Parameters
        ----------
        family : {"cartesian", "spherical", "elements"}
            Native reconstruction family to cache.

        Raises
        ------
        ValueError
            If ``family`` is unknown or the native reconstruction rejects a
            row or returns values outside the documented finite-value contract.
        RuntimeError
            If the native result has an invalid shape.
        """
        if family in self._views:
            return
        from .. import _core

        rows, statuses = _core.reconstruct_canonical_family_batch(
            self._spin, np.ascontiguousarray(self._canonical), family
        )
        rows = np.asarray(rows, dtype=np.float64)
        statuses = np.asarray(statuses)
        expected = 6 if family == "elements" else 7
        if rows.shape != (len(self._canonical), expected) or statuses.shape != (len(self._canonical),):
            raise RuntimeError("native state reconstruction returned invalid shape")
        failed = np.flatnonzero(statuses != 0)
        if failed.size:
            index = int(failed[0])
            raise ValueError(
                f"state row {index} cannot be reconstructed "
                f"({native.describe_reconstruction_status(statuses[index])})"
            )
        if self._scalar:
            rows = rows[0]
        if family == "cartesian":
            if not np.all(np.isfinite(rows)):
                raise ValueError("native Cartesian reconstruction returned non-finite values")
            view = CartesianCoordinates(
                self._time(rows[..., 0]), self._length(rows[..., 1]),
                self._length(rows[..., 2]), self._length(rows[..., 3]),
            )
            velocity = self._speed(rows[..., 4:7])
            ut = self._canonical[0, 4] if self._scalar else self._canonical[:, 4, None]
            spatial_four_velocity = self._speed(rows[..., 4:7] * ut)
            self._views[family] = (view, velocity, spatial_four_velocity,
                                   CartesianStateVector(xyz=view.xyz, vxyz=velocity))
        elif family == "spherical":
            if not np.all(np.isfinite(rows)):
                raise ValueError("native spherical reconstruction returned non-finite values")
            ut = self._canonical[0, 4] if self._scalar else self._canonical[:, 4]
            self._views[family] = (
                SphericalCoordinates(self._time(rows[..., 0]), self._length(rows[..., 1]),
                                     rows[..., 2] * u.rad, rows[..., 3] * u.rad),
                SphericalVelocity(self._speed(rows[..., 4]), self._rate(rows[..., 5]),
                                  self._rate(rows[..., 6])),
                SphericalFourVelocity(self._speed(rows[..., 4] * ut),
                                      self._rate(rows[..., 5] * ut),
                                      self._rate(rows[..., 6] * ut)),
            )
        elif family == "elements":
            # The native contract permits +inf semimajor axis and NaN angles.
            axis = rows[..., 0]
            angles = rows[..., 2:6]
            if not (
                np.all(np.isfinite(axis) | np.isposinf(axis))
                and np.all(np.isfinite(rows[..., 1]))
                and np.all(np.isfinite(angles) | np.isnan(angles))
            ):
                raise ValueError("native orbital elements returned non-finite values")
            self._views[family] = OrbitalElements(
                a=self._length(rows[..., 0]), e=rows[..., 1], inc=rows[..., 2] * u.rad,
                Omega=rows[..., 3] * u.rad, omega=rows[..., 4] * u.rad,
                f=rows[..., 5] * u.rad,
            )
        else:
            raise ValueError(f"unknown state family {family!r}")

    def _select(self, index: object) -> DeferredState:
        """Return a lazy scalar or series selection of stored samples.

        Parameters
        ----------
        index : object
            Validated NumPy sample-axis index.

        Returns
        -------
        DeferredState
            New lazy state preserving the selected canonical rows and physical
            proper times.
        """
        selected = np.asarray(self._canonical[index], dtype=np.float64).reshape(-1, 8)
        tau = self.tau[index]
        scalar = isinstance(index, (int, np.integer))
        return DeferredState(
            canonical=selected, tau=np.zeros(len(selected)), spin=self._spin,
            length_scale=self._length_scale, time_scale=self._time_scale,
            length_unit=self._length_unit, time_unit=self._time_unit,
            input_family=self._input_family, tau_physical=tau, scalar=scalar,
        )

    def __reduce_ex__(self, protocol: int):
        """Preserve canonical data across copy, deep copy, and pickle."""
        return _restore_deferred_state, (
            type(self), np.array(self._canonical, copy=True),
            np.asarray((self.tau / self._time_scale).to_value(u.one)).reshape(-1),
            self._spin, self._length_scale, self._time_scale,
            self._length_unit, self._time_unit, self._input_family,
            self.tau, self._scalar,
        )

    @property
    def tRQP(self) -> BoyerLindquistCoordinates:
        """Return cached Boyer--Lindquist coordinates in stored units."""
        return self._views["bl"][0]

    @property
    def vRQP(self) -> BoyerLindquistVelocity:
        """Return cached Boyer--Lindquist coordinate velocity."""
        return self._views["bl"][1]

    @property
    def uRQP(self) -> BoyerLindquistFourVelocity:
        """Return cached Boyer--Lindquist spatial four-velocity."""
        return self._views["bl"][2]

    @property
    def ut(self) -> u.Quantity:
        """Return read-only dimensionless temporal four-velocity values."""
        return self._views["bl"][3]

    @property
    def t(self) -> u.Quantity:
        """Return read-only scalar or shape ``(n,)`` coordinate time."""
        return self.tRQP.t

    @property
    def txyz(self) -> CartesianCoordinates:
        """Return cached or reconstructed Cartesian coordinates."""
        self._ensure("cartesian")
        return self._views["cartesian"][0]

    @property
    def vxyz(self) -> u.Quantity:
        """Return Cartesian velocity with shape ``(3,)`` or ``(n, 3)``."""
        self._ensure("cartesian")
        return self._views["cartesian"][1]

    @property
    def uxyz(self) -> u.Quantity:
        """Return Cartesian four-velocity with shape ``(3,)`` or ``(n, 3)``."""
        self._ensure("cartesian")
        return self._views["cartesian"][2]

    @property
    def state_vector(self) -> CartesianStateVector:
        """Return the cached Cartesian position and coordinate-velocity pair."""
        self._ensure("cartesian")
        return self._views["cartesian"][3]

    @property
    def trqp(self) -> SphericalCoordinates:
        """Return cached or reconstructed spherical coordinates."""
        self._ensure("spherical")
        return self._views["spherical"][0]

    @property
    def vrqp(self) -> SphericalVelocity:
        """Return cached or reconstructed spherical coordinate velocity."""
        self._ensure("spherical")
        return self._views["spherical"][1]

    @property
    def urqp(self) -> SphericalFourVelocity:
        """Return cached or reconstructed spherical spatial four-velocity."""
        self._ensure("spherical")
        return self._views["spherical"][2]

    def orbital_elements(self) -> OrbitalElements:
        """Return cached or reconstructed instantaneous osculating elements."""
        self._ensure("elements")
        return self._views["elements"]


class DeferredInitialState(DeferredState, InitialState):
    """Lazy representation of one saved initial condition."""
