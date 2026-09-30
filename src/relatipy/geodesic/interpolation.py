"""Private solution interpolation and assembly of reconstructed states."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from astropy import units as u
from astropy.constants import c
from scipy.interpolate import CubicSpline, PchipInterpolator
from scipy.optimize import brentq

from . import _native as native
from .state import State

if TYPE_CHECKING:
    from ..metrics.kerr import Kerr

# State accessor holding the stored value of each reconstructed column.
_STORED_COLUMN = {
    native.BL_R: "R", native.BL_THETA: "Theta", native.BL_PHI: "Phi",
    native.BL_VR: "vR", native.BL_VTHETA: "vTheta", native.BL_VPHI: "vPhi",
    native.UT: "ut",
    native.BL_UR: "uR", native.BL_UTHETA: "uTheta", native.BL_UPHI: "uPhi",
    native.SPH_R: "r", native.SPH_THETA: "theta", native.SPH_PHI: "phi",
    native.SPH_VR: "vr", native.SPH_VTHETA: "vtheta", native.SPH_VPHI: "vphi",
    native.UX: "ux", native.UY: "uy", native.UZ: "uz",
    native.SPH_UR: "ur", native.SPH_UTHETA: "utheta", native.SPH_UPHI: "uphi",
}
_STORED_ELEMENT = {
    native.SEMIMAJOR: "a", native.ECCENTRICITY: "e",
    native.INCLINATION: "inc", native.ASCENDING_NODE: "Omega",
    native.PERIAPSIS_ARGUMENT: "omega", native.TRUE_ANOMALY: "f",
}


@dataclass(frozen=True, slots=True, eq=False)
class InterpolatedCartesian:
    """Cartesian query values in stored units and their source-row matches.

    Attributes
    ----------
    tau, t : numpy.ndarray
        Resolved proper and coordinate times, each with shape ``(m,)`` in the
        source state's presentation units.
    xyz, vxyz : numpy.ndarray
        Cartesian position and coordinate velocity with shapes ``(m, 3)``.
        Their numeric units are carried by the source state, not these arrays.
    exact_indices : numpy.ndarray
        Integer array of shape ``(m,)``. A non-negative entry selects an exact
        stored row; ``-1`` identifies a reconstructed interpolation row.
    scalar : bool
        Whether the original query was scalar.
    """

    tau: np.ndarray
    t: np.ndarray
    xyz: np.ndarray
    vxyz: np.ndarray
    exact_indices: np.ndarray
    scalar: bool


def _time_query(
    value: u.Quantity, unit: u.UnitBase, name: str
) -> tuple[np.ndarray, bool]:
    """Validate a scalar or one-dimensional time query in one unit.

    Returns normalized values with shape ``(m,)`` and whether the input was
    scalar. Raises ``TypeError`` for non-quantities and ``ValueError`` for an
    invalid dimension or non-finite value.
    """
    if not isinstance(value, u.Quantity):
        raise TypeError(f"{name} must be an Astropy time quantity")
    if value.ndim not in (0, 1):
        raise ValueError(f"{name} must be scalar or one-dimensional")
    values = np.atleast_1d(np.asarray(value.to_value(unit), dtype=np.float64))
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{name} must contain only finite values")
    return values, value.ndim == 0


def _exact_indices(samples: np.ndarray, query: np.ndarray) -> np.ndarray:
    """Locate exact sorted-sample matches, marking other queries with ``-1``."""
    positions = np.searchsorted(samples, query)
    result = np.full(query.shape, -1, dtype=np.intp)
    valid = positions < samples.size
    matched = np.zeros(query.shape, dtype=bool)
    matched[valid] = samples[positions[valid]] == query[valid]
    result[matched] = positions[matched]
    return result


def _resolve_tau_from_t(
    tau_samples: np.ndarray,
    t_samples: np.ndarray,
    requested_t: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Invert monotone PCHIP t(tau), preserving exact stored times."""
    if t_samples.size == 1:
        exact = np.where(requested_t == t_samples[0], 0, -1).astype(np.intp)
        if np.any(exact < 0):
            raise ValueError("t lies outside the stored solution domain")
        return np.full(requested_t.shape, tau_samples[0]), exact

    differences = np.diff(t_samples)
    if not (np.all(differences > 0) or np.all(differences < 0)):
        raise ValueError("t(tau) is not strictly monotone; inversion is not unique")

    low, high = min(t_samples[0], t_samples[-1]), max(t_samples[0], t_samples[-1])
    if np.any((requested_t < low) | (requested_t > high)):
        raise ValueError("t lies outside the stored solution domain")

    increasing = differences[0] > 0
    ordered_t = t_samples if increasing else t_samples[::-1]
    ordered_exact = _exact_indices(ordered_t, requested_t)
    exact = ordered_exact if increasing else np.where(
        ordered_exact >= 0, t_samples.size - 1 - ordered_exact, -1
    )
    resolved = np.empty_like(requested_t)
    matched = exact >= 0
    resolved[matched] = tau_samples[exact[matched]]
    if np.any(~matched):
        spline = PchipInterpolator(tau_samples, t_samples, extrapolate=False)
        for row in np.flatnonzero(~matched):
            value = requested_t[row]
            upper_ordered = int(np.searchsorted(ordered_t, value))
            left = upper_ordered - 1
            right = upper_ordered
            if increasing:
                lower_tau, upper_tau = tau_samples[left], tau_samples[right]
            else:
                lower_tau = tau_samples[t_samples.size - 1 - right]
                upper_tau = tau_samples[t_samples.size - 1 - left]
            resolved[row] = brentq(
                lambda tau: float(spline(tau) - value),
                lower_tau,
                upper_tau,
                xtol=np.nextafter(0.0, 1.0),
            )
    return resolved, exact


def interpolate_cartesian(
    state: State,
    *,
    tau: u.Quantity | None = None,
    t: u.Quantity | None = None,
) -> InterpolatedCartesian:
    """Interpolate stored Cartesian position and velocity independently.

    Parameters
    ----------
    state : State
        Stored scalar-series state whose sample time is strictly increasing.
    tau, t : astropy.units.Quantity or None, optional
        Exactly one scalar or one-dimensional query. ``tau`` is proper time;
        ``t`` is coordinate time.

    Returns
    -------
    InterpolatedCartesian
        Query values in the source presentation units. Numeric time arrays
        have shape ``(m,)``; position and velocity arrays have shape ``(m, 3)``.

    Raises
    ------
    TypeError
        If a supplied query is not an Astropy quantity.
    ValueError
        If zero or both time arguments are supplied, a query is non-finite or
        has invalid shape, lies outside the stored domain, cannot uniquely
        invert coordinate time, or a single sample requires interpolation.
    RuntimeError
        If native Cartesian reconstruction returns an invalid result shape.

    Notes
    -----
    Inputs and outputs use the state's presentation units. Exact queries copy
    stored values; only points between samples pass through the splines.
    """
    if (tau is None) == (t is None):
        raise ValueError("exactly one of tau or t is required")

    tau_samples = np.asarray(state.tau.to_value(state.tau.unit), dtype=np.float64)
    t_samples = np.asarray(state.t.to_value(state.t.unit), dtype=np.float64)
    if tau is not None:
        resolved_tau, scalar = _time_query(tau, state.tau.unit, "tau")
        if np.any((resolved_tau < tau_samples[0]) | (resolved_tau > tau_samples[-1])):
            raise ValueError("tau lies outside the stored solution domain")
        exact = _exact_indices(tau_samples, resolved_tau)
    else:
        requested_t, scalar = _time_query(t, state.t.unit, "t")
        resolved_tau, exact = _resolve_tau_from_t(
            tau_samples, t_samples, requested_t
        )

    from ._deferred import DeferredState

    if isinstance(state, DeferredState) and "cartesian" not in state._views:
        from .. import _core

        rows, statuses = _core.reconstruct_canonical_family_batch(
            state._spin, np.ascontiguousarray(state._canonical), "cartesian"
        )
        rows = np.asarray(rows, dtype=np.float64)
        statuses = np.asarray(statuses)
        if rows.shape != (len(tau_samples), 7) or statuses.shape != (len(tau_samples),):
            raise RuntimeError("native Cartesian reconstruction returned invalid shape")
        failed = np.flatnonzero(statuses != 0)
        if failed.size:
            row = int(failed[0])
            raise ValueError(
                f"state row {row} cannot be reconstructed "
                f"({native.describe_reconstruction_status(statuses[row])})"
            )
        xyz_samples = (rows[:, 1:4] * state._length_scale).to_value(state._length_unit)
        vxyz_samples = (rows[:, 4:7] * c).to_value(
            state._length_unit / state._time_unit
        )
    else:
        xyz_samples = np.asarray(state.xyz.to_value(state.x.unit), dtype=np.float64)
        vxyz_samples = np.asarray(state.vxyz.to_value(state.vxyz.unit), dtype=np.float64)
    result_t = np.empty_like(resolved_tau)
    result_xyz = np.empty((resolved_tau.size, 3), dtype=np.float64)
    result_vxyz = np.empty((resolved_tau.size, 3), dtype=np.float64)
    matched = exact >= 0
    result_t[matched] = t_samples[exact[matched]]
    result_xyz[matched] = xyz_samples[exact[matched]]
    result_vxyz[matched] = vxyz_samples[exact[matched]]
    if np.any(~matched):
        if tau_samples.size < 2:
            raise ValueError("a single stored sample only supports exact selection")
        unknown = resolved_tau[~matched]
        result_t[~matched] = PchipInterpolator(
            tau_samples, t_samples, extrapolate=False
        )(unknown)
        result_xyz[~matched] = CubicSpline(
            tau_samples, xyz_samples, axis=0, extrapolate=False
        )(unknown)
        result_vxyz[~matched] = CubicSpline(
            tau_samples, vxyz_samples, axis=0, extrapolate=False
        )(unknown)
    return InterpolatedCartesian(
        tau=resolved_tau,
        t=result_t,
        xyz=result_xyz,
        vxyz=result_vxyz,
        exact_indices=exact,
        scalar=scalar,
    )


def exact_sample_indices(
    state: State, *, tau: u.Quantity | None = None, t: u.Quantity | None = None
) -> np.ndarray | int | None:
    """Return indices when every time query exactly matches stored samples.

    Parameters
    ----------
    state : State
        Source state containing stored sample times.
    tau, t : astropy.units.Quantity or None, optional
        Exactly one scalar or one-dimensional proper- or coordinate-time query.

    Returns
    -------
    int, numpy.ndarray, or None
        Scalar index, one-dimensional integer indices, or ``None`` when at
        least one valid query needs interpolation.

    Raises
    ------
    TypeError
        If a supplied query is not an Astropy quantity.
    ValueError
        If the time arguments are invalid or coordinate time is not uniquely
        monotone.

    Notes
    -----
    The full interpolation path retains responsibility for domain diagnostics.
    """
    if (tau is None) == (t is None):
        raise ValueError("exactly one of tau or t is required")
    if tau is not None:
        query, scalar = _time_query(tau, state.tau.unit, "tau")
        samples = np.asarray(state.tau.to_value(state.tau.unit), dtype=np.float64)
    else:
        query, scalar = _time_query(t, state.t.unit, "t")
        samples = np.asarray(state.t.to_value(state.t.unit), dtype=np.float64)
        if len(samples) > 1 and not (
            np.all(np.diff(samples) > 0) or np.all(np.diff(samples) < 0)
        ):
            raise ValueError("t(tau) is not strictly monotone; inversion is not unique")
    if len(samples) > 1 and samples[0] > samples[-1]:
        positions = _exact_indices(samples[::-1], query)
        indices = np.where(positions >= 0, len(samples) - 1 - positions, -1)
    else:
        indices = _exact_indices(samples, query)
    if np.any(indices < 0):
        return None
    return int(indices[0]) if scalar else indices


def reconstruct_interpolated_deferred(
    source: State, values: InterpolatedCartesian, metric: Kerr
) -> State:
    """Reconstruct only canonical BL rows for interpolated queries.

    Exact query rows are copied from the stored canonical series. Cartesian
    interpolation remains the established interpolation rule; other output
    coordinate families stay deferred.
    """
    from .. import _core
    from ._deferred import DeferredState

    assert isinstance(source, DeferredState)
    count = len(values.tau)
    canonical = np.empty((count, 8), dtype=np.float64)
    exact = values.exact_indices
    matched = exact >= 0
    canonical[matched] = source._canonical[exact[matched]]
    length_scale = metric._length_scale
    time_scale = metric._time_scale
    speed_scale = c.to(u.m / u.s)
    if np.any(~matched):
        cartesian = np.empty((np.count_nonzero(~matched), 7), dtype=np.float64)
        cartesian[:, 0] = (
            values.t[~matched] * source.t.unit / time_scale
        ).to_value(u.one)
        cartesian[:, 1:4] = (
            values.xyz[~matched] * source._length_unit / length_scale
        ).to_value(u.one)
        cartesian[:, 4:7] = (
            values.vxyz[~matched] * (source._length_unit / source._time_unit)
            / speed_scale
        ).to_value(u.one)
        converted, statuses = _core.initial_cartesian_batch(metric.spin, cartesian)
        failed = np.flatnonzero(statuses != 0)
        if failed.size:
            row = int(np.flatnonzero(~matched)[failed[0]])
            raise ValueError(
                f"interpolated state at query row {row} cannot be reconstructed "
                f"({native.describe_reconstruction_status(statuses[failed[0]])})"
            )
        canonical[~matched] = converted

    # Keep the stored azimuth branch across revolutions.
    if np.any(~matched):
        reference = np.interp(
            values.tau[~matched],
            source.tau.to_value(source.tau.unit),
            np.unwrap(source.Phi.to_value(u.rad)),
        )
        canonical[~matched, 3] += (2 * np.pi) * np.rint(
            (reference - canonical[~matched, 3]) / (2 * np.pi)
        )
    return DeferredState(
        canonical=canonical,
        tau=np.asarray((values.tau * source.tau.unit / time_scale).to_value(u.one)),
        spin=metric.spin,
        length_scale=length_scale,
        time_scale=time_scale,
        length_unit=source._length_unit,
        time_unit=source._time_unit,
        input_family=source._input_family,
        tau_physical=values.tau * source.tau.unit,
    )


def reconstruct_interpolated_state(
    source: State,
    values: InterpolatedCartesian,
    metric: Kerr,
) -> State:
    """Assemble an interpolated query, preserving exact samples and units.

    Reconstruct only nonsample rows in one native call. Merge those rows with
    untouched stored components in query order, preserving each source unit.
    The caller handles queries containing only exact samples and scalar shape.
    """
    from .. import _core

    exact = values.exact_indices
    interpolated = exact < 0
    count = int(np.count_nonzero(interpolated))
    length_scale = metric._length_scale
    time_scale = metric._time_scale
    speed_scale = c.to(u.m / u.s)
    scales = {
        "length": length_scale,
        "angle": u.rad,
        "speed": speed_scale,
        "rate": u.rad / time_scale,
        "one": u.one,
    }

    cartesian = np.empty((count, 7), dtype=np.float64)
    cartesian[:, 0] = (
        values.t[interpolated] * source.t.unit / time_scale
    ).to_value(u.one)
    cartesian[:, 1:4] = (
        values.xyz[interpolated] * source.x.unit / length_scale
    ).to_value(u.one)
    cartesian[:, 4:7] = (
        values.vxyz[interpolated] * source.vxyz.unit / speed_scale
    ).to_value(u.one)
    reconstructed, status = _core.reconstruct_batch(
        metric.spin, np.ascontiguousarray(cartesian)
    )
    reconstructed = np.asarray(reconstructed, dtype=np.float64)
    status = np.asarray(status)
    if reconstructed.shape != (count, native.ROW_WIDTH) or status.shape != (count,):
        raise RuntimeError("native state reconstruction returned invalid shape")
    failed = np.flatnonzero(status != 0)
    if failed.size:
        row = int(np.flatnonzero(interpolated)[failed[0]])
        raise ValueError(
            f"interpolated state at query row {row} cannot be reconstructed "
            f"({native.describe_reconstruction_status(status[failed[0]])})"
        )
    native.validate_reconstructed_rows(reconstructed)

    # C returns the principal atan2 azimuth. The stored phase only chooses
    # its 2π branch; the angle within that branch remains the native result.
    stored_tau = source.tau.to_value(source.tau.unit)
    for column, stored_angle in (
        (native.BL_PHI, source.Phi), (native.SPH_PHI, source.phi)
    ):
        reference = np.interp(
            values.tau[interpolated],
            stored_tau,
            np.unwrap(stored_angle.to_value(u.rad)),
        )
        reconstructed[:, column] += (2 * np.pi) * np.rint(
            (reference - reconstructed[:, column]) / (2 * np.pi)
        )

    source_elements = source.orbital_elements()

    def stored(column: int) -> u.Quantity:
        """Return the stored public quantity for one reconstructed column."""
        if column in _STORED_COLUMN:
            return getattr(source, _STORED_COLUMN[column])
        return u.Quantity(getattr(source_elements, _STORED_ELEMENT[column]))

    def merge(column: int) -> u.Quantity:
        """Combine stored exact rows with reconstructed interpolated rows."""
        stored_value = stored(column)
        unit = stored_value.unit
        output = np.empty(values.tau.size, dtype=np.float64)
        output[~interpolated] = stored_value.to_value(unit)[exact[~interpolated]]
        output[interpolated] = (
            reconstructed[:, column] * scales[native.COLUMN_KIND[column]]
        ).to_value(unit)
        return output * unit

    def merge_cartesian(stored_value: u.Quantity, axis: int) -> u.Quantity:
        """Merge one Cartesian component with reconstructed query rows."""
        unit = stored_value.unit
        output = (values.xyz[:, axis] * source.x.unit).to_value(unit)
        output[~interpolated] = stored_value.to_value(unit)[exact[~interpolated]]
        return output * unit

    return native.assemble_state(
        merge,
        tau=values.tau * source.tau.unit,
        t=values.t * source.t.unit,
        x=merge_cartesian(source.x, 0),
        y=merge_cartesian(source.y, 1),
        z=merge_cartesian(source.z, 2),
        vxyz=values.vxyz * source.vxyz.unit,
    )
