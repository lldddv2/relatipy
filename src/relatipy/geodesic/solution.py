"""Define immutable sampled solutions and temporal queries.

A :class:`Solution` wraps a vectorized :class:`~relatipy.geodesic.State`, the
effective integration metadata, and structured termination information.
Indexing copies selected stored samples into a new immutable state; it never
interpolates or reevaluates physical quantities. ``at`` interpolates between
stored samples and defers unrequested coordinate views.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np
from astropy import units as u

from .._validation import immutable_array
from ..coordinates import (
    BoyerLindquistCoordinates,
    BoyerLindquistFourVelocity,
    BoyerLindquistVelocity,
    CartesianCoordinates,
    OrbitalElements,
    SphericalCoordinates,
    SphericalFourVelocity,
    SphericalVelocity,
)
from ._delegation import delegate_state_properties
from ._deferred import DeferredState
from .integration import IntegrationInfo, Termination
from .state import State

if TYPE_CHECKING:
    from ..metrics.kerr import Kerr


def select_state(state: State, index: Any) -> State:
    """Copy selected stored samples into a new immutable state.

    Parameters
    ----------
    state : State
        Vectorized source state.
    index : object
        A validated integer, slice, integer array, or boolean mask applied to
        the sample axis.

    Returns
    -------
    State
        A scalar state for an integer index or a vectorized state for the
        other supported index forms.

    Raises
    ------
    IndexError
        If NumPy reports an invalid or out-of-bounds index.
    """
    if isinstance(state, DeferredState):
        return state._select(index)
    elements = state.orbital_elements()
    return State(
        tau=state.tau[index],
        txyz=CartesianCoordinates(
            t=state.t[index],
            x=state.x[index],
            y=state.y[index],
            z=state.z[index],
        ),
        trqp=SphericalCoordinates(
            t=state.t[index],
            r=state.r[index],
            theta=state.theta[index],
            phi=state.phi[index],
        ),
        tRQP=BoyerLindquistCoordinates(
            t=state.t[index],
            R=state.R[index],
            Theta=state.Theta[index],
            Phi=state.Phi[index],
        ),
        vxyz=state.vxyz[index],
        vrqp=SphericalVelocity(
            vr=state.vr[index],
            vtheta=state.vtheta[index],
            vphi=state.vphi[index],
        ),
        vRQP=BoyerLindquistVelocity(
            vR=state.vR[index],
            vTheta=state.vTheta[index],
            vPhi=state.vPhi[index],
        ),
        ut=state.ut[index],
        uxyz=state.uxyz[index],
        urqp=SphericalFourVelocity(
            ur=state.ur[index],
            utheta=state.utheta[index],
            uphi=state.uphi[index],
        ),
        uRQP=BoyerLindquistFourVelocity(
            uR=state.uR[index],
            uTheta=state.uTheta[index],
            uPhi=state.uPhi[index],
        ),
        orbital_elements=OrbitalElements(
            a=elements.a[index],
            e=np.asarray(elements.e)[index],
            inc=elements.inc[index],
            Omega=elements.Omega[index],
            omega=elements.omega[index],
            f=elements.f[index],
        ),
    )


def constants_of_motion(state: State) -> np.ndarray:
    """Return native normalized ``(E, Lz, Q)`` rows of a canonical state.

    Parameters
    ----------
    state : State
        Scalar or vectorized state built from native canonical samples.

    Returns
    -------
    numpy.ndarray
        Read-only array of shape ``(n, 3)``; a scalar state has ``n == 1``.

    Raises
    ------
    ValueError
        If the state carries no native canonical samples or a row is rejected
        by the native kernel.
    """
    from .. import _core
    from . import _native as native

    if not isinstance(state, DeferredState):
        raise ValueError(
            "constants of motion require native Kerr samples from an orbit"
        )
    rows, statuses = _core.constants_of_motion_batch(
        state._spin, np.ascontiguousarray(state._canonical)
    )
    rows = np.asarray(rows, dtype=np.float64)
    statuses = np.asarray(statuses)
    if rows.shape != (len(state._canonical), 3) or statuses.shape != (len(rows),):
        raise RuntimeError("native constants of motion returned invalid shape")
    failed = np.flatnonzero(statuses != 0)
    if failed.size:
        row = int(failed[0])
        raise ValueError(
            f"constants of motion of state row {row} are undefined "
            f"({native.describe_reconstruction_status(statuses[row])})"
        )
    rows.flags.writeable = False
    return rows


def _validate_index(index: Any, size: int) -> Any:
    """Validate the supported subset of NumPy sample-axis indexing.

    Parameters
    ----------
    index : object
        Integer, slice, one-dimensional integer array, or one-dimensional
        boolean mask.
    size : int
        Number of stored samples; boolean masks must have this length.

    Returns
    -------
    object
        The original integer or slice, or a one-dimensional NumPy array.

    Raises
    ------
    IndexError
        If ``index`` is a tuple, has more than one dimension, is an invalid
        mask shape, or has neither an integer nor boolean dtype.

    Examples
    --------
    >>> _validate_index([True, False], 2)
    array([ True, False])
    """
    if isinstance(index, tuple):
        raise IndexError("Solution accepts one sample-axis index")
    if isinstance(index, (int, np.integer, slice)):
        return index

    array = np.asarray(index)
    if array.ndim != 1:
        raise IndexError("index arrays and masks must be one-dimensional")
    if np.issubdtype(array.dtype, np.bool_):
        if array.shape != (size,):
            raise IndexError(
                f"boolean mask has shape {array.shape}, expected ({size},)"
            )
        return array
    if not np.issubdtype(array.dtype, np.integer):
        raise IndexError("index arrays must contain integers or booleans")
    return array


@dataclass(frozen=True, slots=True, init=False, eq=False)
class Solution:
    """A read-only sampled result from one integration.

    Parameters
    ----------
    state : State
        Vectorized state containing one row per stored integration sample.
        Proper times must be finite and strictly increasing.
    integration : IntegrationInfo
        Stored numerical configuration and diagnostics.
    status : int
        ``-1`` for numerical failure, ``0`` for reaching the requested end,
        or ``1`` for an internal terminal event.
    message : str
        Human-readable stored summary.
    termination : Termination or None
        Structured terminal-event information.  It is required exactly when
        ``status == 1``.

    Attributes
    ----------
    integration : IntegrationInfo
        Effective settings and diagnostics.
    status : int
        Integration outcome code.
    message : str
        Human-readable outcome summary.
    termination : Termination or None
        Terminal-event information when ``status == 1``.
    success : bool
        Whether ``status`` is ``0`` or ``1``.

    Notes
    -----
    Indexing selects stored samples. :meth:`at` interpolates in the stored
    domain. Non-sample queries require the Kerr metric privately supplied
    by the orbit that created this solution.

    Raises
    ------
    TypeError
        If ``state``, ``integration``, ``message``, or ``_metric`` has the
        wrong type.
    ValueError
        If the state is not a non-empty one-dimensional sample series, proper
        times are non-finite or not strictly increasing, ``status`` is not one
        of ``-1``, ``0``, and ``1``, or termination information is
        inconsistent with ``status``.
    """

    _state: State
    integration: IntegrationInfo
    status: int
    message: str
    termination: Termination | None
    _metric: Kerr | None

    def __init__(
        self,
        *,
        state: State,
        integration: IntegrationInfo,
        status: int,
        message: str,
        termination: Termination | None = None,
        _metric: Kerr | None = None,
    ) -> None:
        """Validate and store one sampled integration result."""
        if not isinstance(state, State):
            raise TypeError("state must be a State")
        if state.tau.ndim != 1:
            raise ValueError("Solution state must be a one-dimensional state series")
        if len(state.tau) == 0:
            raise ValueError("Solution must contain at least one stored sample")
        tau_values = state.tau.to_value(state.tau.unit)
        if not np.all(np.isfinite(tau_values)):
            raise ValueError("Solution proper times must be finite")
        if np.any(np.diff(tau_values) <= 0):
            raise ValueError("Solution proper times must be strictly increasing")
        if not isinstance(integration, IntegrationInfo):
            raise TypeError("integration must be an IntegrationInfo")
        if isinstance(status, bool) or status not in (-1, 0, 1):
            raise ValueError("status must be one of -1, 0 or 1")
        if not isinstance(message, str):
            raise TypeError("message must be a string")
        if status == 1 and not isinstance(termination, Termination):
            raise ValueError("termination is required when status is 1")
        if status != 1 and termination is not None:
            raise ValueError("termination must be None unless status is 1")
        if _metric is not None:
            from ..metrics.kerr import Kerr

            if not isinstance(_metric, Kerr):
                raise TypeError("_metric must be a Kerr metric")

        object.__setattr__(self, "_state", state)
        object.__setattr__(self, "integration", integration)
        object.__setattr__(self, "status", int(status))
        object.__setattr__(self, "message", message)
        object.__setattr__(self, "termination", termination)
        object.__setattr__(self, "_metric", _metric)

    def __len__(self) -> int:
        """Return the number of stored samples.

        Returns
        -------
        int
            Number of rows in the vectorized state.
        """
        return len(self._state.tau)

    def __getitem__(self, index: Any) -> State:
        """Select stored samples without interpolation.

        An integer returns a scalar :class:`State`; a slice, one-dimensional
        integer array, or boolean mask returns a vectorized :class:`State`.
        NumPy negative-index and bounds semantics are preserved.

        Parameters
        ----------
        index : int, slice, array-like of int, or array-like of bool
            Sample-axis selection.  Boolean masks must match the solution
            length.

        Returns
        -------
        State
            Immutable scalar or vectorized selection of stored samples.

        Raises
        ------
        IndexError
            If the index form, dimensionality, dtype, mask shape, or bounds
            are invalid.
        """
        return select_state(self._state, _validate_index(index, len(self)))

    def at(
        self,
        *,
        tau: u.Quantity | None = None,
        t: u.Quantity | None = None,
    ) -> State:
        """Select or interpolate states at proper or coordinate times.

        Parameters
        ----------
        tau, t : astropy.units.Quantity or None
            Exactly one time quantity, scalar or one-dimensional. ``tau`` is
            proper time; ``t`` is coordinate time and requires a unique,
            strictly monotone stored ``t(tau)`` relation.

        Returns
        -------
        State
            Immutable scalar or vectorized state matching the query shape.

        Raises
        ------
        TypeError
            If a time is not an Astropy quantity.
        ValueError
            If the query is invalid, outside the stored domain, coordinate
            time is not invertible, or an interpolated state cannot be
            reconstructed in the Kerr chart.

        Notes
        -----
        Stored samples are copied exactly. Between samples, Cartesian
        position and coordinate velocity use separate cubic splines, while
        ``t(tau)`` uses PCHIP. The interpolated canonical BL state is
        reconstructed without geodesic integration. Other coordinate views
        are computed when first accessed and then cached.
        Each reconstructed state reports ``NaN`` for undefined osculating
        angles and positive infinite ``a`` only when its computed Kepler
        energy is exactly zero. Elements describe the instantaneous Cartesian
        position and coordinate velocity; interpolation can change their
        classification. Inspect ``state.orbital_elements().defined`` for
        element availability.
        """
        from .interpolation import (
            exact_sample_indices,
            interpolate_cartesian,
            reconstruct_interpolated_deferred,
            reconstruct_interpolated_state,
        )

        indices = exact_sample_indices(self._state, tau=tau, t=t)
        if indices is not None:
            return select_state(self._state, indices)
        values = interpolate_cartesian(self._state, tau=tau, t=t)
        if np.all(values.exact_indices >= 0):
            indices = (
                int(values.exact_indices[0])
                if values.scalar
                else values.exact_indices
            )
            return select_state(self._state, indices)
        if self._metric is None:
            raise ValueError(
                "interpolation requires the Kerr metric of the creating orbit"
            )
        if isinstance(self._state, DeferredState):
            interpolated = reconstruct_interpolated_deferred(
                self._state, values, self._metric
            )
        else:
            interpolated = reconstruct_interpolated_state(
                self._state, values, self._metric
            )
        return select_state(interpolated, 0) if values.scalar else interpolated

    def plot(self, *, projection="views", interactive=None, style=None, fig=None, ax=None):
        """Plot stored Cartesian samples without integration or interpolation.

        Parameters
        ----------
        projection : {"3d", "xy", "xz", "yz", "views"}, optional
            Cartesian view, default ``"views"``. The ``"views"`` layout
            draws three orthogonal planes at one shared spatial scale.
        interactive : bool or None, optional
            ``True`` returns an interactive Plotly figure; ``False`` returns
            a static Matplotlib figure. ``None`` (default) selects Plotly for
            ``"3d"`` and Matplotlib for plane and ``"views"`` projections.
        style : relatipy.plotting.Style or None, optional
            Appearance settings. ``None`` uses the package default and stored
            length unit; ``length_unit="r_g"`` displays gravitational radii.
        fig : matplotlib.figure.Figure or plotly.graph_objects.Figure, optional
            Existing compatible figure. The ``"views"`` layout creates its
            own figure and does not accept it.
        ax : matplotlib.axes.Axes, optional
            Existing static axes belonging to ``fig``; valid only for static
            figures other than ``"views"``.

        Returns
        -------
        plotly.graph_objects.Figure or tuple
            Interactive Plotly figure, static ``(fig, ax)``, or
            ``(fig, (main, top, right))`` for ``"views"``.

        Raises
        ------
        ImportError
            If the required Plotly or Matplotlib dependency is unavailable.
        TypeError, ValueError
            If plotting arguments or stored positions are invalid.

        Notes
        -----
        Lines connect stored samples only. The method neither generates dense
        output nor changes the solution.
        """
        from ..plotting import plot_solution

        return plot_solution(self, projection=projection, interactive=interactive, style=style,
                             fig=fig, ax=ax)

    def plot_evol(self, *, time="t", coords="cartesian", style=None, **components):
        """Plot stored coordinate components against coordinate or proper time.

        Parameters
        ----------
        time : {"t", "tau"}, optional
            Horizontal coordinate, default coordinate time ``"t"``.
        coords : {"cartesian", "cartesial", "spherical", "bl", "elements"}, optional
            Coordinate or osculating-element family enabled by default.
            ``"elements"`` selects ``a``, ``e`` and ``f``. ``"cartesial"`` aliases
            ``"cartesian"``. Cartesian is the default family.
        style : relatipy.plotting.Style or None, optional
            Figure format.
        **components : bool
            Optional switches named ``x``, ``y``, ``z``, ``r``, ``theta``,
            ``phi``, ``R``, ``Theta``, ``Phi``, ``a``, ``e`` or ``f``. An explicit ``False``
            disables even a default component; ``True`` can enable one from
            another family.

        Returns
        -------
        fig : matplotlib.figure.Figure
            Figure without implicit display or saving.
        axes : numpy.ndarray
            One-dimensional array with one panel per selected component.

        Raises
        ------
        ValueError
            If time, family or component is unknown, or none is selected.
        TypeError
            If a component switch is not boolean or ``style`` is invalid.

        Notes
        -----
        Stored sample values are used without integration or interpolation.
        Lengths and time use stored units; angles are displayed in degrees.
        Large tick values use a power-of-ten label without changing data.
        """
        from ..plotting import plot_evolution

        return plot_evolution(self, time=time, coords=coords, style=style, **components)

    def plot_constants(self, *, time="t", drift=False, style=None, **constants):
        """Plot stored E, Lz and Carter Q against coordinate or proper time.

        Parameters
        ----------
        time : {"t", "tau"}, optional
            Horizontal coordinate, default coordinate time ``"t"``.
        drift : bool, optional
            If ``True``, plot ``(X - X_0) / |X_0|`` relative to the first
            sample, or ``X - X_0`` when ``|X_0| < 1e-12``. Default ``False``.
        style : relatipy.plotting.Style or None, optional
            Figure format.
        **constants : bool
            Switches named ``E``, ``Lz`` or ``Q``; all are enabled by default
            and an explicit ``False`` removes that panel.

        Returns
        -------
        fig : matplotlib.figure.Figure
            Figure without implicit display or saving.
        axes : numpy.ndarray
            One-dimensional array with one panel per selected constant.

        Raises
        ------
        ValueError
            If time or a constant name is unknown, none is selected, or the
            solution carries no native Kerr samples.
        TypeError
            If a switch is not boolean or ``style`` is invalid.

        Notes
        -----
        Values are normalized with ``G = c = M = mu = 1``, as returned by
        :meth:`get_E`, :meth:`get_Lz` and :meth:`get_Q`. Stored samples are
        used without integration or interpolation.
        """
        from ..plotting import plot_constants

        return plot_constants(self, time=time, drift=drift, style=style, **constants)

    @property
    def success(self) -> bool:
        """Report whether integration ended without numerical failure.

        Returns
        -------
        bool
            ``True`` when :attr:`status` is ``0`` or ``1``.
        """
        return self.status in (0, 1)

    def orbital_elements(self) -> OrbitalElements:
        """Return cached or newly computed osculating elements for every sample.

        Returns
        -------
        OrbitalElements
            Immutable scalar-series elements aligned with the stored proper
            times. The ``defined`` mask identifies available elements;
            positive infinite ``a`` represents an exactly parabolic conic.
        """
        return self._state.orbital_elements()


    def get_E(self) -> np.ndarray:
        """Return the specific energy at every stored sample.

        Returns
        -------
        numpy.ndarray
            Read-only float array of shape ``(n,)`` aligned with
            :attr:`tau`. Values are ``E = -u_t``, normalized with
            ``G = c = M = mu = 1``, i.e. ``E / (mu c^2)``.

        Raises
        ------
        ValueError
            If the solution was not produced by an orbit or a sample lies
            on a Boyer--Lindquist coordinate singularity.

        Notes
        -----
        Stored four-velocities are not renormalized, so the series exposes
        the numerical drift of the integration.
        """
        return immutable_array(constants_of_motion(self._state)[:, 0])

    def get_Lz(self) -> np.ndarray:
        """Return the specific axial angular momentum at every stored sample.

        Returns
        -------
        numpy.ndarray
            Read-only float array of shape ``(n,)`` aligned with
            :attr:`tau`. Values are ``Lz = u_phi`` in units of ``mu G M / c``
            (normalized with ``G = c = M = mu = 1``).

        Raises
        ------
        ValueError
            If the solution was not produced by an orbit or a sample lies
            on a Boyer--Lindquist coordinate singularity.
        """
        return immutable_array(constants_of_motion(self._state)[:, 1])

    def get_Q(self) -> np.ndarray:
        """Return the Carter constant at every stored sample.

        Returns
        -------
        numpy.ndarray
            Read-only float array of shape ``(n,)`` aligned with
            :attr:`tau`. Values are
            ``Q = u_theta**2 + cos(theta)**2 * (a**2 (1 - E**2) + Lz**2 / sin(theta)**2)``
            in units of ``(mu G M / c)**2``, with fixed rest mass ``mu = 1``.
            It is not ``K = Q + (Lz - a E)**2``.

        Raises
        ------
        ValueError
            If the solution was not produced by an orbit or a sample lies
            on a Boyer--Lindquist coordinate singularity.
        """
        return immutable_array(constants_of_motion(self._state)[:, 2])

delegate_state_properties(Solution, "_state")
