"""Evaluate Kerr stellar observables directly from MCMC parameters.

The model uses the same osculating Kepler initial state as Kerr.orbit, evolves a
timelike Kerr geodesic, and applies a straight-line Rømer observation model.
It does not trace photons or include Shapiro delay and gravitational lensing.
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
from astropy import constants as const
from astropy import units as u

from .._validation import (
    REAL_TYPES,
    dimensionless_scalar,
    immutable_array,
    non_negative_integer,
    non_negative_real,
    readonly_quantity,
)
from ..geodesic import IntegrationError
from ..geodesic.orbit import (
    _BOUND_KEYS,
    _FAMILIES,
    _classify_inputs,
    _initial_components_from_scales,
)


_RG_AU_PER_MILLION_SOLAR_MASS = (
    const.G * 1e6 * const.M_sun / const.c**2
).to_value(u.au)
_AU_PER_KPC = (1 * u.kpc).to_value(u.au)
_C_AU_PER_YEAR = const.c.to_value(u.au / u.yr)
_C_KM_PER_SECOND = const.c.to_value(u.km / u.s)
_ARCSEC_TO_RAD = (1 * u.arcsec).to_value(u.rad)
_RAD_TO_ARCSEC = (1 * u.rad).to_value(u.arcsec)
_ORBIT_KEYS = frozenset(
    ("elements", "tau", "t", *_BOUND_KEYS, *(
        key for keys in _FAMILIES.values() for key in keys
    ))
)
_SOLVER_METHODS = {"dop853": 0, "radau": 1}
_TIME_KINDS = {"t_eval": 0, "tau_eval": 1, "t_obs": 2}


def _epochs(value: object, name: str) -> np.ndarray:
    """Copy and validate one series of observation epochs in Julian years."""
    if isinstance(value, u.Quantity):
        array = np.asarray(value.to_value(u.yr), dtype=np.float64)
    else:
        array = np.asarray(value, dtype=np.float64)
    if array.ndim != 1 or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be a finite one-dimensional time series")
    return immutable_array(array)


def _epoch(value: object, name: str) -> float:
    """Validate one epoch in Julian years, as a number or time quantity."""
    if isinstance(value, u.Quantity):
        if value.ndim != 0:
            raise ValueError(f"{name} must be a scalar")
        number = float(value.to_value(u.yr))
    elif isinstance(value, (bool, np.bool_)) or not isinstance(value, REAL_TYPES):
        raise TypeError(f"{name} must be a real number of Julian years")
    else:
        number = float(value)
    if not np.isfinite(number):
        raise ValueError(f"{name} must be finite")
    return number


def _positive(value: float, name: str) -> float:
    """Return a finite positive scalar or raise with its field name."""
    if value <= 0:
        raise ValueError(f"{name} must be positive and finite")
    return value


def _body_to_observer(spin_vector: np.ndarray) -> np.ndarray:
    """Construct a proper rotation with the body z axis along the spin."""
    magnitude = float(np.linalg.norm(spin_vector))
    if magnitude == 0.0:
        return np.eye(3, dtype=np.float64)
    axis = spin_vector / magnitude
    reference = np.array([1.0, 0.0, 0.0])
    if abs(float(axis[0])) > 0.9:
        reference = np.array([0.0, 1.0, 0.0])
    body_x = reference - np.dot(reference, axis) * axis
    body_x /= np.linalg.norm(body_x)
    body_y = np.cross(axis, body_x)
    return np.ascontiguousarray(np.column_stack((body_x, body_y, axis)))


def _spin_frame(value: object) -> tuple[float, np.ndarray]:
    """Validate one sky-frame Kerr spin and construct its body rotation."""
    vector = np.asarray(value, dtype=np.float64)
    if vector.shape != (3,) or not np.all(np.isfinite(vector)):
        raise ValueError("spin_vector must contain three finite components")
    spin = float(np.linalg.norm(vector))
    if not np.isfinite(spin) or spin > 1.0:
        raise ValueError("spin_vector magnitude must be at most one")
    return spin, immutable_array(_body_to_observer(vector))


def _mapping(value: object, name: str, allowed: frozenset[str]) -> dict[str, object]:
    """Copy a mapping and reject unknown keys before native evaluation."""
    if not isinstance(value, Mapping):
        raise TypeError(f"{name} must be a mapping")
    result = dict(value)
    unknown = result.keys() - allowed
    if unknown:
        raise ValueError(f"{name} has unknown keys: {', '.join(sorted(unknown))}")
    return result


def _scalar_physical(value: object, unit: u.UnitBase, name: str) -> u.Quantity:
    """Require one finite Astropy scalar with a compatible physical unit."""
    quantity = readonly_quantity(value, unit, name, ndim=(0,))
    if not np.isfinite(quantity.to_value(unit)):
        raise ValueError(f"{name} must be finite")
    return quantity


def _general_spin(kerr: Mapping[str, object]) -> tuple[float, np.ndarray]:
    """Validate mass-independent spin and observer-frame axis direction."""
    spin = dimensionless_scalar(kerr["spin"], "spin")
    if not np.isfinite(spin) or not 0 <= spin <= 1:
        raise ValueError("spin must be finite and in [0, 1]")
    axis = np.asarray(kerr["vec"], dtype=np.float64)
    if axis.shape != (3,) or not np.all(np.isfinite(axis)):
        raise ValueError("vec must contain three finite components")
    norm = float(np.linalg.norm(axis))
    if not np.isfinite(norm) or (norm == 0 and spin != 0):
        raise ValueError("vec must have nonzero finite magnitude when spin is nonzero")
    if norm == 0:
        return spin, immutable_array(np.eye(3, dtype=np.float64))
    return spin, immutable_array(_body_to_observer(axis / norm))


class KerrMcmcModel:
    """Predict Kerr astrometry and redshift without public orbit objects.

    Parameters
    ----------
    astrometry_times, spectroscopy_times : array-like or astropy.units.Quantity, optional
        Legacy observation epochs in Julian years. Omit both for the general
        :meth:`get_ra_dec_vr` API, where epochs are supplied per call.
    reference_epoch : float or astropy.units.Quantity, optional
        Required only with legacy observation epochs; defines the zero point
        of the linear astrometric frame drift.
    spin_vector : array-like, optional
        Three finite components of the dimensionless Kerr spin vector in
        observer axes ``(Dec, RA, away)``, default ``(0.0, 0.0, 0.0)``.
        Its magnitude must not exceed one. This is the default for
        evaluations; individual proposals may override it. The zero vector
        selects Schwarzschild.

    Raises
    ------
    TypeError
        If only one of ``astrometry_times`` and ``spectroscopy_times`` is
        given, ``reference_epoch`` is missing with epochs or given without
        them, or ``reference_epoch`` is not a real number or time quantity.
    ValueError
        If an epoch series is not one-dimensional and finite, both series
        are empty, ``reference_epoch`` is not a finite scalar, or
        ``spin_vector`` does not have three finite components with
        magnitude at most one.
    astropy.units.UnitConversionError
        If an epoch quantity is not convertible to years.

    Notes
    -----
    The general API accepts physical Kerr, orbit, observation-time and distance
    inputs per call. :meth:`set_solver` optionally supplies its solver defaults.
    The legacy ``__call__`` API requires :meth:`set_solver` and takes ``params``
    as the numeric
    vector ``(D, M, t_p, a, e, inc, Omega, omega, xS0, yS0, vxS0, vyS0,
    v_LSR)``. Units are kpc, million solar masses, Julian years, arcsec,
    dimensionless, radians (three angles), arcsec (two offsets), arcsec/year
    (two drifts), and km/s, respectively. This model uses the periapsis state
    as an osculating Kepler convention; ``a`` and ``e`` are not exact Kerr
    turning-point parameters. The small-angle and straight-line light model
    includes Rømer delay but no Shapiro delay or lensing.
    """

    __slots__ = (
        "_astrometry_times",
        "_spectroscopy_times",
        "_sorted_times",
        "_sort_order",
        "_spin",
        "_rotation",
        "_reference_epoch",
        "_solver",
    )

    def __init__(
        self,
        astrometry_times: object | None = None,
        spectroscopy_times: object | None = None,
        *,
        reference_epoch: float | None = None,
        spin_vector: object = (0.0, 0.0, 0.0),
    ) -> None:
        """Initialize general mode or the fixed-epoch legacy evaluation mode."""
        if astrometry_times is None and spectroscopy_times is None:
            if reference_epoch is not None:
                raise TypeError("reference_epoch requires observation epochs")
            spin, rotation = _spin_frame(spin_vector)
            self._astrometry_times = None
            self._spectroscopy_times = None
            self._sorted_times = None
            self._sort_order = None
            self._spin = spin
            self._rotation = rotation
            self._reference_epoch = None
            self._solver = None
            return
        if astrometry_times is None or spectroscopy_times is None:
            raise TypeError("both astrometry_times and spectroscopy_times are required")
        if reference_epoch is None:
            raise TypeError("reference_epoch is required with observation epochs")
        astrometry = _epochs(astrometry_times, "astrometry_times")
        spectroscopy = _epochs(spectroscopy_times, "spectroscopy_times")
        if astrometry.size + spectroscopy.size == 0:
            raise ValueError("at least one observation epoch is required")
        spin, rotation = _spin_frame(spin_vector)
        reference = _epoch(reference_epoch, "reference_epoch")

        joined = np.concatenate((astrometry, spectroscopy))
        order = immutable_array(np.argsort(joined, kind="stable"))
        sorted_times = immutable_array(joined[order])
        self._astrometry_times = astrometry
        self._spectroscopy_times = spectroscopy
        self._sorted_times = sorted_times
        self._sort_order = order
        self._spin = spin
        self._rotation = rotation
        self._reference_epoch = reference
        self._solver: tuple[int, float, float, int] | None = None

    def set_solver(
        self,
        *,
        method: str,
        rtol: float,
        atol: float,
        max_steps: int = 100000,
    ) -> None:
        """Set the native dop853 or radau solver and scalar tolerances.

        Parameters
        ----------
        method : {'dop853', 'radau'}
            Native integration method.
        rtol, atol : float
            Positive scalar relative and absolute tolerances in the internal
            geometric state. Their observational error must be checked for a
            particular fit.
        max_steps : int, optional
            Maximum accepted and rejected steps per endpoint call, default
            ``100000``.

        Raises
        ------
        TypeError
            If a tolerance is not a real number or ``max_steps`` is not an
            integer.
        ValueError
            If ``method`` is unsupported or a value is not positive and finite.
        """
        if method not in _SOLVER_METHODS:
            raise ValueError("method must be 'dop853' or 'radau'")
        rtol_value = _positive(non_negative_real(rtol, "rtol"), "rtol")
        atol_value = _positive(non_negative_real(atol, "atol"), "atol")
        steps = non_negative_integer(max_steps, "max_steps")
        if steps == 0:
            raise ValueError("max_steps must be a positive integer")
        self._solver = (_SOLVER_METHODS[method], rtol_value, atol_value, steps)

    def get_ra_dec_vr(
        self,
        *,
        kerr: Mapping[str, object],
        orbit: Mapping[str, object],
        sol: Mapping[str, object],
        distance: u.Quantity,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Evaluate RA, Dec and redshift velocity for one Kerr proposal.

        Parameters
        ----------
        kerr : mapping
            ``mass`` is a positive mass quantity, ``spin`` a dimensionless
            scalar in ``[0, 1]``, and ``vec`` a three-component spin direction
            in observer axes ``(Dec, RA, away)``. ``vec`` is normalized here;
            zero is allowed only for zero spin.
        orbit : mapping
            Initial conditions accepted by :meth:`relatipy.metrics.Kerr.orbit`,
            including optional physical ``t`` and ``tau`` origin quantities.
            Elements, Cartesian and spherical entries use the observer frame;
            Boyer--Lindquist and bound entries use the spin-aligned frame.
        sol : mapping
            Exactly one nonempty one-dimensional quantity: ``t_eval``,
            ``tau_eval`` or ``t_obs``. Times are absolute physical epochs and
            may be unordered. Optional ``method`` (``'dop853'`` or
            ``'radau'``), ``rtol``, ``atol`` and ``max_steps`` configure the
            native integrator.
        distance : astropy.units.Quantity
            Positive observer distance from the central body.

        Returns
        -------
        alpha, delta : numpy.ndarray
            Read-only angular offsets in arcsec, in the input time order.
        v_los : numpy.ndarray
            Read-only redshift-equivalent velocity in km/s, in that order.

        Raises
        ------
        TypeError, ValueError, astropy.units.UnitConversionError
            If an input mapping, unit, shape or value is invalid.
        relatipy.geodesic.IntegrationError
            If native integration or observed-time matching fails.
        """
        kerr_values = _mapping(kerr, "kerr", frozenset(("mass", "spin", "vec")))
        missing = {"mass", "spin", "vec"} - kerr_values.keys()
        if missing:
            raise ValueError(f"kerr requires {', '.join(sorted(missing))}")
        mass = _scalar_physical(kerr_values["mass"], u.kg, "mass")
        if mass <= 0 * u.kg:
            raise ValueError("mass must be positive")
        spin, rotation = _general_spin(kerr_values)
        observer_distance = _scalar_physical(distance, u.m, "distance")
        if observer_distance <= 0 * u.m:
            raise ValueError("distance must be positive")
        length_scale = (const.G * mass / const.c**2).to(u.m)
        time_scale = (length_scale / const.c).to(u.s)
        if (
            not np.isfinite(length_scale.to_value(u.m))
            or not np.isfinite(time_scale.to_value(u.s))
            or length_scale <= 0 * u.m or time_scale <= 0 * u.s
        ):
            raise ValueError("mass produces nonrepresentable geometric scales")

        provided_orbit = _mapping(orbit, "orbit", _ORBIT_KEYS)
        values = {key: None for key in _ORBIT_KEYS}
        values.update(provided_orbit)
        values["t"] = 0 * u.s if values["t"] is None else values["t"]
        values["tau"] = 0 * u.s if values["tau"] is None else values["tau"]
        t_origin = _scalar_physical(values["t"], u.s, "t")
        tau_origin = _scalar_physical(values["tau"], u.s, "tau")
        family = _classify_inputs(values)
        components, _ = _initial_components_from_scales(
            length_scale, time_scale, family, values, 0 * u.s
        )
        components = np.ascontiguousarray(components)

        supplied_sol = _mapping(
            sol, "sol", frozenset((*_TIME_KINDS, "method", "rtol", "atol", "max_steps"))
        )
        supplied_times = [name for name in _TIME_KINDS if name in supplied_sol]
        if len(supplied_times) != 1:
            raise ValueError("sol requires exactly one of t_eval, tau_eval or t_obs")
        time_name = supplied_times[0]
        raw_times = readonly_quantity(
            supplied_sol[time_name], u.s, time_name, ndim=(1,)
        )
        if raw_times.size == 0:
            raise ValueError(f"{time_name} must be nonempty")
        if not np.all(np.isfinite(raw_times.to_value(u.s))):
            raise ValueError(f"{time_name} must contain finite times")
        origin = tau_origin if time_name == "tau_eval" else t_origin
        normalized = np.asarray(
            ((raw_times - origin) / time_scale).to_value(u.one), dtype=np.float64
        )
        if not np.all(np.isfinite(normalized)):
            raise ValueError(f"{time_name} cannot be represented in geometric units")
        order = np.argsort(normalized, kind="stable")
        sorted_times = np.ascontiguousarray(normalized[order])

        if self._solver is None:
            default_method, default_rtol, default_atol, default_steps = (
                0, 1e-9, 1e-12, 100000
            )
        else:
            default_method, default_rtol, default_atol, default_steps = self._solver
        default_name = next(
            name for name, code in _SOLVER_METHODS.items() if code == default_method
        )
        method_name = supplied_sol.get("method", default_name)
        if method_name not in _SOLVER_METHODS:
            raise ValueError("method must be 'dop853' or 'radau'")
        rtol = _positive(
            non_negative_real(supplied_sol.get("rtol", default_rtol), "rtol"), "rtol"
        )
        atol = _positive(
            non_negative_real(supplied_sol.get("atol", default_atol), "atol"), "atol"
        )
        max_steps = non_negative_integer(
            supplied_sol.get("max_steps", default_steps), "max_steps"
        )
        if max_steps == 0:
            raise ValueError("max_steps must be a positive integer")

        try:
            from .. import _mcmc_core
        except ImportError as exc:
            raise ImportError(
                "native Kerr MCMC extension is unavailable; build relatipy first"
            ) from exc
        try:
            geometric = _mcmc_core.evaluate_general(
                spin, rotation, family, components, sorted_times, 0.0,
                _TIME_KINDS[time_name],
                float((length_scale / observer_distance).to_value(u.one)),
                _SOLVER_METHODS[method_name], rtol, atol, max_steps,
            )
        except RuntimeError as exc:
            raise IntegrationError(str(exc)) from exc
        geometric = np.asarray(geometric, dtype=np.float64)
        if geometric.shape != (raw_times.size, 3) or not np.all(np.isfinite(geometric)):
            raise IntegrationError("native Kerr MCMC evaluation returned invalid observables")
        unsorted = np.empty_like(geometric)
        unsorted[order] = geometric
        alpha = immutable_array(unsorted[:, 0] * _RAD_TO_ARCSEC)
        delta = immutable_array(unsorted[:, 1] * _RAD_TO_ARCSEC)
        v_los = immutable_array(unsorted[:, 2] * _C_KM_PER_SECOND)
        return alpha, delta, v_los

    def __call__(
        self,
        params: object,
        *,
        spin_vector: object | None = None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return ``(alpha, delta, v_los)`` for one numeric MCMC proposal.

        Parameters
        ----------
        params : array-like
            Thirteen finite orbital and observation parameters in the order
            and units specified by the class notes.
        spin_vector : array-like, optional
            Three-component spin for this proposal in observer axes
            ``(Dec, RA, away)``. Its magnitude is the dimensionless Kerr spin
            in ``[0, 1]``. If omitted, use the spin supplied when the model
            was constructed. The model's default remains unchanged after an
            override.

        Returns
        -------
        alpha, delta : numpy.ndarray
            Read-only astrometric offsets in arcsec at ``astrometry_times``.
        v_los : numpy.ndarray
            Read-only redshift-equivalent velocity in km/s at
            ``spectroscopy_times``.

        Raises
        ------
        ValueError
            If the model has no legacy epochs, :meth:`set_solver` was not
            called, ``params`` or ``spin_vector`` is invalid, or the native
            evaluator cannot construct a valid Kerr initial state from the
            proposal.
        relatipy.geodesic.IntegrationError
            If geodesic integration or observed-time matching fails.
        ImportError
            If the compiled Kerr MCMC extension is unavailable.
        """
        if self._astrometry_times is None:
            raise ValueError("legacy evaluation requires observation epochs at construction")
        if self._solver is None:
            raise ValueError("call set_solver before evaluating the model")
        values = np.asarray(params, dtype=np.float64)
        if values.shape != (13,) or not np.all(np.isfinite(values)):
            raise ValueError("params must contain thirteen finite numbers")
        (
            distance_kpc, mass_million_solar, periapsis_year,
            semi_major_arcsec, eccentricity, inclination,
            ascending_node, periapsis_argument, alpha_offset, delta_offset,
            alpha_drift, delta_drift, velocity_offset,
        ) = values
        if distance_kpc <= 0 or mass_million_solar <= 0 or semi_major_arcsec <= 0:
            raise ValueError("D, M and a must be positive")
        if eccentricity < 0 or eccentricity >= 1:
            raise ValueError("e must be in [0, 1)")
        if spin_vector is None:
            spin, rotation = self._spin, self._rotation
        else:
            spin, rotation = _spin_frame(spin_vector)

        gravitational_radius_au = (
            mass_million_solar * _RG_AU_PER_MILLION_SOLAR_MASS
        )
        distance_au = distance_kpc * _AU_PER_KPC
        semi_major_geometric = (
            semi_major_arcsec * _ARCSEC_TO_RAD
            * distance_au / gravitational_radius_au
        )
        time_scale_year = gravitational_radius_au / _C_AU_PER_YEAR
        arrival_times = np.ascontiguousarray(
            (self._sorted_times - periapsis_year) / time_scale_year
        )
        angular_scale = gravitational_radius_au / distance_au
        reference_time = (self._reference_epoch - periapsis_year) / time_scale_year
        method, rtol, atol, max_steps = self._solver
        try:
            from .. import _mcmc_core
        except ImportError as exc:
            raise ImportError(
                "native Kerr MCMC extension is unavailable; build relatipy first"
            ) from exc
        try:
            geometric = _mcmc_core.evaluate(
                spin,
                rotation,
                semi_major_geometric,
                eccentricity,
                inclination,
                ascending_node,
                periapsis_argument,
                arrival_times,
                angular_scale,
                alpha_offset * _ARCSEC_TO_RAD,
                delta_offset * _ARCSEC_TO_RAD,
                alpha_drift * _ARCSEC_TO_RAD * time_scale_year,
                delta_drift * _ARCSEC_TO_RAD * time_scale_year,
                velocity_offset / _C_KM_PER_SECOND,
                reference_time,
                method,
                rtol,
                atol,
                max_steps,
            )
        except RuntimeError as exc:
            raise IntegrationError(str(exc)) from exc

        unsorted = np.empty_like(geometric)
        unsorted[self._sort_order] = geometric
        astrometry_count = self._astrometry_times.size
        alpha = unsorted[:astrometry_count, 0] * _RAD_TO_ARCSEC
        delta = unsorted[:astrometry_count, 1] * _RAD_TO_ARCSEC
        v_los = unsorted[astrometry_count:, 2] * _C_KM_PER_SECOND
        for result in (alpha, delta, v_los):
            result.flags.writeable = False
        return alpha, delta, v_los
