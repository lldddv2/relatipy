"""Evaluate stellar observables of a Kerr orbit perturbed by an external force.

:class:`PerturbedKerrMcmcModel` mirrors :class:`~relatipy.observables.KerrMcmcModel`
(same initial states, light model, units and errors) and integrates the
orbit through the forced native layer with one external force.
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
from astropy import constants as const
from astropy import units as u

from .._validation import (
    immutable_array,
    non_negative_integer,
    non_negative_real,
    readonly_quantity,
)
from ..forces import ExternalForce
from ..forces._bound import BoundForce
from ..geodesic import IntegrationError
from ..geodesic.orbit import _classify_inputs, _initial_components_from_scales
from .kerr_mcmc import (
    _ARCSEC_TO_RAD,
    _AU_PER_KPC,
    _C_AU_PER_YEAR,
    _C_KM_PER_SECOND,
    _ORBIT_KEYS,
    _RAD_TO_ARCSEC,
    _RG_AU_PER_MILLION_SOLAR_MASS,
    _SOLVER_METHODS,
    _TIME_KINDS,
    KerrMcmcModel,
    _general_spin,
    _mapping,
    _positive,
    _scalar_physical,
    _spin_frame,
)

_KG_PER_MILLION_SOLAR_MASS = float((1e6 * u.M_sun).to_value(u.kg))
_LEGACY_SIZE = 13


def _force_core():
    """Import the forced native extension with the error of the other models."""
    try:
        from .. import _force_core
    except ImportError as exc:
        raise ImportError(
            "native forced Kerr extension is unavailable; build relatipy first"
        ) from exc
    return _force_core


def _require_force(force: object) -> ExternalForce:
    if not isinstance(force, ExternalForce):
        raise TypeError("force must be a relatipy.forces.ExternalForce")
    return force


class PerturbedKerrMcmcModel:
    """Predict astrometry and redshift of a Kerr orbit under an external force.

    Parameters
    ----------
    astrometry_times, spectroscopy_times : array-like or astropy.units.Quantity, optional
        Legacy observation epochs in Julian years, as for
        :class:`~relatipy.observables.KerrMcmcModel`.
    reference_epoch : float or astropy.units.Quantity, optional
        Zero point of the linear astrometric drift; required with epochs.
    spin_vector : array-like, optional
        Default dimensionless spin vector in observer axes
        ``(Dec, RA, away)``, as for ``KerrMcmcModel``.
    force : relatipy.forces.ExternalForce
        Default force of :meth:`get_ra_dec_vr`, and the force whose
        parameters extend the legacy vector. Directions of the force, such as
        a disk normal, are read in observer axes ``(Dec, RA, away)``.

    Raises
    ------
    TypeError
        If ``force`` is not an external force, or as ``KerrMcmcModel``.
    ValueError, astropy.units.UnitConversionError
        As ``KerrMcmcModel``.

    Notes
    -----
    Orbits, light travel, units, solver settings and errors are those of
    ``KerrMcmcModel``; only the equation of motion adds the force (see
    :class:`~relatipy.metrics.PerturbedKerr`). With a force that is
    identically zero, such as a thin disk of zero mass, the forced
    right-hand side equals that of ``KerrMcmcModel`` and the predictions
    agree with it to rounding error; the tests bound the difference by
    ``1e-10`` relative.

    The legacy ``__call__`` vector is the thirteen ``KerrMcmcModel``
    parameters followed by the force block, the force's
    :attr:`~relatipy.forces.ExternalForce.parameter_names` in its
    :attr:`~relatipy.forces.ExternalForce.parameter_units`. The force is
    bound once and its parameters are rewritten in place for each call with
    the proposal mass, so a call neither rebuilds tables nor recompiles.
    Instances are picklable; a model is not safe to call from several
    threads at once.
    """

    __slots__ = ("_base", "_force", "_bound")

    def __init__(
        self,
        astrometry_times: object | None = None,
        spectroscopy_times: object | None = None,
        *,
        reference_epoch: float | None = None,
        spin_vector: object = (0.0, 0.0, 0.0),
        force: ExternalForce,
    ) -> None:
        """Validate epochs, spin and force as ``KerrMcmcModel`` does."""
        self._force = _require_force(force)
        self._base = KerrMcmcModel(
            astrometry_times, spectroscopy_times,
            reference_epoch=reference_epoch, spin_vector=spin_vector,
        )
        self._bound: BoundForce | None = None

    def __getstate__(self) -> dict:
        """Return picklable state; the native force handle is rebuilt later."""
        return {"base": self._base, "force": self._force}

    def __setstate__(self, state: dict) -> None:
        """Restore a pickled model without a bound force."""
        self._base = state["base"]
        self._force = state["force"]
        self._bound = None

    @property
    def force(self) -> ExternalForce:
        """Default force and the layout of the legacy force block."""
        return self._force

    @property
    def n_force_parameters(self) -> int:
        """Number of force parameters appended to the legacy vector."""
        return len(self._force.parameter_names)

    def set_solver(
        self,
        *,
        method: str,
        rtol: float,
        atol: float,
        max_steps: int = 100000,
    ) -> None:
        """Set the native solver, as :meth:`KerrMcmcModel.set_solver`.

        Parameters
        ----------
        method : {'dop853', 'radau'}
            Native integration method; Radau uses a finite-difference
            Jacobian, as in ``KerrMcmcModel``.
        rtol, atol : float
            Positive scalar relative and absolute tolerances.
        max_steps : int, optional
            Maximum steps per endpoint call, default ``100000``.

        Raises
        ------
        TypeError, ValueError
            As :meth:`KerrMcmcModel.set_solver`.
        """
        self._base.set_solver(method=method, rtol=rtol, atol=atol, max_steps=max_steps)

    def _binding(self, force: ExternalForce) -> BoundForce:
        """Return a handle for ``force``, reusing the cached one when possible."""
        bound = self._bound
        if bound is None or bound.force is not force:
            bound = force._handle()
            self._bound = bound
        return bound

    def get_ra_dec_vr(
        self,
        *,
        kerr: Mapping[str, object],
        orbit: Mapping[str, object],
        sol: Mapping[str, object],
        distance: u.Quantity,
        force: ExternalForce | None = None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Evaluate RA, Dec and redshift velocity for one perturbed proposal.

        Parameters
        ----------
        kerr, orbit, sol, distance
            As :meth:`KerrMcmcModel.get_ra_dec_vr`.
        force : relatipy.forces.ExternalForce, optional
            Force of this evaluation; the model's force by default. It is
            bound to ``kerr["mass"]`` with its own parameter values, and its
            directions are in observer axes ``(Dec, RA, away)``.

        Returns
        -------
        alpha, delta : numpy.ndarray
            Read-only angular offsets in arcsec, in the input time order.
        v_los : numpy.ndarray
            Read-only redshift-equivalent velocity in km/s, in that order.

        Raises
        ------
        TypeError, ValueError, astropy.units.UnitConversionError
            If an input mapping, unit, shape or value is invalid, or
            ``force`` is not an external force.
        relatipy.geodesic.IntegrationError
            If native integration, the force callback or observed-time
            matching fails.
        ImportError, RuntimeError
            If a run-time force needs SymPy or a C compiler that is missing.
        """
        selected = self._force if force is None else _require_force(force)
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

        solver = self._base._solver
        if solver is None:
            default_method, default_rtol, default_atol, default_steps = (
                0, 1e-9, 1e-12, 100000
            )
        else:
            default_method, default_rtol, default_atol, default_steps = solver
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

        native = _force_core()
        bound = self._binding(selected)
        bound.update(
            float(mass.to_value(u.kg)), selected._parameter_vector(), rotation
        )
        try:
            geometric = native.evaluate_general(
                spin, bound.callback, bound.model, rotation, family, components,
                sorted_times, 0.0, _TIME_KINDS[time_name],
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

    def _check_legacy(self, params: object) -> np.ndarray:
        """Validate a legacy vector of thirteen values plus the force block."""
        base = self._base
        if base._astrometry_times is None:
            raise ValueError("legacy evaluation requires observation epochs at construction")
        if base._solver is None:
            raise ValueError("call set_solver before evaluating the model")
        size = _LEGACY_SIZE + self.n_force_parameters
        values = np.asarray(params, dtype=np.float64)
        if values.shape != (size,) or not np.all(np.isfinite(values)):
            raise ValueError(f"params must contain {size} finite numbers")
        if values[0] <= 0 or values[1] <= 0 or values[3] <= 0:
            raise ValueError("D, M and a must be positive")
        if values[4] < 0 or values[4] >= 1:
            raise ValueError("e must be in [0, 1)")
        self._force._check_vector(values[_LEGACY_SIZE:])
        return values

    def _evaluate_legacy(
        self,
        values: np.ndarray,
        spin: float,
        rotation: np.ndarray,
        bound: BoundForce,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Evaluate thirteen validated legacy values with an updated handle."""
        (
            distance_kpc, mass_million_solar, periapsis_year,
            semi_major_arcsec, eccentricity, inclination,
            ascending_node, periapsis_argument, alpha_offset, delta_offset,
            alpha_drift, delta_drift, velocity_offset,
        ) = values[:_LEGACY_SIZE]
        base = self._base
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
            (base._sorted_times - periapsis_year) / time_scale_year
        )
        angular_scale = gravitational_radius_au / distance_au
        reference_time = (base._reference_epoch - periapsis_year) / time_scale_year
        method, rtol, atol, max_steps = base._solver
        native = _force_core()
        try:
            geometric = native.evaluate(
                spin,
                bound.callback,
                bound.model,
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
        unsorted[base._sort_order] = geometric
        astrometry_count = base._astrometry_times.size
        alpha = unsorted[:astrometry_count, 0] * _RAD_TO_ARCSEC
        delta = unsorted[:astrometry_count, 1] * _RAD_TO_ARCSEC
        v_los = unsorted[astrometry_count:, 2] * _C_KM_PER_SECOND
        for result in (alpha, delta, v_los):
            result.flags.writeable = False
        return alpha, delta, v_los

    def __call__(
        self,
        params: object,
        *,
        spin_vector: object | None = None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return ``(alpha, delta, v_los)`` for one numeric proposal.

        Parameters
        ----------
        params : array-like
            The thirteen ``KerrMcmcModel`` parameters followed by the force
            block, in the order and units of the class notes.
        spin_vector : array-like, optional
            Spin for this proposal in observer axes, as for
            ``KerrMcmcModel``; the model default when omitted.

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
            called, ``params``, a force parameter or ``spin_vector`` is
            invalid, or the native evaluator rejects the initial state.
        relatipy.geodesic.IntegrationError
            If integration, the force callback or arrival matching fails.
        ImportError, RuntimeError
            If the extension is unavailable, or a run-time force needs SymPy
            or a C compiler that is missing.
        """
        values = self._check_legacy(params)
        if spin_vector is None:
            spin, rotation = self._base._spin, self._base._rotation
        else:
            spin, rotation = _spin_frame(spin_vector)
        bound = self._binding(self._force)
        bound.update(
            values[1] * _KG_PER_MILLION_SOLAR_MASS, values[_LEGACY_SIZE:], rotation
        )
        return self._evaluate_legacy(values, spin, rotation, bound)

