# cython: language_level=3, boundscheck=False, wraparound=False
"""Private binding of the forced Kerr layer and its built-in force models.

This module is an adapter for :mod:`relatipy.forces`,
:class:`relatipy.metrics.PerturbedKerr` and the perturbed observable models;
it is not public API. It contains no Kerr or force formulas and passes no
Python callback to an integrator: each function validates arrays, calls one
native entry point with the GIL released, and copies results into
NumPy-owned arrays.

A native force is a pair of integers, ``callback`` (the address of an
``rp_external_acceleration_fn``) and ``model`` (the address of its read-only
model data, or zero). ``callback == 0`` selects the inactive force, which is
bitwise the homogeneous Kerr path. The caller keeps the function and the
model memory alive, and unchanged, for the duration of each call.

Error mapping mirrors ``_core`` and ``_mcmc_core``: invalid arguments and
native rejection of an input state raise ``ValueError``; failed allocations
raise ``MemoryError``; sizes beyond Python limits raise ``OverflowError``.
``integrate_kerr`` returns the native status (-1, 0 or 1) with the samples
reached; a failing force callback is an RHS failure (status -1). The
observable evaluators raise ``RuntimeError`` for integration failure, which
includes a failing force callback, and for unresolved arrival times. Direct
force evaluations raise ``RuntimeError`` when the callback reports failure
(``RP_FORCE_STATUS_CALLBACK_FAILURE`` for the four-acceleration). All values
use normalized geometric units with M = 1.
"""

import numpy as np
cimport numpy as cnp

from libc.stddef cimport size_t
from libc.stdint cimport uintptr_t
from libc.string cimport memset


cdef extern from "relatipy/kerr_geometry.h":
    ctypedef enum rp_kerr_status:
        RP_KERR_STATUS_OK

cdef extern from "force/external_force.h":
    ctypedef int (*rp_external_acceleration_fn)(
        const double *position,
        const double *velocity,
        const void *model,
        double *acceleration
    ) nogil
    ctypedef struct rp_external_force:
        rp_external_acceleration_fn evaluate
        const void *model

cdef extern from "force/external.h":
    int RP_FORCE_STATUS_CALLBACK_FAILURE
    rp_kerr_status rp_force_bl_to_cartesian(
        double spin, const double coordinates[4],
        double position[3], double (*jacobian)[3]
    ) nogil
    int rp_force_four_acceleration(
        double spin, const rp_external_force *force,
        const double coordinates[4], const double four_velocity[4],
        double four_acceleration[4]
    ) nogil

cdef extern from "force/models/plummer.h":
    int RP_FORCE_PLUMMER_PARAMETERS
    int rp_force_plummer_acceleration(
        const double *position, const double *velocity,
        const void *model, double *acceleration
    ) nogil
    int rp_force_plummer_potential(
        const double *position, const void *model, double *potential
    ) nogil

cdef extern from "geodesic/integrators/integrator.h":
    ctypedef enum rp_integrator_status:
        RP_INTEGRATOR_STATUS_OK
        RP_INTEGRATOR_STATUS_OBSERVER_STOPPED
    ctypedef enum rp_integrator_method:
        RP_INTEGRATOR_METHOD_DOP853
        RP_INTEGRATOR_METHOD_RADAU
        RP_INTEGRATOR_METHOD_DP45
    ctypedef struct rp_integrator_config:
        rp_integrator_method method
        size_t dimension
        double relative_tolerance
        double absolute_tolerance
        double initial_step
        double maximum_step
        size_t maximum_steps
        const double *absolute_tolerances
        void *step_observer
        void *step_observer_context
        void *jacobian
        void *projector
    ctypedef struct rp_integrator_stats:
        size_t accepted_steps
        size_t rejected_steps
        size_t rhs_evaluations
        size_t jacobian_evaluations
        size_t linear_solves
        double final_independent_variable

cdef extern from "geodesic/solution/solve.h":
    ctypedef struct rp_kerr_trajectory:
        double *taus
        double *states
        size_t count
        size_t capacity
        rp_integrator_stats stats
        rp_integrator_status integrator_status
        int crossed_outer_horizon
        int allocation_failed
    void rp_kerr_trajectory_free(rp_kerr_trajectory *trajectory) nogil

cdef extern from "force/solution/solve.h":
    rp_integrator_status rp_kerr_forced_trajectory_solve(
        double spin, const rp_external_force *force,
        const double initial_state[8],
        double tau_initial, double tau_final, rp_integrator_method method,
        double rtol, double scalar_atol, const double *vector_atol,
        double first_step, double max_step, int store_steps,
        rp_kerr_trajectory *trajectory
    ) nogil

cdef extern from "geodesic/kerr_mcmc.h":
    ctypedef enum rp_kerr_mcmc_status:
        RP_KERR_MCMC_OK
        RP_KERR_MCMC_INVALID_INPUT
        RP_KERR_MCMC_INITIAL_STATE_FAILURE
        RP_KERR_MCMC_INTEGRATION_FAILURE
        RP_KERR_MCMC_ARRIVAL_FAILURE
    ctypedef enum rp_kerr_mcmc_time_kind:
        RP_KERR_MCMC_TIME_COORDINATE
        RP_KERR_MCMC_TIME_PROPER
        RP_KERR_MCMC_TIME_ARRIVAL
    ctypedef struct rp_kerr_mcmc_observation:
        double angular_scale
        double alpha_offset
        double delta_offset
        double alpha_drift
        double delta_drift
        double velocity_offset
        double reference_time

cdef extern from "force/kerr_mcmc.h":
    rp_kerr_mcmc_status rp_kerr_forced_mcmc_evaluate(
        const rp_integrator_config *solver, double spin,
        const rp_external_force *force, const double rotation[3][3],
        double semi_major_axis, double eccentricity, double inclination,
        double ascending_node, double periapsis_argument,
        const double *arrival_times, size_t sample_count,
        const rp_kerr_mcmc_observation *observation,
        double *output, rp_integrator_stats *statistics
    ) nogil
    rp_kerr_mcmc_status rp_kerr_forced_mcmc_evaluate_state(
        const rp_integrator_config *solver, double spin,
        const rp_external_force *force, const double rotation[3][3],
        const double initial[8], double tau_initial,
        rp_kerr_mcmc_time_kind time_kind, const double *sorted_times,
        size_t sample_count, double angular_scale,
        double *output, rp_integrator_stats *statistics
    ) nogil

cdef extern from "geodesic/initial/convert.h":
    rp_kerr_status rp_initial_elements_to_cartesian(
        const double elements[7], double cartesian[7]
    ) nogil
    rp_kerr_status rp_initial_spherical_to_cartesian(
        const double spherical[7], double cartesian[7]
    ) nogil
    rp_kerr_status rp_initial_observer_cartesian_to_canonical(
        double spin, const double rotation[3][3],
        const double observer[7], double canonical[8]
    ) nogil
    rp_kerr_status rp_initial_from_bl(
        double spin, const double bl[7], double canonical[8]
    ) nogil

cdef extern from "geodesic/initial/bound.h":
    rp_kerr_status rp_initial_from_bound(
        double spin, const double bound[7], double canonical[8]
    ) nogil


#: Status of a failing force callback in ``four_acceleration``.
CALLBACK_FAILURE = RP_FORCE_STATUS_CALLBACK_FAILURE
#: Number of doubles in the Plummer model data ``(mass, scale)``.
PLUMMER_PARAMETERS = RP_FORCE_PLUMMER_PARAMETERS


cdef rp_external_force _force(uintptr_t callback, uintptr_t model) noexcept:
    """Build the native force pair; a zero callback is the inactive force."""
    cdef rp_external_force force
    force.evaluate = <rp_external_acceleration_fn><void *>callback
    force.model = <const void *>model
    return force


def plummer_address():
    """Return the address of the built-in Plummer acceleration callback."""
    return <uintptr_t><void *>rp_force_plummer_acceleration


def plummer_potential(
    cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] model,
    cnp.ndarray[cnp.float64_t, ndim=2, mode="c"] positions,
):
    """Evaluate the Plummer potential at normalized Cartesian positions."""
    cdef Py_ssize_t count = positions.shape[0]
    cdef Py_ssize_t i
    cdef int status = 0
    cdef Py_ssize_t failed = -1
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] output
    if model.shape[0] != RP_FORCE_PLUMMER_PARAMETERS:
        raise ValueError("Plummer model data must contain (mass, scale)")
    if positions.shape[1] != 3:
        raise ValueError("positions must have shape (n, 3)")
    output = np.zeros(count, dtype=np.float64)
    with nogil:
        for i in range(count):
            status = rp_force_plummer_potential(
                &positions[i, 0], &model[0], &output[i]
            )
            if status != 0:
                failed = i
                break
    if failed >= 0:
        raise ValueError(f"Plummer potential failed at row {failed}")
    return output


def acceleration(
    uintptr_t callback,
    uintptr_t model,
    cnp.ndarray[cnp.float64_t, ndim=2, mode="c"] positions,
    cnp.ndarray[cnp.float64_t, ndim=2, mode="c"] velocities,
):
    """Call one force callback on normalized Cartesian rows.

    Returns the ``(n, 3)`` Newtonian accelerations. A nonzero callback value
    raises ``RuntimeError`` naming the first failing row.
    """
    cdef rp_external_acceleration_fn evaluate
    cdef const void *data = <const void *>model
    cdef Py_ssize_t count = positions.shape[0]
    cdef Py_ssize_t i
    cdef int status = 0
    cdef Py_ssize_t failed = -1
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="c"] output
    if callback == 0:
        raise ValueError("callback must be a nonzero function address")
    if positions.shape[1] != 3 or velocities.shape[1] != 3:
        raise ValueError("positions and velocities must have shape (n, 3)")
    if velocities.shape[0] != count:
        raise ValueError("positions and velocities must have the same rows")
    evaluate = <rp_external_acceleration_fn><void *>callback
    output = np.zeros((count, 3), dtype=np.float64)
    with nogil:
        for i in range(count):
            status = evaluate(
                &positions[i, 0], &velocities[i, 0], data, &output[i, 0]
            )
            if status != 0:
                failed = i
                break
    if failed >= 0:
        raise RuntimeError(
            f"external force callback failed at row {failed} (status {status})"
        )
    return output


def bl_to_cartesian(
    double spin,
    cnp.ndarray[cnp.float64_t, ndim=2, mode="c"] coordinates,
):
    """Map normalized ``(t, r, theta, phi)`` rows to the black-hole frame."""
    cdef Py_ssize_t count = coordinates.shape[0]
    cdef Py_ssize_t i
    cdef rp_kerr_status status = RP_KERR_STATUS_OK
    cdef Py_ssize_t failed = -1
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="c"] output
    if coordinates.shape[1] != 4:
        raise ValueError("coordinates must have shape (n, 4)")
    if not np.isfinite(spin) or spin < 0.0 or spin > 1.0:
        raise ValueError("spin must be finite and in [0, 1]")
    output = np.zeros((count, 3), dtype=np.float64)
    with nogil:
        for i in range(count):
            status = rp_force_bl_to_cartesian(
                spin, &coordinates[i, 0], &output[i, 0], NULL
            )
            if status != RP_KERR_STATUS_OK:
                failed = i
                break
    if failed >= 0:
        raise ValueError(
            f"Boyer-Lindquist row {failed} is invalid (native status {<int>status})"
        )
    return output


def four_acceleration(
    double spin,
    uintptr_t callback,
    uintptr_t model,
    coordinates,
    four_velocity,
):
    """Return the projected four-acceleration ``f^mu`` at one state."""
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] x
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] velocity
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] output
    cdef rp_external_force force = _force(callback, model)
    cdef int status
    x = np.ascontiguousarray(coordinates, dtype=np.float64)
    velocity = np.ascontiguousarray(four_velocity, dtype=np.float64)
    if x.ndim != 1 or x.shape[0] != 4 or velocity.ndim != 1 or velocity.shape[0] != 4:
        raise ValueError("coordinates and four_velocity must have shape (4,)")
    output = np.zeros(4, dtype=np.float64)
    with nogil:
        status = rp_force_four_acceleration(
            spin, &force, &x[0], &velocity[0], &output[0]
        )
    if status == RP_FORCE_STATUS_CALLBACK_FAILURE:
        raise RuntimeError("external force callback failed")
    if status != 0:
        raise ValueError(f"forced Kerr state is invalid (native status {status})")
    return output


cdef rp_integrator_method _trajectory_method(method) except *:
    if method == "radau":
        return RP_INTEGRATOR_METHOD_RADAU
    if method == "dop853":
        return RP_INTEGRATOR_METHOD_DOP853
    if method == "dp45":
        return RP_INTEGRATOR_METHOD_DP45
    if method == "projection_radau":
        raise ValueError(
            "projection_radau is not available with an external force: "
            "its projector restores Carter's constant, which a force does "
            "not conserve"
        )
    raise ValueError("method must be one of radau, dop853, dp45")


def integrate_kerr(
    double spin, uintptr_t callback, uintptr_t model, y0, double tau0,
    double tau1, method="radau", double rtol=1e-3, atol=1e-6,
    first_step=None, max_step=None, bint store_steps=True,
):
    """Integrate a canonical Kerr state under one force; mirrors ``_core``.

    Returns ``(taus, states, stats, status)`` exactly as
    ``relatipy._core.integrate_kerr``; ``stats`` also records whether the
    force was active.
    """
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] initial
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] vector_atol
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] taus
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="c"] states
    cdef rp_integrator_method selected_method
    cdef rp_integrator_status native_status
    cdef rp_kerr_trajectory trajectory
    cdef rp_external_force force = _force(callback, model)
    cdef const double *atol_pointer = NULL
    cdef double scalar_atol = 0.0
    cdef double initial_step = 0.0
    cdef double maximum_step = 0.0
    cdef Py_ssize_t i, j, count
    cdef int status
    selected_method = _trajectory_method(method)
    if not np.isfinite(spin) or spin < 0.0 or spin > 1.0:
        raise ValueError("spin must be finite and in [0, 1]")
    initial = np.ascontiguousarray(y0, dtype=np.float64)
    if initial.ndim != 1 or initial.shape[0] != 8 or not np.all(np.isfinite(initial)):
        raise ValueError("y0 must be a finite state of shape (8,)")
    if not np.isfinite(tau0) or not np.isfinite(tau1) or tau1 < tau0:
        raise ValueError("tau interval must be finite and nondecreasing")
    if not np.isfinite(rtol) or rtol < 0.0:
        raise ValueError("rtol must be finite and nonnegative")
    if np.ndim(atol) == 0:
        scalar_atol = float(atol)
        if not np.isfinite(scalar_atol) or scalar_atol < 0.0:
            raise ValueError("atol must be finite and nonnegative")
    else:
        vector_atol = np.ascontiguousarray(atol, dtype=np.float64)
        if vector_atol.ndim != 1 or vector_atol.shape[0] != 8:
            raise ValueError("atol must be scalar or have shape (8,)")
        if not np.all(np.isfinite(vector_atol)) or np.any(vector_atol < 0.0):
            raise ValueError("atol must be finite and nonnegative")
        atol_pointer = &vector_atol[0]
    if first_step is not None:
        initial_step = float(first_step)
        if not np.isfinite(initial_step) or initial_step <= 0.0:
            raise ValueError("first_step must be finite and positive")
    if max_step is not None:
        maximum_step = float(max_step)
        if maximum_step == np.inf:
            maximum_step = 0.0
        elif not np.isfinite(maximum_step) or maximum_step <= 0.0:
            raise ValueError("max_step must be positive")
    memset(&trajectory, 0, sizeof(trajectory))
    with nogil:
        native_status = rp_kerr_forced_trajectory_solve(
            spin, &force, &initial[0], tau0, tau1, selected_method,
            rtol, scalar_atol, atol_pointer, initial_step, maximum_step,
            <int>store_steps, &trajectory
        )
    try:
        if trajectory.allocation_failed:
            raise MemoryError("native trajectory allocation failed")
        if trajectory.count == 0:
            raise ValueError(
                f"native forced Kerr integration rejected input (status {<int>native_status})"
            )
        if trajectory.count > ((<size_t>-1) >> 1):
            raise OverflowError("trajectory exceeds Python array limits")
        count = <Py_ssize_t>trajectory.count
        taus = np.empty(count, dtype=np.float64)
        states = np.empty((count, 8), dtype=np.float64)
        for i in range(count):
            taus[i] = trajectory.taus[i]
            for j in range(8):
                states[i, j] = trajectory.states[i * 8 + j]
        stats = {
            "n_steps": int(trajectory.stats.accepted_steps),
            "nfev": int(trajectory.stats.rhs_evaluations),
            "rejected_steps": int(trajectory.stats.rejected_steps),
            "jacobian_evaluations": int(trajectory.stats.jacobian_evaluations),
            "linear_solves": int(trajectory.stats.linear_solves),
            "native_status": int(native_status),
            "crossed_outer_horizon": bool(trajectory.crossed_outer_horizon),
            "force_active": callback != 0,
        }
        if native_status == RP_INTEGRATOR_STATUS_OK:
            status = 0
        elif native_status == RP_INTEGRATOR_STATUS_OBSERVER_STOPPED and trajectory.crossed_outer_horizon:
            status = 1
        else:
            status = -1
        return taus, states, stats, status
    finally:
        with nogil:
            rp_kerr_trajectory_free(&trajectory)


cdef void _solver(
    rp_integrator_config *solver, int method, double rtol, double atol,
    size_t maximum_steps,
) noexcept:
    solver.method = <rp_integrator_method>method
    solver.dimension = 8
    solver.relative_tolerance = rtol
    solver.absolute_tolerance = atol
    solver.initial_step = 0.0
    solver.maximum_step = 0.0
    solver.maximum_steps = maximum_steps
    solver.absolute_tolerances = NULL
    solver.step_observer = NULL
    solver.step_observer_context = NULL
    solver.jacobian = NULL
    solver.projector = NULL


cdef int _check_mcmc_method(int method) except -1:
    if method == 3:
        raise ValueError("projection_radau is not available with an external force")
    if method not in (0, 1, 2):
        raise ValueError("method must be 0 (dop853), 1 (radau) or 2 (dp45)")
    return 0


def evaluate(
    double spin,
    uintptr_t callback,
    uintptr_t model,
    cnp.ndarray[cnp.float64_t, ndim=2, mode="c"] rotation,
    double semi_major_axis,
    double eccentricity,
    double inclination,
    double ascending_node,
    double periapsis_argument,
    cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] arrival_times,
    double angular_scale,
    double alpha_offset,
    double delta_offset,
    double alpha_drift,
    double delta_drift,
    double velocity_offset,
    double reference_time,
    int method,
    double rtol,
    double atol,
    size_t maximum_steps,
):
    """Forced counterpart of ``_mcmc_core.evaluate`` (fixed epochs)."""
    cdef rp_integrator_config solver
    cdef rp_integrator_stats statistics
    cdef rp_kerr_mcmc_observation observation
    cdef rp_kerr_mcmc_status status
    cdef rp_external_force force = _force(callback, model)
    cdef Py_ssize_t count = arrival_times.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="c"] output

    if rotation.shape[0] != 3 or rotation.shape[1] != 3:
        raise ValueError("rotation must have shape (3, 3)")
    _check_mcmc_method(method)
    if maximum_steps == 0:
        raise ValueError("maximum_steps must be positive")

    output = np.empty((count, 3), dtype=np.float64)
    if count == 0:
        return output

    _solver(&solver, method, rtol, atol, maximum_steps)
    observation.angular_scale = angular_scale
    observation.alpha_offset = alpha_offset
    observation.delta_offset = delta_offset
    observation.alpha_drift = alpha_drift
    observation.delta_drift = delta_drift
    observation.velocity_offset = velocity_offset
    observation.reference_time = reference_time

    with nogil:
        status = rp_kerr_forced_mcmc_evaluate(
            &solver, spin, &force,
            <const double (*)[3]>&rotation[0, 0],
            semi_major_axis, eccentricity, inclination,
            ascending_node, periapsis_argument,
            &arrival_times[0], <size_t>count, &observation,
            &output[0, 0], &statistics,
        )

    if status == RP_KERR_MCMC_INVALID_INPUT:
        raise ValueError("invalid Kerr MCMC model input")
    if status == RP_KERR_MCMC_INITIAL_STATE_FAILURE:
        raise ValueError("orbital elements do not define a valid Kerr initial state")
    if status == RP_KERR_MCMC_INTEGRATION_FAILURE:
        raise RuntimeError(
            "forced Kerr geodesic integration failed (numerical failure or "
            "external force callback failure)"
        )
    if status == RP_KERR_MCMC_ARRIVAL_FAILURE:
        raise RuntimeError("observed arrival time could not be resolved")
    if status != RP_KERR_MCMC_OK:
        raise RuntimeError("unknown forced Kerr MCMC backend status")
    return output


def evaluate_general(
    double spin,
    uintptr_t callback,
    uintptr_t model,
    cnp.ndarray[cnp.float64_t, ndim=2, mode="c"] rotation,
    family,
    cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] components,
    cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] times,
    double tau_initial,
    int time_kind,
    double angular_scale,
    int method,
    double rtol,
    double atol,
    size_t maximum_steps,
):
    """Forced counterpart of ``_mcmc_core.evaluate_general``."""
    cdef rp_integrator_config solver
    cdef rp_integrator_stats statistics
    cdef rp_kerr_status initial_status
    cdef rp_kerr_mcmc_status status = RP_KERR_MCMC_OK
    cdef rp_external_force force = _force(callback, model)
    cdef double observer[7]
    cdef double initial[8]
    cdef Py_ssize_t count = times.shape[0]
    cdef int family_code
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="c"] output

    if rotation.shape[0] != 3 or rotation.shape[1] != 3:
        raise ValueError("rotation must have shape (3, 3)")
    if components.shape[0] != 7 or not np.all(np.isfinite(components)):
        raise ValueError("components must contain seven finite values")
    if count == 0:
        raise ValueError("at least one evaluation time is required")
    _check_mcmc_method(method)
    if time_kind not in (0, 1, 2):
        raise ValueError("time_kind must be coordinate, proper, or arrival")
    if maximum_steps == 0:
        raise ValueError("maximum_steps must be positive")
    if family == "elements":
        family_code = 0
    elif family == "cartesian":
        family_code = 1
    elif family == "spherical":
        family_code = 2
    elif family == "bl":
        family_code = 3
    elif family == "bound":
        family_code = 4
    else:
        raise ValueError("unsupported orbit input family")

    _solver(&solver, method, rtol, atol, maximum_steps)
    output = np.empty((count, 3), dtype=np.float64)

    with nogil:
        if family_code == 0:
            initial_status = rp_initial_elements_to_cartesian(
                &components[0], observer
            )
            if initial_status == RP_KERR_STATUS_OK:
                initial_status = rp_initial_observer_cartesian_to_canonical(
                    spin, <const double (*)[3]>&rotation[0, 0],
                    observer, initial
                )
        elif family_code == 1:
            initial_status = rp_initial_observer_cartesian_to_canonical(
                spin, <const double (*)[3]>&rotation[0, 0],
                &components[0], initial
            )
        elif family_code == 2:
            initial_status = rp_initial_spherical_to_cartesian(
                &components[0], observer
            )
            if initial_status == RP_KERR_STATUS_OK:
                initial_status = rp_initial_observer_cartesian_to_canonical(
                    spin, <const double (*)[3]>&rotation[0, 0],
                    observer, initial
                )
        elif family_code == 3:
            initial_status = rp_initial_from_bl(
                spin, &components[0], initial
            )
        else:
            initial_status = rp_initial_from_bound(
                spin, &components[0], initial
            )

        if initial_status == RP_KERR_STATUS_OK:
            status = rp_kerr_forced_mcmc_evaluate_state(
                &solver, spin, &force, <const double (*)[3]>&rotation[0, 0],
                initial, tau_initial,
                <rp_kerr_mcmc_time_kind>time_kind,
                &times[0], <size_t>count, angular_scale,
                &output[0, 0], &statistics,
            )

    if initial_status != RP_KERR_STATUS_OK:
        raise ValueError(
            f"orbit does not define a valid Kerr initial state "
            f"(native status {<int>initial_status})"
        )
    if status == RP_KERR_MCMC_INVALID_INPUT:
        raise ValueError("invalid Kerr observable evaluation input")
    if status == RP_KERR_MCMC_INITIAL_STATE_FAILURE:
        raise ValueError("initial Kerr state is invalid")
    if status == RP_KERR_MCMC_INTEGRATION_FAILURE:
        raise RuntimeError(
            "forced Kerr geodesic integration failed (numerical failure or "
            "external force callback failure)"
        )
    if status == RP_KERR_MCMC_ARRIVAL_FAILURE:
        raise RuntimeError("observed arrival time could not be resolved")
    if status != RP_KERR_MCMC_OK:
        raise RuntimeError("unknown Kerr observable backend status")
    return output
