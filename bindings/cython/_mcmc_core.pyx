# cython: language_level=3, boundscheck=False, wraparound=False
"""Private one-call binding for the native Kerr MCMC evaluator.

This module is an adapter for :class:`relatipy.observables.KerrMcmcModel`
and is not public API. It contains no Kerr formulas and passes no Python
callback to the integrator: each call validates shapes, releases the GIL for
one native evaluation, and returns a new NumPy array of shape ``(n, 3)`` with
``(alpha_rad, delta_rad, v_los_over_c)`` per requested time.

Error mapping: invalid arguments, invalid model input and a proposal that
does not define a valid Kerr initial state raise ``ValueError``; native
integration failure, unresolved observed arrival time and unknown native
status raise ``RuntimeError``, which the Python frontend converts to
``IntegrationError``. ``evaluate`` accepts zero requested times and returns
an empty array; ``evaluate_general`` rejects them with ``ValueError``.
"""

import numpy as np
cimport numpy as cnp

from libc.stddef cimport size_t

cdef extern from "geodesic/integrators/integrator.h":
    ctypedef enum rp_integrator_method:
        RP_INTEGRATOR_METHOD_DOP853
        RP_INTEGRATOR_METHOD_RADAU

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

    rp_kerr_mcmc_status rp_kerr_mcmc_evaluate(
        const rp_integrator_config *solver,
        double spin,
        const double rotation[3][3],
        double semi_major_axis,
        double eccentricity,
        double inclination,
        double ascending_node,
        double periapsis_argument,
        const double *arrival_times,
        size_t sample_count,
        const rp_kerr_mcmc_observation *observation,
        double *output,
        rp_integrator_stats *statistics
    ) nogil

    rp_kerr_mcmc_status rp_kerr_mcmc_evaluate_state(
        const rp_integrator_config *solver,
        double spin,
        const double rotation[3][3],
        const double initial[8],
        double tau_initial,
        rp_kerr_mcmc_time_kind time_kind,
        const double *evaluation_times,
        size_t sample_count,
        double angular_scale,
        double *output,
        rp_integrator_stats *statistics
    ) nogil

cdef extern from "relatipy/kerr_geometry.h":
    ctypedef enum rp_kerr_status:
        RP_KERR_STATUS_OK

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


def evaluate(
    double spin,
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
    """Evaluate angular observables and normalized v_los in native code."""
    cdef rp_integrator_config solver
    cdef rp_integrator_stats statistics
    cdef rp_kerr_mcmc_observation observation
    cdef rp_kerr_mcmc_status status
    cdef Py_ssize_t count = arrival_times.shape[0]
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="c"] output

    if rotation.shape[0] != 3 or rotation.shape[1] != 3:
        raise ValueError("rotation must have shape (3, 3)")
    if method not in (0, 1):
        raise ValueError("method must be 0 (dop853) or 1 (radau)")
    if maximum_steps == 0:
        raise ValueError("maximum_steps must be positive")

    output = np.empty((count, 3), dtype=np.float64)
    if count == 0:
        return output

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
    observation.angular_scale = angular_scale
    observation.alpha_offset = alpha_offset
    observation.delta_offset = delta_offset
    observation.alpha_drift = alpha_drift
    observation.delta_drift = delta_drift
    observation.velocity_offset = velocity_offset
    observation.reference_time = reference_time

    with nogil:
        status = rp_kerr_mcmc_evaluate(
            &solver,
            spin,
            <const double (*)[3]>&rotation[0, 0],
            semi_major_axis,
            eccentricity,
            inclination,
            ascending_node,
            periapsis_argument,
            &arrival_times[0],
            <size_t>count,
            &observation,
            &output[0, 0],
            &statistics,
        )

    if status == RP_KERR_MCMC_INVALID_INPUT:
        raise ValueError("invalid Kerr MCMC model input")
    if status == RP_KERR_MCMC_INITIAL_STATE_FAILURE:
        raise ValueError("orbital elements do not define a valid Kerr initial state")
    if status == RP_KERR_MCMC_INTEGRATION_FAILURE:
        raise RuntimeError("Kerr geodesic integration failed")
    if status == RP_KERR_MCMC_ARRIVAL_FAILURE:
        raise RuntimeError("observed arrival time could not be resolved")
    if status != RP_KERR_MCMC_OK:
        raise RuntimeError("unknown Kerr MCMC backend status")
    return output


def evaluate_general(
    double spin,
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
    """Build one Kerr state and evaluate three observables in one native call."""
    cdef rp_integrator_config solver
    cdef rp_integrator_stats statistics
    cdef rp_kerr_status initial_status
    cdef rp_kerr_mcmc_status status
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
    if method not in (0, 1):
        raise ValueError("method must be 0 (dop853) or 1 (radau)")
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
            status = rp_kerr_mcmc_evaluate_state(
                &solver, spin, <const double (*)[3]>&rotation[0, 0],
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
        raise RuntimeError("Kerr geodesic integration failed")
    if status == RP_KERR_MCMC_ARRIVAL_FAILURE:
        raise RuntimeError("observed arrival time could not be resolved")
    if status != RP_KERR_MCMC_OK:
        raise RuntimeError("unknown Kerr observable backend status")
    return output
