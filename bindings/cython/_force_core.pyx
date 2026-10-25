# cython: language_level=3, boundscheck=False, wraparound=False
"""Private binding of the forced Kerr MCMC evaluators and the thin disk.

This module is an adapter for :mod:`relatipy.forces` and the perturbed
observable models; it is not public API. It contains no Kerr or force
formulas and passes no Python callback to an integrator: each function
validates arrays, calls one native entry point with the GIL released, and
copies results into NumPy-owned arrays.

A native force is a pair of integers, ``callback`` (the address of an
``rp_external_acceleration_fn``) and ``model`` (the address of its read-only
model data, or zero). ``callback == 0`` selects the inactive force, which is
bitwise the homogeneous Kerr path. The caller keeps the function and the
model memory alive, and unchanged, for the duration of each call.

Error mapping mirrors ``_mcmc_core``: invalid arguments and native rejection
of an input state raise ``ValueError``; failed allocations raise
``MemoryError``. The observable evaluators raise ``RuntimeError`` for
integration failure, which includes a failing force callback, and for
unresolved arrival times; :func:`acceleration` raises ``RuntimeError`` when
the callback reports failure. All values use normalized geometric units with
M = 1.
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

cdef extern from "geodesic/integrators/integrator.h":
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


cdef rp_external_force _force(uintptr_t callback, uintptr_t model) noexcept:
    """Build the native force pair; a zero callback is the inactive force."""
    cdef rp_external_force force
    force.evaluate = <rp_external_acceleration_fn><void *>callback
    force.model = <const void *>model
    return force


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


cdef extern from "force/models/axisymmetric_table.h":
    ctypedef enum rp_axisymmetric_status:
        RP_AXISYMMETRIC_STATUS_OK
        RP_AXISYMMETRIC_STATUS_ALLOCATION_FAILURE
    ctypedef struct rp_axisymmetric_table:
        size_t bytes
    ctypedef struct rp_axisymmetric_model:
        const rp_axisymmetric_table *table
        double mass
        double length_scale
        double inverse_length_scale
        double rotation[3][3]
    rp_axisymmetric_status rp_axisymmetric_model_init(
        rp_axisymmetric_model *model, const rp_axisymmetric_table *table,
        double mass, double length_scale, const double normal[3]
    ) nogil

cdef extern from "force/models/thin_disk.h":
    ctypedef struct rp_thin_disk_params:
        double alpha
        double inner_ratio
        double softening
    void rp_thin_disk_default_params(rp_thin_disk_params *params) nogil
    rp_axisymmetric_status rp_thin_disk_table_build(
        const rp_thin_disk_params *params, rp_axisymmetric_table **table
    ) nogil
    void rp_thin_disk_table_free(rp_axisymmetric_table *table) nogil
    int rp_thin_disk_acceleration(
        const double *position, const double *velocity,
        const void *model, double *acceleration
    ) nogil


cdef class AxisymmetricTable:
    """Owner of one built thin-disk table (read-only, shareable).

    Build with :meth:`thin_disk`; the GIL is released
    while the table is built. The table is freed when this object is
    collected, so every placed model keeps a reference to it.
    """

    cdef rp_axisymmetric_table *table
    cdef readonly str kind

    def __cinit__(self):
        self.table = NULL
        self.kind = ""

    def __dealloc__(self):
        if self.table != NULL:
            rp_thin_disk_table_free(self.table)
            self.table = NULL

    @staticmethod
    def thin_disk(double alpha, double inner_ratio, double softening):
        """Build the unit thin-disk table (``R_out = 1``, unit mass)."""
        cdef AxisymmetricTable owner = AxisymmetricTable()
        cdef rp_thin_disk_params params
        cdef rp_axisymmetric_status status
        rp_thin_disk_default_params(&params)
        params.alpha = alpha
        params.inner_ratio = inner_ratio
        params.softening = softening
        with nogil:
            status = rp_thin_disk_table_build(&params, &owner.table)
        owner.kind = "thin_disk"
        _check_table(status, owner.table)
        return owner

    @property
    def nbytes(self):
        """Heap bytes owned by the table."""
        return int(self.table.bytes) if self.table != NULL else 0

    @property
    def callback(self):
        """Address of the thin-disk acceleration callback."""
        return <uintptr_t><void *>rp_thin_disk_acceleration


cdef int _check_table(rp_axisymmetric_status status, rp_axisymmetric_table *table) except -1:
    if status == RP_AXISYMMETRIC_STATUS_ALLOCATION_FAILURE:
        raise MemoryError("native force table allocation failed")
    if status != RP_AXISYMMETRIC_STATUS_OK or table == NULL:
        raise RuntimeError(f"native force table build failed (status {<int>status})")
    return 0


cdef class AxisymmetricModel:
    """One placed thin disk: borrowed table, mass, scale, axis.

    :attr:`address` is the ``rp_axisymmetric_model`` passed to the callback;
    :meth:`set` rewrites it in place.
    """

    cdef rp_axisymmetric_model model
    cdef readonly AxisymmetricTable table

    def __cinit__(self, AxisymmetricTable table not None):
        if table.table == NULL:
            raise ValueError("table is not built")
        self.table = table
        memset(&self.model, 0, sizeof(self.model))
        self.model.table = table.table
        self.model.length_scale = 1.0
        self.model.inverse_length_scale = 1.0
        self.model.rotation[0][0] = 1.0
        self.model.rotation[1][1] = 1.0
        self.model.rotation[2][2] = 1.0

    @property
    def address(self):
        """Address of the native model data."""
        return <uintptr_t><void *>&self.model

    def set(self, double mass, double length_scale, normal):
        """Place the model: geometric mass and length, unit black-hole axis."""
        cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] axis
        cdef rp_axisymmetric_status status
        axis = np.ascontiguousarray(normal, dtype=np.float64)
        if axis.ndim != 1 or axis.shape[0] != 3:
            raise ValueError("normal must have shape (3,)")
        status = rp_axisymmetric_model_init(
            &self.model, self.table.table, mass, length_scale, &axis[0]
        )
        if status != RP_AXISYMMETRIC_STATUS_OK:
            raise ValueError(
                f"invalid placed force model (native status {<int>status})"
            )

    @property
    def parameters(self):
        """``(mass, length_scale, n_x, n_y, n_z)`` of the placed model."""
        return np.array([
            self.model.mass, self.model.length_scale,
            self.model.rotation[2][0], self.model.rotation[2][1],
            self.model.rotation[2][2],
        ], dtype=np.float64)
