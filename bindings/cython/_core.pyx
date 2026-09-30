# cython: language_level=3, boundscheck=False, wraparound=False
"""Private thin binding for coarse-grained Kerr state reconstruction.

This module is an adapter between :mod:`relatipy.geodesic`/:mod:`relatipy.metrics`
and the native C core; it is not public API. It contains no Kerr formulas and
passes no Python callback to the integrator: each function validates array
shapes and finiteness, calls one native entry point with the GIL released, and
copies the result into NumPy-owned arrays. Native buffers allocated here are
freed before the function returns.

Error mapping: invalid arguments and native rejection of an input state raise
``ValueError``; failed allocations raise ``MemoryError``; sizes beyond Python
or native limits raise ``OverflowError``. Numerical integration failure and
terminal events are not raised here: ``integrate_kerr`` returns the native
status (-1, 0 or 1) with the samples reached, and the Python frontend maps it
to ``Solution.status`` or to ``IntegrationError``/``IntegrationTerminated``.
All values use normalized geometric units with M = 1.
"""

import numpy as np
cimport numpy as cnp

from libc.stddef cimport size_t
from libc.stdlib cimport malloc, free
from libc.string cimport memset

cdef extern from "geodesic/solution/reconstruct.h":
    ctypedef enum rp_kerr_status:
        RP_KERR_STATUS_OK

    rp_kerr_status rp_solution_reconstruct_batch(
        double spin,
        const double *cartesian,
        size_t count,
        double *reconstructed,
        rp_kerr_status *row_status
    ) nogil
    rp_kerr_status rp_solution_reconstruct_canonical_batch(
        double spin, const double *canonical, size_t count,
        double *reconstructed, double *cartesian,
        rp_kerr_status *row_status
    ) nogil
    ctypedef enum rp_solution_family:
        RP_SOLUTION_FAMILY_CARTESIAN
        RP_SOLUTION_FAMILY_SPHERICAL
        RP_SOLUTION_FAMILY_ELEMENTS
    rp_kerr_status rp_solution_reconstruct_canonical_family_batch(
        double spin, const double *canonical, size_t count,
        rp_solution_family family, double *output,
        rp_kerr_status *row_status
    ) nogil

cdef extern from "metric/kerr_properties.h":
    rp_kerr_status rp_kerr_characteristic_radii(
        double spin, double radii[6]
    ) nogil
    rp_kerr_status rp_kerr_ergosurface_radii(
        double spin, const double *theta, size_t count, double *radii
    ) nogil
    ctypedef enum rp_kerr_surface:
        RP_KERR_SURFACE_OUTER_HORIZON
        RP_KERR_SURFACE_ERGOSURFACE
    rp_kerr_status rp_kerr_surface_profile(
        double spin, rp_kerr_surface surface, const double *theta,
        size_t count, double *rho, double *z
    ) nogil


def kerr_properties(double spin):
    """Return normalized horizon, ISCO and photon radii from native Kerr."""
    cdef double radii[6]
    cdef rp_kerr_status status
    if not np.isfinite(spin) or spin < 0.0 or spin > 1.0:
        raise ValueError("spin must be finite and in [0, 1]")
    with nogil:
        status = rp_kerr_characteristic_radii(spin, radii)
    if status != RP_KERR_STATUS_OK:
        raise ValueError(f"native Kerr radii failed (status {<int>status})")
    return np.array([radii[i] for i in range(6)], dtype=np.float64)


def kerr_ergosurface(
    double spin,
    cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] theta,
):
    """Return normalized outer stationary-limit radii for polar angles."""
    cdef Py_ssize_t count = theta.shape[0]
    cdef rp_kerr_status status
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] output
    if not np.isfinite(spin) or spin < 0.0 or spin > 1.0:
        raise ValueError("spin must be finite and in [0, 1]")
    if not np.all(np.isfinite(theta)) or np.any(theta < 0.0) or np.any(theta > np.pi):
        raise ValueError("theta must be finite and in [0, pi]")
    output = np.empty(count, dtype=np.float64)
    if count == 0:
        return output
    with nogil:
        status = rp_kerr_ergosurface_radii(
            spin, &theta[0], <size_t>count, &output[0]
        )
    if status != RP_KERR_STATUS_OK:
        raise ValueError(f"native Kerr ergosurface failed (status {<int>status})")
    return output


def kerr_surface_profile(
    double spin,
    int surface,
    cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] theta,
):
    """Return normalized Cartesian ``(rho, z)`` of a Kerr reference surface.

    ``surface`` is 0 for the outer horizon and 1 for the outer ergosurface.
    """
    cdef Py_ssize_t count = theta.shape[0]
    cdef rp_kerr_status status
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] rho
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] z
    if not np.isfinite(spin) or spin < 0.0 or spin > 1.0:
        raise ValueError("spin must be finite and in [0, 1]")
    if surface != RP_KERR_SURFACE_OUTER_HORIZON and surface != RP_KERR_SURFACE_ERGOSURFACE:
        raise ValueError("surface must be 0 (outer horizon) or 1 (ergosurface)")
    if not np.all(np.isfinite(theta)) or np.any(theta < 0.0) or np.any(theta > np.pi):
        raise ValueError("theta must be finite and in [0, pi]")
    rho = np.empty(count, dtype=np.float64)
    z = np.empty(count, dtype=np.float64)
    if count == 0:
        return rho, z
    with nogil:
        status = rp_kerr_surface_profile(
            spin, <rp_kerr_surface>surface, &theta[0], <size_t>count,
            &rho[0], &z[0]
        )
    if status != RP_KERR_STATUS_OK:
        raise ValueError(f"native Kerr surface profile failed (status {<int>status})")
    return rho, z


def reconstruct_batch(
    double spin,
    cnp.ndarray[cnp.float64_t, ndim=2, mode="c"] cartesian,
):
    """Return normalized reconstructed rows and their native status codes."""
    cdef Py_ssize_t count = cartesian.shape[0]
    cdef Py_ssize_t i
    cdef rp_kerr_status result
    cdef rp_kerr_status *row_status = NULL
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="c"] output
    cdef cnp.ndarray[cnp.int32_t, ndim=1, mode="c"] status

    if cartesian.shape[1] != 7:
        raise ValueError("cartesian must have shape (n, 7)")
    if not np.isfinite(spin) or spin < 0.0 or spin > 1.0:
        raise ValueError("spin must be finite and in [0, 1]")
    if not np.all(np.isfinite(cartesian)):
        raise ValueError("cartesian must contain only finite values")

    output = np.empty((count, 29), dtype=np.float64)
    status = np.empty(count, dtype=np.int32)
    if count == 0:
        return output, status
    if <size_t>count > (<size_t>-1) // sizeof(rp_kerr_status):
        raise OverflowError("too many rows for native status buffer")
    row_status = <rp_kerr_status *>malloc(<size_t>count * sizeof(rp_kerr_status))
    if row_status == NULL:
        raise MemoryError("could not allocate native status buffer")
    try:
        with nogil:
            result = rp_solution_reconstruct_batch(
                spin,
                &cartesian[0, 0],
                <size_t>count,
                &output[0, 0],
                row_status,
            )
        for i in range(count):
            status[i] = <int>row_status[i]
    finally:
        free(row_status)
    return output, status


cdef extern from "geodesic/initial/convert.h":
    rp_kerr_status rp_initial_cartesian_to_canonical(
        double spin, const double cartesian[7], double canonical[8]
    ) nogil
    rp_kerr_status rp_initial_from_spherical(
        double spin, const double spherical[7], double canonical[8]
    ) nogil
    rp_kerr_status rp_initial_from_bl(
        double spin, const double bl[7], double canonical[8]
    ) nogil
    rp_kerr_status rp_initial_from_elements(
        double spin, const double elements[7], double canonical[8]
    ) nogil
    rp_kerr_status rp_initial_canonical_to_cartesian(
        double spin, const double canonical[8], double cartesian[7]
    ) nogil

cdef extern from "geodesic/integrators/integrator.h":
    ctypedef enum rp_integrator_status:
        RP_INTEGRATOR_STATUS_OK
        RP_INTEGRATOR_STATUS_OBSERVER_STOPPED
        RP_INTEGRATOR_STATUS_OBSERVER_FAILURE
        RP_INTEGRATOR_STATUS_PROJECTION_FAILURE
    ctypedef enum rp_integrator_method:
        RP_INTEGRATOR_METHOD_DOP853
        RP_INTEGRATOR_METHOD_RADAU
        RP_INTEGRATOR_METHOD_DP45
        RP_INTEGRATOR_METHOD_PROJECTION_RADAU
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
    rp_integrator_status rp_kerr_trajectory_solve(
        double spin, const double initial_state[8],
        double tau_initial, double tau_final, rp_integrator_method method,
        double rtol, double scalar_atol, const double *vector_atol,
        double first_step, double max_step, int store_steps,
        rp_kerr_trajectory *trajectory
    ) nogil
    void rp_kerr_trajectory_free(rp_kerr_trajectory *trajectory) nogil


cdef extern from "geodesic/initial/bound.h":
    rp_kerr_status rp_initial_from_bound(
        double spin, const double *elements, double *canonical
    ) nogil


def initial_kerr(double spin, family, components):
    """Convert one normalized position and coordinate velocity in native C."""
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] source
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] output
    cdef rp_kerr_status status
    if not np.isfinite(spin) or spin < 0.0 or spin > 1.0:
        raise ValueError("spin must be finite and in [0, 1]")
    source = np.ascontiguousarray(components, dtype=np.float64)
    if source.ndim != 1 or source.shape[0] != 7:
        raise ValueError("components must have shape (7,)")
    if not np.all(np.isfinite(source)):
        raise ValueError("components must be finite")
    output = np.empty(8, dtype=np.float64)
    if family == "cartesian":
        with nogil:
            status = rp_initial_cartesian_to_canonical(spin, &source[0], &output[0])
    elif family == "spherical":
        with nogil:
            status = rp_initial_from_spherical(spin, &source[0], &output[0])
    elif family == "bl":
        with nogil:
            status = rp_initial_from_bl(spin, &source[0], &output[0])
    elif family == "elements":
        with nogil:
            status = rp_initial_from_elements(spin, &source[0], &output[0])
    elif family == "bound":
        with nogil:
            status = rp_initial_from_bound(spin, &source[0], &output[0])
    else:
        raise ValueError("family must be cartesian, spherical, bl, elements or bound")
    if status != RP_KERR_STATUS_OK:
        if family == "bound":
            raise ValueError(
                "bound Kerr parameters must define a stable exterior orbit "
                "within the Boyer-Lindquist chart "
                f"(native status {<int>status})"
            )
        raise ValueError(f"initial Kerr state is invalid (native status {<int>status})")
    return output


def initial_cartesian_batch(double spin, cartesian):
    """Convert normalized Cartesian rows to canonical BL states in native C."""
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="c"] source
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="c"] output
    cdef cnp.ndarray[cnp.int32_t, ndim=1, mode="c"] statuses
    cdef rp_kerr_status status
    cdef Py_ssize_t count, i
    if not np.isfinite(spin) or spin < 0.0 or spin > 1.0:
        raise ValueError("spin must be finite and in [0, 1]")
    source = np.ascontiguousarray(cartesian, dtype=np.float64)
    if source.ndim != 2 or source.shape[1] != 7:
        raise ValueError("cartesian must have shape (n, 7)")
    count = source.shape[0]
    output = np.zeros((count, 8), dtype=np.float64)
    statuses = np.zeros(count, dtype=np.int32)
    if count:
        with nogil:
            for i in range(count):
                status = rp_initial_cartesian_to_canonical(
                    spin, &source[i, 0], &output[i, 0]
                )
                statuses[i] = <int>status
    output.flags.writeable = False
    statuses.flags.writeable = False
    return output, statuses


def reconstruct_canonical_batch(double spin, states):
    """Convert canonical BL states in one C call, preserving supplied x/u."""
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="c"] source
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="c"] cartesian
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="c"] output
    cdef cnp.ndarray[cnp.int32_t, ndim=1, mode="c"] statuses
    cdef rp_kerr_status *row_status = NULL
    cdef rp_kerr_status native_status
    cdef Py_ssize_t count, i
    source = np.ascontiguousarray(states, dtype=np.float64)
    if source.ndim != 2 or source.shape[1] != 8:
        raise ValueError("states must have shape (n, 8)")
    count = source.shape[0]
    cartesian = np.zeros((count, 7), dtype=np.float64)
    output = np.zeros((count, 29), dtype=np.float64)
    statuses = np.zeros(count, dtype=np.int32)
    if count == 0:
        return output, cartesian, statuses
    if <size_t>count > (<size_t>-1) // sizeof(rp_kerr_status):
        raise OverflowError("too many states")
    row_status = <rp_kerr_status *>malloc(<size_t>count * sizeof(rp_kerr_status))
    if row_status == NULL:
        raise MemoryError("could not allocate native status buffer")
    try:
        with nogil:
            native_status = rp_solution_reconstruct_canonical_batch(
                spin, &source[0, 0], <size_t>count,
                &output[0, 0], &cartesian[0, 0], row_status
            )
        for i in range(count):
            statuses[i] = <int>row_status[i]
    finally:
        free(row_status)
    return output, cartesian, statuses


def reconstruct_canonical_family_batch(double spin, states, family):
    """Return one selected geometric-unit family and native row statuses."""
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="c"] source
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="c"] output
    cdef cnp.ndarray[cnp.int32_t, ndim=1, mode="c"] statuses
    cdef rp_kerr_status *row_status = NULL
    cdef rp_kerr_status native_status
    cdef rp_solution_family selected
    cdef Py_ssize_t count, i, width

    if family == "cartesian":
        selected = RP_SOLUTION_FAMILY_CARTESIAN
        width = 7
    elif family == "spherical":
        selected = RP_SOLUTION_FAMILY_SPHERICAL
        width = 7
    elif family == "elements":
        selected = RP_SOLUTION_FAMILY_ELEMENTS
        width = 6
    else:
        raise ValueError("family must be cartesian, spherical or elements")
    source = np.ascontiguousarray(states, dtype=np.float64)
    if source.ndim != 2 or source.shape[1] != 8:
        raise ValueError("states must have shape (n, 8)")
    count = source.shape[0]
    output = np.zeros((count, width), dtype=np.float64)
    statuses = np.zeros(count, dtype=np.int32)
    if count == 0:
        output.flags.writeable = False
        statuses.flags.writeable = False
        return output, statuses
    if <size_t>count > (<size_t>-1) // sizeof(rp_kerr_status):
        raise OverflowError("too many states")
    row_status = <rp_kerr_status *>malloc(<size_t>count * sizeof(rp_kerr_status))
    if row_status == NULL:
        raise MemoryError("could not allocate native status buffer")
    try:
        with nogil:
            native_status = rp_solution_reconstruct_canonical_family_batch(
                spin, &source[0, 0], <size_t>count, selected,
                &output[0, 0], row_status
            )
        for i in range(count):
            statuses[i] = <int>row_status[i]
    finally:
        free(row_status)
    output.flags.writeable = False
    statuses.flags.writeable = False
    return output, statuses


def integrate_kerr(
    double spin, y0, double tau0, double tau1, method="radau",
    double rtol=1e-3, atol=1e-6, first_step=None, max_step=None,
    bint store_steps=True,
):
    """Integrate a canonical Kerr state in C and return accepted samples."""
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] initial
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] vector_atol
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] taus
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="c"] states
    cdef rp_integrator_method selected_method
    cdef rp_integrator_status native_status
    cdef rp_kerr_trajectory trajectory
    cdef const double *atol_pointer = NULL
    cdef double scalar_atol = 0.0
    cdef double initial_step = 0.0
    cdef double maximum_step = 0.0
    cdef Py_ssize_t i, j, count
    cdef int status
    if method == "radau":
        selected_method = RP_INTEGRATOR_METHOD_RADAU
    elif method == "dop853":
        selected_method = RP_INTEGRATOR_METHOD_DOP853
    elif method == "dp45":
        selected_method = RP_INTEGRATOR_METHOD_DP45
    elif method == "projection_radau":
        selected_method = RP_INTEGRATOR_METHOD_PROJECTION_RADAU
    else:
        raise ValueError("method must be one of radau, dop853, dp45, projection_radau")
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
        native_status = rp_kerr_trajectory_solve(
            spin, &initial[0], tau0, tau1, selected_method,
            rtol, scalar_atol, atol_pointer, initial_step, maximum_step,
            <int>store_steps, &trajectory
        )
    try:
        if trajectory.allocation_failed:
            raise MemoryError("native trajectory allocation failed")
        if trajectory.count == 0:
            raise ValueError(f"native Kerr integration rejected input (status {<int>native_status})")
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


cdef extern from "geodesic/solution/preview.h":
    rp_kerr_status rp_solution_osculating_preview(
        double spin, const double canonical[8], size_t count,
        double *xyz, double references[3]
    ) nogil


def osculating_preview(double spin, canonical, int samples=256):
    """Return normalized Newtonian osculating positions and reference radii."""
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] state
    cdef cnp.ndarray[cnp.float64_t, ndim=2, mode="c"] points
    cdef cnp.ndarray[cnp.float64_t, ndim=1, mode="c"] references
    cdef rp_kerr_status status
    if samples < 2 or samples > 1000000:
        raise ValueError("samples must be in [2, 1000000]")
    state = np.ascontiguousarray(canonical, dtype=np.float64)
    if state.ndim != 1 or state.shape[0] != 8:
        raise ValueError("canonical must have shape (8,)")
    points = np.empty((samples, 3), dtype=np.float64)
    references = np.empty(3, dtype=np.float64)
    with nogil:
        status = rp_solution_osculating_preview(
            spin, &state[0], <size_t>samples, &points[0, 0], &references[0]
        )
    if status != RP_KERR_STATUS_OK:
        raise ValueError(f"osculating preview is undefined (native status {<int>status})")
    return points, references
