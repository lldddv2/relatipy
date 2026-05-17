"""
Radau IIA 3-stage (order 5) geodesic integrator for the Kerr metric.

This module wraps a compiled C implementation that integrates geodesics in
covariant phase space :math:`(q^\\mu, p_\\mu)` to avoid the Boyer-Lindquist
coordinate singularity at the polar axis. Stored output uses
:math:`(q^\\mu, u^\\mu)` with constraint projection (normalization, conserved
energy :math:`E`, :math:`L_z`, and Carter constant :math:`Q`) applied via
Newton iteration at every stored output point.

The shared library ``radau_core`` (e.g. ``radau_core.so`` / ``radau_core.dll``)
is loaded from this package directory at
import time; if loading fails, public functions raise ``RuntimeError`` when
called.

Notes
-----
Integration uses the Radau IIA method (3 stages, classical order 5) as in
Hairer & Wanner. The C backend exposes ``kerr_radau2`` for stepping and
``kerr_project_trajectory`` / ``kerr_project_constraints`` for constraint
projection.

References
----------
.. [1] Hairer, E., & Wanner, G. (1996). *Solving Ordinary Differential
       Equations II: Stiff and Differential-Algebraic Problems*.
       Springer, §IV.5 (Radau IIA).

Examples
--------
The following assumes ``radau_core.so`` is built and import succeeds:

>>> from relatipy.numeric.geodesic.integrators import radau
>>> radau._C_LIB is not None  # doctest: +SKIP
True
"""

import ctypes
import os
import sys

import numpy as np


def _find_radau_core_path() -> str | None:
    """Path to the ctypes shared library next to this module, if present."""
    base = os.path.dirname(os.path.abspath(__file__))
    if sys.platform == "win32":
        names = ("radau_core.dll", "radau_core.so")
    else:
        names = ("radau_core.so", "radau_core.dylib")
    for name in names:
        p = os.path.join(base, name)
        if os.path.isfile(p):
            return p
    return None


# ---- Load C backend (radau_core.{so,dll,...}) ------------------------------
_C_LIB = None
_RADAU_LOAD_ERROR: str | None = None
try:
    _SO_PATH = _find_radau_core_path()
    if _SO_PATH is None:
        raise FileNotFoundError("radau_core shared library not found next to radau.py")
    _lib = ctypes.CDLL(_SO_PATH)
    _dbl_p = ctypes.POINTER(ctypes.c_double)
    _int_p = ctypes.POINTER(ctypes.c_int)

    _lib.kerr_project_constraints.restype = None
    _lib.kerr_project_constraints.argtypes = [
        ctypes.c_double,  # Rs
        ctypes.c_double,  # a
        _dbl_p,  # q[4]  fixed position
        _dbl_p,  # u[4]  velocity, corrected in-place
        ctypes.c_double,  # E0
        ctypes.c_double,  # Lz0
        ctypes.c_double,  # Q0
        ctypes.c_double,  # tol
        ctypes.c_int,  # max_iter
    ]
    _lib.kerr_project_trajectory.restype = None
    _lib.kerr_project_trajectory.argtypes = [
        ctypes.c_double,  # Rs
        ctypes.c_double,  # a
        _dbl_p,  # q_arr (N,4) row-major, fixed
        _dbl_p,  # u_arr (N,4) row-major, corrected in-place
        ctypes.c_double,  # E0
        ctypes.c_double,  # Lz0
        ctypes.c_double,  # Q0
        ctypes.c_double,  # tol
        ctypes.c_int,  # max_iter
        ctypes.c_int,  # N
    ]
    _lib.kerr_radau2.restype = None
    _lib.kerr_radau2.argtypes = [
        ctypes.c_double,  # Rs
        ctypes.c_double,  # a
        _dbl_p,  # state0 [8]
        ctypes.c_double,  # tau0
        ctypes.c_double,  # tau_end
        ctypes.c_int,  # n_steps
        ctypes.c_int,  # stride
        ctypes.c_int,  # n_fix_iter
        ctypes.c_double,  # E0
        ctypes.c_double,  # Lz0
        ctypes.c_double,  # Q0
        _dbl_p,  # out_t
        _dbl_p,  # out_y
        _int_p,  # n_out
    ]
    _lib.kerr_radau2_adaptive.restype = ctypes.c_int
    _lib.kerr_radau2_adaptive.argtypes = [
        ctypes.c_double,   # Rs
        ctypes.c_double,   # a
        _dbl_p,            # y0 [8]
        ctypes.c_double,   # t0
        ctypes.c_double,   # t_bound
        ctypes.c_double,   # rtol
        ctypes.c_double,   # atol
        ctypes.c_double,   # max_step
        ctypes.c_double,   # E0
        ctypes.c_double,   # Lz0
        ctypes.c_double,   # Q0
        ctypes.c_double,   # tol_proj
        ctypes.c_int,      # max_iter_proj
        _dbl_p,            # t_eval (NULL if not used)
        ctypes.c_int,      # n_eval
        _dbl_p,            # t_out
        _dbl_p,            # y_out
        ctypes.c_int,      # step_cap
        _int_p,            # n_out
        _int_p,            # n_steps_taken
        _int_p,            # n_jac
        _int_p,            # n_lu
    ]
    _C_LIB = _lib
except Exception as _e:
    _C_LIB = None
    _RADAU_LOAD_ERROR = str(_e)


def _radau_backend_error(what: str) -> RuntimeError:
    detail = f" ({_RADAU_LOAD_ERROR})" if _RADAU_LOAD_ERROR else ""
    return RuntimeError(
        f"C backend not loaded; cannot call {what}.{detail} "
        "If the shared library is for another OS (e.g. Linux ELF inside a py3-none-any wheel on macOS), "
        "reinstall from source on this machine: pip install --no-cache-dir /path/to/relatipy"
    )


def project_kerr_trajectory(Rs, a, q_arr, u_arr, E0, Lz0, Q0, tol=1e-12, max_iter=20):
    """
    Project Kerr four-velocities onto the constraint surface in place.

    For each row of ``q_arr``, ``u_arr`` is corrected by Newton iteration so
    that the four-velocity is normalized (:math:`g_{\\mu\\nu} u^\\mu u^\\nu = 1`)
    and the conserved quantities match ``E0``, ``Lz0``, and ``Q0``.

    Parameters
    ----------
    Rs : float
        Schwarzschild radius parameter for the Kerr metric.
    a : float
        Spin parameter :math:`a` (Boyer-Lindquist).
    q_arr : ndarray, shape (N, 4), dtype float64
        Positions :math:`q^\\mu = (t, r, \\theta, \\phi)` in Boyer-Lindquist
        coordinates. Must be C-contiguous; not modified.
    u_arr : ndarray, shape (N, 4), dtype float64
        Four-velocities :math:`u^\\mu`. Must be C-contiguous; updated in place.
    E0 : float
        Target conserved energy (per unit mass) :math:`E`.
    Lz0 : float
        Target conserved :math:`z`-component of angular momentum :math:`L_z`.
    Q0 : float
        Target Carter constant :math:`Q`.
    tol : float, optional
        Convergence tolerance for the Newton iteration (default ``1e-12``).
    max_iter : int, optional
        Maximum Newton iterations per point (default ``20``).

    Raises
    ------
    RuntimeError
        If the C shared library was not loaded at import time.

    Examples
    --------
    >>> import numpy as np
    >>> from relatipy.numeric.geodesic.integrators.radau import project_kerr_trajectory
    >>> q = np.zeros((1, 4), dtype=np.float64)
    >>> u = np.zeros((1, 4), dtype=np.float64)
    >>> project_kerr_trajectory(2.0, 0.5, q, u, 0.9, 0.1, 1.0)  # doctest: +SKIP

    See Also
    --------
    _integrate_kerr_radau2 : Full Radau geodesic integration with projection.
    """
    if _C_LIB is None:
        raise _radau_backend_error("project_kerr_trajectory")
    N = len(q_arr)
    _ptr = lambda arr: arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    _C_LIB.kerr_project_trajectory(
        ctypes.c_double(Rs),
        ctypes.c_double(a),
        _ptr(q_arr),
        _ptr(u_arr),
        ctypes.c_double(E0),
        ctypes.c_double(Lz0),
        ctypes.c_double(Q0),
        ctypes.c_double(tol),
        ctypes.c_int(max_iter),
        ctypes.c_int(N),
    )


def _integrate_kerr_radau2(
    Rs, a, state0, tau_span, n_steps, E0, Lz0, Q0, stride=1, n_fix_iter=5
):
    """
    Fixed-step Radau IIA integration for Kerr geodesics.

    Integrates from proper time ``tau_span[0]`` to ``tau_span[1]`` in ``n_steps``
    steps of the Radau IIA method. Constraint projection is applied at every
    stored output point inside the C implementation.

    Parameters
    ----------
    Rs : float
        Schwarzschild radius parameter for the Kerr metric.
    a : float
        Spin parameter :math:`a` (Boyer-Lindquist).
    state0 : array_like, shape (8,)
        Initial state :math:`[q^0, q^1, q^2, q^3, u^0, u^1, u^2, u^3]` as
        contiguous ``float64``.
    tau_span : tuple of float
        Pair ``(tau0, tau_end)`` giving the proper-time interval.
    n_steps : int
        Number of Radau steps between ``tau0`` and ``tau_end``.
    E0 : float
        Conserved energy (per unit mass) used in projection.
    Lz0 : float
        Conserved :math:`L_z` used in projection.
    Q0 : float
        Carter constant used in projection.
    stride : int, optional
        Store output every ``stride`` internal steps (default ``1``).
    n_fix_iter : int, optional
        Number of fixed-point constraint iterations per step in the C code
        (default ``5``).

    Returns
    -------
    result : object
        Namespace with attributes:

        ``t`` : ndarray, shape (N_out,)
            Proper times at stored points. The last entry is set to
            ``tau_span[1]`` to correct floating-point accumulation error.
        ``y`` : ndarray, shape (8, N_out)
            States :math:`[q^\\mu; u^\\mu]` at each stored time (columns are
            snapshots).

    Raises
    ------
    RuntimeError
        If the C shared library was not loaded at import time.

    Notes
    -----
    Preallocated buffers in C may be larger than ``N_out``; only the first
    ``N_out`` samples are returned.

    Examples
    --------
    >>> import numpy as np
    >>> from relatipy.numeric.geodesic.integrators.radau import _integrate_kerr_radau2
    >>> s0 = np.zeros(8, dtype=np.float64)
    >>> out = _integrate_kerr_radau2(2.0, 0.0, s0, (0.0, 1.0), 4, 1.0, 0.0, 0.0)  # doctest: +SKIP
    >>> out.t.shape == out.y.shape[1:]  # doctest: +SKIP
    True
    """
    if _C_LIB is None:
        raise _radau_backend_error("_integrate_kerr_radau2")

    N_alloc = n_steps + 2
    out_t = np.zeros(N_alloc, dtype=np.float64)
    out_y = np.zeros(N_alloc * 8, dtype=np.float64)
    n_out = ctypes.c_int(0)

    s0 = np.ascontiguousarray(state0, dtype=np.float64)
    _ptr = lambda arr: arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    _C_LIB.kerr_radau2(
        ctypes.c_double(Rs),
        ctypes.c_double(a),
        _ptr(s0),
        ctypes.c_double(tau_span[0]),
        ctypes.c_double(tau_span[1]),
        ctypes.c_int(n_steps),
        ctypes.c_int(stride),
        ctypes.c_int(n_fix_iter),
        ctypes.c_double(E0),
        ctypes.c_double(Lz0),
        ctypes.c_double(Q0),
        _ptr(out_t),
        _ptr(out_y),
        ctypes.byref(n_out),
    )

    N = n_out.value

    class _Result:
        """
        Lightweight holder for Radau integration output.

        Attributes
        ----------
        t : ndarray
            Proper times at output points, shape ``(N_out,)``.
        y : ndarray
            State vectors, shape ``(8, N_out)``.
        """

    result = _Result()
    result.t = out_t[:N].copy()
    result.t[-1] = float(tau_span[1])  # correct float-accumulation error
    result.y = out_y[: N * 8].reshape(N, 8).T.copy()  # (8, N)
    return result


def _integrate_kerr_radau2_adaptive(
    Rs, a, state0, tau_span, E0, Lz0, Q0,
    rtol=1e-9, atol=1e-12, t_eval=None, max_step=0.0,
    tol_proj=1e-12, max_iter_proj=20,
):
    """Adaptive Radau IIA integration for Kerr geodesics (C backend).

    Parameters
    ----------
    Rs, a : float
        Metric parameters.
    state0 : array_like, shape (8,)
        Initial state [q, u] (contravariant).
    tau_span : (float, float)
        Integration interval.
    E0, Lz0, Q0 : float
        Conserved quantities for projection.
    rtol, atol : float
        Tolerances.
    t_eval : array_like or None
        If given, output is sampled at these proper times via dense interpolation.
        Otherwise every accepted step is stored.
    max_step : float
        Maximum step size (0 = unbounded).
    tol_proj, max_iter_proj : float, int
        Newton projection parameters.

    Returns
    -------
    result : namespace with `.t` (N,) and `.y` (8, N).
    """
    if _C_LIB is None:
        raise _radau_backend_error("_integrate_kerr_radau2_adaptive")

    _ptr = lambda arr: arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    s0 = np.ascontiguousarray(state0, dtype=np.float64)

    if t_eval is not None:
        t_eval_np = np.ascontiguousarray(t_eval, dtype=np.float64)
        n_eval = len(t_eval_np)
        N_alloc = n_eval
        t_out = np.zeros(N_alloc, dtype=np.float64)
        y_out = np.zeros(N_alloc * 8, dtype=np.float64)
        t_eval_ptr = _ptr(t_eval_np)
        step_cap = N_alloc
    else:
        n_eval = 0
        N_alloc = 100_000
        t_out = np.zeros(N_alloc, dtype=np.float64)
        y_out = np.zeros(N_alloc * 8, dtype=np.float64)
        t_eval_ptr = ctypes.cast(None, ctypes.POINTER(ctypes.c_double))
        step_cap = N_alloc

    n_out = ctypes.c_int(0)
    nst   = ctypes.c_int(0)
    njac  = ctypes.c_int(0)
    nlu   = ctypes.c_int(0)

    ret = _C_LIB.kerr_radau2_adaptive(
        ctypes.c_double(Rs),
        ctypes.c_double(a),
        _ptr(s0),
        ctypes.c_double(tau_span[0]),
        ctypes.c_double(tau_span[1]),
        ctypes.c_double(rtol),
        ctypes.c_double(atol),
        ctypes.c_double(max_step),
        ctypes.c_double(E0),
        ctypes.c_double(Lz0),
        ctypes.c_double(Q0),
        ctypes.c_double(tol_proj),
        ctypes.c_int(max_iter_proj),
        t_eval_ptr,
        ctypes.c_int(n_eval),
        _ptr(t_out),
        _ptr(y_out),
        ctypes.c_int(step_cap),
        ctypes.byref(n_out),
        ctypes.byref(nst),
        ctypes.byref(njac),
        ctypes.byref(nlu),
    )

    if ret not in (0,):
        msg = {-1: "max steps exceeded", -2: "Newton/LU failure",
               -4: "step underflow"}.get(ret, f"error {ret}")
        raise RuntimeError(f"kerr_radau2_adaptive: {msg}")

    N = n_out.value

    class _Result:
        pass

    result = _Result()
    result.t = t_out[:N].copy()
    result.y = y_out[:N * 8].reshape(N, 8).T.copy()  # (8, N)
    result.nsteps = nst.value
    result.njac   = njac.value
    result.nlu    = nlu.value
    return result
