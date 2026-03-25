"""
Radau IIA 3-stage (order 5) geodesic integrator for the Kerr metric.

Integrates in covariant phase space (q^mu, p_mu) to avoid the
Boyer-Lindquist coordinate singularity at the polar axis.
Constraint projection (normalization, E, Lz, Carter Q) is applied
via Newton iteration at every stored output point.

The compiled C backend (radau_core.so) is loaded automatically.
"""

import ctypes
import os
import numpy as np

# ---- Load C backend (radau_core.so) ----------------------------------------
_C_LIB = None
try:
    _SO_PATH = os.path.join(os.path.dirname(__file__), "radau_core.so")
    _lib = ctypes.CDLL(_SO_PATH)
    _dbl_p = ctypes.POINTER(ctypes.c_double)
    _int_p = ctypes.POINTER(ctypes.c_int)

    _lib.kerr_project_constraints.restype = None
    _lib.kerr_project_constraints.argtypes = [
        ctypes.c_double,   # Rs
        ctypes.c_double,   # a
        _dbl_p,            # q[4]  fixed position
        _dbl_p,            # u[4]  velocity, corrected in-place
        ctypes.c_double,   # E0
        ctypes.c_double,   # Lz0
        ctypes.c_double,   # Q0
        ctypes.c_double,   # tol
        ctypes.c_int,      # max_iter
    ]
    _lib.kerr_project_trajectory.restype = None
    _lib.kerr_project_trajectory.argtypes = [
        ctypes.c_double,   # Rs
        ctypes.c_double,   # a
        _dbl_p,            # q_arr (N,4) row-major, fixed
        _dbl_p,            # u_arr (N,4) row-major, corrected in-place
        ctypes.c_double,   # E0
        ctypes.c_double,   # Lz0
        ctypes.c_double,   # Q0
        ctypes.c_double,   # tol
        ctypes.c_int,      # max_iter
        ctypes.c_int,      # N
    ]
    _lib.kerr_radau2.restype = None
    _lib.kerr_radau2.argtypes = [
        ctypes.c_double,   # Rs
        ctypes.c_double,   # a
        _dbl_p,            # state0 [8]
        ctypes.c_double,   # tau0
        ctypes.c_double,   # tau_end
        ctypes.c_int,      # n_steps
        ctypes.c_int,      # stride
        ctypes.c_int,      # n_fix_iter
        ctypes.c_double,   # E0
        ctypes.c_double,   # Lz0
        ctypes.c_double,   # Q0
        _dbl_p,            # out_t  [n_steps+2]
        _dbl_p,            # out_y  [(n_steps+2)*8], row-major (N,8)
        _int_p,            # n_out
    ]
    _C_LIB = _lib
except Exception:
    pass  # C backend unavailable


def project_kerr_trajectory(Rs, a, q_arr, u_arr, E0, Lz0, Q0,
                             tol=1e-12, max_iter=20):
    """
    In-place Newton projection of Kerr 4-velocities onto the constraint surface
    {normalization=1, E=E0, Lz=Lz0, Carter Q=Q0}.

    Parameters
    ----------
    Rs, a   : Kerr metric parameters
    q_arr   : (N, 4) C-contiguous float64 — positions (fixed)
    u_arr   : (N, 4) C-contiguous float64 — velocities, corrected in-place
    E0,Lz0,Q0 : conserved quantities
    """
    if _C_LIB is None:
        raise RuntimeError("C backend not loaded; cannot call project_kerr_trajectory")
    N = len(q_arr)
    _ptr = lambda arr: arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    _C_LIB.kerr_project_trajectory(
        ctypes.c_double(Rs), ctypes.c_double(a),
        _ptr(q_arr), _ptr(u_arr),
        ctypes.c_double(E0), ctypes.c_double(Lz0), ctypes.c_double(Q0),
        ctypes.c_double(tol), ctypes.c_int(max_iter), ctypes.c_int(N),
    )


def _integrate_kerr_radau2(Rs, a, state0, tau_span, n_steps,
                            E0, Lz0, Q0, stride=1, n_fix_iter=5):
    """
    Fixed-step Radau IIA (order 5) integration for Kerr geodesics.

    Constraint projection is applied at every stored output point inside C.

    Parameters
    ----------
    Rs, a       : Kerr metric parameters
    state0      : array_like (8,) — [q^mu, u^mu]
    tau_span    : (tau0, tau_end)
    n_steps     : number of Radau steps
    E0,Lz0,Q0  : conserved quantities (energy, z-angular momentum, Carter Q)
    stride      : store output every stride steps (default 1)
    n_fix_iter  : fixed-point iterations per step (default 5)

    Returns
    -------
    result.t : ndarray (N_out,)
    result.y : ndarray (8, N_out)
    """
    if _C_LIB is None:
        raise RuntimeError("C backend not loaded; cannot call _integrate_kerr_radau2")

    N_alloc = n_steps + 2
    out_t = np.zeros(N_alloc, dtype=np.float64)
    out_y = np.zeros(N_alloc * 8, dtype=np.float64)
    n_out = ctypes.c_int(0)

    s0 = np.ascontiguousarray(state0, dtype=np.float64)
    _ptr = lambda arr: arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    _C_LIB.kerr_radau2(
        ctypes.c_double(Rs), ctypes.c_double(a),
        _ptr(s0),
        ctypes.c_double(tau_span[0]), ctypes.c_double(tau_span[1]),
        ctypes.c_int(n_steps), ctypes.c_int(stride), ctypes.c_int(n_fix_iter),
        ctypes.c_double(E0), ctypes.c_double(Lz0), ctypes.c_double(Q0),
        _ptr(out_t), _ptr(out_y),
        ctypes.byref(n_out),
    )

    N = n_out.value

    class _Result:
        pass

    result = _Result()
    result.t = out_t[:N].copy()
    result.t[-1] = float(tau_span[1])           # correct float-accumulation error
    result.y = out_y[:N * 8].reshape(N, 8).T.copy()  # (8, N)
    return result
