"""
Yoshida 6th-order symplectic integrator for geodesic equations.

Works in Hamiltonian phase space (q^μ, p_μ) where p_μ = g_{μν} u^ν.
Uses Störmer-Verlet (DKD) as the 2nd-order base and a 7-stage Yoshida
composition for 6th-order accuracy.

The kick force is computed using the identity

    F_μ = g_{μρ} a^ρ,   a^ρ = -Γ^ρ_{αβ} u^α u^β

which avoids finite differences entirely and reuses the analytical
Christoffel symbols already implemented for each metric.

Reference: Yoshida, H. (1990), Phys. Lett. A 150, 262-268.
"""

import ctypes
import os
import numpy as np

# ---- Load C backend (yoshida6_core.so) if available ----
_C_LIB = None
try:
    _SO_PATH = os.path.join(os.path.dirname(__file__), "yoshida6_core.so")
    _lib = ctypes.CDLL(_SO_PATH)
    _dbl_p = ctypes.POINTER(ctypes.c_double)
    _int_p = ctypes.POINTER(ctypes.c_int)

    _lib.yoshida6_schwarzschild.restype = None
    _lib.yoshida6_schwarzschild.argtypes = [
        ctypes.c_double,   # Rs
        _dbl_p,            # state0
        ctypes.c_double,   # tau0
        ctypes.c_double,   # tau_end
        ctypes.c_int,      # n_steps
        ctypes.c_int,      # stride
        _dbl_p,            # out_t
        _dbl_p,            # out_y  (N*8, row-major)
        _int_p,            # n_out
    ]
    _lib.yoshida6_kerr.restype = None
    _lib.yoshida6_kerr.argtypes = [
        ctypes.c_double,   # Rs
        ctypes.c_double,   # a
        _dbl_p,            # state0
        ctypes.c_double,   # tau0
        ctypes.c_double,   # tau_end
        ctypes.c_int,      # n_steps
        ctypes.c_int,      # stride
        _dbl_p,            # out_t
        _dbl_p,            # out_y  (N*8, row-major)
        _int_p,            # n_out
    ]
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
    pass  # fall back to pure Python


def project_kerr_trajectory(Rs, a, q_arr, u_arr, E0, Lz0, Q0,
                             tol=1e-12, max_iter=20):
    """
    In-place Newton projection of Kerr 4-velocities onto constraint surface.

    Parameters
    ----------
    Rs, a   : Kerr metric parameters (self.metric.R_s, self.metric.a)
    q_arr   : (N, 4) C-contiguous float64 — positions (fixed)
    u_arr   : (N, 4) C-contiguous float64 — velocities, corrected in-place
    E0, Lz0, Q0 : conserved quantities from initial conditions
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
    Fixed-step Radau IIA (order 5) integration for Kerr geodesics in C.

    Constraint projection (normalization, E, Lz, Carter Q) is applied
    at every stored output point inside the C loop.

    Parameters
    ----------
    Rs, a       : Kerr metric parameters
    state0      : array_like (8,) — [q^mu, u^mu]
    tau_span    : (tau0, tau_end)
    n_steps     : number of Radau steps
    E0,Lz0,Q0  : conserved quantities (computed from initial state)
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
    result.t[-1] = float(tau_span[1])          # correct float-accumulation error
    result.y = out_y[:N * 8].reshape(N, 8).T.copy()  # (8, N)
    return result


# Yoshida 6th-order composition weights (7 stages, palindromic)
_W1 =  0.784513610477560
_W2 =  0.235573213359357
_W3 = -1.17767998417887
_W0 =  1.0 - 2.0 * (_W1 + _W2 + _W3)   # ≈ 1.31518632068391
_YOSHIDA6_WEIGHTS = [_W1, _W2, _W3, _W0, _W3, _W2, _W1]


class Yoshida6Integrator:
    """
    Yoshida 6th-order symplectic integrator for geodesic equations.

    State input/output uses the same convention as all other integrators:
        state = [q^0, q^1, q^2, q^3, u^0, u^1, u^2, u^3]

    Internally works in (q, p) phase space with p_μ = g_{μν} u^ν.
    """

    def __init__(self, metric):
        self.metric = metric

    # ------------------------------------------------------------------
    # Core building blocks
    # ------------------------------------------------------------------

    def _kick_force(self, q, p, g, g_inv):
        """
        Kick force in covariant phase space.

        From the geodesic Hamiltonian H = ½ g^{μν} p_μ p_ν the canonical
        equation gives:

            dp_μ/dτ = Γ_{σρμ} u^σ u^ρ = g_{σα} Γ^α_{ρμ} u^σ u^ρ

        where Γ^α_{ρμ} are the Christoffel symbols (already implemented
        analytically for each metric).

        Parameters
        ----------
        q      : ndarray(4,) — position
        p      : ndarray(4,) — covariant momenta p_μ
        g      : ndarray(4,4) — metric g_{μν}(q)
        g_inv  : ndarray(4,4) — inverse metric g^{μν}(q)

        Returns
        -------
        ndarray(4,) — kick force F_μ
        """
        u     = g_inv @ p                                   # u^α = g^{αβ} p_β
        chris = self.metric.get_christoffel_symbols(q)      # Γ^α_{ρμ}, shape (4,4,4)
        # F_m = g_{s,a} Γ^a_{r,m} u^s u^r
        return np.einsum('sa,arm,s,r->m', g, chris, u, u)

    def _verlet_step(self, q, p, h):
        """
        One Störmer-Verlet DKD step.

        1. half-drift: q_{1/2} = q       + (h/2) g^{μν}(q)       p_ν
        2. kick:       p'      = p       + h     F(q_{1/2}, p)
        3. half-drift: q'      = q_{1/2} + (h/2) g^{μν}(q_{1/2}) p'_ν
        """
        # --- first half-drift ---
        g0     = self.metric.metric(q)
        g0_inv = np.linalg.inv(g0)
        q_half = q + (0.5 * h) * (g0_inv @ p)

        # --- kick at q_half (reuse g_half for second half-drift) ---
        g_half     = self.metric.metric(q_half)
        g_half_inv = np.linalg.inv(g_half)
        F          = self._kick_force(q_half, p, g_half, g_half_inv)
        p_new      = p + h * F

        # --- second half-drift (metric at q_half, updated p) ---
        q_new = q_half + (0.5 * h) * (g_half_inv @ p_new)

        return q_new, p_new

    def _yoshida_step(self, q, p, h):
        """One Yoshida 6th-order step = 7 Verlet sub-steps."""
        for w in _YOSHIDA6_WEIGHTS:
            q, p = self._verlet_step(q, p, w * h)
        return q, p

    # ------------------------------------------------------------------
    # C backend dispatch
    # ------------------------------------------------------------------

    def _integrate_c(self, state0, tau0, tau_end, n_steps, stride):
        """Call the compiled C integration function (Schwarzschild or Kerr)."""
        metric_name = type(self.metric).__name__
        N_alloc = n_steps + 2
        out_t = np.zeros(N_alloc, dtype=np.float64)
        out_y = np.zeros(N_alloc * 8, dtype=np.float64)
        n_out = ctypes.c_int(0)

        s0 = np.ascontiguousarray(state0, dtype=np.float64)
        _dbl_p = ctypes.POINTER(ctypes.c_double)

        def _ptr(arr):
            return arr.ctypes.data_as(_dbl_p)

        if metric_name == "Schwarzschild":
            _C_LIB.yoshida6_schwarzschild(
                ctypes.c_double(self.metric.R_s),
                _ptr(s0),
                ctypes.c_double(tau0), ctypes.c_double(tau_end),
                ctypes.c_int(n_steps), ctypes.c_int(stride),
                _ptr(out_t), _ptr(out_y),
                ctypes.byref(n_out),
            )
        elif metric_name == "Kerr":
            _C_LIB.yoshida6_kerr(
                ctypes.c_double(self.metric.R_s),
                ctypes.c_double(self.metric.a),
                _ptr(s0),
                ctypes.c_double(tau0), ctypes.c_double(tau_end),
                ctypes.c_int(n_steps), ctypes.c_int(stride),
                _ptr(out_t), _ptr(out_y),
                ctypes.byref(n_out),
            )
        else:
            raise ValueError(f"C backend not available for metric '{metric_name}'")

        N = n_out.value
        class _Result:
            pass
        result   = _Result()
        result.t = out_t[:N].copy()
        result.y = out_y[: N * 8].reshape(N, 8).T.copy()  # (8, N)
        return result

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def integrate(self, state0, tau_span, n_steps, stride=1):
        """
        Integrate geodesic equations with fixed step size.

        Parameters
        ----------
        state0 : array_like, shape (8,)
            Initial state [q^0..q^3, u^0..u^3].
        tau_span : (tau0, tau_end)
            Proper-time integration interval.
        n_steps : int
            Total number of Yoshida steps (h = Δτ / n_steps).
        stride : int
            Store output every `stride` steps (default 1).
            The initial state is always stored.

        Returns
        -------
        result : namespace
            .t — ndarray(N_out,)   proper time at output points
            .y — ndarray(8, N_out) state [q^μ, u^μ] at output points
        """
        state0 = np.asarray(state0, dtype=float)
        if state0.shape != (8,):
            raise ValueError("state0 must have shape (8,)")

        tau0, tau_end = float(tau_span[0]), float(tau_span[1])
        n_steps = int(n_steps)
        stride  = max(1, int(stride))

        # Try the compiled C backend first
        metric_name = type(self.metric).__name__
        if _C_LIB is not None and metric_name in ("Schwarzschild", "Kerr"):
            return self._integrate_c(state0, tau0, tau_end, n_steps, stride)

        # ---- Pure Python fallback ----
        h  = (tau_end - tau0) / n_steps
        q  = state0[:4].copy()
        u0 = state0[4:].copy()

        g0     = self.metric.metric(q)
        g0_inv = np.linalg.inv(g0)
        p      = g0 @ u0

        tau        = tau0
        tau_list   = [tau0]
        state_list = [np.concatenate([q, g0_inv @ p])]

        for i in range(n_steps):
            q, p = self._yoshida_step(q, p, h)
            tau += h

            if (i + 1) % stride == 0 or (i + 1) == n_steps:
                g     = self.metric.metric(q)
                g_inv = np.linalg.inv(g)
                tau_list.append(tau)
                state_list.append(np.concatenate([q, g_inv @ p]))

        tau_arr   = np.array(tau_list)
        state_arr = np.array(state_list).T   # (8, N_out)

        class _Result:
            pass
        result   = _Result()
        result.t = tau_arr
        result.y = state_arr
        return result


def yoshida6_integrate_geodesic(
    metric,
    state0,
    tau0,
    tau_end,
    n_steps,
    adaptive_output=True,
    tau_eval=None,
):
    """
    Convenience wrapper for Yoshida 6th-order geodesic integration.

    Returns
    -------
    tau   : ndarray, shape (N,)
    state : ndarray, shape (N, 8) — [q^μ, u^μ] at each output point
    """
    state0 = np.asarray(state0, dtype=float)

    if adaptive_output or tau_eval is None:
        stride = 1
    else:
        tau_eval = np.asarray(tau_eval, dtype=float)
        stride   = max(1, n_steps // len(tau_eval))

    result = Yoshida6Integrator(metric).integrate(
        state0, (tau0, tau_end), n_steps, stride=stride
    )
    return result.t, result.y.T
