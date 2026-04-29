"""
Yoshida 6th-order symplectic integrator for geodesic equations.

This module provides a fixed-step Yoshida composition integrator in Hamiltonian
phase space :math:`(q^\\mu, p_\\mu)` with :math:`p_\\mu = g_{\\mu\\nu} u^\\nu`.
It uses Störmer–Verlet (drift–kick–drift) as the second-order base scheme and a
seven-stage palindromic Yoshida composition for sixth-order accuracy.

The kick force uses the analytical Christoffel symbols (no finite differences):

.. math::

    F_\\mu = g_{\\mu\\rho} a^\\rho, \\quad
    a^\\rho = -\\Gamma^\\rho_{\\alpha\\beta} u^\\alpha u^\\beta

When the compiled extension ``yoshida6_core.so`` is present, Schwarzschild and
Kerr geodesics can be integrated in C; otherwise a pure Python path is used.

Notes
-----
The composition weights are fixed constants (palindromic seven-stage scheme).

References
----------
Yoshida, H. (1990). Construction of higher order symplectic integrators.
*Physics Letters A*, 150(5–7), 262–268.

Examples
--------
>>> import numpy as np
>>> # Integration requires a metric instance; see relatipy.numeric.metrics.
>>> # from relatipy.numeric.metrics import SchwarzschildMetric
>>> # m = SchwarzschildMetric(R_s=2.0)
>>> # y = Yoshida6Integrator(m)
>>> # out = y.integrate(np.zeros(8), (0.0, 1.0), n_steps=10)
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
    Project Kerr four-velocities onto the constraint surface (in-place).

    Uses a Newton iteration in the C backend to correct ``u_arr`` so that the
    normalization and conserved quantities (energy :math:`E`, axial angular
    momentum :math:`L_z`, Carter constant :math:`Q`) match the target values at
    each fixed position row in ``q_arr``.

    Parameters
    ----------
    Rs : float
        Schwarzschild radius parameter for the Kerr metric.
    a : float
        Spin parameter :math:`a` of the Kerr metric.
    q_arr : ndarray, shape (N, 4), dtype float64
        C-contiguous position coordinates :math:`q^\\mu` (read-only).
    u_arr : ndarray, shape (N, 4), dtype float64
        C-contiguous four-velocities :math:`u^\\mu`; updated in-place.
    E0 : float
        Target conserved energy (from initial conditions).
    Lz0 : float
        Target conserved axial angular momentum.
    Q0 : float
        Target Carter constant :math:`Q`.
    tol : float, optional
        Convergence tolerance for the projection (default ``1e-12``).
    max_iter : int, optional
        Maximum Newton iterations per point (default ``20``).

    Raises
    ------
    RuntimeError
        If the C extension is not loaded.

    Returns
    -------
    None
        Velocities are written in-place to ``u_arr``.

    Examples
    --------
    >>> import numpy as np
    >>> q = np.zeros((2, 4), dtype=np.float64)
    >>> u = np.zeros((2, 4), dtype=np.float64)
    >>> # project_kerr_trajectory(2.0, 0.0, q, u, E0, Lz0, Q0)  # needs C backend
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

    Constraint projection (normalization, :math:`E`, :math:`L_z`, Carter
    :math:`Q`) is applied at every stored output point inside the C loop.

    Parameters
    ----------
    Rs : float
        Schwarzschild radius parameter.
    a : float
        Kerr spin parameter.
    state0 : array_like, shape (8,)
        Initial phase-space vector :math:`[q^\\mu, u^\\mu]`.
    tau_span : tuple of float
        Proper-time interval ``(tau0, tau_end)``.
    n_steps : int
        Number of Radau steps.
    E0 : float
        Conserved energy (from the initial state).
    Lz0 : float
        Conserved axial angular momentum.
    Q0 : float
        Carter constant :math:`Q`.
    stride : int, optional
        Store output every ``stride`` steps (default ``1``).
    n_fix_iter : int, optional
        Fixed-point iterations per step in the C integrator (default ``5``).

    Returns
    -------
    result : simple namespace
        Object with attributes:

        - ``t`` : ndarray, shape (N_out,)
            Proper time at each output point.
        - ``y`` : ndarray, shape (8, N_out)
            States :math:`[q^\\mu, u^\\mu]` as columns.

    Raises
    ------
    RuntimeError
        If the C extension is not loaded.

    Examples
    --------
    >>> import numpy as np
    >>> s0 = np.zeros(8, dtype=np.float64)
    >>> # _integrate_kerr_radau2(Rs, a, s0, (0.0, 1.0), 4, E0, Lz0, Q0)
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
    Sixth-order Yoshida symplectic integrator for Hamiltonian geodesic flow.

    State vectors follow the project convention ``state = [q^0, q^1, q^2, q^3,
    u^0, u^1, u^2, u^3]``. Internally the integrator advances :math:`(q, p)` with
    :math:`p_\\mu = g_{\\mu\\nu} u^\\nu` using a seven-stage Yoshida composition
    of Störmer–Verlet steps.

    Parameters
    ----------
    metric : object
        Metric instance providing ``metric(q)`` (and Christoffel symbols via
        ``get_christoffel_symbols`` for the Python path). For the C backend,
        the class name must be ``Schwarzschild`` or ``Kerr``.

    Attributes
    ----------
    metric : object
        The spacetime metric used for integration.

    Examples
    --------
    >>> import numpy as np
    >>> # Yoshida6Integrator(metric).integrate(state0, (0.0, 1.0), n_steps=8)
    """

    def __init__(self, metric):
        """
        Initialize the integrator with a metric instance.

        Parameters
        ----------
        metric : object
            Metric providing ``metric(q)`` and ``get_christoffel_symbols(q)``
            for the pure Python integrator.

        Examples
        --------
        >>> # Yoshida6Integrator(schwarzschild_metric_instance)
        """
        self.metric = metric

    # ------------------------------------------------------------------
    # Core building blocks
    # ------------------------------------------------------------------

    def _kick_force(self, q, p, g, g_inv):
        """
        Covariant kick force :math:`F_\\mu` for the geodesic Hamiltonian flow.

        From :math:`H = \\tfrac{1}{2} g^{\\mu\\nu} p_\\mu p_\\nu`, the canonical
        equation gives

        .. math::

            \\frac{\\mathrm{d} p_\\mu}{\\mathrm{d}\\tau}
            = \\Gamma_{\\sigma\\rho\\mu} u^\\sigma u^\\rho
            = g_{\\sigma\\alpha} \\Gamma^\\alpha_{\\rho\\mu} u^\\sigma u^\\rho

        with Christoffel symbols :math:`\\Gamma^\\alpha_{\\rho\\mu}` supplied by
        ``self.metric.get_christoffel_symbols``.

        Parameters
        ----------
        q : ndarray, shape (4,)
            Position :math:`q^\\mu`.
        p : ndarray, shape (4,)
            Covariant momentum :math:`p_\\mu`.
        g : ndarray, shape (4, 4)
            Metric :math:`g_{\\mu\\nu}(q)`.
        g_inv : ndarray, shape (4, 4)
            Inverse metric :math:`g^{\\mu\\nu}(q)`.

        Returns
        -------
        ndarray, shape (4,)
            Kick force components :math:`F_\\mu`.

        Examples
        --------
        >>> import numpy as np
        >>> # F = integrator._kick_force(q, p, g, g_inv)  # inside Yoshida6Integrator
        """
        u     = g_inv @ p                                   # u^α = g^{αβ} p_β
        chris = self.metric.get_christoffel_symbols(q)      # Γ^α_{ρμ}, shape (4,4,4)
        # F_m = g_{s,a} Γ^a_{r,m} u^s u^r
        return np.einsum('sa,arm,s,r->m', g, chris, u, u)

    def _verlet_step(self, q, p, h):
        """
        One Störmer–Verlet drift–kick–drift sub-step of step size ``h``.

        The scheme applies: (1) half drift with :math:`g^{\\mu\\nu}(q)`;
        (2) kick using :meth:`_kick_force` at the midpoint; (3) second half
        drift using the metric at the midpoint and updated :math:`p`.

        Parameters
        ----------
        q : ndarray, shape (4,)
            Position before the step.
        p : ndarray, shape (4,)
            Covariant momentum before the step.
        h : float
            Sub-step size in proper time.

        Returns
        -------
        q_new : ndarray, shape (4,)
            Position after the sub-step.
        p_new : ndarray, shape (4,)
            Covariant momentum after the sub-step.

        Examples
        --------
        >>> import numpy as np
        >>> # q_new, p_new = integrator._verlet_step(q, p, h)
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
        """
        One full Yoshida sixth-order step as seven Verlet sub-steps.

        Parameters
        ----------
        q : ndarray, shape (4,)
            Position at the start of the step.
        p : ndarray, shape (4,)
            Covariant momentum at the start of the step.
        h : float
            Total step size (partitioned by internal composition weights).

        Returns
        -------
        q : ndarray, shape (4,)
            Position at the end of the step.
        p : ndarray, shape (4,)
            Covariant momentum at the end of the step.

        Examples
        --------
        >>> import numpy as np
        >>> # q, p = integrator._yoshida_step(q, p, h)
        """
        for w in _YOSHIDA6_WEIGHTS:
            q, p = self._verlet_step(q, p, w * h)
        return q, p

    # ------------------------------------------------------------------
    # C backend dispatch
    # ------------------------------------------------------------------

    def _integrate_c(self, state0, tau0, tau_end, n_steps, stride):
        """
        Integrate using the compiled C routine (Schwarzschild or Kerr).

        Parameters
        ----------
        state0 : ndarray, shape (8,)
            Initial state :math:`[q^\\mu, u^\\mu]`.
        tau0 : float
            Start of proper time.
        tau_end : float
            End of proper time.
        n_steps : int
            Number of Yoshida steps.
        stride : int
            Store output every ``stride`` steps.

        Returns
        -------
        result : simple namespace
            ``result.t`` — proper times; ``result.y`` — states with shape
            ``(8, N_out)``.

        Raises
        ------
        ValueError
            If the metric class is not supported by the C backend.

        Examples
        --------
        >>> import numpy as np
        >>> # integrator._integrate_c(s0, 0.0, 1.0, n_steps=4, stride=1)
        """
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
        Integrate the geodesic equations with a fixed step size.

        Uses the C backend when available for ``Schwarzschild`` or ``Kerr``
        metrics; otherwise integrates in pure Python via :meth:`_yoshida_step`.

        Parameters
        ----------
        state0 : array_like, shape (8,)
            Initial state :math:`[q^0, \\ldots, q^3, u^0, \\ldots, u^3]`.
        tau_span : tuple of float
            Proper-time interval ``(tau0, tau_end)``.
        n_steps : int
            Number of Yoshida steps; step size is
            :math:`h = (\\tau_{\\mathrm{end}} - \\tau_0) / n_{\\mathrm{steps}}`.
        stride : int, optional
            Store output every ``stride`` steps (default ``1``). The initial
            state is always stored.

        Returns
        -------
        result : simple namespace
            Object with attributes:

            - ``t`` : ndarray, shape (N_out,)
                Proper time at output points.
            - ``y`` : ndarray, shape (8, N_out)
                State :math:`[q^\\mu, u^\\mu]` at each output (columns).

        Raises
        ------
        ValueError
            If ``state0`` does not have shape ``(8,)``.

        Examples
        --------
        >>> import numpy as np
        >>> state0 = np.zeros(8)
        >>> # Yoshida6Integrator(metric).integrate(state0, (0.0, 1.0), n_steps=10)
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
    Convenience wrapper around :class:`Yoshida6Integrator` for geodesic output.

    When ``adaptive_output`` is False and ``tau_eval`` is provided, the stride
    for stored steps is chosen from ``n_steps`` and the length of ``tau_eval``.

    Parameters
    ----------
    metric : object
        Metric instance passed to :class:`Yoshida6Integrator`.
    state0 : array_like, shape (8,)
        Initial state :math:`[q^\\mu, u^\\mu]`.
    tau0 : float
        Start of proper time.
    tau_end : float
        End of proper time.
    n_steps : int
        Number of Yoshida steps in ``[tau0, tau_end]``.
    adaptive_output : bool, optional
        If True (default) or if ``tau_eval`` is None, output is stored every
        step (``stride=1``). If False and ``tau_eval`` is set, stride is
        ``max(1, n_steps // len(tau_eval))``.
    tau_eval : array_like or None, optional
        Optional array of evaluation times; only affects ``stride`` when
        ``adaptive_output`` is False.

    Returns
    -------
    tau : ndarray, shape (N,)
        Proper times at output points.
    state : ndarray, shape (N, 8)
        States :math:`[q^\\mu, u^\\mu]` at each output row.

    Examples
    --------
    >>> import numpy as np
    >>> state0 = np.zeros(8)
    >>> # tau, state = yoshida6_integrate_geodesic(m, state0, 0.0, 1.0, 10)
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
