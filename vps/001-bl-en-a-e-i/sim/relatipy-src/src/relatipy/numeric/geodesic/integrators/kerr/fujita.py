"""
Kerr geodesic integrator using Mino time parametrization (Fujita method).

This module wraps a compiled C implementation that integrates Kerr geodesics
in Boyer-Lindquist coordinates using the Mino time parameter :math:`\\lambda`,
defined by

.. math::

    \\frac{\\mathrm{d}\\tau}{\\mathrm{d}\\lambda} = \\Sigma(r,\\theta)
    = r^2 + a^2\\cos^2\\theta.

In Mino time the radial and polar equations of motion decouple:

.. math::

    \\left(\\frac{\\mathrm{d}r}{\\mathrm{d}\\lambda}\\right)^2 = R(r), \\qquad
    \\left(\\frac{\\mathrm{d}\\theta}{\\mathrm{d}\\lambda}\\right)^2 = \\Theta(\\theta),

where :math:`R` and :math:`\\Theta` are the standard Kerr radial and polar
potentials. The azimuthal and time equations are integrated as quadratures.

The shared library ``fujita_core`` (``fujita_core.so`` / ``fujita_core.dll``) is
loaded from this package directory at import time; if loading fails, the
public functions raise ``RuntimeError`` when called.

Notes
-----
Integration uses 4th-order Runge-Kutta in :math:`\\lambda`. The step size
:math:`\\mathrm{d}\\lambda` is estimated from the requested proper-time span
and the initial :math:`\\Sigma`.

References
----------
.. [1] Mino, Y. (2003). *Perturbative approach to an orbital evolution around
       a supermassive black hole*. Phys. Rev. D, 67, 084027.
.. [2] Drasco, S. & Hughes, S. A. (2004). Phys. Rev. D, 69, 044708.
.. [3] Fujita, R. & Hikida, W. (2009). *Analytical solutions of bound timelike
       geodesic orbits in Kerr spacetime*. Class. Quantum Grav., 26, 135002.

Examples
--------
>>> from relatipy.numeric.geodesic.integrators.kerr import fujita
>>> fujita._C_LIB is not None  # doctest: +SKIP
True
"""

import ctypes
import os
import sys

import numpy as np


def _find_fujita_core_path() -> str | None:
    """Path to the ctypes shared library next to this module, if present."""
    base = os.path.dirname(os.path.abspath(__file__))
    if sys.platform == "win32":
        names = ("fujita_core.dll", "fujita_core.so")
    else:
        names = ("fujita_core.so", "fujita_core.dylib")
    for name in names:
        p = os.path.join(base, name)
        if os.path.isfile(p):
            return p
    return None


# ---- Load C backend (fujita_core.{so,dll,...}) --------------------------------
_C_LIB = None
_FUJITA_LOAD_ERROR: str | None = None
try:
    _SO_PATH = _find_fujita_core_path()
    if _SO_PATH is None:
        raise FileNotFoundError("fujita_core shared library not found next to fujita.py")
    _lib = ctypes.CDLL(_SO_PATH)
    _dbl_p = ctypes.POINTER(ctypes.c_double)
    _int_p = ctypes.POINTER(ctypes.c_int)

    _lib.kerr_mino.restype = None
    _lib.kerr_mino.argtypes = [
        ctypes.c_double,  # Rs
        ctypes.c_double,  # a
        ctypes.c_double,  # E
        ctypes.c_double,  # Lz
        ctypes.c_double,  # Q
        _dbl_p,           # state0[8]
        ctypes.c_double,  # tau0
        ctypes.c_double,  # tau_end
        ctypes.c_int,     # n_steps
        ctypes.c_int,     # stride
        _dbl_p,           # out_tau
        _dbl_p,           # out_y  (N_out * 8, row-major)
        _int_p,           # n_out
    ]
    _C_LIB = _lib
except Exception as _e:
    _C_LIB = None
    _FUJITA_LOAD_ERROR = str(_e)


def _fujita_backend_error(what: str) -> RuntimeError:
    detail = f" ({_FUJITA_LOAD_ERROR})" if _FUJITA_LOAD_ERROR else ""
    return RuntimeError(
        f"C backend not loaded; cannot call {what}.{detail} "
        "Rebuild with: cc -O3 -march=native -shared -fPIC "
        "-o fujita_core.so fujita_core.c -lm"
    )


def _integrate_kerr_fujita(
    Rs, a, state0, tau_span, n_steps, E0, Lz0, Q0, stride=1
):
    """
    Integrate a Kerr geodesic using Mino-time RK4 (Fujita parametrization).

    The integration runs in Mino time :math:`\\lambda` and stores output in
    proper time :math:`\\tau`. The output state format matches all other Kerr
    integrators: :math:`[q^\\mu, u^\\mu]` in Boyer-Lindquist coordinates.

    Parameters
    ----------
    Rs : float
        Schwarzschild radius :math:`R_s = 2M` for the Kerr metric.
    a : float
        Spin parameter :math:`a` (Boyer-Lindquist).
    state0 : array_like, shape (8,)
        Initial state :math:`[t, r, \\theta, \\phi, u^t, u^r, u^\\theta, u^\\phi]`
        as contiguous ``float64``.
    tau_span : tuple of float
        Pair ``(tau0, tau_end)`` giving the proper-time interval.
    n_steps : int
        Number of RK4 steps in Mino time :math:`\\lambda`.
    E0 : float
        Conserved energy (per unit mass).
    Lz0 : float
        Conserved :math:`z`-component of angular momentum.
    Q0 : float
        Carter constant.
    stride : int, optional
        Store output every ``stride`` internal steps (default ``1``).

    Returns
    -------
    result : object
        Namespace with attributes:

        ``t`` : ndarray, shape (N_out,)
            Proper times at stored output points.
        ``y`` : ndarray, shape (8, N_out)
            State vectors :math:`[q^\\mu; u^\\mu]` at each stored time.

    Raises
    ------
    RuntimeError
        If the C shared library was not loaded at import time.

    Examples
    --------
    >>> import numpy as np
    >>> from relatipy.numeric.geodesic.integrators.kerr.fujita import _integrate_kerr_fujita
    >>> s0 = np.array([0., 10., 1.5708, 0., 1.0, 0., 0., 0.1], dtype=np.float64)
    >>> out = _integrate_kerr_fujita(2.0, 0.5, s0, (0.0, 100.0), 1000, 0.95, 2.0, 1.0)  # doctest: +SKIP
    >>> out.t.shape == out.y.shape[1:]  # doctest: +SKIP
    True
    """
    if _C_LIB is None:
        raise _fujita_backend_error("_integrate_kerr_fujita")

    N_alloc = n_steps // max(stride, 1) + 4
    out_tau = np.zeros(N_alloc, dtype=np.float64)
    out_y   = np.zeros(N_alloc * 8, dtype=np.float64)
    n_out   = ctypes.c_int(0)

    s0   = np.ascontiguousarray(state0, dtype=np.float64)
    _ptr = lambda arr: arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    _C_LIB.kerr_mino(
        ctypes.c_double(Rs),
        ctypes.c_double(a),
        ctypes.c_double(E0),
        ctypes.c_double(Lz0),
        ctypes.c_double(Q0),
        _ptr(s0),
        ctypes.c_double(tau_span[0]),
        ctypes.c_double(tau_span[1]),
        ctypes.c_int(n_steps),
        ctypes.c_int(stride),
        _ptr(out_tau),
        _ptr(out_y),
        ctypes.byref(n_out),
    )

    N = n_out.value

    class _Result:
        """
        Lightweight holder for Fujita integration output.

        Attributes
        ----------
        t : ndarray
            Proper times at output points, shape ``(N_out,)``.
        y : ndarray
            State vectors, shape ``(8, N_out)``.
        """

    result   = _Result()
    result.t = out_tau[:N].copy()
    result.y = out_y[: N * 8].reshape(N, 8).T.copy()  # (8, N)
    return result
