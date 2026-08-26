"""
Latin-Hypercube sampling of Kerr orbital elements for a spin-determinability study.

The sampler draws :math:`N` triplets :math:`(a, e, i)` in geometric units
(:math:`G = c = M = 1`) from

- :math:`a \\in [a_g/10,\\; 1.3\\, a_g]`, uniform in :math:`\\ln a`;
- :math:`e \\in [0, 0.999]`, uniform in :math:`\\log_{10}(1-e)`;
- :math:`i \\in [0, \\pi/2]`, uniform in :math:`i` (parametric sweep, *not* an
  isotropic population).

Points whose periapsis :math:`r_p = a(1-e)` falls below the Kerr separatrix for
the given spin and inclination are rejected and re-drawn, so the returned sample
always has exactly ``N`` dynamically stable orbits.

Notes
-----
The separatrix is the locus of homoclinic orbits: the periapsis of a separatrix
orbit coincides with an *unstable* spherical orbit of the same
:math:`(E, L_z, Q)`. Writing the radial potential as

.. math::

    R(r) = (E^2-1)r^4 + 2r^3
           + \\left[a_\\bullet^2 (E^2-1) - L_z^2 - Q\\right] r^2
           + 2\\left[(L_z - a_\\bullet E)^2 + Q\\right] r
           - a_\\bullet^2 Q,

the separatrix is parametrised by the unstable spherical-orbit radius
:math:`r_u`: :math:`R` has a double root at :math:`r_u`, and its remaining real
root is the apoapsis :math:`r_a`, giving
:math:`e = (r_a - r_u)/(r_a + r_u)` and :math:`r_p^{\\rm sep}(e) = r_u`.

Equatorial-prograde (:math:`x \\equiv \\cos i = 1`) and polar (:math:`x = 0`)
branches are computed exactly from Bardeen's circular-orbit constants and from a
:math:`2\\times2` linear system in :math:`(E^2, Q)` respectively. Intermediate
inclinations are interpolated linearly in :math:`x`, which reproduces the exact
result at both ends and stays within a few percent in between (the separatrix is
smooth and nearly linear in :math:`x`; see Stein & Warburton 2020).

Examples
--------
>>> import numpy as np
>>> pts = sample_orbits_lhs(64, spin=0.9, seed=42)
>>> pts.shape
(64, 3)
>>> bool(np.all(pts[:, 0] * (1 - pts[:, 1]) > 0))
True
"""

import numpy

A_G = 25240.142762436135
"""float: Reference semi-major axis (geometric units) of the target orbit."""

A_MIN = A_G / 10.0
A_MAX = 1.3 * A_G
E_MIN = 0.0
E_MAX = 0.999
INC_MIN = 0.0
INC_MAX = numpy.pi / 2.0

_ONE_MINUS_E_MIN = 1.0 - E_MAX      # 1e-3
_ONE_MINUS_E_MAX = 1.0 - E_MIN      # 1.0


# ----------------------------------------------------------------------------
# Kerr separatrix
# ----------------------------------------------------------------------------
def _radial_quartic_roots(spin, E, Lz, Q):
    """
    Real roots of the Kerr radial potential :math:`R(r)`.

    Parameters
    ----------
    spin : float
        Dimensionless Kerr spin :math:`a_\\bullet \\in [0, 1)`.
    E, Lz, Q : float
        Specific energy, axial angular momentum, and Carter constant.

    Returns
    -------
    ndarray
        Real parts of the roots with negligible imaginary component, sorted
        ascending.

    Examples
    --------
    >>> r = _radial_quartic_roots(0.9, 0.97, 3.5, 0.0)
    >>> r.ndim
    1
    """
    coeffs = [
        E**2 - 1.0,
        2.0,
        spin**2 * (E**2 - 1.0) - Lz**2 - Q,
        2.0 * ((Lz - spin * E) ** 2 + Q),
        -spin**2 * Q,
    ]
    roots = numpy.roots(coeffs)
    real = roots[numpy.abs(roots.imag) < 1e-8 * numpy.maximum(1.0, numpy.abs(roots.real))]
    return numpy.sort(real.real)


def _eccentricity_from_double_root(spin, E, Lz, Q, r_u):
    """
    Eccentricity of the separatrix orbit whose periapsis is the double root ``r_u``.

    Parameters
    ----------
    spin : float
        Dimensionless Kerr spin.
    E, Lz, Q : float
        Constants of motion of the unstable spherical orbit at ``r_u``.
    r_u : float
        Radius of the unstable spherical orbit (double root of :math:`R`).

    Returns
    -------
    float
        Eccentricity :math:`(r_a - r_u)/(r_a + r_u)`, or ``nan`` when no apoapsis
        root :math:`r_a > r_u` exists.

    Examples
    --------
    >>> import numpy as np
    >>> np.isnan(_eccentricity_from_double_root(0.9, 1.5, 0.0, 0.0, 5.0))
    True
    """
    roots = _radial_quartic_roots(spin, E, Lz, Q)
    outer = roots[roots > r_u * (1.0 + 1e-6)]
    if outer.size == 0:
        return numpy.nan
    r_a = outer.max()
    return (r_a - r_u) / (r_a + r_u)


def _equatorial_branch(spin, n_grid=4000):
    """
    Separatrix branch for prograde equatorial orbits (:math:`x = \\cos i = 1`).

    Parameters
    ----------
    spin : float
        Dimensionless Kerr spin.
    n_grid : int, optional
        Number of unstable-circular-orbit radii scanned. Default ``4000``.

    Returns
    -------
    e_grid, rp_grid : ndarray
        Eccentricity and separatrix periapsis, sorted by increasing ``e_grid``.

    Notes
    -----
    Uses Bardeen's equatorial circular-orbit constants between the prograde
    photon orbit (:math:`e \\to 1`) and the ISCO (:math:`e = 0`).

    Examples
    --------
    >>> e, rp = _equatorial_branch(0.0, n_grid=200)
    >>> bool(abs(rp[0] - 6.0) < 1e-2)          # Schwarzschild: r_p(e=0) = 6M
    True
    """
    r_ph = 2.0 * (1.0 + numpy.cos((2.0 / 3.0) * numpy.arccos(-spin)))
    z1 = 1.0 + (1.0 - spin**2) ** (1.0 / 3.0) * (
        (1.0 + spin) ** (1.0 / 3.0) + (1.0 - spin) ** (1.0 / 3.0)
    )
    z2 = numpy.sqrt(3.0 * spin**2 + z1**2)
    r_isco = 3.0 + z2 - numpy.sqrt((3.0 - z1) * (3.0 + z1 + 2.0 * z2))

    radii = numpy.linspace(r_ph * (1.0 + 1e-9), r_isco, n_grid)
    es, rps = [], []
    for r in radii:
        root = r**1.5 - 3.0 * numpy.sqrt(r) + 2.0 * spin
        if root <= 0.0:
            continue
        den = r**0.75 * numpy.sqrt(root)
        E = (r**1.5 - 2.0 * numpy.sqrt(r) + spin) / den
        Lz = (r**2 - 2.0 * spin * numpy.sqrt(r) + spin**2) / den
        if not (0.0 < E < 1.0):
            continue
        e = _eccentricity_from_double_root(spin, E, Lz, 0.0, r)
        if numpy.isfinite(e) and 0.0 <= e < 1.0:
            es.append(e)
            rps.append(r)

    es = numpy.asarray(es)
    rps = numpy.asarray(rps)
    order = numpy.argsort(es)
    return es[order], rps[order]


def _polar_branch(spin, n_grid=4000, r_max=25.0):
    """
    Separatrix branch for polar orbits (:math:`x = \\cos i = 0`, :math:`L_z = 0`).

    Parameters
    ----------
    spin : float
        Dimensionless Kerr spin.
    n_grid : int, optional
        Number of spherical-orbit radii scanned. Default ``4000``.
    r_max : float, optional
        Upper end of the scan. Default ``25.0``.

    Returns
    -------
    e_grid, rp_grid : ndarray
        Eccentricity and separatrix periapsis, sorted by increasing ``e_grid``.

    Notes
    -----
    With :math:`L_z = 0` both :math:`R(r) = 0` and :math:`R'(r) = 0` are linear
    in :math:`(E^2, Q)`, so the unstable spherical-orbit constants follow from a
    :math:`2\\times2` solve at each radius.

    Examples
    --------
    >>> e, rp = _polar_branch(0.0, n_grid=200)
    >>> bool(abs(rp[0] - 6.0) < 1e-2)          # spin 0: polar == equatorial
    True
    """
    r_h = 1.0 + numpy.sqrt(max(1.0 - spin**2, 0.0))
    radii = numpy.linspace(r_h * (1.0 + 1e-6), r_max, n_grid)
    a2 = spin**2
    es, rps = [], []
    for r in radii:
        # R  = u (r^4 + a^2 r^2 + 2 a^2 r) + (-r^4 + 2 r^3 - a^2 r^2) + q (-r^2 + 2r - a^2)
        # R' = u (4r^3 + 2 a^2 r + 2 a^2) + (-4r^3 + 6 r^2 - 2 a^2 r) + q (-2r + 2)
        M = numpy.array(
            [
                [r**4 + a2 * r**2 + 2.0 * a2 * r, -(r**2) + 2.0 * r - a2],
                [4.0 * r**3 + 2.0 * a2 * r + 2.0 * a2, -2.0 * r + 2.0],
            ]
        )
        rhs = numpy.array(
            [
                r**4 - 2.0 * r**3 + a2 * r**2,
                4.0 * r**3 - 6.0 * r**2 + 2.0 * a2 * r,
            ]
        )
        if abs(numpy.linalg.det(M)) < 1e-12:
            continue
        u, q = numpy.linalg.solve(M, rhs)
        if not (0.0 < u < 1.0) or q <= 0.0:
            continue
        e = _eccentricity_from_double_root(spin, numpy.sqrt(u), 0.0, q, r)
        if numpy.isfinite(e) and 0.0 <= e < 1.0:
            es.append(e)
            rps.append(r)

    es = numpy.asarray(es)
    rps = numpy.asarray(rps)
    order = numpy.argsort(es)
    return es[order], rps[order]


class KerrSeparatrix:
    """
    Interpolator for the Kerr separatrix periapsis :math:`r_p^{\\rm sep}(e, i)`.

    Parameters
    ----------
    spin : float, optional
        Dimensionless Kerr spin :math:`a_\\bullet`. Default ``0.9``.
    n_grid : int, optional
        Radii sampled per branch when tabulating. Default ``4000``.

    Attributes
    ----------
    spin : float
        Spin used to build the tables.
    e_eq, rp_eq : ndarray
        Prograde-equatorial branch (:math:`x = 1`).
    e_po, rp_po : ndarray
        Polar branch (:math:`x = 0`).

    Examples
    --------
    >>> sep = KerrSeparatrix(spin=0.0, n_grid=400)
    >>> bool(abs(sep(0.0, 0.0) - 6.0) < 5e-2)      # Schwarzschild circular ISCO
    True
    >>> bool(abs(sep(0.5, 0.0) - 7.0) < 5e-2)      # r_p = (6 + 2e)/(1 + e) * ...
    True
    """

    def __init__(self, spin=0.9, n_grid=4000):
        self.spin = float(spin)
        self.e_eq, self.rp_eq = _equatorial_branch(self.spin, n_grid=n_grid)
        self.e_po, self.rp_po = _polar_branch(self.spin, n_grid=n_grid)

    def __call__(self, e, inc):
        """
        Separatrix periapsis for eccentricity ``e`` and inclination ``inc``.

        Parameters
        ----------
        e : array_like
            Eccentricity in :math:`[0, 1)`.
        inc : array_like
            Inclination in radians; :math:`x = \\cos(\\mathrm{inc})`.

        Returns
        -------
        ndarray or float
            Separatrix periapsis in units of :math:`M`, linearly interpolated in
            :math:`x` between the polar and prograde-equatorial branches.

        Examples
        --------
        >>> sep = KerrSeparatrix(spin=0.9, n_grid=400)
        >>> float(sep(0.9, 0.0)) > 1.0
        True
        """
        e = numpy.asarray(e, dtype=float)
        inc = numpy.asarray(inc, dtype=float)
        x = numpy.clip(numpy.cos(inc), 0.0, 1.0)
        rp_eq = numpy.interp(e, self.e_eq, self.rp_eq)
        rp_po = numpy.interp(e, self.e_po, self.rp_po)
        return rp_po + (rp_eq - rp_po) * x


# ----------------------------------------------------------------------------
# Sampling
# ----------------------------------------------------------------------------
def _unit_to_elements(unit):
    """
    Map unit-cube samples to :math:`(a, e, i)` through the sampling transforms.

    Parameters
    ----------
    unit : ndarray, shape (n, 3)
        Points in :math:`[0, 1)^3` ordered as
        :math:`(\\ln a, \\log_{10}(1-e), i)`.

    Returns
    -------
    ndarray, shape (n, 3)
        Columns ``[a, e, i]`` with ``a`` in geometric units and ``i`` in radians.

    Examples
    --------
    >>> import numpy as np
    >>> out = _unit_to_elements(np.array([[0.0, 1.0, 0.0]]))
    >>> bool(abs(out[0, 0] - A_MIN) < 1e-6 and abs(out[0, 1]) < 1e-12)
    True
    """
    unit = numpy.atleast_2d(numpy.asarray(unit, dtype=float))

    ln_a = numpy.log(A_MIN) + unit[:, 0] * (numpy.log(A_MAX) - numpy.log(A_MIN))
    a = numpy.exp(ln_a)

    log_ome = numpy.log10(_ONE_MINUS_E_MIN) + unit[:, 1] * (
        numpy.log10(_ONE_MINUS_E_MAX) - numpy.log10(_ONE_MINUS_E_MIN)
    )
    e = 1.0 - 10.0**log_ome

    inc = INC_MIN + unit[:, 2] * (INC_MAX - INC_MIN)

    return numpy.column_stack([a, e, inc])


def sample_orbits_lhs(
    n,
    spin=0.9,
    seed=20260731,
    safety=1.0,
    rp_floor=10.0,
    use_separatrix=True,
    max_batches=200,
    return_stats=False,
):
    """
    Latin-Hypercube sample of ``n`` stable Kerr orbits :math:`(a, e, i)`.

    Sampling is performed with :class:`scipy.stats.qmc.LatinHypercube` in the
    transformed coordinates :math:`(\\ln a, \\log_{10}(1-e), i)`; the transform
    is then inverted. Points failing the stability cut are discarded and the
    remaining slots refilled with fresh LHS batches until ``n`` valid orbits are
    collected.

    Parameters
    ----------
    n : int
        Number of valid orbits requested.
    spin : float, optional
        Dimensionless Kerr spin used by the separatrix. Default ``0.9``.
    seed : int, optional
        Seed for reproducibility. Default ``20260731``.
    safety : float, optional
        Multiplicative margin on the separatrix; a point is kept when
        :math:`r_p > \\mathrm{safety} \\times r_p^{\\rm sep}`. Default ``1.0``.
    rp_floor : float, optional
        Conservative absolute floor on the periapsis, in units of :math:`M`.
        Applied together with (or instead of) the separatrix. Default ``10.0``.
    use_separatrix : bool, optional
        If ``False``, only the ``rp_floor`` cut is applied. Default ``True``.
    max_batches : int, optional
        Maximum number of LHS batches drawn before giving up. Default ``200``.
    return_stats : bool, optional
        If ``True``, also return a dictionary with acceptance diagnostics.

    Returns
    -------
    ndarray, shape (n, 3)
        Columns ``[a, e, i]``; ``a`` in geometric units, ``i`` in radians.
    dict, optional
        Returned only when ``return_stats`` is ``True``: keys ``n_drawn``,
        ``n_accepted``, ``acceptance``, ``n_batches``, ``spin``, ``seed``.

    Raises
    ------
    RuntimeError
        If ``n`` valid orbits are not collected within ``max_batches`` batches.

    Notes
    -----
    Rejection breaks the strict one-point-per-stratum property of a single LHS
    design; the result is a stratified sample restricted to the stable region,
    which is the intended behaviour here.

    Examples
    --------
    >>> pts, stats = sample_orbits_lhs(32, seed=1, return_stats=True)
    >>> pts.shape, stats["n_accepted"]
    ((32, 3), 32)
    """
    from scipy.stats import qmc

    sep = KerrSeparatrix(spin) if use_separatrix else None
    rng = numpy.random.default_rng(seed)

    kept = []
    n_kept = 0
    n_drawn = 0
    n_batches = 0

    while n_kept < n and n_batches < max_batches:
        missing = n - n_kept
        # Over-draw on refills so the loop converges in a couple of batches.
        n_batch = missing if n_batches == 0 else max(int(numpy.ceil(missing * 1.5)), 16)
        sampler = qmc.LatinHypercube(d=3, seed=rng.integers(0, 2**32 - 1))
        cand = _unit_to_elements(sampler.random(n_batch))
        n_drawn += n_batch
        n_batches += 1

        a, e, inc = cand[:, 0], cand[:, 1], cand[:, 2]
        r_p = a * (1.0 - e)

        mask = r_p > rp_floor
        if sep is not None:
            mask &= r_p > safety * sep(e, inc)

        good = cand[mask]
        if good.shape[0] > missing:
            good = good[:missing]
        if good.size:
            kept.append(good)
            n_kept += good.shape[0]

    if n_kept < n:
        raise RuntimeError(
            f"only {n_kept}/{n} valid orbits after {n_batches} batches; "
            "relax rp_floor/safety or raise max_batches"
        )

    points = numpy.vstack(kept)[:n]
    if not return_stats:
        return points

    stats = {
        "n_drawn": n_drawn,
        "n_accepted": int(points.shape[0]),
        "acceptance": points.shape[0] / n_drawn,
        "n_batches": n_batches,
        "spin": spin,
        "seed": seed,
    }
    return points, stats


def save_sample(points, path="orbit_sample", fmt="npy"):
    """
    Persist an ``(N, 3)`` sample to ``.npy`` and/or ``.csv``.

    Parameters
    ----------
    points : ndarray, shape (N, 3)
        Columns ``[a, e, i]``.
    path : str, optional
        Path without extension. Default ``"orbit_sample"``.
    fmt : {'npy', 'csv', 'both'}, optional
        Output format. Default ``'npy'``.

    Returns
    -------
    list of str
        Paths actually written.

    Examples
    --------
    >>> import numpy as np, tempfile, os
    >>> p = np.zeros((2, 3))
    >>> files = save_sample(p, os.path.join(tempfile.mkdtemp(), "s"), fmt="both")
    >>> len(files)
    2
    """
    points = numpy.asarray(points, dtype=float)
    written = []
    if fmt in ("npy", "both"):
        numpy.save(f"{path}.npy", points)
        written.append(f"{path}.npy")
    if fmt in ("csv", "both"):
        numpy.savetxt(
            f"{path}.csv", points, delimiter=",", header="a,e,inc_rad", comments=""
        )
        written.append(f"{path}.csv")
    return written


def plot_diagnostic(points, spin=0.9, rp_floor=10.0, ax=None):
    """
    Diagnostic scatter of the sample in the :math:`(a, 1-e)` plane.

    Both axes are logarithmic, points are coloured by inclination, and the
    constant-periapsis threshold :math:`a(1-e) = r_p^{\\rm floor}` appears as a
    straight line, together with the equatorial and polar separatrix limits.

    Parameters
    ----------
    points : ndarray, shape (N, 3)
        Columns ``[a, e, i]`` with ``i`` in radians.
    spin : float, optional
        Spin used to draw the separatrix reference curves. Default ``0.9``.
    rp_floor : float, optional
        Periapsis threshold drawn as a solid line. Default ``10.0``.
    ax : matplotlib.axes.Axes, optional
        Target axes; created if omitted.

    Returns
    -------
    matplotlib.axes.Axes
        The axes containing the plot.

    Examples
    --------
    >>> import matplotlib
    >>> matplotlib.use("Agg")
    >>> pts = sample_orbits_lhs(16, seed=3)
    >>> ax = plot_diagnostic(pts)
    >>> ax.get_xscale()
    'log'
    """
    import matplotlib.pyplot as plt

    points = numpy.asarray(points, dtype=float)
    a, e, inc = points[:, 0], points[:, 1], points[:, 2]

    if ax is None:
        _, ax = plt.subplots(figsize=(7.0, 5.0))

    sc = ax.scatter(
        a, 1.0 - e, c=numpy.degrees(inc), cmap="viridis", s=14, alpha=0.85,
        edgecolors="none",
    )
    cb = ax.figure.colorbar(sc, ax=ax)
    cb.set_label(r"$i$ [deg]")

    a_line = numpy.logspace(numpy.log10(A_MIN), numpy.log10(A_MAX), 200)
    ax.plot(a_line, rp_floor / a_line, "k-", lw=1.4,
            label=rf"$r_p = {rp_floor:g}\,M$")

    # Separatrix curves: parametrise by e, then a = r_sep(e, i) / (1 - e).
    sep = KerrSeparatrix(spin)
    ome_ref = numpy.logspace(numpy.log10(_ONE_MINUS_E_MIN), 0.0, 400)
    e_ref = 1.0 - ome_ref
    for inc_ref, style, lbl in ((0.0, "--", "separatrix (equatorial)"),
                                (numpy.pi / 2, ":", "separatrix (polar)")):
        rp_sep = sep(e_ref, numpy.full_like(e_ref, inc_ref))
        a_sep = rp_sep / ome_ref
        ax.plot(a_sep, ome_ref, style, color="crimson", lw=1.2, label=lbl)
    ax.set_xlim(A_MIN * 0.8, A_MAX * 1.25)
    ax.set_ylim(_ONE_MINUS_E_MIN * 0.7, 1.4)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$a$ [$M$]")
    ax.set_ylabel(r"$1 - e$")
    ax.set_title(rf"LHS sample, $N={len(a)}$, $a_\bullet={spin}$")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax.tick_params(direction="in", which="both", top=True, right=True)
    return ax
