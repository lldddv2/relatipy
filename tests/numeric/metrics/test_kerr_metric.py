"""
Tests for the numeric Kerr metric and related geodesic quantities.

This module compares ``relatipy``'s Kerr implementation with EinsteinPy's
covariant metric in Boyer–Lindquist coordinates, checks Christoffel symbol
symmetry and dimensionless versus SI consistency, and validates timelike
normalization and conserved quantities along numerically integrated geodesics.

Three independent initial-condition sets (``xs_i``, ``vs_i``, masses ``M_i``,
and dimensionless spins ``a_i``) are built at import time using helpers from
``initial_conditions``.

Notes
-----
Full trajectory comparison against EinsteinPy geodesics is not performed here
because of convention differences in initial data and four-velocity between
packages.

Examples
--------
Run this test module with pytest from the repository root:

.. code-block:: console

    python -m pytest tests/numeric/metrics/test_kerr_metric.py -v
"""

import numpy as np
import astropy.units as u
from einsteinpy.metric import Kerr as ep_Kerr
from einsteinpy.coordinates import BoyerLindquistDifferential
from relatipy.numeric.metrics import Kerr as rp_Kerr
from relatipy.numeric.coordinates import BoyerLindquist
from einsteinpy.geodesic import Timelike
from initial_conditions import (
    M_1,
    M_2,
    M_3,
    xs_1,
    vs_1,
    xs_2,
    vs_2,
    xs_3,
    vs_3,
    position_ep_1,
    position_ep_2,
    position_ep_3,
    momentum_ep_1,
    momentum_ep_2,
    momentum_ep_3,
)

# CI 1
a_1 = 0.5  # dimensionless spin parameter

initial_conditions_1 = BoyerLindquistDifferential(
    t=xs_1[0],
    r=xs_1[1],
    theta=xs_1[2],
    phi=xs_1[3],
    v_r=vs_1[0],
    v_th=vs_1[1],
    v_p=vs_1[2],
)

initial_conditions_1_rp = BoyerLindquist(xs_1, vs_1, a=a_1)
x_vec_1 = np.array(initial_conditions_1.position())

# CI 2
a_2 = 0.9

initial_conditions_2 = BoyerLindquistDifferential(
    t=xs_2[0],
    r=xs_2[1],
    theta=xs_2[2],
    phi=xs_2[3],
    v_r=vs_2[0],
    v_th=vs_2[1],
    v_p=vs_2[2],
)

initial_conditions_2_rp = BoyerLindquist(xs_2, vs_2, a=a_2)
x_vec_2 = np.array(initial_conditions_2.position())

# CI 3
a_3 = 0.7

initial_conditions_3 = BoyerLindquistDifferential(
    t=xs_3[0],
    r=xs_3[1],
    theta=xs_3[2],
    phi=xs_3[3],
    v_r=vs_3[0],
    v_th=vs_3[1],
    v_p=vs_3[2],
)

initial_conditions_3_rp = BoyerLindquist(xs_3, vs_3, a=a_3)
x_vec_3 = np.array(initial_conditions_3.position())


class TestKerrMetric:
    """
    Regression tests for the Kerr metric and Kerr geodesic invariants.

    Uses three fixed initial-condition sets defined at module level together
    with EinsteinPy ``BoyerLindquistDifferential`` and relatipy
    ``BoyerLindquist`` state objects.

    Examples
    --------
    .. code-block:: console

        python -m pytest tests/numeric/metrics/test_kerr_metric.py::TestKerrMetric -v
    """

    def test_kerr_metric(self):
        """
        Match covariant Kerr metric components to EinsteinPy.

        For each initial-condition set, compares the full :math:`4 \\times 4`
        covariant metric from EinsteinPy's ``metric_covariant`` at the position
        vector with ``relatipy.numeric.metrics.Kerr.metric(..., dimensionless=False)``.

        Notes
        -----
        Uses SI units for the relatipy call (``dimensionless=False``) to align
        with EinsteinPy's coordinate scaling for the same masses and spins.

        Examples
        --------
        .. code-block:: console

            python -m pytest tests/numeric/metrics/test_kerr_metric.py::TestKerrMetric::test_kerr_metric -v
        """
        # CI 1
        kerr = ep_Kerr(coords=initial_conditions_1, M=M_1, a=a_1 * u.one)
        g = kerr.metric_covariant(x_vec_1)
        g_rp = rp_Kerr(M_1, a_1).metric(xs_1, dimensionless=False)
        assert np.isclose(g, g_rp).all()

        # CI 2
        kerr = ep_Kerr(coords=initial_conditions_2, M=M_2, a=a_2 * u.one)
        g = kerr.metric_covariant(x_vec_2)
        g_rp = rp_Kerr(M_2, a_2).metric(xs_2, dimensionless=False)
        assert np.isclose(g, g_rp).all()

        # CI 3
        kerr = ep_Kerr(coords=initial_conditions_3, M=M_3, a=a_3 * u.one)
        g = kerr.metric_covariant(x_vec_3)
        g_rp = rp_Kerr(M_3, a_3).metric(xs_3, dimensionless=False)
        assert np.isclose(g, g_rp).all()

    def test_kerr_christoffel_symbols(self):
        """
        Check Christoffel symmetry in the last two indices.

        Verifies that the Levi-Civita connection satisfies
        :math:`\\Gamma^{\\rho}{}_{\\mu\\nu} = \\Gamma^{\\rho}{}_{\\nu\\mu}` by
        comparing each array with its transpose over the last two axes, for
        both ``dimensionless=True`` and ``dimensionless=False``.

        Notes
        -----
        The relatipy Kerr implementation returns Christoffel components as a
        three-index array; symmetry is checked independently for dimensionless
        and SI-normalized outputs.

        Examples
        --------
        .. code-block:: console

            python -m pytest tests/numeric/metrics/test_kerr_metric.py::TestKerrMetric::test_kerr_christoffel_symbols -v
        """
        # CI 1
        kerr = rp_Kerr(M_1, a_1)
        G_dimless = kerr.get_christoffel_symbols(xs_1, dimensionless=True)
        G_si = kerr.get_christoffel_symbols(xs_1, dimensionless=False)
        assert np.isclose(G_dimless, G_dimless.transpose(0, 2, 1)).all(), (
            "CI1: Christoffel is not symmetric in the lower pair of indices (dimensionless)"
        )
        assert np.isclose(G_si, G_si.transpose(0, 2, 1)).all(), (
            "CI1: Christoffel is not symmetric in the lower pair of indices (SI)"
        )

        # CI 2
        kerr = rp_Kerr(M_2, a_2)
        G_dimless = kerr.get_christoffel_symbols(xs_2, dimensionless=True)
        G_si = kerr.get_christoffel_symbols(xs_2, dimensionless=False)
        assert np.isclose(G_dimless, G_dimless.transpose(0, 2, 1)).all(), (
            "CI2: Christoffel is not symmetric in the lower pair of indices (dimensionless)"
        )
        assert np.isclose(G_si, G_si.transpose(0, 2, 1)).all(), (
            "CI2: Christoffel is not symmetric in the lower pair of indices (SI)"
        )

        # CI 3
        kerr = rp_Kerr(M_3, a_3)
        G_dimless = kerr.get_christoffel_symbols(xs_3, dimensionless=True)
        G_si = kerr.get_christoffel_symbols(xs_3, dimensionless=False)
        assert np.isclose(G_dimless, G_dimless.transpose(0, 2, 1)).all(), (
            "CI3: Christoffel is not symmetric in the lower pair of indices (dimensionless)"
        )
        assert np.isclose(G_si, G_si.transpose(0, 2, 1)).all(), (
            "CI3: Christoffel is not symmetric in the lower pair of indices (SI)"
        )

    # Full geodesic comparison with EinsteinPy is omitted due to convention
    # differences in initial conditions and four-velocity.

    def test_kerr_ds_dtau(self):
        """
        Assert unit timelike normalization :math:`\\mathrm{d}s/\\mathrm{d}\\tau = 1`.

        Integrates the geodesic for each initial-condition set over a fixed
        proper-time grid and checks that the path's ``_get_ds_dtau`` values are
        unity everywhere (natural units).

        Notes
        -----
        This encodes the timelike constraint after numerical integration.

        Examples
        --------
        .. code-block:: console

            python -m pytest tests/numeric/metrics/test_kerr_metric.py::TestKerrMetric::test_kerr_ds_dtau -v
        """
        # CI 1
        taus_1 = np.linspace(0, 100, 100)
        kerr = rp_Kerr(M_1, a_1)
        path = kerr.geodesic.get_path(initial_conditions_1_rp, taus_1)
        ds_dtau = path._get_ds_dtau(kerr)
        assert np.isclose(ds_dtau, 1).all(), "ds/dtau deviates from 1 along the path"

        # CI 2
        taus_2 = np.linspace(0, 100, 100)
        kerr = rp_Kerr(M_2, a_2)
        path = kerr.geodesic.get_path(initial_conditions_2_rp, taus_2)
        ds_dtau = path._get_ds_dtau(kerr)
        assert np.isclose(ds_dtau, 1).all(), "ds/dtau deviates from 1 along the path"

        # CI 3
        taus_3 = np.linspace(0, 100, 100)
        kerr = rp_Kerr(M_3, a_3)
        path = kerr.geodesic.get_path(initial_conditions_3_rp, taus_3)
        ds_dtau = path._get_ds_dtau(kerr)
        assert np.isclose(ds_dtau, 1).all(), "ds/dtau deviates from 1 along the path"

    def test_kerr_E(self):
        """
        Energy constant :math:`E` is conserved along the integrated geodesic.

        Examples
        --------
        .. code-block:: console

            python -m pytest tests/numeric/metrics/test_kerr_metric.py::TestKerrMetric::test_kerr_E -v
        """
        # CI 1
        taus_1 = np.linspace(0, 100, 100)
        kerr = rp_Kerr(M_1, a_1)
        path = kerr.geodesic.get_path(initial_conditions_1_rp, taus_1)
        E = path._get_E(kerr)
        assert np.isclose(E, E[0]).all(), "E is not constant over the trajectory"

        # CI 2
        taus_2 = np.linspace(0, 100, 100)
        kerr = rp_Kerr(M_2, a_2)
        path = kerr.geodesic.get_path(initial_conditions_2_rp, taus_2)
        E = path._get_E(kerr)
        assert np.isclose(E, E[0]).all(), "E is not constant over the trajectory"

        # CI 3
        taus_3 = np.linspace(0, 100, 100)
        kerr = rp_Kerr(M_3, a_3)
        path = kerr.geodesic.get_path(initial_conditions_3_rp, taus_3)
        E = path._get_E(kerr)
        assert np.isclose(E, E[0]).all(), "E is not constant over the trajectory"

    def test_kerr_Lz(self):
        """
        Azimuthal angular momentum :math:`L_z` is conserved along the path.

        Examples
        --------
        .. code-block:: console

            python -m pytest tests/numeric/metrics/test_kerr_metric.py::TestKerrMetric::test_kerr_Lz -v
        """
        # CI 1
        taus_1 = np.linspace(0, 100, 100)
        kerr = rp_Kerr(M_1, a_1)
        path = kerr.geodesic.get_path(initial_conditions_1_rp, taus_1)
        Lz = path._get_Lz(kerr)
        assert np.isclose(Lz, Lz[0]).all(), "CI1: Lz is not constant along the trajectory"

        # CI 2
        taus_2 = np.linspace(0, 100, 100)
        kerr = rp_Kerr(M_2, a_2)
        path = kerr.geodesic.get_path(initial_conditions_2_rp, taus_2)
        Lz = path._get_Lz(kerr)
        assert np.isclose(Lz, Lz[0]).all(), "CI2: Lz is not constant along the trajectory"

        # CI 3
        taus_3 = np.linspace(0, 100, 100)
        kerr = rp_Kerr(M_3, a_3)
        path = kerr.geodesic.get_path(initial_conditions_3_rp, taus_3)
        Lz = path._get_Lz(kerr)
        assert np.isclose(Lz, Lz[0]).all(), "CI3: Lz is not constant along the trajectory"

    def test_kerr_Q(self):
        """
        Carter constant :math:`Q` is conserved along the integrated geodesic.

        Examples
        --------
        .. code-block:: console

            python -m pytest tests/numeric/metrics/test_kerr_metric.py::TestKerrMetric::test_kerr_Q -v
        """
        # CI 1
        taus_1 = np.linspace(0, 100, 100)
        kerr = rp_Kerr(M_1, a_1)
        path = kerr.geodesic.get_path(initial_conditions_1_rp, taus_1)
        Q = path._get_Q(kerr)
        assert np.isclose(Q, Q[0]).all(), "Q is not constant over the trajectory"

        # CI 2
        taus_2 = np.linspace(0, 100, 100)
        kerr = rp_Kerr(M_2, a_2)
        path = kerr.geodesic.get_path(initial_conditions_2_rp, taus_2)
        Q = path._get_Q(kerr)
        assert np.isclose(Q, Q[0]).all(), "Q is not constant over the trajectory"

        # CI 3
        taus_3 = np.linspace(0, 100, 100)
        kerr = rp_Kerr(M_3, a_3)
        path = kerr.geodesic.get_path(initial_conditions_3_rp, taus_3)
        Q = path._get_Q(kerr)
        assert np.isclose(Q, Q[0]).all(), "Q is not constant over the trajectory"
