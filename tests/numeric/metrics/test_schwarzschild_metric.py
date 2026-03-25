"""
Tests for the numeric Schwarzschild metric and related geodesic consistency.

This module compares ``relatipy.numeric.metrics.Schwarzschild`` against
``einsteinpy`` for the covariant metric tensor and Christoffel symbols at
several initial-condition sets. It also checks timelike geodesic invariants
(:math:`ds/d\\tau`), conserved angular momentum :math:`L_z`, and conserved
energy :math:`E`, and verifies that integrating in spherical, Cartesian, or
orbital-element coordinates yields the same trajectory after conversion.

Notes
-----
Full geodesic comparison against ``einsteinpy`` is intentionally not included
because of differing conventions for initial data and four-velocity.

Examples
--------
Run this module's tests:

.. code-block:: bash

    pytest tests/numeric/metrics/test_schwarzschild_metric.py -v
"""

import numpy as np
import astropy.units as u
from einsteinpy.metric import Schwarzschild
from einsteinpy.coordinates import SphericalDifferential
from relatipy.numeric.metrics import Schwarzschild as rp_Schwarzschild
from relatipy.numeric.coordinates import Spherical
from relatipy.numeric.constants import _L_ref, _T_ref, _c_SI
from einsteinpy.geodesic import Timelike
from initial_conditions import (
    position_ep_1,
    momentum_ep_1,
    xs_1,
    vs_1,
    M_1,
    position_ep_2,
    momentum_ep_2,
    xs_2,
    vs_2,
    M_2,
    position_ep_3,
    momentum_ep_3,
    xs_3,
    vs_3,
    M_3,
)

initial_conditions_1 = SphericalDifferential(
    t=xs_1[0],
    r=xs_1[1],
    theta=xs_1[2],
    phi=xs_1[3],
    v_r=vs_1[0],
    v_th=vs_1[1],
    v_p=vs_1[2],
)
initial_conditions_1_rp = Spherical(xs_1, vs_1)
x_vec_1 = np.array(initial_conditions_1.position())

initial_conditions_2 = SphericalDifferential(
    t=xs_2[0],
    r=xs_2[1],
    theta=xs_2[2],
    phi=xs_2[3],
    v_r=vs_2[0],
    v_th=vs_2[1],
    v_p=vs_2[2],
)
initial_conditions_2_rp = Spherical(xs_2, vs_2)
x_vec_2 = np.array(initial_conditions_2.position())

initial_conditions_3 = SphericalDifferential(
    t=xs_3[0],
    r=xs_3[1],
    theta=xs_3[2],
    phi=xs_3[3],
    v_r=vs_3[0],
    v_th=vs_3[1],
    v_p=vs_3[2],
)
initial_conditions_3_rp = Spherical(xs_3, vs_3)
x_vec_3 = np.array(initial_conditions_3.position())


def _rp_traj_to_si(traj_rp):
    """
    Convert a relatipy trajectory from geometric units to SI length and time.

    Columns are spherical coordinates :math:`(t, r, \\theta, \\phi)`; time and
    radial distance are scaled by ``_T_ref`` and ``_L_ref`` respectively.
    Angular coordinates are dimensionless and unchanged.

    Parameters
    ----------
    traj_rp : ndarray, shape (4, N)
        Trajectory with rows ``[t_geom, r_geom, theta, phi]`` in geometric
        units for :math:`t` and :math:`r`.

    Returns
    -------
    ndarray, shape (4, N)
        Same layout with ``t`` in seconds and ``r`` in meters.

    Examples
    --------
    >>> import numpy as np
    >>> from relatipy.numeric.constants import _L_ref, _T_ref
    >>> traj = np.zeros((4, 2))
    >>> traj[0, :] = 1.0
    >>> traj[1, :] = 2.0
    >>> out = _rp_traj_to_si(traj)
    >>> np.allclose(out[0], _T_ref) and np.allclose(out[1], 2 * _L_ref)
    True
    """
    out = traj_rp.copy().astype(float)
    out[0] *= _T_ref
    out[1] *= _L_ref
    return out


class TestSchwarzschildMetric:
    """
    Regression tests for Schwarzschild geometry and timelike geodesics.

    Uses three independent initial-condition sets (masses ``M_1``, ``M_2``,
    ``M_3``) defined in ``initial_conditions`` to cross-check the covariant
    metric and connection against ``einsteinpy`` and to validate internal
    geodesic invariants and coordinate representations.

    Examples
    --------
    .. code-block:: bash

        pytest tests/numeric/metrics/test_schwarzschild_metric.py::TestSchwarzschildMetric -v
    """

    def test_schwarzschild_metric(self):
        """
        Covariant metric tensor matches ``einsteinpy`` for all IC sets.

        Compares ``Schwarzschild.metric_covariant`` from ``einsteinpy`` with
        ``relatipy.numeric.metrics.Schwarzschild.metric(..., dimensionless=False)``
        at the same spherical coordinates for three initial conditions.

        Parameters
        ----------
        self : TestSchwarzschildMetric
            Test instance (unused).

        Examples
        --------
        .. code-block:: bash

            pytest tests/numeric/metrics/test_schwarzschild_metric.py::TestSchwarzschildMetric::test_schwarzschild_metric -v
        """
        # IC set 1
        sch = Schwarzschild(coords=initial_conditions_1, M=M_1)
        g_ep = sch.metric_covariant(x_vec_1)
        g_rp = rp_Schwarzschild(M_1).metric(xs_1, dimensionless=False)
        assert np.isclose(g_ep, g_rp).all()

        # IC set 2
        sch = Schwarzschild(coords=initial_conditions_2, M=M_2)
        g_ep = sch.metric_covariant(x_vec_2)
        g_rp = rp_Schwarzschild(M_2).metric(xs_2, dimensionless=False)
        assert np.isclose(g_ep, g_rp).all()

        # IC set 3
        sch = Schwarzschild(coords=initial_conditions_3, M=M_3)
        g_ep = sch.metric_covariant(x_vec_3)
        g_rp = rp_Schwarzschild(M_3).metric(xs_3, dimensionless=False)
        assert np.isclose(g_ep, g_rp).all()

    def test_schwarzschild_christoffel_symbols(self):
        """
        Christoffel symbols match ``einsteinpy`` for all IC sets.

        Compares ``einsteinpy`` ``christoffels`` with
        ``get_christoffel_symbols(..., dimensionless=False)`` from relatipy.

        Parameters
        ----------
        self : TestSchwarzschildMetric
            Test instance (unused).

        Examples
        --------
        .. code-block:: bash

            pytest tests/numeric/metrics/test_schwarzschild_metric.py::TestSchwarzschildMetric::test_schwarzschild_christoffel_symbols -v
        """
        # IC set 1
        sch = Schwarzschild(coords=initial_conditions_1, M=M_1)
        christoffel_ep = sch.christoffels(x_vec_1)
        christoffel_rp = rp_Schwarzschild(M_1).get_christoffel_symbols(xs_1, dimensionless=False)
        assert np.isclose(christoffel_ep, christoffel_rp).all()

        # IC set 2
        sch = Schwarzschild(coords=initial_conditions_2, M=M_2)
        christoffel_ep = sch.christoffels(x_vec_2)
        christoffel_rp = rp_Schwarzschild(M_2).get_christoffel_symbols(xs_2, dimensionless=False)
        assert np.isclose(christoffel_ep, christoffel_rp).all()

        # IC set 3
        sch = Schwarzschild(coords=initial_conditions_3, M=M_3)
        christoffel_ep = sch.christoffels(x_vec_3)
        christoffel_rp = rp_Schwarzschild(M_3).get_christoffel_symbols(xs_3, dimensionless=False)
        assert np.isclose(christoffel_ep, christoffel_rp).all()

    def test_schwarzschild_ds_dtau(self):
        """
        Proper-time normalization :math:`ds/d\\tau = 1` along timelike geodesics.

        Integrates the geodesic for each IC set and asserts the internal
        ``_get_ds_dtau`` diagnostic is unity along the sampled affine parameter.

        Parameters
        ----------
        self : TestSchwarzschildMetric
            Test instance (unused).

        Examples
        --------
        .. code-block:: bash

            pytest tests/numeric/metrics/test_schwarzschild_metric.py::TestSchwarzschildMetric::test_schwarzschild_ds_dtau -v
        """
        # IC set 1
        taus_1 = np.linspace(0, 100, 100)
        sch = rp_Schwarzschild(M_1)
        path = sch.geodesic.get_path(initial_conditions_1_rp, taus_1)
        ds_dtau = path._get_ds_dtau(sch)
        assert np.isclose(ds_dtau, 1).all(), "IC1: ds/dtau is not 1"

        # IC set 2
        taus_2 = np.linspace(0, 100, 100)
        sch = rp_Schwarzschild(M_2)
        path = sch.geodesic.get_path(initial_conditions_2_rp, taus_2)
        ds_dtau = path._get_ds_dtau(sch)
        assert np.isclose(ds_dtau, 1).all(), "IC2: ds/dtau is not 1"

        # IC set 3
        taus_3 = np.linspace(0, 100, 100)
        sch = rp_Schwarzschild(M_3)
        path = sch.geodesic.get_path(initial_conditions_3_rp, taus_3)
        ds_dtau = path._get_ds_dtau(sch)
        assert np.isclose(ds_dtau, 1).all(), "IC3: ds/dtau is not 1"

    def test_schwarzschild_Lz(self):
        """
        Conserved azimuthal angular momentum :math:`L_z` along the geodesic.

        Asserts ``_get_Lz()`` is constant in proper time for each IC set.

        Parameters
        ----------
        self : TestSchwarzschildMetric
            Test instance (unused).

        Examples
        --------
        .. code-block:: bash

            pytest tests/numeric/metrics/test_schwarzschild_metric.py::TestSchwarzschildMetric::test_schwarzschild_Lz -v
        """
        # IC set 1
        taus_1 = np.linspace(0, 100, 100)
        sch = rp_Schwarzschild(M_1)
        path = sch.geodesic.get_path(initial_conditions_1_rp, taus_1)
        Lz = path._get_Lz()
        assert np.allclose(Lz, Lz[0]), "IC1: Lz is not constant"

        # IC set 2
        taus_2 = np.linspace(0, 100, 100)
        sch = rp_Schwarzschild(M_2)
        path = sch.geodesic.get_path(initial_conditions_2_rp, taus_2)
        Lz = path._get_Lz()
        assert np.allclose(Lz, Lz[0]), "IC2: Lz is not constant"

        # IC set 3
        taus_3 = np.linspace(0, 100, 100)
        sch = rp_Schwarzschild(M_3)
        path = sch.geodesic.get_path(initial_conditions_3_rp, taus_3)
        Lz = path._get_Lz()
        assert np.allclose(Lz, Lz[0]), "IC3: Lz is not constant"

    def test_schwarzschild_E(self):
        """
        Conserved energy :math:`E` along the timelike geodesic.

        Asserts ``_get_E`` is constant for each IC set when evaluated with the
        same ``Schwarzschild`` instance.

        Parameters
        ----------
        self : TestSchwarzschildMetric
            Test instance (unused).

        Examples
        --------
        .. code-block:: bash

            pytest tests/numeric/metrics/test_schwarzschild_metric.py::TestSchwarzschildMetric::test_schwarzschild_E -v
        """
        # IC set 1
        taus_1 = np.linspace(0, 100, 100)
        sch = rp_Schwarzschild(M_1)
        path = sch.geodesic.get_path(initial_conditions_1_rp, taus_1)
        E = path._get_E(sch)
        assert np.allclose(E, E[0]), "IC1: E is not constant"

        # IC set 2
        taus_2 = np.linspace(0, 100, 100)
        sch = rp_Schwarzschild(M_2)
        path = sch.geodesic.get_path(initial_conditions_2_rp, taus_2)
        E = path._get_E(sch)
        assert np.allclose(E, E[0]), "IC2: E is not constant"

        # IC set 3
        taus_3 = np.linspace(0, 100, 100)
        sch = rp_Schwarzschild(M_3)
        path = sch.geodesic.get_path(initial_conditions_3_rp, taus_3)
        E = path._get_E(sch)
        assert np.allclose(E, E[0]), "IC3: E is not constant"

    def test_same_trajectory_with_different_coordinates(self):
        """
        Geodesic state matches when integrated in spherical vs Cartesian form.

        For each IC set, compares seven trajectory components after converting
        between spherical and Cartesian representations in both directions.

        Parameters
        ----------
        self : TestSchwarzschildMetric
            Test instance (unused).

        Examples
        --------
        .. code-block:: bash

            pytest tests/numeric/metrics/test_schwarzschild_metric.py::TestSchwarzschildMetric::test_same_trajectory_with_different_coordinates -v
        """
        # IC set 1
        initial_conditions_1_cartesian = initial_conditions_1_rp.convert_to("Cartesian")
        taus_1 = np.linspace(0, 100, 100)
        sch = rp_Schwarzschild(M_1)

        path_spherical = sch.geodesic.get_path(initial_conditions_1_rp, taus_1)
        path_cartesian = sch.geodesic.get_path(initial_conditions_1_cartesian, taus_1).convert_to("Spherical")
        for i in range(7):
            assert np.isclose(path_spherical[i], path_cartesian[i]).all(), (
                f"IC1 Sph->Cart: component {i} does not match"
            )

        path_cartesian = sch.geodesic.get_path(initial_conditions_1_cartesian, taus_1)
        path_spherical = sch.geodesic.get_path(initial_conditions_1_rp, taus_1).convert_to("Cartesian")
        for i in range(7):
            assert np.isclose(path_cartesian[i], path_spherical[i]).all(), (
                f"IC1 Cart->Sph: component {i} does not match"
            )

        # IC set 2
        initial_conditions_2_cartesian = initial_conditions_2_rp.convert_to("Cartesian")
        taus_2 = np.linspace(0, 100, 100)
        sch = rp_Schwarzschild(M_2)

        path_spherical = sch.geodesic.get_path(initial_conditions_2_rp, taus_2)
        path_cartesian = sch.geodesic.get_path(initial_conditions_2_cartesian, taus_2).convert_to("Spherical")
        for i in range(7):
            assert np.isclose(path_spherical[i], path_cartesian[i]).all(), (
                f"IC2 Sph->Cart: component {i} does not match"
            )

        path_cartesian = sch.geodesic.get_path(initial_conditions_2_cartesian, taus_2)
        path_spherical = sch.geodesic.get_path(initial_conditions_2_rp, taus_2).convert_to("Cartesian")
        for i in range(7):
            assert np.isclose(path_cartesian[i], path_spherical[i]).all(), (
                f"IC2 Cart->Sph: component {i} does not match"
            )

        # IC set 3
        initial_conditions_3_cartesian = initial_conditions_3_rp.convert_to("Cartesian")
        taus_3 = np.linspace(0, 100, 100)
        sch = rp_Schwarzschild(M_3)

        path_spherical = sch.geodesic.get_path(initial_conditions_3_rp, taus_3)
        path_cartesian = sch.geodesic.get_path(initial_conditions_3_cartesian, taus_3).convert_to("Spherical")
        for i in range(7):
            assert np.isclose(path_spherical[i], path_cartesian[i]).all(), (
                f"IC3 Sph->Cart: component {i} does not match"
            )

        path_cartesian = sch.geodesic.get_path(initial_conditions_3_cartesian, taus_3)
        path_spherical = sch.geodesic.get_path(initial_conditions_3_rp, taus_3).convert_to("Cartesian")
        for i in range(7):
            assert np.isclose(path_cartesian[i], path_spherical[i]).all(), (
                f"IC3 Cart->Sph: component {i} does not match"
            )

    def test_schwarzschild_geodesic_with_orbital_elements(self):
        """
        Geodesic in spherical form matches when starting from orbital elements.

        Converts initial data to orbital elements, integrates, converts back to
        spherical, and compares all seven components to the direct spherical
        integration (relaxed tolerance for set 3).

        Parameters
        ----------
        self : TestSchwarzschildMetric
            Test instance (unused).

        Examples
        --------
        .. code-block:: bash

            pytest tests/numeric/metrics/test_schwarzschild_metric.py::TestSchwarzschildMetric::test_schwarzschild_geodesic_with_orbital_elements -v
        """
        # IC set 1
        initial_conditions_1_orbital_elements = initial_conditions_1_rp.convert_to("OrbitalElements", mass=M_1)
        taus_1 = np.linspace(0, 100, 100)
        sch = rp_Schwarzschild(M_1)

        path_spherical = sch.geodesic.get_path(initial_conditions_1_rp, taus_1).convert_to("Spherical")
        path_orbital_elements = sch.geodesic.get_path(initial_conditions_1_orbital_elements, taus_1).convert_to("Spherical")
        for i in range(7):
            assert np.isclose(path_spherical[i], path_orbital_elements[i]).all(), (
                f"IC1 OE: component {i} does not match"
            )

        # IC set 2
        initial_conditions_2_orbital_elements = initial_conditions_2_rp.convert_to("OrbitalElements", mass=M_2)
        taus_2 = np.linspace(0, 100, 100)
        sch = rp_Schwarzschild(M_2)

        path_spherical = sch.geodesic.get_path(initial_conditions_2_rp, taus_2).convert_to("Spherical")
        path_orbital_elements = sch.geodesic.get_path(initial_conditions_2_orbital_elements, taus_2).convert_to("Spherical")
        for i in range(7):
            assert np.isclose(path_spherical[i], path_orbital_elements[i]).all(), (
                f"IC2 OE: component {i} does not match"
            )

        # IC set 3
        initial_conditions_3_orbital_elements = initial_conditions_3_rp.convert_to("OrbitalElements", mass=M_3)
        taus_3 = np.linspace(0, 100, 100)
        sch = rp_Schwarzschild(M_3)

        path_spherical = sch.geodesic.get_path(initial_conditions_3_rp, taus_3).convert_to("Spherical")
        path_orbital_elements = sch.geodesic.get_path(initial_conditions_3_orbital_elements, taus_3).convert_to("Spherical")
        for i in range(7):
            assert np.isclose(path_spherical[i], path_orbital_elements[i], atol=1e-4).all(), (
                f"IC3 OE: component {i} does not match"
            )
