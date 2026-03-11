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
    Convierte trayectoria de relatipy (unidades geométricas) a SI.
    traj_rp shape: (4, N) → [t_geom, r_geom, theta, phi]
    Devuelve (4, N) con [t_s, r_m, theta, phi]
    """
    out = traj_rp.copy().astype(float)
    out[0] *= _T_ref   # t: unidades geométricas → segundos
    out[1] *= _L_ref   # r: unidades geométricas → metros
    # theta y phi son adimensionales, sin cambio
    return out


class TestSchwarzschildMetric:
    def test_schwarzschild_metric(self):
        # CI 1
        sch = Schwarzschild(coords=initial_conditions_1, M=M_1)
        g_ep = sch.metric_covariant(x_vec_1)
        g_rp = rp_Schwarzschild(M_1).metric(xs_1, dimensionless=False)
        assert np.isclose(g_ep, g_rp).all()

        # CI 2
        sch = Schwarzschild(coords=initial_conditions_2, M=M_2)
        g_ep = sch.metric_covariant(x_vec_2)
        g_rp = rp_Schwarzschild(M_2).metric(xs_2, dimensionless=False)
        assert np.isclose(g_ep, g_rp).all()

        # CI 3
        sch = Schwarzschild(coords=initial_conditions_3, M=M_3)
        g_ep = sch.metric_covariant(x_vec_3)
        g_rp = rp_Schwarzschild(M_3).metric(xs_3, dimensionless=False)
        assert np.isclose(g_ep, g_rp).all()

    def test_schwarzschild_christoffel_symbols(self):
        # CI 1
        sch = Schwarzschild(coords=initial_conditions_1, M=M_1)
        christoffel_ep = sch.christoffels(x_vec_1)
        christoffel_rp = rp_Schwarzschild(M_1).get_christoffel_symbols(xs_1, dimensionless=False)
        assert np.isclose(christoffel_ep, christoffel_rp).all()

        # CI 2
        sch = Schwarzschild(coords=initial_conditions_2, M=M_2)
        christoffel_ep = sch.christoffels(x_vec_2)
        christoffel_rp = rp_Schwarzschild(M_2).get_christoffel_symbols(xs_2, dimensionless=False)
        assert np.isclose(christoffel_ep, christoffel_rp).all()

        # CI 3
        sch = Schwarzschild(coords=initial_conditions_3, M=M_3)
        christoffel_ep = sch.christoffels(x_vec_3)
        christoffel_rp = rp_Schwarzschild(M_3).get_christoffel_symbols(xs_3, dimensionless=False)
        assert np.isclose(christoffel_ep, christoffel_rp).all()

    # Comparar geodésicas completas con einsteinpy se ha descartado por
    # diferencias de convención en las condiciones iniciales y 4-velocidad.

    def test_schwarzschild_ds_dtau(self):
        # CI 1
        taus_1 = np.linspace(0, 100, 100)
        sch = rp_Schwarzschild(M_1)
        path = sch.geodesic.get_path(initial_conditions_1_rp, taus_1)
        ds_dtau = path._get_ds_dtau(sch)
        assert np.isclose(ds_dtau, 1).all(), "CI1: ds/dtau no es 1"

        # CI 2
        taus_2 = np.linspace(0, 100, 100)
        sch = rp_Schwarzschild(M_2)
        path = sch.geodesic.get_path(initial_conditions_2_rp, taus_2)
        ds_dtau = path._get_ds_dtau(sch)
        assert np.isclose(ds_dtau, 1).all(), "CI2: ds/dtau no es 1"

        # CI 3
        taus_3 = np.linspace(0, 100, 100)
        sch = rp_Schwarzschild(M_3)
        path = sch.geodesic.get_path(initial_conditions_3_rp, taus_3)
        ds_dtau = path._get_ds_dtau(sch)
        assert np.isclose(ds_dtau, 1).all(), "CI3: ds/dtau no es 1"

    def test_schwarzschild_Lz(self):
        # CI 1
        taus_1 = np.linspace(0, 100, 100)
        sch = rp_Schwarzschild(M_1)
        path = sch.geodesic.get_path(initial_conditions_1_rp, taus_1)
        Lz = path._get_Lz()
        assert np.allclose(Lz, Lz[0]), "CI1: Lz no es constante"

        # CI 2
        taus_2 = np.linspace(0, 100, 100)
        sch = rp_Schwarzschild(M_2)
        path = sch.geodesic.get_path(initial_conditions_2_rp, taus_2)
        Lz = path._get_Lz()
        assert np.allclose(Lz, Lz[0]), "CI2: Lz no es constante"

        # CI 3
        taus_3 = np.linspace(0, 100, 100)
        sch = rp_Schwarzschild(M_3)
        path = sch.geodesic.get_path(initial_conditions_3_rp, taus_3)
        Lz = path._get_Lz()
        assert np.allclose(Lz, Lz[0]), "CI3: Lz no es constante"

    def test_schwarzschild_E(self):
        # CI 1
        taus_1 = np.linspace(0, 100, 100)
        sch = rp_Schwarzschild(M_1)
        path = sch.geodesic.get_path(initial_conditions_1_rp, taus_1)
        E = path._get_E(sch)
        assert np.allclose(E, E[0]), "CI1: E no es constante"

        # CI 2
        taus_2 = np.linspace(0, 100, 100)
        sch = rp_Schwarzschild(M_2)
        path = sch.geodesic.get_path(initial_conditions_2_rp, taus_2)
        E = path._get_E(sch)
        assert np.allclose(E, E[0]), "CI2: E no es constante"

        # CI 3
        taus_3 = np.linspace(0, 100, 100)
        sch = rp_Schwarzschild(M_3)
        path = sch.geodesic.get_path(initial_conditions_3_rp, taus_3)
        E = path._get_E(sch)
        assert np.allclose(E, E[0]), "CI3: E no es constante"

    def test_same_trajectory_with_different_coordinates(self):
        # CI 1
        initial_conditions_1_cartesian = initial_conditions_1_rp.convert_to("Cartesian")
        taus_1 = np.linspace(0, 100, 100)
        sch = rp_Schwarzschild(M_1)

        path_spherical = sch.geodesic.get_path(initial_conditions_1_rp, taus_1)
        path_cartesian = sch.geodesic.get_path(initial_conditions_1_cartesian, taus_1).convert_to("Spherical")
        for i in range(7):
            assert np.isclose(path_spherical[i], path_cartesian[i]).all(), f"CI1 Sph->Cart: componente {i} no coincide"

        path_cartesian = sch.geodesic.get_path(initial_conditions_1_cartesian, taus_1)
        path_spherical = sch.geodesic.get_path(initial_conditions_1_rp, taus_1).convert_to("Cartesian")
        for i in range(7):
            assert np.isclose(path_cartesian[i], path_spherical[i]).all(), f"CI1 Cart->Sph: componente {i} no coincide"

        # CI 2
        initial_conditions_2_cartesian = initial_conditions_2_rp.convert_to("Cartesian")
        taus_2 = np.linspace(0, 100, 100)
        sch = rp_Schwarzschild(M_2)

        path_spherical = sch.geodesic.get_path(initial_conditions_2_rp, taus_2)
        path_cartesian = sch.geodesic.get_path(initial_conditions_2_cartesian, taus_2).convert_to("Spherical")
        for i in range(7):
            assert np.isclose(path_spherical[i], path_cartesian[i]).all(), f"CI2 Sph->Cart: componente {i} no coincide"

        path_cartesian = sch.geodesic.get_path(initial_conditions_2_cartesian, taus_2)
        path_spherical = sch.geodesic.get_path(initial_conditions_2_rp, taus_2).convert_to("Cartesian")
        for i in range(7):
            assert np.isclose(path_cartesian[i], path_spherical[i]).all(), f"CI2 Cart->Sph: componente {i} no coincide"

        # CI 3
        initial_conditions_3_cartesian = initial_conditions_3_rp.convert_to("Cartesian")
        taus_3 = np.linspace(0, 100, 100)
        sch = rp_Schwarzschild(M_3)

        path_spherical = sch.geodesic.get_path(initial_conditions_3_rp, taus_3)
        path_cartesian = sch.geodesic.get_path(initial_conditions_3_cartesian, taus_3).convert_to("Spherical")
        for i in range(7):
            assert np.isclose(path_spherical[i], path_cartesian[i]).all(), f"CI3 Sph->Cart: componente {i} no coincide"

        path_cartesian = sch.geodesic.get_path(initial_conditions_3_cartesian, taus_3)
        path_spherical = sch.geodesic.get_path(initial_conditions_3_rp, taus_3).convert_to("Cartesian")
        for i in range(7):
            assert np.isclose(path_cartesian[i], path_spherical[i]).all(), f"CI3 Cart->Sph: componente {i} no coincide"

    def test_schwarzschild_geodesic_with_orbital_elements(self):
        # CI 1
        initial_conditions_1_orbital_elements = initial_conditions_1_rp.convert_to("OrbitalElements", mass=M_1)
        taus_1 = np.linspace(0, 100, 100)
        sch = rp_Schwarzschild(M_1)

        path_spherical = sch.geodesic.get_path(initial_conditions_1_rp, taus_1).convert_to("Spherical")
        path_orbital_elements = sch.geodesic.get_path(initial_conditions_1_orbital_elements, taus_1).convert_to("Spherical")
        for i in range(7):
            assert np.isclose(path_spherical[i], path_orbital_elements[i]).all(), f"CI1 OE: componente {i} no coincide"

        # CI 2
        initial_conditions_2_orbital_elements = initial_conditions_2_rp.convert_to("OrbitalElements", mass=M_2)
        taus_2 = np.linspace(0, 100, 100)
        sch = rp_Schwarzschild(M_2)

        path_spherical = sch.geodesic.get_path(initial_conditions_2_rp, taus_2).convert_to("Spherical")
        path_orbital_elements = sch.geodesic.get_path(initial_conditions_2_orbital_elements, taus_2).convert_to("Spherical")
        for i in range(7):
            assert np.isclose(path_spherical[i], path_orbital_elements[i]).all(), f"CI2 OE: componente {i} no coincide"

        # CI 3
        initial_conditions_3_orbital_elements = initial_conditions_3_rp.convert_to("OrbitalElements", mass=M_3)
        taus_3 = np.linspace(0, 100, 100)
        sch = rp_Schwarzschild(M_3)

        path_spherical = sch.geodesic.get_path(initial_conditions_3_rp, taus_3).convert_to("Spherical")
        path_orbital_elements = sch.geodesic.get_path(initial_conditions_3_orbital_elements, taus_3).convert_to("Spherical")
        for i in range(7):
            assert np.isclose(path_spherical[i], path_orbital_elements[i], atol=1e-4).all(), f"CI3 OE: componente {i} no coincide"