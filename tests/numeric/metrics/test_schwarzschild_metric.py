import numpy as np
import astropy.units as u
from einsteinpy.metric import Schwarzschild
from einsteinpy.coordinates import SphericalDifferential
from relatipy.numeric.metrics import Schwarzschild as rp_Schwarzschild
from relatipy.numeric.coordinates import Spherical

# CI 1
M_1 = 5.972e24 * u.kg
xs_1 = [0.0 * u.s, 7000e3 * u.m, np.pi / 2 * u.rad, 0.0 * u.rad]
vs_1 = [0.0 * u.m / u.s, 70 * u.rad / u.s, 0.0 * u.rad / u.s]

initial_conditions_1 = SphericalDifferential(
    t=xs_1[0],
    r=xs_1[1],
    theta=xs_1[2],
    phi=xs_1[3],
    v_r=vs_1[0],
    v_th=vs_1[1],
    v_p=vs_1[2],
)

x_vec_1 = np.array(initial_conditions_1.position())

# CI 2
M_2 = 5.972e24 * u.kg
xs_2 = [0.0 * u.s, 7000e3 * u.m, np.pi / 2 * u.rad, 0.0 * u.rad]
vs_2 = [0.0 * u.m / u.s, 70 * u.rad / u.s, 0.0 * u.rad / u.s]

initial_conditions_2 = SphericalDifferential(
    t=xs_2[0],
    r=xs_2[1],
    theta=xs_2[2],
    phi=xs_2[3],
    v_r=vs_2[0],
    v_th=vs_2[1],
    v_p=vs_2[2],
)

x_vec_2 = np.array(initial_conditions_2.position())

# CI 3
M_3 = 5.972e24 * u.kg
xs_3 = [0.0 * u.s, 7000e3 * u.m, np.pi / 2 * u.rad, 0.0 * u.rad]
vs_3 = [0.0 * u.m / u.s, 70 * u.rad / u.s, 0.0 * u.rad / u.s]

initial_conditions_3 = SphericalDifferential(
    t=xs_3[0],
    r=xs_3[1],
    theta=xs_3[2],
    phi=xs_3[3],
    v_r=vs_3[0],
    v_th=vs_3[1],
    v_p=vs_3[2],
)

x_vec_3 = np.array(initial_conditions_3.position())

class TestSchwarzschildMetric:
    def test_schwarzschild_metric(self):
        # CI 1
        sch = Schwarzschild(coords=initial_conditions_1, M=M_1)
        g = sch.metric_covariant(x_vec_1)
        g_rp = rp_Schwarzschild(M_1).metric(xs_1, dimensionless=False)
        assert np.isclose(g, g_rp).all()

        # CI 2
        sch = Schwarzschild(coords=initial_conditions_2, M=M_2)
        g = sch.metric_covariant(x_vec_2)
        g_rp = rp_Schwarzschild(M_2).metric(xs_2, dimensionless=False)
        assert np.isclose(g, g_rp).all()

        # CI 3
        sch = Schwarzschild(coords=initial_conditions_3, M=M_3)
        g = sch.metric_covariant(x_vec_3)
        g_rp = rp_Schwarzschild(M_3).metric(xs_3, dimensionless=False)
        assert np.isclose(g, g_rp).all()
    
    def test_schwarzschild_christoffel_symbols(self):
        # CI 1
        sch = Schwarzschild(coords=initial_conditions_1, M=M_1)
        christoffel = sch.christoffels(x_vec_1)
        christoffel_rp = rp_Schwarzschild(M_1).get_christoffel_symbols(xs_1, dimensionless=False)
        assert np.isclose(christoffel, christoffel_rp).all()

        # CI 2
        sch = Schwarzschild(coords=initial_conditions_2, M=M_2)
        christoffel = sch.christoffels(x_vec_2)
        christoffel_rp = rp_Schwarzschild(M_2).get_christoffel_symbols(xs_2, dimensionless=False)
        assert np.isclose(christoffel, christoffel_rp).all()

        # CI 3
        sch = Schwarzschild(coords=initial_conditions_3, M=M_3)
        christoffel = sch.christoffels(x_vec_3)
        christoffel_rp = rp_Schwarzschild(M_3).get_christoffel_symbols(xs_3, dimensionless=False)
        assert np.isclose(christoffel, christoffel_rp).all()
