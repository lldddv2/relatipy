import numpy as np
import astropy.units as u
from einsteinpy.metric import Schwarzschild
from einsteinpy.coordinates import SphericalDifferential
from relatipy.numeric.metrics import Schwarzschild as rp_Schwarzschild
from relatipy.numeric.coordinates import Spherical
from einsteinpy.geodesic.geodesic import Geodesic
from einsteinpy.geodesic import Timelike

# CI 1
M_1 = 5.972e24 * u.kg

position_ep_1 = [7000e3, np.pi / 2, 0.0]
momentum_ep_1 = [0.0, 70, 0.0]

xs_1 = [0.0 * u.s, position_ep_1[0] * u.m, position_ep_1[1] * u.rad, position_ep_1[2] * u.rad]
vs_1 = [momentum_ep_1[0] * u.m / u.s, momentum_ep_1[1] * u.rad / u.s, momentum_ep_1[2] * u.rad / u.s]

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

# CI 2
M_2 = 1.989e30 * u.kg

position_ep_2 = [900e3, np.pi / 3, 0.0]
momentum_ep_2 = [0.0, 0, 10]

xs_2 = [0.0 * u.s, position_ep_2[0] * u.m, position_ep_2[1] * u.rad, position_ep_2[2] * u.rad]
vs_2 = [momentum_ep_2[0] * u.m / u.s, momentum_ep_2[1] * u.rad / u.s, momentum_ep_2[2] * u.rad / u.s]

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

# CI 3
M_3 = 2e32 * u.kg

position_ep_3 = [1.5e11, np.pi / 2, 0.0]
momentum_ep_3 = [10, 0, 0]

xs_3 = [0.0 * u.s, position_ep_3[0] * u.m, position_ep_3[1] * u.rad, position_ep_3[2] * u.rad]
vs_3 = [momentum_ep_3[0] * u.m / u.s, momentum_ep_3[1] * u.rad / u.s, momentum_ep_3[2] * u.rad / u.s]

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

    def test_schwarzschild_geodesic(self):
        # CI 1
        taus_1 = np.linspace(0, 100, 100)
        steps_1 = len(taus_1)
        delta_1 = (taus_1[-1] - taus_1[0]) / steps_1

        geod_ep = Timelike(
            metric="Schwarzschild",
            metric_params=(0.0,),
            position=position_ep_1,
            momentum=momentum_ep_1,
            steps=steps_1,
            delta=delta_1,
            suppress_warnings=True,
            return_cartesian=False,
        )

        traj_ep = geod_ep.trajectory[1]  # (100, 8)

        sch = rp_Schwarzschild(M_1)
        traj_rp = sch.geodesic.get_path(initial_conditions_1_rp, taus_1)  # (8, 100)

        pos_ep = traj_ep[:, :4].T           # (4, 100)
        pos_rp = traj_rp[:4, :]        # (4, 100)

        assert np.isclose(pos_ep[1], pos_rp[1]).all(), "The second position is not the same" # check if the second position is the same
        assert np.isclose(pos_ep[2], pos_rp[2]).all(), "The third position is not the same" # check if the third position is the same
        assert np.isclose(pos_ep[3], pos_rp[3]).all(), "The fourth position is not the same" # check if the fourth position is the same
        
        # CI 2
        taus_2 = np.linspace(0, 100, 100)
        steps_2 = len(taus_2)
        delta_2 = (taus_2[-1] - taus_2[0]) / steps_2

        geod_ep = Timelike(
            metric="Schwarzschild",
            metric_params=(0.0,),
            position=position_ep_2,
            momentum=momentum_ep_2,
            steps=steps_2,
            delta=delta_2,
            suppress_warnings=True,
            return_cartesian=False,
        )

        traj_ep = geod_ep.trajectory[1]  # (100, 8)

        sch = rp_Schwarzschild(M_2)
        traj_rp = sch.geodesic.get_path(initial_conditions_2_rp, taus_2)  # (8, 100)

        pos_ep = traj_ep[:, :4].T           # (4, 100)
        pos_rp = traj_rp[:4, :]        # (4, 100)

        assert np.isclose(pos_ep[1], pos_rp[1]).all(), "The second position is not the same" # check if the second position is the same
        assert np.isclose(pos_ep[2], pos_rp[2]).all(), "The third position is not the same" # check if the third position is the same
        assert np.isclose(pos_ep[3], pos_rp[3]).all(), "The fourth position is not the same" # check if the fourth position is the same
        
        # CI 3
        taus_3 = np.linspace(0, 100, 100)
        steps_3 = len(taus_3)
        delta_3 = (taus_3[-1] - taus_3[0]) / steps_3

        geod_ep = Timelike(
            metric="Schwarzschild",
            metric_params=(0.0,),
            position=position_ep_3,
            momentum=momentum_ep_3,
            steps=steps_3,
            delta=delta_3,
            suppress_warnings=True,
            return_cartesian=False,
        )
        traj_ep = geod_ep.trajectory[1]  # (100, 8)

        sch = rp_Schwarzschild(M_3)
        traj_rp = sch.geodesic.get_path(initial_conditions_3_rp, taus_3)  # (8, 100)

        pos_ep = traj_ep[:, :4].T           # (4, 100)
        pos_rp = traj_rp[:4, :]        # (4, 100)

        assert np.isclose(pos_ep[1], pos_rp[1]).all(), "The second position is not the same" # check if the second position is the same
        assert np.isclose(pos_ep[2], pos_rp[2]).all(), "The third position is not the same" # check if the third position is the same
        assert np.isclose(pos_ep[3], pos_rp[3]).all(), "The fourth position is not the same" # check if the fourth position is the same       