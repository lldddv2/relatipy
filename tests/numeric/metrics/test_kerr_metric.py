# test_kerr_metric.py

import numpy as np
import astropy.units as u
from einsteinpy.metric import Kerr as ep_Kerr
from einsteinpy.coordinates import BoyerLindquistDifferential
from relatipy.numeric.metrics import Kerr as rp_Kerr
from relatipy.numeric.coordinates import BoyerLindquist
from einsteinpy.geodesic import Timelike

# CI 1
M_1 = 5.972e24 * u.kg
a_1 = 0.5  # dimensionless spin parameter

position_ep_1 = [700e3, np.pi / 2, 0.0]
momentum_ep_1 = [0.0, 70, 0.0]

xs_1 = [0.0 * u.s, position_ep_1[0] * u.m, position_ep_1[1] * u.rad, position_ep_1[2] * u.rad]
vs_1 = [momentum_ep_1[0] * u.m / u.s, momentum_ep_1[1] * u.rad / u.s, momentum_ep_1[2] * u.rad / u.s]

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
M_2 = 1.989e30 * u.kg
a_2 = 0.9

position_ep_2 = [900e3, np.pi / 3, 0.0]
momentum_ep_2 = [0.0, 0, 10]

xs_2 = [0.0 * u.s, position_ep_2[0] * u.m, position_ep_2[1] * u.rad, position_ep_2[2] * u.rad]
vs_2 = [momentum_ep_2[0] * u.m / u.s, momentum_ep_2[1] * u.rad / u.s, momentum_ep_2[2] * u.rad / u.s]

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
M_3 = 2000 * u.kg
a_3 = 0.7

position_ep_3 = [1.5e11, np.pi / 2, 0.0]
momentum_ep_3 = [10, 0, 0]

xs_3 = [0.0 * u.s, position_ep_3[0] * u.m, position_ep_3[1] * u.rad, position_ep_3[2] * u.rad]
vs_3 = [momentum_ep_3[0] * u.m / u.s, momentum_ep_3[1] * u.rad / u.s, momentum_ep_3[2] * u.rad / u.s]

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
    def test_kerr_metric(self):
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
        # CI 1
        kerr = rp_Kerr(M_1, a_1)
        christoffel = kerr.get_christoffel_symbols(xs_1, dimensionless=False)
        christoffel_rp = kerr.get_christoffel_symbols(xs_1, dimensionless=False)
        assert np.isclose(christoffel, christoffel_rp).all()

        # CI 2
        kerr = rp_Kerr(M_2, a_2)
        christoffel = kerr.get_christoffel_symbols(xs_2, dimensionless=False)
        christoffel_rp = kerr.get_christoffel_symbols(xs_2, dimensionless=False)
        assert np.isclose(christoffel, christoffel_rp).all()

        # CI 3
        kerr = rp_Kerr(M_3, a_3)
        christoffel = kerr.get_christoffel_symbols(xs_3, dimensionless=False)
        christoffel_rp = kerr.get_christoffel_symbols(xs_3, dimensionless=False)
        assert np.isclose(christoffel, christoffel_rp).all()

    def test_kerr_geodesic(self):
        # CI 1
        taus_1 = np.linspace(0, 100, 100)
        steps_1 = len(taus_1)
        delta_1 = (taus_1[-1] - taus_1[0]) / steps_1

        geod_ep = Timelike(
            metric="Kerr",
            metric_params=(a_1,),
            position=position_ep_1,
            momentum=momentum_ep_1,
            steps=steps_1,
            delta=delta_1,
            suppress_warnings=True,
            return_cartesian=False,
        )

        traj_ep = geod_ep.trajectory[1]

        kerr = rp_Kerr(M_1, a_1)
        ys0 = kerr.get_4state_vector(initial_conditions_1_rp)
        traj_rp = kerr.geodesic._get_path_from_4state_vector(ys0, taus_1)

        pos_ep = traj_ep[:, :4].T
        pos_rp = traj_rp[:4, :]

        assert np.isclose(pos_ep[1], pos_rp[1]).all(), "The second position is not the same"
        assert np.isclose(pos_ep[2], pos_rp[2]).all(), "The third position is not the same"
        assert np.isclose(pos_ep[3], pos_rp[3]).all(), "The fourth position is not the same"

        # CI 2
        taus_2 = np.linspace(0, 100, 100)
        steps_2 = len(taus_2)
        delta_2 = (taus_2[-1] - taus_2[0]) / steps_2

        geod_ep = Timelike(
            metric="Kerr",
            metric_params=(a_2,),
            position=position_ep_2,
            momentum=momentum_ep_2,
            steps=steps_2,
            delta=delta_2,
            suppress_warnings=True,
            return_cartesian=False,
        )

        traj_ep = geod_ep.trajectory[1]

        kerr = rp_Kerr(M_2, a_2)
        ys0 = kerr.get_4state_vector(initial_conditions_2_rp)
        traj_rp = kerr.geodesic._get_path_from_4state_vector(ys0, taus_2)

        pos_ep = traj_ep[:, :4].T
        pos_rp = traj_rp[:4, :]

        assert np.isclose(pos_ep[1], pos_rp[1]).all(), "The second position is not the same"
        assert np.isclose(pos_ep[2], pos_rp[2]).all(), "The third position is not the same"
        assert np.isclose(pos_ep[3], pos_rp[3]).all(), "The fourth position is not the same"

        # CI 3
        taus_3 = np.linspace(0, 100, 100)
        steps_3 = len(taus_3)
        delta_3 = (taus_3[-1] - taus_3[0]) / steps_3

        geod_ep = Timelike(
            metric="Kerr",
            metric_params=(a_3,),
            position=position_ep_3,
            momentum=momentum_ep_3,
            steps=steps_3,
            delta=delta_3,
            suppress_warnings=True,
            return_cartesian=False,
        )

        traj_ep = geod_ep.trajectory[1]

        kerr = rp_Kerr(M_3, a_3)
        ys0 = kerr.get_4state_vector(initial_conditions_3_rp)
        traj_rp = kerr.geodesic._get_path_from_4state_vector(ys0, taus_3)

        pos_ep = traj_ep[:, :4].T
        pos_rp = traj_rp[:4, :]

        assert np.isclose(pos_ep[1], pos_rp[1]).all(), "The second position is not the same"
        assert np.isclose(pos_ep[2], pos_rp[2]).all(), "The third position is not the same"
        assert np.isclose(pos_ep[3], pos_rp[3]).all(), "The fourth position is not the same"