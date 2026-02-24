# test_kerr_metric.py

import numpy as np
import astropy.units as u
from einsteinpy.metric import Kerr as ep_Kerr
from einsteinpy.coordinates import BoyerLindquistDifferential
from relatipy.numeric.metrics import Kerr as rp_Kerr
from relatipy.numeric.coordinates import BoyerLindquist
from einsteinpy.geodesic import Timelike
from initial_conditions import M_1, M_2, M_3, xs_1, vs_1, xs_2, vs_2, xs_3, vs_3, position_ep_1, position_ep_2, position_ep_3, momentum_ep_1, momentum_ep_2, momentum_ep_3

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
        traj_rp = kerr.geodesic.get_path(initial_conditions_1_rp, taus_1)

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
        traj_rp = kerr.geodesic.get_path(initial_conditions_2_rp, taus_2)

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
        traj_rp = kerr.geodesic.get_path(initial_conditions_3_rp, taus_3)

        pos_ep = traj_ep[:, :4].T
        pos_rp = traj_rp[:4, :]

        assert np.isclose(pos_ep[1], pos_rp[1]).all(), "The second position is not the same"
        assert np.isclose(pos_ep[2], pos_rp[2]).all(), "The third position is not the same"
        assert np.isclose(pos_ep[3], pos_rp[3]).all(), "The fourth position is not the same"

    def test_kerr_ds_dtau(self):
        # CI 1
        taus_1 = np.linspace(0, 100, 100)
        kerr = rp_Kerr(M_1, a_1)
        path = kerr.geodesic.get_path(initial_conditions_1_rp, taus_1)
        ds_dtau = path._get_ds_dtau(kerr)
        assert np.isclose(ds_dtau, 1).all(), "The ds/dtau is not c**2"

        # CI 2
        taus_2 = np.linspace(0, 100, 100)
        kerr = rp_Kerr(M_2, a_2)
        path = kerr.geodesic.get_path(initial_conditions_2_rp, taus_2)
        ds_dtau = path._get_ds_dtau(kerr)
        assert np.isclose(ds_dtau, 1).all(), "The ds/dtau is not c**2"

        # CI 3
        taus_3 = np.linspace(0, 100, 100)
        kerr = rp_Kerr(M_3, a_3)
        path = kerr.geodesic.get_path(initial_conditions_3_rp, taus_3)
        ds_dtau = path._get_ds_dtau(kerr)
        assert np.isclose(ds_dtau, 1).all(), "The ds/dtau is not c**2"

    def test_kerr_E(self):
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

    def test_kerr_Q(self):
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