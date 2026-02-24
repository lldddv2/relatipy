import numpy as np
import astropy.units as u

M_1 = 5.972e24 * u.kg
position_ep_1 = [7e6, np.pi / 2, 0.0]
momentum_ep_1 = [0, 70, 10]
xs_1 = [0.0 * u.s, position_ep_1[0] * u.m, position_ep_1[1] * u.rad, position_ep_1[2] * u.rad]
vs_1 = [momentum_ep_1[0] * u.m / u.s, momentum_ep_1[1] * u.rad / u.s, momentum_ep_1[2] * u.rad / u.s]


M_2 = 1.989e30 * u.kg
position_ep_2 = [9e8, np.pi / 3, 0.0]
momentum_ep_2 = [0, 6, 10]
xs_2 = [0.0 * u.s, position_ep_2[0] * u.m, position_ep_2[1] * u.rad, position_ep_2[2] * u.rad]
vs_2 = [momentum_ep_2[0] * u.m / u.s, momentum_ep_2[1] * u.rad / u.s, momentum_ep_2[2] * u.rad / u.s]


M_3 = 2e32 * u.kg
position_ep_3 = [1.5e11, np.pi / 2, 0.0]
momentum_ep_3 = [1e4, 0, 10]
xs_3 = [0.0 * u.s, position_ep_3[0] * u.m, position_ep_3[1] * u.rad, position_ep_3[2] * u.rad]
vs_3 = [momentum_ep_3[0] * u.m / u.s, momentum_ep_3[1] * u.rad / u.s, momentum_ep_3[2] * u.rad / u.s]
