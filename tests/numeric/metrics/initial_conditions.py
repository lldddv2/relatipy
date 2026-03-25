"""
Test fixtures: masses, positions, momenta, and coordinate states for metric tests.

This module defines three independent scenarios used by
``test_schwarzschild_metric`` and ``test_kerr_metric``. Each scenario provides a
central mass :math:`M`, spherical-like coordinates
:math:`(r, \\theta, \\phi)`, conjugate momentum components, and Astropy
``Quantity`` arrays ``xs`` (time plus positions) and ``vs`` (velocity-like
components) suitable for building :class:`~relatipy.numeric.coordinates.spherical.Spherical`
or Boyer–Lindquist differential objects.

Notes
-----
For each suffix ``N`` in ``{1, 2, 3}``:

* ``M_N`` — Gravitating mass as an :class:`astropy.units.Quantity` in kilograms.
* ``position_ep_N`` — List ``[r, \\theta, \\phi]`` with :math:`r` in meters and
  angles in radians (unitless floats).
* ``momentum_ep_N`` — List of three momentum components matching the coordinate
  convention used by the tests (see the metric/coordinate modules).
* ``xs_N`` — ``[t, r, \\theta, \\phi]`` with :math:`t = 0` and spatial parts
  taken from ``position_ep_N`` with :mod:`astropy.units` attached.
* ``vs_N`` — Components of the velocity/momentum vector with units consistent
  with ``xs_N``.

Set 1 uses an Earth-like mass and a low-altitude orbit scale; set 2 a solar
mass and a larger radius; set 3 a very large mass with non-zero radial momentum.

Examples
--------
The structure of the first scenario (same construction as ``M_1``, ``xs_1``):

>>> import numpy as np
>>> import astropy.units as u
>>> M = 5.972e24 * u.kg
>>> r, theta, phi = 7e6, np.pi / 2, 0.0
>>> xs = [0.0 * u.s, r * u.m, theta * u.rad, phi * u.rad]
>>> len(xs)
4
>>> float(M.to(u.kg).value) > 0
True
"""

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
