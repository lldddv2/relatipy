"""
Tests for Keplerian :class:`~relatipy.numeric.coordinates.OrbitalElements` conversions.

These tests verify that converting orbital elements to Cartesian state and back
preserves all osculating elements (time, semi-major axis, eccentricity,
inclination, longitude of ascending node, argument of periapsis, true anomaly,
and central mass). A third case checks Cartesian :math:`(t, x, y, z, v_x, v_y,
v_z)` round-trip for a bound orbit (:math:`e < 1`).

Notes
-----
Test data are grouped as *CI 1* (scalar elements), *CI 2* (broadcast lists of
equal length), and *CI 3* (reference orbit from elements, then Cartesian
round-trip with explicit central mass in ``astropy`` units).

References
----------
The implementation uses SPICE ``conics`` / ``oscelt`` via
:class:`~relatipy.numeric.coordinates.OrbitalElements`.

Examples
--------
>>> import numpy as np
>>> from relatipy.numeric.coordinates import OrbitalElements
>>> oe = OrbitalElements(0, 1.0, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0)
>>> oe_back = oe.convert_to("Cartesian").convert_to("OrbitalElements", mass=1.0)
>>> bool(np.isclose(oe.a, oe_back.a).all())
True
"""

import numpy as np
import astropy.units as u_astropy
from relatipy.numeric.coordinates import OrbitalElements, Cartesian

# CI 1: scalar orbital elements
a = 1.0
e = 0.1
inc = 0.2
Omega = 0.3
omega = 0.4
f = 0.5
mass = 1.0


# CI 2: equal-length sequences (broadcast test data)
ts = [0.0, 1.0, 2.0, 3.0]
as_ = [1.0, 2.0, 3.0, 4.0]
es = [0.1, 0.2, 0.3, 0.4]
incs = [0.2, 0.3, 0.4, 0.5]
Omegas = [0.3, 0.4, 0.5, 0.6]
omegas = [0.4, 0.5, 0.6, 0.7]
fs = [0.5, 0.6, 0.7, 0.8]
mass_2 = 1.0

# CI 3: Cartesian state of a bound orbit (e < 1) for a valid round-trip
_mass_3 = 1e30 * u_astropy.kg  # explicit central mass in kg
_oe_ref = OrbitalElements(0.0, 1e7, 0.1, 0.2, 0.3, 0.4, 0.5, _mass_3)
_cart_ref = _oe_ref.convert_to("Cartesian")
xs = _cart_ref.xs
vs = _cart_ref.vs
mass_3 = _mass_3


class TestOrbitalElements:
    """
    Regression tests for :class:`~relatipy.numeric.coordinates.OrbitalElements`.

    Covers scalar inputs, sequence inputs, and Cartesian round-trip with
    physical mass units.

    Examples
    --------
    Run with pytest::

        pytest tests/numeric/coordinates/test_orbital_elements.py -v
    """

    def test_orbital_elements_conservation(self):
        """
        Assert element conservation under Cartesian round-trip conversion.

        Three cases are checked:

        - **CI 1**: scalar ``OrbitalElements`` converted
          ``OrbitalElements -> Cartesian -> OrbitalElements`` with fixed mass.
        - **CI 2**: same pipeline with list-valued elements of equal length.
        - **CI 3**: ``Cartesian -> OrbitalElements -> Cartesian`` for indices
          ``0``–``6`` (time, position, velocity), using ``mass_3`` in kg.

        Parameters
        ----------
        None
            Pytest injects no parameters.

        Returns
        -------
        None
            Assertions raise on failure.

        Examples
        --------
        >>> import numpy as np
        >>> from relatipy.numeric.coordinates import OrbitalElements
        >>> oe = OrbitalElements(0, 1.0, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0)
        >>> oe2 = oe.convert_to("Cartesian").convert_to("OrbitalElements", mass=1.0)
        >>> np.isclose(oe.a, oe2.a).all()
        True
        """
        # CI 1
        orbital_elements = OrbitalElements(0, a, e, inc, Omega, omega, f, mass_2)
        orbital_elements_cartesian = orbital_elements.convert_to("Cartesian").convert_to(
            "OrbitalElements", mass=mass_2
        )

        assert np.isclose(orbital_elements.t, orbital_elements_cartesian.t).all()
        assert np.isclose(orbital_elements.a, orbital_elements_cartesian.a).all()
        assert np.isclose(orbital_elements.e, orbital_elements_cartesian.e).all()
        assert np.isclose(orbital_elements.inc, orbital_elements_cartesian.inc).all()
        assert np.isclose(orbital_elements.Omega, orbital_elements_cartesian.Omega).all()
        assert np.isclose(orbital_elements.omega, orbital_elements_cartesian.omega).all()
        assert np.isclose(orbital_elements.f, orbital_elements_cartesian.f).all()
        assert np.isclose(orbital_elements.mass, orbital_elements_cartesian.mass).all()

        # CI 2
        orbital_elements = OrbitalElements(ts, as_, es, incs, Omegas, omegas, fs, mass)
        orbital_elements_cartesian = orbital_elements.convert_to("Cartesian").convert_to(
            "OrbitalElements", mass=mass_2
        )

        assert np.isclose(orbital_elements.t, orbital_elements_cartesian.t).all()
        assert np.isclose(orbital_elements.a, orbital_elements_cartesian.a).all()
        assert np.isclose(orbital_elements.e, orbital_elements_cartesian.e).all()
        assert np.isclose(orbital_elements.inc, orbital_elements_cartesian.inc).all()
        assert np.isclose(orbital_elements.Omega, orbital_elements_cartesian.Omega).all()
        assert np.isclose(orbital_elements.omega, orbital_elements_cartesian.omega).all()
        assert np.isclose(orbital_elements.f, orbital_elements_cartesian.f).all()
        assert np.isclose(orbital_elements.mass, orbital_elements_cartesian.mass)

        # CI 3: round-trip Cartesian -> OrbitalElements -> Cartesian (indices 0–6: t,x,y,z,vx,vy,vz)
        orbital_elements = Cartesian(xs, vels=vs)
        orbital_elements_orbital_elements = orbital_elements.convert_to(
            "OrbitalElements", mass=mass_3
        ).convert_to("Cartesian")

        assert np.isclose(orbital_elements[0], orbital_elements_orbital_elements[0]).all()
        assert np.isclose(orbital_elements[1], orbital_elements_orbital_elements[1]).all()
        assert np.isclose(orbital_elements[2], orbital_elements_orbital_elements[2]).all()
        assert np.isclose(orbital_elements[3], orbital_elements_orbital_elements[3]).all()
        assert np.isclose(orbital_elements[4], orbital_elements_orbital_elements[4]).all()
        assert np.isclose(orbital_elements[5], orbital_elements_orbital_elements[5]).all()
        assert np.isclose(orbital_elements[6], orbital_elements_orbital_elements[6]).all()
