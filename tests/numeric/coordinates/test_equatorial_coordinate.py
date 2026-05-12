"""
Tests for :class:`~relatipy.numeric.coordinates.EquatorialCoordinate`.

Covers:

- Zero-position → zero RA and Dec offsets.
- Pure east offset (BH frame, on-axis spin) → non-zero RA, zero Dec.
- Pure north offset (BH frame, on-axis spin) → zero RA, non-zero Dec.
- Exact angular values for a known geometry.
- ``from_apparent_orbital_elements`` factory: ``zeta``/``eta`` are forwarded
  and the result is an :class:`EquatorialCoordinate`.
- Input validation: missing ``distance`` Quantity, wrong ``distance`` units,
  plain-float position without ``mass``.
- ``to_sexagesimal`` returns two strings.
- ``__repr__`` and ``__eq__`` behave correctly.
- Registry: class is registered in ``coordinate_systems``.

Examples
--------
Run with pytest::

    pytest tests/numeric/coordinates/test_equatorial_coordinate.py -v
"""

import numpy as np
import pytest
import astropy.units as u

from relatipy.numeric.coordinates import EquatorialCoordinate, coordinate_systems
from relatipy.numeric.coordinates.apparent_orbital_elements import ApparentOrbitalElements


DISTANCE = 8.5 * u.kpc
MASS = 4.0e6  # solar masses (Sgr A* scale)


class TestEquatorialCoordinateBasic:
    """
    Basic construction and attribute checks for :class:`EquatorialCoordinate`.

    Examples
    --------
    >>> import astropy.units as u
    >>> from relatipy.numeric.coordinates import EquatorialCoordinate
    >>> ec = EquatorialCoordinate(0.0, 0.0, 0.0, zeta=0.0, eta=0.0,
    ...                           distance=8.5 * u.kpc, mass=4e6)
    >>> ec.ra, ec.dec
    (0.0, 0.0)
    """

    def test_zero_position_zero_angles(self):
        """Zero BH-frame position → zero RA and Dec offsets."""
        ec = EquatorialCoordinate(0.0, 0.0, 0.0, zeta=0.0, eta=0.0,
                                  distance=DISTANCE, mass=MASS)
        assert ec.ra == 0.0
        assert ec.dec == 0.0

    def test_name_metric(self):
        """name_metric attribute is 'EquatorialCoordinate'."""
        ec = EquatorialCoordinate(0.0, 0.0, 0.0, zeta=0.0, eta=0.0,
                                  distance=DISTANCE, mass=MASS)
        assert ec.name_metric == "EquatorialCoordinate"

    def test_zeta_eta_stored(self):
        """zeta and eta are stored as floats on the instance."""
        ec = EquatorialCoordinate(0.0, 0.0, 0.0, zeta=30.0, eta=45.0,
                                  distance=DISTANCE, mass=MASS)
        assert ec.zeta == 30.0
        assert ec.eta == 45.0

    def test_distance_stored(self):
        """distance attribute equals the Quantity passed at construction."""
        ec = EquatorialCoordinate(0.0, 0.0, 0.0, zeta=0.0, eta=0.0,
                                  distance=DISTANCE, mass=MASS)
        assert ec.distance == DISTANCE

    def test_in_registry(self):
        """EquatorialCoordinate is registered in coordinate_systems."""
        assert "EquatorialCoordinate" in coordinate_systems
        assert coordinate_systems["EquatorialCoordinate"] is EquatorialCoordinate


class TestEquatorialCoordinateAngles:
    """
    Angular projection tests with on-axis BH spin (:math:`\\zeta = 0`).

    When :math:`\\zeta = 0` the rotation :math:`R` is the identity, so
    BH frame and sky frame are the same:

    - ``y_BH`` maps to ``y_east`` → drives RA.
    - ``x_BH`` maps to ``x_north`` → drives Dec.
    """

    def test_pure_east_offset(self):
        """y_BH ≠ 0, x_BH = z_BH = 0 → ra ≠ 0, dec = 0 (zeta=0)."""
        ec = EquatorialCoordinate(0.0, 1000.0, 0.0, zeta=0.0, eta=0.0,
                                  distance=DISTANCE, mass=MASS)
        assert ec.ra != 0.0
        assert np.isclose(ec.dec, 0.0, atol=1e-12)

    def test_pure_north_offset(self):
        """x_BH ≠ 0, y_BH = z_BH = 0 → ra = 0, dec ≠ 0 (zeta=0)."""
        ec = EquatorialCoordinate(1000.0, 0.0, 0.0, zeta=0.0, eta=0.0,
                                  distance=DISTANCE, mass=MASS)
        assert np.isclose(ec.ra, 0.0, atol=1e-12)
        assert ec.dec != 0.0

    def test_ra_positive_east(self):
        """Positive y_BH (east) → positive RA offset."""
        ec = EquatorialCoordinate(0.0, 1000.0, 0.0, zeta=0.0, eta=0.0,
                                  distance=DISTANCE, mass=MASS)
        assert ec.ra > 0.0

    def test_dec_positive_north(self):
        """Positive x_BH (north, for zeta=0) → positive Dec offset."""
        ec = EquatorialCoordinate(1000.0, 0.0, 0.0, zeta=0.0, eta=0.0,
                                  distance=DISTANCE, mass=MASS)
        assert ec.dec > 0.0

    def test_known_angle_ra(self):
        """
        East offset of 1 kpc at distance 57.296 kpc → RA ≈ 1° (arctan(1/57.296)).

        57.296 kpc ≈ 180/π kpc so that arctan(1/57.296) ≈ 1°.
        """
        d = 57.296 * u.kpc
        ec = EquatorialCoordinate(0.0 * u.kpc, 1.0 * u.kpc, 0.0 * u.kpc,
                                  zeta=0.0, eta=0.0, distance=d)
        assert np.isclose(ec.ra, np.degrees(np.arctan(1.0 / 57.296)), rtol=1e-4)

    def test_known_angle_dec(self):
        """
        North offset of 1 kpc at distance 57.296 kpc → Dec ≈ 1°.
        """
        d = 57.296 * u.kpc
        ec = EquatorialCoordinate(1.0 * u.kpc, 0.0 * u.kpc, 0.0 * u.kpc,
                                  zeta=0.0, eta=0.0, distance=d)
        assert np.isclose(ec.dec, np.degrees(np.arctan(1.0 / 57.296)), rtol=1e-4)

    def test_antisymmetry(self):
        """Negating east offset negates RA."""
        d = 10.0 * u.kpc
        ec_pos = EquatorialCoordinate(0.0, 500.0, 0.0, zeta=0.0, eta=0.0,
                                      distance=d, mass=MASS)
        ec_neg = EquatorialCoordinate(0.0, -500.0, 0.0, zeta=0.0, eta=0.0,
                                      distance=d, mass=MASS)
        assert np.isclose(ec_pos.ra, -ec_neg.ra)


class TestEquatorialCoordinateUnits:
    """
    Unit-handling tests: astropy Quantities vs plain floats with mass.
    """

    def test_quantity_position_no_mass(self):
        """Quantity positions do not need mass."""
        ec = EquatorialCoordinate(0.0 * u.au, 1.0 * u.au, 0.0 * u.au,
                                  zeta=0.0, eta=0.0, distance=DISTANCE)
        assert isinstance(ec.ra, float)

    def test_float_position_with_mass(self):
        """Float positions in geometric units convert correctly when mass given."""
        ec = EquatorialCoordinate(0.0, 1000.0, 0.0, zeta=0.0, eta=0.0,
                                  distance=DISTANCE, mass=MASS)
        assert isinstance(ec.ra, float)

    def test_quantity_equals_float_with_mass(self):
        """
        Quantity position and equivalent float + mass produce the same RA/Dec.
        """
        from relatipy.numeric.constants import _L_ref
        offset_geom = 1000.0  # in GM_sun/c^2
        offset_m = offset_geom * MASS * _L_ref * u.m

        ec_qty = EquatorialCoordinate(0.0 * u.m, offset_m, 0.0 * u.m,
                                      zeta=0.0, eta=0.0, distance=DISTANCE)
        ec_float = EquatorialCoordinate(0.0, offset_geom, 0.0, zeta=0.0, eta=0.0,
                                        distance=DISTANCE, mass=MASS)
        assert np.isclose(ec_qty.ra, ec_float.ra, rtol=1e-10)
        assert np.isclose(ec_qty.dec, ec_float.dec, rtol=1e-10)


class TestEquatorialCoordinateFromAOE:
    """
    Tests for :meth:`EquatorialCoordinate.from_apparent_orbital_elements`.
    """

    def _make_aoe(self, zeta=45.0, eta=90.0):
        return ApparentOrbitalElements(
            0.0, 1e7, 0.1, 30.0, 40.0, 50.0, 60.0,
            zeta=zeta, eta=eta, mass=MASS,
        )

    def test_returns_equatorial_coordinate(self):
        """from_apparent_orbital_elements returns an EquatorialCoordinate."""
        aoe = self._make_aoe()
        ec = EquatorialCoordinate.from_apparent_orbital_elements(aoe, DISTANCE)
        assert isinstance(ec, EquatorialCoordinate)

    def test_zeta_eta_forwarded(self):
        """zeta and eta from the AOE are stored on the result."""
        aoe = self._make_aoe(zeta=30.0, eta=60.0)
        ec = EquatorialCoordinate.from_apparent_orbital_elements(aoe, DISTANCE)
        assert ec.zeta == 30.0
        assert ec.eta == 60.0

    def test_distance_stored(self):
        """distance passed to classmethod is stored on the result."""
        aoe = self._make_aoe()
        ec = EquatorialCoordinate.from_apparent_orbital_elements(aoe, DISTANCE)
        assert ec.distance == DISTANCE

    def test_identity_spin(self):
        """
        On-axis spin (zeta=0) — BH frame = sky frame.

        A particle at (0, r, 0) in BH frame (no x/z) has zero Dec offset.
        """
        aoe = ApparentOrbitalElements(
            0.0, 1e7, 0.0, 0.0, 0.0, 0.0, 0.0,
            zeta=0.0, eta=0.0, mass=MASS,
        )
        ec = EquatorialCoordinate.from_apparent_orbital_elements(aoe, DISTANCE)
        assert isinstance(ec, EquatorialCoordinate)


class TestEquatorialCoordinateValidation:
    """
    Input-validation error-path tests.
    """

    def test_distance_not_quantity_raises_type_error(self):
        """Plain float distance raises TypeError."""
        with pytest.raises(TypeError, match="astropy Quantity"):
            EquatorialCoordinate(0.0, 0.0, 0.0, zeta=0.0, eta=0.0, distance=8.5)

    def test_distance_wrong_units_raises_value_error(self):
        """Distance in non-length units raises ValueError."""
        with pytest.raises(ValueError, match="length units"):
            EquatorialCoordinate(0.0, 0.0, 0.0, zeta=0.0, eta=0.0,
                                 distance=8.5 * u.kg)

    def test_float_position_no_mass_raises_value_error(self):
        """Float position without mass raises ValueError."""
        with pytest.raises(ValueError, match="mass"):
            EquatorialCoordinate(1000.0, 0.0, 0.0, zeta=0.0, eta=0.0,
                                 distance=DISTANCE)

    def test_position_wrong_units_raises_type_error(self):
        """Quantity position in non-length units raises TypeError."""
        with pytest.raises(TypeError, match="length units"):
            EquatorialCoordinate(1.0 * u.kg, 0.0 * u.m, 0.0 * u.m,
                                 zeta=0.0, eta=0.0, distance=DISTANCE)


class TestEquatorialCoordinateDisplay:
    """
    Tests for display helpers.
    """

    def test_to_sexagesimal_returns_strings(self):
        """to_sexagesimal returns a tuple of two strings."""
        ec = EquatorialCoordinate(0.0, 1000.0, 0.0, zeta=0.0, eta=0.0,
                                  distance=DISTANCE, mass=MASS)
        ra_str, dec_str = ec.to_sexagesimal()
        assert isinstance(ra_str, str)
        assert isinstance(dec_str, str)

    def test_repr_contains_ra_dec(self):
        """__repr__ contains 'EquatorialCoordinate', 'ra', and 'dec'."""
        ec = EquatorialCoordinate(0.0, 0.0, 0.0, zeta=0.0, eta=0.0,
                                  distance=DISTANCE, mass=MASS)
        r = repr(ec)
        assert "EquatorialCoordinate" in r
        assert "ra=" in r
        assert "dec=" in r

    def test_eq_same_instance(self):
        """Two instances with identical inputs compare equal."""
        ec1 = EquatorialCoordinate(0.0, 1000.0, 0.0, zeta=0.0, eta=0.0,
                                   distance=DISTANCE, mass=MASS)
        ec2 = EquatorialCoordinate(0.0, 1000.0, 0.0, zeta=0.0, eta=0.0,
                                   distance=DISTANCE, mass=MASS)
        assert ec1 == ec2

    def test_eq_different_returns_false(self):
        """Two instances with different offsets are not equal."""
        ec1 = EquatorialCoordinate(0.0, 1000.0, 0.0, zeta=0.0, eta=0.0,
                                   distance=DISTANCE, mass=MASS)
        ec2 = EquatorialCoordinate(0.0, 2000.0, 0.0, zeta=0.0, eta=0.0,
                                   distance=DISTANCE, mass=MASS)
        assert ec1 != ec2
