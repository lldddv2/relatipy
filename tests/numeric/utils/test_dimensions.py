"""
Tests for :mod:`relatipy.numeric.utils.dimensions`.

These tests cover the default :obj:`~relatipy.numeric.utils.dimensions.validator`
instance and verify that physical :class:`astropy.units.Quantity` inputs are
converted to dimensionless geometric units, while plain floats are treated as
already normalized (except where validation applies).

Notes
-----
Reference scales :math:`M_{\\mathrm{ref}}`, :math:`L_{\\mathrm{ref}}`, and
:math:`T_{\\mathrm{ref}}` are imported from :mod:`relatipy.numeric.constants`.
"""

import astropy.units as u
import pytest

from relatipy.numeric.utils.dimensions import validator
from relatipy.numeric.constants import _c_SI, _M_ref_kg, _L_ref, _T_ref


class TestValidator:
    """
    Tests for :class:`~relatipy.numeric.utils.dimensions.Validator` methods.

    The default module-level ``validator`` singleton is used throughout.

    Examples
    --------
    Run this test module:

    >>> # pytest tests/numeric/utils/test_dimensions.py -q
    """

    def test_validate_mass(self):
        """
        Check mass normalization and rejection of negative values.

        Verifies that a plain float passes through unchanged, that one solar
        mass in kilograms maps to dimensionless mass 1, and that negative
        masses raise ``ValueError``.

        Raises
        ------
        AssertionError
            If any assertion fails.

        Examples
        --------
        >>> # pytest tests/numeric/utils/test_dimensions.py::TestValidator::test_validate_mass -q
        """
        assert validator.validate_mass(1) == 1

        assert validator.validate_mass(_M_ref_kg * u.kg) == 1

        with pytest.raises(ValueError):
            validator.validate_mass(-1)

    def test_validate_length(self):
        """
        Check length normalization from meters to geometric units.

        A plain float is unchanged. One reference length :math:`L_{\\mathrm{ref}}`
        in meters maps to dimensionless length 1.

        Raises
        ------
        AssertionError
            If any assertion fails.

        Examples
        --------
        >>> # pytest tests/numeric/utils/test_dimensions.py::TestValidator::test_validate_length -q
        """
        assert validator.validate_length(1) == 1

        assert pytest.approx(validator.validate_length(_L_ref * u.m), rel=1e-12) == 1

    def test_validate_time(self):
        """
        Check time normalization from seconds to geometric units.

        A plain float is unchanged. One reference time :math:`T_{\\mathrm{ref}}`
        in seconds maps to dimensionless time 1.

        Raises
        ------
        AssertionError
            If any assertion fails.

        Examples
        --------
        >>> # pytest tests/numeric/utils/test_dimensions.py::TestValidator::test_validate_time -q
        """
        assert validator.validate_time(1) == 1

        assert pytest.approx(validator.validate_time(_T_ref * u.s), rel=1e-12) == 1

    def test_validate_angle(self):
        """
        Check angle handling for floats and radian quantities.

        Both a dimensionless radian value and ``1 * u.rad`` normalize to 1.

        Raises
        ------
        AssertionError
            If any assertion fails.

        Examples
        --------
        >>> # pytest tests/numeric/utils/test_dimensions.py::TestValidator::test_validate_angle -q
        """
        assert validator.validate_angle(1) == 1
        assert validator.validate_angle(1 * u.rad) == 1

    def test_validate_velocity(self):
        """
        Check velocity normalization to :math:`v/c` and superluminal guard.

        A plain float in units of :math:`c` passes through. One SI speed of
        light maps to 1. Values above 1 in geometric units raise ``ValueError``.

        Raises
        ------
        AssertionError
            If any assertion fails.

        Examples
        --------
        >>> # pytest tests/numeric/utils/test_dimensions.py::TestValidator::test_validate_velocity -q
        """
        assert validator.validate_velocity(1) == 1

        assert validator.validate_velocity(_c_SI * u.m / u.s) == 1

        with pytest.raises(ValueError):
            validator.validate_velocity(1.1)

    def test_validate_angular_velocity(self):
        """
        Check angular velocity scaling by :math:`T_{\\mathrm{ref}}`.

        A plain float is unchanged. An angular rate of :math:`1/T_{\\mathrm{ref}}`
        in rad/s maps to dimensionless angular velocity 1.

        Raises
        ------
        AssertionError
            If any assertion fails.

        Examples
        --------
        >>> # pytest tests/numeric/utils/test_dimensions.py::TestValidator::test_validate_angular_velocity -q
        """
        assert validator.validate_angular_velocity(1) == 1

        assert pytest.approx(
            validator.validate_angular_velocity((1 / _T_ref) * u.rad / u.s),
            rel=1e-12,
        ) == 1
