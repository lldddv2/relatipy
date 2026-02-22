import astropy.units as u
from relatipy.numeric.utils.dimensions import validator
from relatipy.numeric.constants import _G_SI, _c_SI
import pytest

class TestValidator:
    def test_validate_mass(self):
        assert validator.validate_mass(1) == 1 # sin unidades
        assert validator.validate_mass(_c_SI**2 / _G_SI * u.kg) == 1 # en unidades SI
        with pytest.raises(ValueError):
            validator.validate_mass(-1)

    def test_validate_length(self):
        assert validator.validate_length(1) == 1
        assert validator.validate_length(1 * u.m) == 1

    def test_validate_time(self):
        assert validator.validate_time(1) == 1
        assert validator.validate_time(_c_SI * u.s) == 1

    def test_validate_angle(self):
        assert validator.validate_angle(1) == 1
        assert validator.validate_angle(1 * u.rad) == 1

    def test_validate_velocity(self):
        assert validator.validate_velocity(1) == 1
        assert validator.validate_velocity(_c_SI * u.m / u.s) == 1
        # outrange, velocity must be less than the speed of light
        with pytest.raises(ValueError):
            validator.validate_velocity(1.1)

    def test_validate_angular_velocity(self):
        assert validator.validate_angular_velocity(1) == 1
        assert validator.validate_angular_velocity(_c_SI * u.rad / u.s) == 1