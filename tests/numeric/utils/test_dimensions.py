import astropy.units as u
import pytest

from relatipy.numeric.utils.dimensions import validator
from relatipy.numeric.constants import _G_SI, _c_SI, _M_ref_kg, _L_ref, _T_ref


class TestValidator:
    def test_validate_mass(self):
        # sin unidades → se deja igual
        assert validator.validate_mass(1) == 1

        # 1 masa solar en kg → 1 en unidades geométricas
        assert validator.validate_mass(_M_ref_kg * u.kg) == 1

        with pytest.raises(ValueError):
            validator.validate_mass(-1)

    def test_validate_length(self):
        # sin unidades → se deja igual
        assert validator.validate_length(1) == 1

        # 1 unidad geométrica de longitud L_ref = G M_ref / c^2 en metros → 1
        assert pytest.approx(validator.validate_length(_L_ref * u.m), rel=1e-12) == 1

    def test_validate_time(self):
        # sin unidades → se deja igual
        assert validator.validate_time(1) == 1

        # 1 unidad geométrica de tiempo T_ref = G M_ref / c^3 en segundos → 1
        assert pytest.approx(validator.validate_time(_T_ref * u.s), rel=1e-12) == 1

    def test_validate_angle(self):
        assert validator.validate_angle(1) == 1
        assert validator.validate_angle(1 * u.rad) == 1

    def test_validate_velocity(self):
        # adimensional ya en unidades de c
        assert validator.validate_velocity(1) == 1

        # 1 * c_SI → 1 en unidades geométricas
        assert validator.validate_velocity(_c_SI * u.m / u.s) == 1

        # outrange, velocity must be less than the speed of light
        with pytest.raises(ValueError):
            validator.validate_velocity(1.1)

    def test_validate_angular_velocity(self):
        # adimensional ya en unidades geométricas
        assert validator.validate_angular_velocity(1) == 1

        # rad/s → adimensional mediante factor T_ref
        assert pytest.approx(
            validator.validate_angular_velocity((1 / _T_ref) * u.rad / u.s),
            rel=1e-12,
        ) == 1