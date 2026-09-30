"""Internal validation helpers check types and freeze quantity data."""

import numpy as np
import pytest
from astropy import units as u

from relatipy import coordinates
from relatipy._validation import immutable_array, readonly_quantity, require_instance


def test_immutable_array_copies_into_read_only_contiguous_storage() -> None:
    source = np.arange(6.0).reshape(2, 3).T
    result = immutable_array(source)

    assert not result.flags.writeable
    assert result.flags.c_contiguous
    np.testing.assert_array_equal(result, source)
    source[0, 0] = 99.0
    assert result[0, 0] == 0.0
    with pytest.raises(ValueError):
        result.setflags(write=True)


@pytest.mark.parametrize("value", [["a", "b"], [True, False], np.array([None])])
def test_immutable_array_rejects_non_numeric(value: object) -> None:
    with pytest.raises(TypeError, match="numeric"):
        immutable_array(value)


def test_immutable_array_keeps_scalar_and_empty_shapes() -> None:
    assert immutable_array(1.5).shape == ()
    assert immutable_array([]).shape == (0,)


def test_readonly_quantity_validates_type_unit_and_ndim() -> None:
    with pytest.raises(TypeError, match="Quantity"):
        readonly_quantity(1.0, u.m, "d")
    with pytest.raises(u.UnitConversionError):
        readonly_quantity(1.0 * u.s, u.m, "d")
    with pytest.raises(ValueError, match="scalar or a one-dimensional"):
        readonly_quantity(np.ones((2, 2)) * u.m, u.m, "d")
    assert readonly_quantity(np.ones((2, 3)) * u.km, u.m, "d", ndim=(2,)).unit == u.km


def test_require_instance_names_field_and_expected_type() -> None:
    require_instance(3, int, "count")
    with pytest.raises(TypeError, match="count must be an instance of int"):
        require_instance("3", int, "count")


def test_coordinates_package_exposes_no_private_helpers() -> None:
    assert not [name for name in vars(coordinates) if name.startswith("_readonly")]
    assert all(not name.startswith("_") for name in coordinates.__all__)


def test_top_level_package_has_no_legacy_module_aliases() -> None:
    import relatipy

    for name in ("state", "solution", "exceptions", "orbital_elements", "numerical"):
        assert name not in relatipy.__all__
