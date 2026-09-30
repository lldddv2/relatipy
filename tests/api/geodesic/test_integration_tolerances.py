"""Absolute-tolerance storage at the public diagnostic boundary."""

from dataclasses import FrozenInstanceError

import numpy as np
import pytest
from astropy import units as u

from relatipy.geodesic import IntegrationInfo


def make_info(atol):
    return IntegrationInfo("TEST", 1e-8, atol, None, None, 1, 2)


@pytest.mark.parametrize("value", [0, 2, 0.25, np.float64(0.5), np.array(0.75)])
def test_scalar_absolute_tolerance_is_a_float(value) -> None:
    info = make_info(value)
    assert isinstance(info.atol, float)
    assert info.atol == float(value)


def test_vector_absolute_tolerance_is_copied_and_strictly_read_only() -> None:
    source = np.array([0.0, 1e-6, 2e-6])
    info = make_info(source)
    source[:] = 9.0

    np.testing.assert_array_equal(info.atol, [0.0, 1e-6, 2e-6])
    assert info.atol.dtype == np.dtype(float)
    assert not info.atol.flags.writeable
    with pytest.raises(ValueError):
        info.atol[0] = 3.0
    with pytest.raises(ValueError):
        info.atol.flags.writeable = True
    with pytest.raises((AttributeError, FrozenInstanceError)):
        info.atol = np.array([3.0])


@pytest.mark.parametrize("value", [[0], (1, 2), np.array([1, 2])])
def test_nonempty_vectors_accept_multiple_lengths_without_state_assumptions(
    value,
) -> None:
    info = make_info(value)
    np.testing.assert_array_equal(info.atol, value)
    assert info.atol.ndim == 1


@pytest.mark.parametrize("value", [-1, np.nan, np.inf, -np.inf])
def test_invalid_scalar_value_is_rejected(value) -> None:
    with pytest.raises(ValueError, match="finite and non-negative"):
        make_info(value)


@pytest.mark.parametrize("value", [[1, -1], [0, np.nan], [np.inf], [-np.inf]])
def test_invalid_vector_component_is_rejected(value) -> None:
    with pytest.raises(ValueError, match="finite and non-negative"):
        make_info(value)


@pytest.mark.parametrize("value", [[], np.array([]), [[1, 2]], np.ones((2, 2))])
def test_invalid_vector_shape_is_rejected(value) -> None:
    with pytest.raises(ValueError, match="nonempty one-dimensional"):
        make_info(value)


@pytest.mark.parametrize(
    "value",
    ["1", ["1"], 1 + 2j, [1 + 2j], True, [True], [True, 1],
     1 * u.one, [1, 2] * u.one, [1 * u.one],
     (np.array(1) * u.one,), object()],
)
def test_non_real_or_unit_bearing_tolerance_is_rejected(value) -> None:
    with pytest.raises(TypeError, match="atol"):
        make_info(value)
