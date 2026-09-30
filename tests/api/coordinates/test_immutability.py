"""Mutation-resistance tests for public domain objects."""

from dataclasses import FrozenInstanceError

import numpy as np
import pytest
from astropy import units as u

from relatipy.coordinates import CartesianCoordinates, OrbitalElements
from relatipy.geodesic import IntegrationInfo

from support.builders import make_state


def test_frozen_attributes_reject_assignment() -> None:
    state = make_state()
    with pytest.raises((AttributeError, FrozenInstanceError)):
        state.tau = 3 * u.s
    with pytest.raises((AttributeError, FrozenInstanceError)):
        state.txyz.x = 3 * u.km
    with pytest.raises((AttributeError, FrozenInstanceError)):
        state.orbital_elements().e = 0.8


@pytest.mark.parametrize(
    "extract",
    [
        lambda state: state.tau,
        lambda state: state.txyz.x,
        lambda state: state.txyz.xyz,
        lambda state: state.vxyz,
        lambda state: state.uxyz,
        lambda state: state.state_vector.xyz,
        lambda state: state.orbital_elements().a,
        lambda state: state.orbital_elements().e,
    ],
)
def test_public_arrays_reject_item_assignment_and_writeable_flag(extract) -> None:
    value = extract(make_state(3))

    with pytest.raises(ValueError):
        value[0] = 0
    array = value.value if isinstance(value, u.Quantity) else value
    with pytest.raises(ValueError):
        array.flags.writeable = True


def test_constructor_copies_mutable_inputs() -> None:
    source = np.array([1.0, 2.0])
    coordinates = CartesianCoordinates(
        source * u.s,
        source * u.km,
        source * u.km,
        source * u.km,
    )
    elements = OrbitalElements(
        source * u.km,
        source,
        source * u.rad,
        source * u.rad,
        source * u.rad,
        source * u.rad,
    )

    source[:] = 99.0
    np.testing.assert_array_equal(coordinates.x.value, [1.0, 2.0])
    np.testing.assert_array_equal(elements.e, [1.0, 2.0])


def test_integration_info_is_frozen_and_owns_step_quantities() -> None:
    source = np.array(2.0) * u.s
    info = IntegrationInfo("TEST", 1e-8, 1e-10, source, None, 1, 2)

    with pytest.raises((AttributeError, FrozenInstanceError)):
        info.n_steps = 4
    source[...] = 7.0 * u.s
    assert info.max_step == 2.0 * u.s
