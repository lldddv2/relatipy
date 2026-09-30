"""Array-backed records compare and hash by identity instead of raising."""

import numpy as np
from astropy import units as u

from relatipy.coordinates import CartesianCoordinates
from relatipy.geodesic import IntegrationInfo

from support.builders import make_state


def test_array_records_compare_by_identity_and_are_hashable() -> None:
    coordinates = CartesianCoordinates([0, 1] * u.s, [1, 2] * u.m, [1, 2] * u.m, [1, 2] * u.m)
    info = IntegrationInfo("Radau", 1e-3, np.ones(8), None, None, 1, 1)
    state = make_state()
    for record in (coordinates, info, state):
        assert record == record
        assert {record: 1}[record] == 1
    assert info != IntegrationInfo("Radau", 1e-3, np.ones(8), None, None, 1, 1)
