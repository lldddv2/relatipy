"""Orbit and Solution expose the same documented State properties."""

from relatipy.geodesic import Orbit, Solution, State
from relatipy.geodesic._delegation import STATE_PROPERTIES


def test_orbit_and_solution_delegate_every_state_property() -> None:
    for name in STATE_PROPERTIES:
        assert hasattr(State, name)
        for owner in (Orbit, Solution):
            descriptor = getattr(owner, name)
            assert isinstance(descriptor, property)
            assert descriptor.fset is None
            assert descriptor.__doc__
