"""Unit validation tests for the public domain model."""

import numpy as np
import pytest
from astropy import units as u

from relatipy.coordinates import (
    BoyerLindquistFourVelocity,
    CartesianCoordinates,
    OrbitalElements,
    SphericalFourVelocity,
    SphericalVelocity,
)
from relatipy.geodesic import State

from support.builders import make_state


def test_compatible_units_are_preserved_without_silent_reformatting() -> None:
    coordinates = CartesianCoordinates(
        np.arange(2) * u.day,
        np.arange(2) * u.au,
        np.arange(2) * u.km,
        np.arange(2) * u.pc,
    )

    assert coordinates.t.unit == u.day
    assert coordinates.x.unit == u.au
    assert coordinates.y.unit == u.km
    assert coordinates.z.unit == u.pc
    assert coordinates.xyz.unit == u.au


def test_heterogeneous_velocity_container_preserves_linear_and_angular_units() -> None:
    velocity = SphericalVelocity(
        3 * u.km / u.s,
        4 * u.deg / u.s,
        5 * u.rad / u.day,
    )

    assert velocity.vr.unit == u.km / u.s
    assert velocity.vtheta.unit == u.deg / u.s
    assert velocity.vphi.unit == u.rad / u.day


def test_eccentricity_accepts_only_dimensionless_quantities() -> None:
    elements = OrbitalElements(
        10 * u.km,
        20 * u.percent,
        1 * u.deg,
        2 * u.deg,
        3 * u.deg,
        4 * u.deg,
    )
    assert elements.e == pytest.approx(0.2)

    with pytest.raises(u.UnitConversionError):
        OrbitalElements(
            10 * u.km,
            0.2 * u.km,
            1 * u.deg,
            2 * u.deg,
            3 * u.deg,
            4 * u.deg,
        )


def test_four_velocity_units_follow_proper_time_contract() -> None:
    spherical = SphericalFourVelocity(
        3 * u.km / u.s,
        4 * u.deg / u.s,
        5 * u.rad / u.day,
    )
    boyer_lindquist = BoyerLindquistFourVelocity(
        6 * u.au / u.yr,
        7 * u.deg / u.s,
        8 * u.rad / u.day,
    )

    assert spherical.ur.unit == u.km / u.s
    assert spherical.utheta.unit == u.deg / u.s
    assert boyer_lindquist.uR.unit == u.au / u.yr
    assert boyer_lindquist.uPhi.unit == u.rad / u.day

    scalar = make_state()
    with pytest.raises(u.UnitConversionError):
        State(
            tau=scalar.tau,
            txyz=scalar.txyz,
            trqp=scalar.trqp,
            tRQP=scalar.tRQP,
            vxyz=scalar.vxyz,
            vrqp=scalar.vrqp,
            vRQP=scalar.vRQP,
            ut=1 * u.s,
            uxyz=scalar.uxyz,
            urqp=scalar.urqp,
            uRQP=scalar.uRQP,
            orbital_elements=scalar.orbital_elements(),
        )
