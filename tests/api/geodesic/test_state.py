"""Contract tests for immutable scalar and vectorized states."""

import numpy as np
import pytest
from astropy import units as u

from relatipy.coordinates import (
    CartesianCoordinates,
    OrbitalElements,
    SphericalCoordinates,
)
from relatipy.geodesic import (
    InitialState,
    IntegrationError,
    IntegrationInfo,
    IntegrationTerminated,
    State,
    Termination,
)

from support.builders import make_state


def test_scalar_state_has_contract_shapes_units_and_delegates() -> None:
    state = make_state()

    assert state.tau.shape == ()
    assert state.txyz.xyz.shape == (3,)
    assert state.vxyz.shape == (3,)
    assert state.uxyz.shape == (3,)
    assert state.state_vector.xyz.shape == (3,)
    assert state.state_vector.vxyz.shape == (3,)
    assert state.x is state.txyz.x
    assert state.vr is state.vrqp.vr
    assert state.uPhi is state.uRQP.uPhi
    assert state.t.unit == u.s
    assert state.x.unit == u.km
    assert state.vx.unit == u.km / u.s
    assert state.orbital_elements().a.unit == u.km
    assert isinstance(state.orbital_elements().e, float)


def test_vectorized_state_has_contract_shapes() -> None:
    state = make_state(4)

    for component in (state.tau, state.t, state.x, state.vx, state.ut, state.uR):
        assert component.shape == (4,)
    assert state.xyz.shape == (4, 3)
    assert state.vxyz.shape == (4, 3)
    assert state.uxyz.shape == (4, 3)
    assert state.state_vector.xyz.shape == (4, 3)
    assert state.orbital_elements().a.shape == (4,)
    assert np.shape(state.orbital_elements().e) == (4,)


def test_initial_state_uses_the_same_validated_contract() -> None:
    initial = make_state(cls=InitialState)

    assert isinstance(initial, State)
    assert initial.tau == 0 * u.s


@pytest.mark.parametrize(
    ("factory", "error"),
    [
        (
            lambda: CartesianCoordinates(0 * u.s, 1 * u.km, 2 * u.km, 3 * u.s),
            u.UnitConversionError,
        ),
        (
            lambda: CartesianCoordinates(
                np.arange(2) * u.s,
                np.arange(3) * u.km,
                np.arange(2) * u.km,
                np.arange(2) * u.km,
            ),
            ValueError,
        ),
        (
            lambda: OrbitalElements(
                1 * u.km,
                0.1,
                1 * u.km,
                0 * u.rad,
                0 * u.rad,
                0 * u.rad,
            ),
            u.UnitConversionError,
        ),
    ],
)
def test_containers_reject_invalid_units_and_shapes(factory, error) -> None:
    with pytest.raises(error):
        factory()


def test_state_rejects_inconsistent_cardinality() -> None:
    scalar = make_state()
    with pytest.raises(ValueError, match="vxyz has state shape"):
        State(
            tau=scalar.tau,
            txyz=scalar.txyz,
            trqp=scalar.trqp,
            tRQP=scalar.tRQP,
            vxyz=np.ones((2, 3)) * u.km / u.s,
            vrqp=scalar.vrqp,
            vRQP=scalar.vRQP,
            ut=scalar.ut,
            uxyz=scalar.uxyz,
            urqp=scalar.urqp,
            uRQP=scalar.uRQP,
            orbital_elements=scalar.orbital_elements(),
        )


def test_state_rejects_disagreeing_coordinate_times() -> None:
    scalar = make_state()
    with pytest.raises(ValueError, match="same coordinate time"):
        State(
            tau=scalar.tau,
            txyz=scalar.txyz,
            trqp=SphericalCoordinates(
                1 * u.s,
                scalar.r,
                scalar.theta,
                scalar.phi,
            ),
            tRQP=scalar.tRQP,
            vxyz=scalar.vxyz,
            vrqp=scalar.vrqp,
            vRQP=scalar.vRQP,
            ut=scalar.ut,
            uxyz=scalar.uxyz,
            urqp=scalar.urqp,
            uRQP=scalar.uRQP,
            orbital_elements=scalar.orbital_elements(),
        )


def test_integration_info_validates_configuration_and_diagnostics() -> None:
    info = IntegrationInfo(
        method="TEST",
        rtol=1e-8,
        atol=1e-10,
        max_step=2 * u.s,
        first_step=0.1 * u.s,
        n_steps=8,
        nfev=32,
    )

    assert info.method == "TEST"
    assert info.max_step == 2 * u.s
    with pytest.raises(u.UnitConversionError):
        IntegrationInfo("TEST", 1e-8, 1e-10, 2 * u.km, None, 8, 32)
    with pytest.raises(ValueError, match="non-negative"):
        IntegrationInfo("TEST", 1e-8, 1e-10, None, None, -1, 32)
    for invalid_step in (np.nan * u.s, np.inf * u.s):
        with pytest.raises(ValueError, match="finite and positive"):
            IntegrationInfo("TEST", 1e-8, 1e-10, invalid_step, None, 8, 32)
        with pytest.raises(ValueError, match="finite and positive"):
            IntegrationInfo("TEST", 1e-8, 1e-10, None, invalid_step, 8, 32)


def test_termination_requires_a_scalar_state() -> None:
    termination = Termination("event", 0 * u.s, make_state())
    assert termination.reason == "event"
    with pytest.raises(ValueError, match="scalar"):
        Termination("event", 2 * u.s, make_state(2))
    with pytest.raises(ValueError, match="match state.tau"):
        Termination("event", 2 * u.s, make_state())


def test_integration_exceptions_have_the_closed_public_shape() -> None:
    state = make_state()
    error = IntegrationError("numerical failure")
    terminated = IntegrationTerminated("event", 0 * u.s, state)

    assert str(error) == "numerical failure"
    assert terminated.reason == "event"
    assert terminated.tau == 0 * u.s
    assert terminated.state is state
    with pytest.raises(AttributeError):
        terminated.reason = "different"
    with pytest.raises(ValueError, match="match state.tau"):
        IntegrationTerminated("event", 2 * u.s, state)
    assert isinstance(terminated.termination, Termination)
    assert terminated.termination.state is state
    assert str(IntegrationTerminated("event", 0 * u.s, state, "detail")) == "detail"
