"""Solution contracts: read-only delegation, indexing and stored-series validation."""

from __future__ import annotations

import numpy as np
import pytest
from astropy import units as u

from relatipy.geodesic import Solution, State, Termination

from support.builders import make_series_state, make_solution


def test_solution_delegates_read_only_series() -> None:
    solution = make_solution()
    assert len(solution) == 4
    assert solution.success is True
    assert solution.xyz.shape == (4, 3)
    assert solution.state_vector.vxyz.shape == (4, 3)
    assert solution.orbital_elements().a.shape == (4,)
    with pytest.raises(ValueError):
        solution.xyz[0, 0] = 0 * u.km


def test_integer_index_returns_scalar_state() -> None:
    selected = make_solution()[-1]
    assert isinstance(selected, State)
    assert selected.tau.shape == ()
    assert selected.xyz.shape == (3,)
    assert selected.orbital_elements().a.shape == ()


@pytest.mark.parametrize(
    "index, expected",
    [
        (slice(1, 3), [1.0, 2.0]),
        (np.array([3, 1, 1]), [3.0, 1.0, 1.0]),
        (np.array([True, False, True, False]), [0.0, 2.0]),
    ],
)
def test_vector_selection_preserves_numpy_order(index, expected) -> None:
    selected = make_solution()[index]
    assert selected.tau.shape == (len(expected),)
    np.testing.assert_allclose(selected.tau.to_value(u.s), expected)


def test_invalid_indices_are_rejected_clearly() -> None:
    solution = make_solution()
    with pytest.raises(IndexError, match="one-dimensional"):
        solution[np.array([[0, 1]])]
    with pytest.raises(IndexError, match="boolean mask"):
        solution[np.array([True, False])]
    with pytest.raises(IndexError, match="integers or booleans"):
        solution[np.array([0.5])]


def test_status_and_termination_are_consistent() -> None:
    with pytest.raises(ValueError, match="termination is required"):
        make_solution(status=1)
    with pytest.raises(ValueError, match="must be None"):
        make_solution(
            status=0,
            termination=Termination(
                reason="event",
                tau=0 * u.s,
                state=make_solution()[0],
            ),
        )
    failed = make_solution(status=-1)
    assert failed.success is False
    termination = Termination(
        reason="event",
        tau=0 * u.s,
        state=make_solution()[0],
    )
    terminated = make_solution(status=1, termination=termination)
    assert terminated.success is True
    assert terminated.termination is termination


def test_solution_requires_increasing_nonempty_series() -> None:
    scalar = make_solution()[0]
    with pytest.raises(ValueError, match="one-dimensional"):
        Solution(
            state=scalar,
            integration=make_solution().integration,
            status=0,
            message="scalar",
        )


@pytest.mark.parametrize("bad_tau", [[0.0, np.nan, 2.0, 3.0], [0.0, 1.0, 2.0, np.inf]])
def test_solution_rejects_non_finite_proper_times(bad_tau) -> None:
    state = make_series_state()
    state_with_bad_tau = State(
        tau=np.asarray(bad_tau) * u.s,
        txyz=state.txyz,
        trqp=state.trqp,
        tRQP=state.tRQP,
        vxyz=state.vxyz,
        vrqp=state.vrqp,
        vRQP=state.vRQP,
        ut=state.ut,
        uxyz=state.uxyz,
        urqp=state.urqp,
        uRQP=state.uRQP,
        orbital_elements=state.orbital_elements(),
    )
    with pytest.raises(ValueError, match="finite"):
        Solution(
            state=state_with_bad_tau,
            integration=make_solution().integration,
            status=0,
            message="invalid domain",
        )


@pytest.mark.parametrize(
    "index",
    [slice(1, 3), np.array([3, 1, 1]), np.array([True, False, True, False])],
)
def test_selected_states_remain_strictly_read_only(index) -> None:
    selected = make_solution()[index]
    with pytest.raises(ValueError):
        selected.xyz[0, 0] = 0 * u.km
    with pytest.raises(ValueError):
        selected.xyz.value.flags.writeable = True
