"""Public coordinate-time integration, termination and tolerance contracts."""

import warnings

import numpy as np
import pytest
from astropy import units as u
from astropy.constants import c

import relatipy as rp


def _radial(sign=1, *, t=0 * u.s):
    metric = rp.Kerr(mass=1 * u.Msun, spin=0)
    photon = metric.null(
        R=10 * metric.r_g, Theta=np.pi / 2 * u.rad, Phi=0 * u.rad,
        vR=sign * c, t=t,
    )
    return metric, photon, (metric.r_g / c).to(u.s)


@pytest.mark.parametrize("method", ["radau", "dop853", "dp45"])
def test_integrate_reaches_exact_future_coordinate_time(method) -> None:
    _, photon, scale = _radial(t=2 * u.ms)
    target = photon.t + 0.1 * scale
    assert photon.integrate(target, method=method) is None
    assert photon.t == target


@pytest.mark.parametrize("offset", [0, -1])
def test_integrate_rejects_present_and_past(offset) -> None:
    _, photon, scale = _radial(t=2 * u.ms)
    with pytest.raises(ValueError, match="future"):
        photon.integrate(photon.t + offset * scale, method="dp45")


@pytest.mark.parametrize(
    "target, error",
    [(1.0, TypeError), (np.nan * u.s, ValueError),
     (np.inf * u.s, ValueError), ([1, 2] * u.s, ValueError)],
)
def test_integrate_requires_finite_scalar_time_quantity(target, error) -> None:
    _, photon, _ = _radial()
    with pytest.raises(error):
        photon.integrate(target, method="dp45")


@pytest.mark.parametrize("operation", ["integrate", "solve"])
@pytest.mark.parametrize("method", ["projection_radau", "unknown"])
def test_unsupported_methods_are_rejected(operation, method) -> None:
    _, photon, scale = _radial()
    with pytest.raises(ValueError):
        if operation == "integrate":
            photon.integrate(scale, method=method)
        else:
            photon.solve(t_span=(photon.t, scale), method=method)


@pytest.mark.parametrize("operation", ["integrate", "solve"])
@pytest.mark.parametrize("option", ["first_step", "max_step"])
def test_affine_step_controls_are_not_public(operation, option) -> None:
    _, photon, scale = _radial()
    with pytest.raises(TypeError):
        if operation == "integrate":
            photon.integrate(scale, method="dp45", **{option: scale / 10})
        else:
            photon.solve(t_span=(photon.t, scale), method="dp45",
                         **{option: scale / 10})


@pytest.mark.parametrize("operation", ["integrate", "solve"])
@pytest.mark.parametrize("radius, error", [(9, ValueError), (10, ValueError),
                                          (None, TypeError)])
def test_escape_radius_must_be_a_length_above_current_radius(
    operation, radius, error,
) -> None:
    metric, photon, scale = _radial()
    escape = 50.0 if radius is None else radius * metric.r_g
    with pytest.raises(error):
        if operation == "integrate":
            photon.integrate(scale, method="dp45", r_escape=escape)
        else:
            photon.solve(t_span=(photon.t, scale), method="dp45", r_escape=escape)


def test_integrate_validates_escape_against_current_point() -> None:
    metric, photon, scale = _radial()
    photon.integrate(2 * scale, method="dp45")
    assert photon.R > 10 * metric.r_g
    with pytest.raises(ValueError):
        photon.integrate(3 * scale, method="dp45", r_escape=10.5 * metric.r_g)


@pytest.mark.parametrize("reason, sign", [("horizon", -1), ("escape", 1)])
def test_integrate_terminal_events_preserve_last_valid_state(reason, sign) -> None:
    metric, photon, scale = _radial(sign)
    target = 100 * scale
    options = {"r_escape": 50 * metric.r_g} if reason == "escape" else {}
    with pytest.raises(rp.IntegrationTerminated) as caught:
        photon.integrate(target, method="dp45", **options)
    event = caught.value
    assert event.reason == reason
    assert event.t == photon.t == event.state.t
    assert event.termination.reason == reason
    assert event.termination.t == event.t
    assert event.termination.state.t == event.state.t
    np.testing.assert_array_equal(event.state.xyz.value, photon.xyz.value)
    assert photon.t < target
    assert photon.R > metric.horizons.event
    if reason == "escape":
        assert photon.R >= 50 * metric.r_g
    assert not hasattr(event, "tau")


def test_copy_reset_initial_and_invariants_are_independent_of_evolution() -> None:
    metric = rp.Kerr(mass=1 * u.Msun, spin=0.7)
    scale = (metric.r_g / c).to(u.s)
    photon = metric.null(
        R=10 * metric.r_g, Theta=1.1 * u.rad, Phi=0.2 * u.rad,
        b=3 * metric.r_g, eta=5 * metric.r_g**2,
        radial_sign=1, polar_sign=-1,
    )
    initial_t = photon.initial.t.copy()
    initial_xyz = photon.initial.xyz.copy()
    initial_velocity = photon.initial.vxyz.copy()
    b, eta = photon.b.copy(), photon.eta.copy()
    photon.integrate(0.1 * scale, method="dp45")
    assert photon.b == b
    assert photon.eta == eta
    assert photon.initial.t == initial_t
    np.testing.assert_array_equal(photon.initial.xyz.value, initial_xyz.value)
    np.testing.assert_array_equal(photon.initial.vxyz.value, initial_velocity.value)
    copied = photon.copy()
    assert copied is not photon
    np.testing.assert_array_equal(copied.xyz.value, photon.xyz.value)
    assert copied.t == photon.t
    copied.integrate(0.2 * scale, method="dp45")
    assert photon.t == 0.1 * scale
    copied.reset()
    assert copied.t == initial_t
    np.testing.assert_array_equal(copied.xyz.value, initial_xyz.value)
    np.testing.assert_array_equal(copied.vxyz.value, initial_velocity.value)
    assert photon.t == 0.1 * scale
    photon.reset()
    assert photon.t == initial_t
    np.testing.assert_array_equal(photon.xyz.value, initial_xyz.value)
    np.testing.assert_array_equal(photon.vxyz.value, initial_velocity.value)


@pytest.mark.parametrize("method", ["radau", "dop853", "dp45"])
def test_solve_records_effective_options_and_exact_endpoint(method) -> None:
    _, photon, scale = _radial()
    target = 0.1 * scale
    solution = photon.solve(t_span=(photon.t, target), method=method)
    assert isinstance(solution, rp.NullSolution)
    assert solution.status == 0
    assert solution.success is True
    assert isinstance(solution.message, str)
    assert solution.termination is None
    assert solution.t[-1] == target
    assert solution.integration.method == method
    assert solution.integration.rtol == 1e-10
    assert np.shape(solution.integration.atol) == (8,)
    assert solution.integration.first_step is None
    assert solution.integration.max_step is None
    assert photon.t == photon.initial.t


def test_default_solve_method_is_radau() -> None:
    _, photon, scale = _radial()
    assert photon.solve(t_span=(photon.t, 0.1 * scale)).integration.method == "radau"


def test_solve_uses_initial_state_without_mutating_current_point() -> None:
    _, photon, scale = _radial()
    photon.integrate(0.1 * scale, method="dp45")
    current_t = photon.t.copy()
    current_xyz = photon.xyz.copy()
    current_velocity = photon.vxyz.copy()
    solution = photon.solve(t_span=(photon.initial.t, 0.2 * scale), method="dp45")
    assert solution.t[0] == photon.initial.t
    np.testing.assert_array_equal(solution.xyz[0].value, photon.initial.xyz.value)
    assert photon.t == current_t
    np.testing.assert_array_equal(photon.xyz.value, current_xyz.value)
    np.testing.assert_array_equal(photon.vxyz.value, current_velocity.value)


@pytest.mark.parametrize("reason, sign", [("horizon", -1), ("escape", 1)])
def test_solve_terminal_events_return_structured_success(reason, sign) -> None:
    metric, photon, scale = _radial(sign)
    options = {"r_escape": 50 * metric.r_g} if reason == "escape" else {}
    solution = photon.solve(t_span=(photon.t, 100 * scale), method="dp45", **options)
    assert solution.status == 1
    assert solution.success is True
    assert isinstance(solution.message, str)
    event = solution.termination
    assert event.reason == reason
    assert event.t == event.state.t == solution.t[-1]
    assert event.state.R > metric.horizons.event
    assert event.t < 100 * scale
    if reason == "escape":
        assert event.state.R >= 50 * metric.r_g
    assert photon.t == photon.initial.t


def test_solve_requires_span_or_samples() -> None:
    _, photon, _ = _radial()
    with pytest.raises(ValueError):
        photon.solve(method="dp45")


@pytest.mark.parametrize("start, stop, match", [(1, 2, None),
                                             (0, 0, "future"), (0, -1, "future")])
def test_solve_span_starts_at_initial_and_ends_in_future(start, stop, match) -> None:
    _, photon, scale = _radial()
    with pytest.raises(ValueError, match=match):
        photon.solve(t_span=(start * scale, stop * scale), method="dp45")


@pytest.mark.parametrize(
    "samples",
    [np.array(0.1), np.array([]), np.array([0, np.nan]),
     np.array([0, np.inf]), np.array([0, 0]), np.array([0.2, 0.1]),
     np.array([[0, 0.1]]), np.array([-0.1, 0.1]), np.array([0, 0.3])],
)
def test_solve_rejects_invalid_sample_shapes_values_and_domain(samples) -> None:
    _, photon, scale = _radial()
    with pytest.raises(ValueError):
        photon.solve(t_span=(photon.t, 0.2 * scale), t_eval=samples * scale,
                     method="dp45")


def test_solve_samples_require_units() -> None:
    _, photon, _ = _radial()
    with pytest.raises(TypeError):
        photon.solve(t_eval=[0, 1], method="dp45")


@pytest.mark.parametrize("with_span", [False, True])
def test_solve_returns_exact_requested_sample_times(with_span) -> None:
    _, photon, scale = _radial()
    times = np.array([0, 0.03, 0.07, 0.1]) * scale
    options = {"t_span": (photon.t, times[-1])} if with_span else {}
    solution = photon.solve(t_eval=times, method="dp45", **options)
    assert solution.status == 0
    np.testing.assert_array_equal(solution.t.to_value(times.unit), times.value)


def test_single_future_sample_is_valid_but_initial_only_is_not() -> None:
    _, photon, scale = _radial()
    times = np.array([0.1]) * scale
    solution = photon.solve(t_eval=times, method="dp45")
    np.testing.assert_array_equal(solution.t.to_value(times.unit), times.value)
    with pytest.raises(ValueError):
        photon.solve(t_eval=np.array([0.0]) * scale, method="dp45")


def test_horizon_keeps_only_reached_samples_and_terminal_state() -> None:
    metric, photon, scale = _radial(-1)
    times = np.array([0, 1, 100]) * scale
    solution = photon.solve(t_eval=times, method="dp45")
    assert solution.status == 1
    np.testing.assert_array_equal(solution.t.to_value(times.unit), times[:2].value)
    event = solution.termination
    assert event.reason == "horizon"
    assert event.t == event.state.t
    assert event.t > solution.t[-1]
    assert event.state.R > metric.horizons.event
    assert event.t < times[-1]


def test_horizon_before_first_sample_raises() -> None:
    _, photon, scale = _radial(-1)
    with pytest.raises(rp.IntegrationTerminated) as caught:
        photon.solve(t_eval=np.array([100, 101]) * scale, method="dp45")
    assert caught.value.reason == "horizon"
    assert caught.value.state.t == caught.value.t
    assert photon.t == photon.initial.t


@pytest.mark.parametrize("operation", ["integrate", "solve"])
def test_explicit_loose_tolerance_warns_without_kepler_warning(operation) -> None:
    _, photon, scale = _radial()
    with pytest.warns(rp.IntegrationWarning) as caught:
        if operation == "integrate":
            photon.integrate(0.1 * scale, method="dp45", rtol=1e-3)
        else:
            solution = photon.solve(t_span=(photon.t, 0.1 * scale),
                                    method="dp45", rtol=1e-3)
            assert solution.integration.rtol == 1e-3
    assert all("Kepler" not in str(item.message) for item in caught)


@pytest.mark.parametrize("operation", ["integrate", "solve"])
def test_automatic_tolerances_are_silent(operation) -> None:
    _, photon, scale = _radial()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        if operation == "integrate":
            photon.integrate(0.1 * scale, method="dp45")
        else:
            photon.solve(t_eval=np.array([0, 0.1]) * scale, method="dp45")
    assert not caught


@pytest.mark.parametrize("operation", ["integrate", "solve"])
def test_absolute_tolerance_vector_has_eight_components(operation) -> None:
    _, photon, scale = _radial()
    with pytest.raises(ValueError):
        if operation == "integrate":
            photon.integrate(scale, method="dp45", atol=np.full(7, 1e-12))
        else:
            photon.solve(t_span=(photon.t, scale), method="dp45",
                         atol=np.full(7, 1e-12))


def _tangential_failure_case():
    """Use finite validated inputs whose Radau error estimate overflows."""
    metric = rp.Kerr(mass=1 * u.Msun, spin=0)
    photon = metric.null(
        R=10 * metric.r_g, Theta=np.pi / 2 * u.rad, Phi=0 * u.rad,
        vPhi=1 * u.rad / u.s,
    )
    return photon, 0.001 * (metric.r_g / c).to(u.s)


@pytest.mark.parametrize("sampling", ["steps", "initial", "future"])
def test_nonfinite_numerical_failure_keeps_partial_solve(sampling) -> None:
    photon, target = _tangential_failure_case()
    options = {"t_span": (photon.initial.t, target)}
    if sampling == "initial":
        options["t_eval"] = np.array([0, 1]) * target
    elif sampling == "future":
        options["t_eval"] = np.array([0.5, 1]) * target
    with pytest.warns(rp.IntegrationWarning):
        if sampling == "future":
            with pytest.raises(rp.IntegrationError, match="non-finite numerical value; native status 3"):
                photon.solve(method="radau", rtol=1e-170, atol=1e-170, **options)
        else:
            solution = photon.solve(method="radau", rtol=1e-170, atol=1e-170, **options)
            assert solution.status == -1
            assert not solution.success
            assert len(solution) >= 1
            assert solution.t[0] == photon.initial.t
            np.testing.assert_array_equal(solution[0]._canonical, photon.initial._canonical)
            assert "non-finite numerical value; native status 3" in solution.message
    np.testing.assert_array_equal(photon._current_y, photon._initial_y)


def test_nonfinite_numerical_failure_integrate_retains_last_valid_point() -> None:
    photon, target = _tangential_failure_case()
    initial = photon._initial_y.copy()
    with pytest.warns(rp.IntegrationWarning):
        with pytest.raises(rp.IntegrationError, match="non-finite numerical value; native status 3"):
            photon.integrate(target, method="radau", rtol=1e-170, atol=1e-170)
    np.testing.assert_array_equal(photon._current_y, initial)
    assert photon.t == photon.initial.t


def test_nonfinite_native_failure_returns_owned_readonly_last_state() -> None:
    from relatipy import _core

    photon, _ = _tangential_failure_case()
    initial = photon._initial_y.copy()
    states, final, stats, status, reason = _core.integrate_kerr_null(
        photon._metric.spin, initial, 0.001, None, 0, "radau", 1e-170, 1e-170,
    )
    assert status == -1
    assert stats["native_status"] == 3
    assert reason is None
    assert len(states) >= 1
    np.testing.assert_array_equal(states[0], initial)
    np.testing.assert_array_equal(final, initial)
    assert not np.shares_memory(final, initial)
    assert not np.shares_memory(states, initial)
    assert not np.shares_memory(states, final)
    assert not states.flags.writeable
    assert not final.flags.writeable


def test_native_initial_invalid_argument_remains_internal_error() -> None:
    from relatipy import _core

    photon, _ = _tangential_failure_case()
    with pytest.raises(ValueError, match="internal error.*native status 2"):
        _core.integrate_kerr_null(photon._metric.spin, photon._initial_y, 0.0)


@pytest.mark.parametrize("sampling", ["steps", "samples"])
def test_native_empty_zeroed_failure_falls_back_to_independent_initial_state(sampling) -> None:
    from relatipy import _core

    photon, _ = _tangential_failure_case()
    initial = photon._initial_y.copy()
    evaluation = None if sampling == "steps" else np.array([0.0005, 0.001])
    # Bypass public tolerance validation to exercise a native failure that
    # returns no allocated rows and an untouched, all-zero final_state.
    states, final, stats, status, reason = _core.integrate_kerr_null(
        photon._metric.spin, initial, 0.001, evaluation, rtol=np.nan,
    )
    assert status == -1
    assert stats["native_status"] == 3
    assert reason is None
    np.testing.assert_array_equal(final, initial)
    assert not np.shares_memory(final, initial)
    assert not final.flags.writeable
    assert not states.flags.writeable
    if evaluation is None:
        np.testing.assert_array_equal(states, initial.reshape(1, 8))
        assert not np.shares_memory(states, initial)
        assert not np.shares_memory(states, final)
    else:
        assert states.shape == (0, 8)
