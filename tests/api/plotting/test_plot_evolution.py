"""Temporal plots use stored public coordinates and explicit component choices."""

from __future__ import annotations

import numpy as np
import pytest
from astropy import units as u

from relatipy.geodesic.orbit import Orbit
from relatipy.geodesic.solution import Solution
from relatipy.metrics.kerr import Kerr


@pytest.fixture
def pyplot(monkeypatch):
    """Use Matplotlib without a display and close figures after each test."""
    monkeypatch.setenv("MPLBACKEND", "Agg")
    matplotlib = pytest.importorskip("matplotlib")
    with matplotlib.rc_context():
        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt

        yield plt
        plt.close("all")


@pytest.fixture
def solution():
    metric = Kerr(mass=1 * u.Msun, spin=0.5)
    orbit = metric.orbit(
        a=30 * metric.r_g, e=0.2, inc=0.7 * u.rad,
        Omega=0.4 * u.rad, omega=0.3 * u.rad, f=0.2 * u.rad,
    )
    return orbit.solve(
        tau_eval=np.linspace(0, 2, 5) * metric._time_scale,
        method="dop853", rtol=1e-10, atol=1e-13,
    )


def _line_data(ax):
    """Return the plotted sample series, excluding decorative artists."""
    return [(np.asarray(line.get_xdata()), np.asarray(line.get_ydata()))
            for line in ax.lines if len(line.get_xdata())]


def _assert_panel(ax, time, values):
    expected_x = time.to_value(time.unit)
    expected_y = values.to_value(u.deg if values.unit.is_equivalent(u.rad) else values.unit)
    assert any(np.allclose(x, expected_x, rtol=0, atol=1e-12)
               and np.allclose(y, expected_y, rtol=0, atol=1e-12)
               for x, y in _line_data(ax))


def test_default_plots_cartesian_samples_against_coordinate_time(solution, pyplot):
    fig, axes = solution.plot_evol()

    assert isinstance(axes, np.ndarray) and axes.shape == (3,)
    assert r"$t$" in axes[-1].get_xlabel()
    assert all(ax.figure is fig for ax in axes)
    assert fig.get_size_inches()[1] < fig.get_size_inches()[0]
    for upper, lower in zip(axes[:-1], axes[1:]):
        assert upper.get_position().y0 == pytest.approx(lower.get_position().y1)
        assert lower.spines["top"].get_visible()
    assert not np.allclose(solution.t.to_value(u.s), solution.tau.to_value(u.s))
    for ax, name in zip(axes, ("x", "y", "z")):
        _assert_panel(ax, solution.t, getattr(solution, name))


@pytest.mark.parametrize("coords,names", (
    ("cartesial", ("x", "y", "z")),
    ("spherical", ("r", "theta", "phi")),
    ("bl", ("R", "Theta", "Phi")),
))
def test_tau_and_each_coordinate_family_use_public_samples(solution, pyplot, coords, names):
    fig, axes = solution.plot_evol(time="tau", coords=coords)

    assert isinstance(axes, np.ndarray) and axes.shape == (3,)
    assert r"$\tau$" in axes[-1].get_xlabel()
    for ax, name in zip(axes, names):
        assert ax.figure is fig
        _assert_panel(ax, solution.tau, getattr(solution, name))


@pytest.mark.parametrize("coords,labels", (
    ("spherical", (r"$r$", r"$\theta$", r"$\phi$")),
    ("bl", (r"$R$", r"$\Theta$", r"$\Phi$")),
))
def test_angular_axis_names_use_latex(solution, pyplot, coords, labels):
    fig, axes = solution.plot_evol(coords=coords)

    for ax, label in zip(axes, labels):
        assert label in ax.get_ylabel()
    for ax in axes[1:]:
        assert ax.get_ylabel().endswith("[deg]")
    fig.canvas.draw()


def test_elements_plot_stored_a_e_f_against_proper_time(solution, pyplot):
    elements = solution.orbital_elements()
    fig, axes = solution.plot_evol(time="tau", coords="elements")

    assert isinstance(axes, np.ndarray) and axes.shape == (3,)
    expected = (
        elements.a.to_value(elements.a.unit),
        np.asarray(elements.e),
        elements.f.to_value(u.deg),
    )
    for ax, name, values in zip(axes, ("a", "e", "f"), expected):
        assert ax.figure is fig
        assert name in ax.get_ylabel()
        assert len(_line_data(ax)) == 1
        x, y = _line_data(ax)[0]
        np.testing.assert_array_equal(x, solution.tau.to_value(solution.tau.unit))
        np.testing.assert_array_equal(y, values)
    assert elements.a.unit.to_string() in axes[0].get_ylabel()
    assert "deg" in axes[2].get_ylabel()
    fig.canvas.draw()
    assert not axes[1].yaxis.get_offset_text().get_text()
    assert not axes[2].yaxis.get_offset_text().get_text()


def test_elements_explicit_false_excludes_a_default(solution, pyplot):
    elements = solution.orbital_elements()
    _, axes = solution.plot_evol(coords="elements", a=False)

    assert axes.shape == (2,)
    assert "e" in axes[0].get_ylabel()
    assert "f" in axes[1].get_ylabel()
    np.testing.assert_array_equal(_line_data(axes[0])[0][1], elements.e)
    np.testing.assert_array_equal(_line_data(axes[1])[0][1], elements.f.to_value(u.deg))


def test_explicit_false_excludes_default_components(solution, pyplot):
    fig, axes = solution.plot_evol(x=False, y=False)

    assert isinstance(axes, np.ndarray) and axes.shape == (1,)
    assert axes[0].figure is fig
    _assert_panel(axes[0], solution.t, solution.z)


def test_true_flag_adds_component_from_another_family(solution, pyplot):
    fig, axes = solution.plot_evol(coords="cartesian", r=True)

    assert isinstance(axes, np.ndarray) and axes.shape == (4,)
    assert all(ax.figure is fig for ax in axes)
    for name in ("x", "y", "z", "r"):
        assert any(
            any(np.allclose(y, getattr(solution, name).to_value(getattr(solution, name).unit),
                            rtol=0, atol=1e-12) for _, y in _line_data(ax))
            for ax in axes
        )


@pytest.mark.parametrize("options,error", (
    ({"time": "proper"}, ValueError),
    ({"coords": "cylindrical"}, ValueError),
    ({"unknown": True}, ValueError),
    ({"x": 1}, TypeError),
    ({"x": None}, TypeError),
    ({"x": False, "y": False, "z": False}, ValueError),
    ({"coords": "elements", "a": False, "e": False, "f": False}, ValueError),
    ({"coords": "elements", "a": None}, TypeError),
))
def test_invalid_options_fail_before_creating_a_figure(solution, pyplot, options, error):
    before = pyplot.get_fignums()
    with pytest.raises(error):
        solution.plot_evol(**options)
    assert pyplot.get_fignums() == before


def test_single_sample_is_visible(solution, pyplot):
    single = Solution(
        state=solution[:1], integration=solution.integration,
        status=0, message="one stored sample",
    )
    _, axes = single.plot_evol()

    for ax in axes:
        lines = [line for line in ax.lines if len(line.get_xdata()) == 1]
        visible_line = any(line.get_marker() not in (None, "", " ", "None", "none")
                           for line in lines)
        visible_scatter = any(len(collection.get_offsets()) == 1
                              for collection in ax.collections)
        assert visible_line or visible_scatter


def test_failed_solution_marks_partial_output(solution, pyplot):
    failed = Solution(
        state=solution._state, integration=solution.integration,
        status=-1, message="numerical failure after stored samples",
    )
    fig, _ = failed.plot_evol()

    text = " ".join(artist.get_text() for artist in fig.findobj(
        lambda artist: hasattr(artist, "get_text"))).lower()
    assert "partial" in text


def test_plot_does_not_integrate_or_mutate_solution(monkeypatch, solution, pyplot):
    before = {name: getattr(solution, name).copy()
              for name in ("tau", "t", "x", "y", "z", "r", "theta", "phi")}

    def unexpected_integration(*args, **kwargs):
        raise AssertionError("plotting invoked geodesic integration")

    monkeypatch.setattr(Orbit, "integrate", unexpected_integration)
    monkeypatch.setattr(Orbit, "solve", unexpected_integration)
    solution.plot_evol(coords="spherical", x=True)

    for name, values in before.items():
        np.testing.assert_array_equal(getattr(solution, name), values)


def test_large_position_uses_visible_power_of_ten_without_changing_data(pyplot):
    metric = Kerr(mass=1e5 * u.Msun, spin=0.5)
    orbit = metric.orbit(
        a=30 * metric.r_g, e=0.2, inc=0.7 * u.rad,
        Omega=0.4 * u.rad, omega=0.3 * u.rad, f=0.2 * u.rad,
    )
    solution = orbit.solve(
        tau_eval=np.linspace(0, 2, 5) * metric._time_scale,
        method="dop853", rtol=1e-10, atol=1e-13,
    )
    fig, axes = solution.plot_evol()
    fig.canvas.draw()

    for ax, name in zip(axes, ("x", "y", "z")):
        values = getattr(solution, name)
        np.testing.assert_array_equal(
            _line_data(ax)[0][1], values.to_value(values.unit),
        )
        factor = ax.get_ylabel()
        assert "10" in factor and r"\times" in factor
        assert not ax.yaxis.get_offset_text().get_text()
        labels = [label.get_text() for label in ax.get_yticklabels()]
        assert all("e+" not in label and "e-" not in label for label in labels)
