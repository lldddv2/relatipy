"""Several stored trajectories in one frame: same data, one black hole.

Agg is selected only by this test fixture. Assertions inspect plotted data,
legend names and colours that differ, not a fixed palette or pixel output.
"""

from __future__ import annotations

import numpy as np
import pytest
from astropy import units as u

from relatipy.geodesic.solution import Solution
from relatipy.metrics.kerr import Kerr
from relatipy.plotting import Style


@pytest.fixture
def pyplot(monkeypatch):
    """Provide a headless Matplotlib context and close test figures."""
    monkeypatch.setenv("MPLBACKEND", "Agg")
    matplotlib = pytest.importorskip("matplotlib")
    with matplotlib.rc_context():
        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt

        yield plt
        plt.close("all")


@pytest.fixture
def metric():
    return Kerr(mass=1 * u.Msun, spin=0.5)


def _solve(metric, a):
    orbit = metric.orbit(a=(a * metric.r_g).to(u.km), e=0.2, inc=0.7 * u.rad,
                         Omega=0.4 * u.rad, omega=0.3 * u.rad, f=0.2 * u.rad)
    return orbit.solve(tau_eval=np.linspace(0, 2, 5) * metric._time_scale,
                       method="dop853", rtol=1e-10, atol=1e-13)


@pytest.fixture
def solutions(metric):
    return _solve(metric, 30), _solve(metric, 45)


def _trajectory_lines(ax):
    return [line for line in ax.lines if line.get_gid() == "trajectory"]


def test_top_level_alias_is_the_plotting_function():
    import relatipy
    from relatipy.plotting import plot_sols

    assert relatipy.plot_sols is plot_sols and "plot_sols" in relatipy.__all__


@pytest.mark.parametrize("projection", ("xy", "xz", "3d"))
def test_each_solution_is_drawn_once_with_its_own_colour(solutions, projection, pyplot):
    from relatipy import plot_sols

    fig, ax = plot_sols(*solutions, projection=projection, interactive=False)
    lines = _trajectory_lines(ax)
    assert len(lines) == len(solutions)
    columns = (0, 1, 2) if projection == "3d" else tuple("xyz".index(a) for a in projection)
    for line, solution in zip(lines, solutions):
        data = (np.column_stack(line.get_data_3d()) if projection == "3d"
                else np.column_stack(line.get_data()))
        np.testing.assert_allclose(data, solution.xyz.value[:, columns], rtol=0, atol=1e-12)
    assert lines[0].get_color() != lines[1].get_color()
    centers = [c for c in ax.collections if c.get_gid() == "center"]
    assert len(centers) == 1


def test_legend_lists_default_and_custom_labels(solutions, pyplot):
    from relatipy import plot_sols

    fig, ax = plot_sols(*solutions, projection="xy")
    names = [text.get_text() for text in ax.get_legend().texts]
    assert "Trajectory 1" in names and "Trajectory 2" in names
    fig, ax = plot_sols(*solutions, projection="xy", labels=("a = 30", "a = 45"))
    names = [text.get_text() for text in ax.get_legend().texts]
    assert names.count("a = 30") == 1 and names.count("a = 45") == 1


def test_views_share_extent_of_all_solutions(solutions, pyplot):
    from relatipy import plot_sols

    fig, (main, top, right) = plot_sols(*solutions)
    for ax in (main, top, right):
        assert len(_trajectory_lines(ax)) == 2
    names = [text.get_text() for text in fig.axes[3].get_legend().texts]
    assert "Trajectory 1" in names and "Trajectory 2" in names


def test_lengths_use_the_unit_of_the_first_solution(metric, pyplot):
    from relatipy import plot_sols

    first, second = _solve(metric, 30), _solve(metric, 45)
    assert first.x.unit == u.km
    fig, ax = plot_sols(first, second, projection="xy")
    line = _trajectory_lines(ax)[1]
    np.testing.assert_allclose(np.column_stack(line.get_data()),
                               second.xyz.to_value(u.km)[:, :2], rtol=1e-12)


def test_one_unlabelled_solution_matches_plot_solution(solutions, pyplot):
    from relatipy import plot_sols
    from relatipy.plotting import plot_solution

    _, ax = plot_sols(solutions[0], projection="xy")
    _, reference = plot_solution(solutions[0], projection="xy")
    assert ([line.get_color() for line in _trajectory_lines(ax)]
            == [line.get_color() for line in _trajectory_lines(reference)])
    legend = ax.get_legend()
    expected = reference.get_legend()
    assert (legend is None) == (expected is None)


def test_early_end_is_reported_with_the_trajectory_name(solutions, pyplot):
    from relatipy import plot_sols

    failed = Solution(state=solutions[1]._state, integration=solutions[1].integration,
                      status=-1, message="numerical failure after stored samples")
    fig, ax = plot_sols(solutions[0], failed, projection="xy")
    text = " ".join(t.get_text() for t in fig.findobj(lambda a: hasattr(a, "get_text")))
    assert "Trajectory 2: partial" in text and "Trajectory 1: partial" not in text


def test_interactive_adds_one_trace_per_solution(solutions):
    pytest.importorskip("plotly")
    from relatipy import plot_sols

    figure = plot_sols(*solutions, projection="3d", labels=["inner", "outer"])
    names = [trace.name for trace in figure.data if trace.type == "scatter3d"
             and trace.mode == "lines"]
    assert "inner" in names and "outer" in names
    colors = {trace.name: trace.line.color for trace in figure.data
              if trace.name in ("inner", "outer")}
    assert colors["inner"] != colors["outer"]


@pytest.mark.parametrize("arguments,error", (
    ((), ValueError),
    (("not a solution",), TypeError),
))
def test_invalid_solutions_are_rejected(arguments, error, pyplot):
    from relatipy import plot_sols

    existing = pyplot.get_fignums()
    with pytest.raises(error):
        plot_sols(*arguments, projection="xy")
    assert pyplot.get_fignums() == existing


def test_labels_and_black_holes_must_match(solutions, pyplot):
    from relatipy import plot_sols

    with pytest.raises(ValueError, match="labels"):
        plot_sols(*solutions, labels=["only one"])
    other = _solve(Kerr(mass=1 * u.Msun, spin=0.9), 30)
    with pytest.raises(ValueError, match="same black hole"):
        plot_sols(solutions[0], other, projection="xy")


def test_trajectory_colours_come_from_the_style(solutions, pyplot):
    from relatipy import plot_sols

    style = Style(trajectory_colors=("red", "blue"))
    _, ax = plot_sols(*solutions, *solutions[:1], projection="xy", style=style)
    colors = [line.get_color() for line in _trajectory_lines(ax)]
    assert colors == ["red", "blue", "red"]
    with pytest.raises(ValueError):
        Style(trajectory_colors=())
    with pytest.raises(ValueError):
        Style(trajectory_colors="red")
