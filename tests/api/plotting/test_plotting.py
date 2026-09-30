"""Plotting contracts: physical coordinates, honest scales, and no evolution.

Agg is selected only by this test fixture. Assertions inspect displayed data
and scale, rather than freezing colors, artist order, or rendering pixels.
Matplotlib keeps data in the solution unit and shows lengths through tick
labels; Plotly shows them in the data (hover labels), so its data are scaled.
"""

from __future__ import annotations

import subprocess
import sys

import numpy as np
import pytest
from astropy import units as u

from relatipy.geodesic.orbit import Orbit
from relatipy.geodesic.solution import Solution
from relatipy.metrics.kerr import Kerr
from relatipy.plotting import Style

GRAVITATIONAL_STYLE = Style(length_unit="r_g")


@pytest.fixture
def pyplot(monkeypatch):
    """Provide an optional headless Matplotlib context and close test figures."""
    monkeypatch.setenv("MPLBACKEND", "Agg")
    matplotlib = pytest.importorskip("matplotlib")
    with matplotlib.rc_context():
        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt

        yield plt
        plt.close("all")


@pytest.fixture
def orbit():
    metric = Kerr(mass=1 * u.Msun, spin=0.5)
    return metric.orbit(
        a=(30 * metric.r_g).to(u.km), e=0.2,
        inc=0.7 * u.rad, Omega=0.4 * u.rad, omega=0.3 * u.rad, f=0.2 * u.rad,
    )


@pytest.fixture
def solution(orbit):
    return orbit.solve(
        tau_eval=np.linspace(0, 2, 5) * orbit._metric._time_scale,
        method="dop853", rtol=1e-10, atol=1e-13,
    )


@pytest.fixture
def metre_orbit():
    metric = Kerr(mass=1 * u.Msun, spin=0.9)
    return metric.orbit(x=1e11 * u.m, vy=30 * u.km / u.s)


@pytest.fixture
def metre_solution(metre_orbit):
    return metre_orbit.solve(tau_eval=[0, 1e-9] * u.yr, method="dop853")


def _line_coordinates(line, projection):
    if projection == "3d":
        return np.column_stack(line.get_data_3d())
    return np.column_stack(line.get_data())


def _find_data_line(ax, expected, projection):
    for line in ax.lines:
        values = _line_coordinates(line, projection)
        if values.shape == expected.shape and np.allclose(values, expected, rtol=0, atol=1e-12):
            return line
    raise AssertionError("no plotted line matches the public physical samples")


def _axes_snapshot(fig, ax):
    """Snapshot caller-owned drawing state before a rejected operation."""
    return (
        tuple(fig.axes), tuple(ax.lines), tuple(ax.collections), tuple(ax.texts),
        ax.get_title(), ax.get_xlabel(), ax.get_ylabel(),
        ax.get_xlim(), ax.get_ylim(),
    )


@pytest.mark.parametrize("projection", ("xy", "xz", "yz", "3d"))
@pytest.mark.parametrize("style", (None, GRAVITATIONAL_STYLE), ids=("unit", "r_g"))
def test_plot_uses_stored_public_xyz_with_matching_labels(solution, projection, style,
                                                          pyplot):
    from relatipy.plotting import plot_solution

    expected = solution.xyz.to_value(solution.xyz.unit)
    if projection != "3d":
        expected = expected[:, ["xyz".index(axis) for axis in projection]]
    fig, ax = plot_solution(solution, projection=projection, interactive=False,
                            style=style)
    assert ax.figure is fig
    _find_data_line(ax, expected, projection)
    unit = solution.xyz.unit.to_string() if style is None else r"r_\mathrm{g}"
    axes = "xyz" if projection == "3d" else projection
    for index, name in enumerate(axes):
        label = getattr(ax, f"get_{'xyz'[index]}label")()
        assert f"${name}$" in label and unit in label
    if projection != "3d":
        assert ax.get_aspect() == 1.0
    fig.canvas.draw()


@pytest.mark.parametrize("projection", ("xy", "xz"))
@pytest.mark.parametrize("style", (None, GRAVITATIONAL_STYLE), ids=("unit", "r_g"))
def test_tick_labels_show_selected_unit(solution, projection, style, pyplot):
    from relatipy.plotting import plot_solution

    fig, ax = plot_solution(solution, projection=projection, style=style)
    fig.canvas.draw()
    scale = 1 if style is None else solution._metric.r_g.to_value(solution.xyz.unit)
    ticks = ax.get_xticks()
    shown = [label.get_text() for label in ax.get_xticklabels()]
    for value, text in zip(ticks, shown):
        if text:
            assert float(text.replace("\N{MINUS SIGN}", "-")) == pytest.approx(value / scale)
    assert all("-" not in text for text in shown)  # typographic minus only


@pytest.mark.parametrize("kind", ("solution", "preview"))
def test_default_views_show_input_metres(metre_solution, metre_orbit, kind, pyplot):
    source = metre_solution if kind == "solution" else metre_orbit
    fig, axes = source.plot() if kind == "solution" else source.preview()
    assert source.xyz.unit == u.m
    assert all("[m]" in label for ax in axes
               for label in (ax.get_xlabel(), ax.get_ylabel()) if label)
    fig.canvas.draw()


@pytest.mark.parametrize("interactive", (False, True), ids=("matplotlib", "plotly"))
@pytest.mark.parametrize("style", (None, GRAVITATIONAL_STYLE), ids=("metres", "r_g"))
def test_metre_input_uses_selected_display_unit(
    metre_solution, interactive, style, pyplot,
):
    if interactive:
        pytest.importorskip("plotly")
    expected = metre_solution.xyz.to_value(u.m)[:, :2].copy()
    unit = "m"
    if style is not None:
        expected /= metre_solution._metric.r_g.to_value(u.m)
        unit = "<i>r</i><sub>g</sub>"
    result = metre_solution.plot(projection="xy", interactive=interactive, style=style)
    if interactive:
        traces = [trace for trace in result.data if trace.name == "Trajectory"]
        assert len(traces) == 1
        np.testing.assert_allclose(np.column_stack((traces[0].x, traces[0].y)), expected)
        assert unit in result.layout.xaxis.title.text
        assert unit in result.layout.yaxis.title.text
    else:
        fig, ax = result
        _find_data_line(ax, metre_solution.xyz.to_value(u.m)[:, :2], "xy")
        assert ("[m]" if style is None else r"r_\mathrm{g}") in ax.get_xlabel()
        assert ("[m]" if style is None else r"r_\mathrm{g}") in ax.get_ylabel()
        fig.canvas.draw()


@pytest.mark.parametrize("kind", ("solution", "preview"))
def test_three_dimensional_plot_has_equal_physical_scale(solution, orbit, kind, pyplot):
    from relatipy.plotting import plot_solution, preview_orbit

    function, source = (plot_solution, solution) if kind == "solution" else (preview_orbit, orbit)
    fig, ax = function(source, projection="3d", interactive=False)
    fig.canvas.draw()
    spans = np.diff(np.asarray(ax.get_w_lims()).reshape(3, 2), axis=1).ravel()
    physical_scale = np.asarray(ax.get_box_aspect()) / spans
    assert np.all(spans > 0)
    np.testing.assert_allclose(physical_scale, [physical_scale[0]] * 3,
                               rtol=1e-10, atol=1e-12)


@pytest.mark.parametrize("projection", ("xz", "yz", "3d"))
def test_existing_compatible_axes_are_reused(solution, projection, pyplot):
    from relatipy.plotting import plot_solution

    fig = pyplot.figure()
    ax = fig.add_subplot(111, projection="3d" if projection == "3d" else None)
    original = ax.plot([0, 1], [0, 1], label="caller data")[0]
    returned_fig, returned_ax = plot_solution(solution, projection=projection,
                                              interactive=False, ax=ax)
    assert returned_fig is fig and returned_ax is ax
    assert fig.axes == [ax] and original in ax.lines


@pytest.mark.parametrize("kind", ("solution", "preview"))
@pytest.mark.parametrize("incompatibility", ("axes_dimension", "figure_owner"))
def test_incompatible_axes_rejection_does_not_mutate_figures(
    solution, orbit, kind, incompatibility, pyplot,
):
    from relatipy.plotting import plot_solution, preview_orbit

    fig, ax = pyplot.subplots()
    ax.plot([1, 2], [3, 4], label="caller data")
    ax.set_title("caller title")
    other = pyplot.figure()
    before = _axes_snapshot(fig, ax)
    other_axes = tuple(other.axes)
    projection = "3d" if incompatibility == "axes_dimension" else "xy"
    supplied_fig = fig if incompatibility == "axes_dimension" else other
    function, source = (plot_solution, solution) if kind == "solution" else (preview_orbit, orbit)
    with pytest.raises(ValueError):
        function(source, projection=projection, interactive=False,
                 fig=supplied_fig, ax=ax)
    assert _axes_snapshot(fig, ax) == before
    assert tuple(other.axes) == other_axes


@pytest.mark.parametrize("kind", ("solution", "preview"))
@pytest.mark.parametrize("keyword,value", (
    ("fig", object()), ("ax", object()), ("projection", "polar"), ("interactive", "yes"),
    ("style", object()),
))
@pytest.mark.parametrize("interactive", (False, None))
def test_invalid_plot_arguments_are_rejected(solution, orbit, kind, keyword, value,
                                             interactive, pyplot):
    from relatipy.plotting import plot_solution, preview_orbit

    function, source = (plot_solution, solution) if kind == "solution" else (preview_orbit, orbit)
    existing = pyplot.get_fignums()
    arguments = {"interactive": interactive, keyword: value}
    with pytest.raises((TypeError, ValueError)):
        function(source, **arguments)
    assert pyplot.get_fignums() == existing


@pytest.mark.parametrize("projection", ("xy", "3d"))
def test_one_sample_is_visible_as_a_marker(solution, projection, pyplot):
    from relatipy.plotting import plot_solution

    single = Solution(
        state=solution[:1], integration=solution.integration,
        status=0, message="one stored sample",
    )
    fig, ax = plot_solution(single, projection=projection, interactive=False)
    expected = single.xyz.to_value(single.xyz.unit)
    if projection == "xy":
        expected = expected[:, :2]
    matching_line = _find_data_line(ax, expected, projection)
    visible_line_marker = matching_line.get_marker() not in (None, "", " ", "None", "none")
    visible_scatter = False
    for collection in ax.collections:
        # Matplotlib exposes 2D offsets publicly; 3D marker coordinates are
        # retained in _offsets3d until projection by the canvas.
        values = np.column_stack(collection._offsets3d) if projection == "3d" else collection.get_offsets()
        if np.shape(values) == expected.shape and np.allclose(values, expected, rtol=0, atol=1e-12):
            visible_scatter = bool(len(collection.get_paths()) and np.any(collection.get_sizes() > 0))
    assert visible_line_marker or visible_scatter
    fig.canvas.draw()


def test_failed_solution_is_visibly_labeled_as_partial(solution, pyplot):
    from relatipy.plotting import plot_solution

    failed = Solution(
        state=solution._state, integration=solution.integration,
        status=-1, message="numerical failure after stored samples",
    )
    for projection in ("xz", "views", "3d"):
        fig, axes = plot_solution(failed, projection=projection, interactive=False)
        # Shown even though titles are off by default.
        text = " ".join(artist.get_text() for artist in fig.findobj(
            lambda artist: hasattr(artist, "get_text"))).lower()
        assert "partial" in text and "fail" in text
        fig.canvas.draw()


def test_preview_marker_uses_current_public_position_and_projection(orbit, pyplot):
    from relatipy.plotting import preview_orbit

    orbit.integrate(0.1 * orbit._metric._time_scale, method="dp45")
    current = orbit.xyz.to_value(orbit.xyz.unit)
    for projection, columns in (("xz", (0, 2)), ("yz", (1, 2))):
        fig, ax = preview_orbit(orbit, projection=projection, show_horizon=False,
                                show_isco=False)
        expected = current[list(columns)]
        assert any(np.allclose(collection.get_offsets(), [expected], rtol=0, atol=1e-12)
                   for collection in ax.collections)
        names = " ".join(text.get_text() for text in ax.get_legend().texts).lower()
        assert "preview" in names and "trajectory" not in names
        fig.canvas.draw()


def test_rendering_never_integrates_or_mutates_sources(monkeypatch, solution, orbit, pyplot):
    from relatipy import _core
    from relatipy.plotting import plot_solution, preview_orbit

    def unexpected_integration(*args, **kwargs):
        raise AssertionError("rendering invoked geodesic integration")

    solution_xyz = solution.xyz.copy()
    solution_tangent = solution.uxyz.copy()
    current_xyz = orbit.xyz.copy()
    current_tau = orbit.tau.copy()
    initial_xyz = orbit.initial.xyz.copy()
    monkeypatch.setattr(_core, "integrate_kerr", unexpected_integration)
    monkeypatch.setattr(Orbit, "solve", unexpected_integration)
    monkeypatch.setattr(Orbit, "integrate", unexpected_integration)
    for function, source in ((plot_solution, solution), (preview_orbit, orbit)):
        fig, ax = function(source, projection="3d", interactive=False)
        assert ax.figure is fig
        fig.canvas.draw()
    for source, method in ((solution, "plot"), (orbit, "preview")):
        fig, ax = getattr(source, method)(projection="xz")
        assert ax.figure is fig
        fig, axes = getattr(source, method)(projection="views")
        assert len(axes) == 3
    np.testing.assert_array_equal(solution.xyz, solution_xyz)
    np.testing.assert_array_equal(solution.uxyz, solution_tangent)
    np.testing.assert_array_equal(orbit.xyz, current_xyz)
    np.testing.assert_array_equal(orbit.initial.xyz, initial_xyz)
    assert orbit.tau == current_tau
    assert not solution.xyz.value.flags.writeable


def test_matplotlib_is_optional_and_imported_only_when_rendering():
    """Block Matplotlib in a fresh process so prior tests cannot mask eagerness."""
    script = r"""
import importlib.abc
import sys
class BlockMatplotlib(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == 'matplotlib' or fullname.startswith('matplotlib.'):
            raise ModuleNotFoundError('Matplotlib deliberately unavailable')
sys.meta_path.insert(0, BlockMatplotlib())
import relatipy
from relatipy.plotting import plot_solution, preview_orbit
assert not any(name == 'matplotlib' or name.startswith('matplotlib.') for name in sys.modules)
from astropy import units as u
metric = relatipy.Kerr(mass=1*u.Msun, spin=0)
orbit = metric.orbit(a=30*metric.r_g, e=.2, inc=.7*u.rad, Omega=.4*u.rad, omega=.3*u.rad)
solution = orbit.solve(tau_eval=[0, .001]*metric._time_scale, method='dop853')
for function, source in ((plot_solution, solution), (preview_orbit, orbit)):
    try:
        function(source, projection='xy')
    except ImportError as error:
        assert 'matplotlib' in str(error).lower()
    else:
        raise AssertionError('rendering succeeded with Matplotlib blocked')
"""
    result = subprocess.run([sys.executable, "-c", script], capture_output=True,
                            text=True, check=False, timeout=30)
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("projection", ("xy", "xz", "yz", "3d"))
@pytest.mark.parametrize("style", (None, GRAVITATIONAL_STYLE), ids=("unit", "r_g"))
def test_interactive_plot_preserves_public_coordinates_and_units(
    solution, projection, style,
):
    pytest.importorskip("plotly")
    from relatipy.plotting import plot_solution_interactive

    figure = plot_solution_interactive(solution, projection=projection, style=style)
    columns = (0, 1, 2) if projection == "3d" else tuple("xyz".index(axis) for axis in projection)
    # Hover labels show the trace data in the selected display unit.
    expected = (solution.xyz.to_value(solution.xyz.unit) if style is None else
                (solution.xyz / solution._metric.r_g).decompose().value)[:, columns]
    coordinate_names = "xyz" if projection == "3d" else "xy"
    candidates = [np.column_stack([getattr(trace, axis) for axis in coordinate_names])
                  for trace in figure.data if trace.type in ("scatter", "scatter3d")]
    assert any(values.shape == expected.shape and np.allclose(values, expected, rtol=1e-12)
               for values in candidates)
    unit = solution.xyz.unit.to_string() if style is None else "<i>r</i><sub>g</sub>"
    axes = "xyz" if projection == "3d" else projection
    layout = figure.layout.scene if projection == "3d" else figure.layout
    for index, coordinate in enumerate(axes):
        title = getattr(layout, f"{'xyz'[index]}axis").title.text
        assert coordinate in title and unit in title
    if projection == "3d":
        assert figure.layout.scene.aspectmode == "data"
    else:
        assert figure.layout.yaxis.scaleanchor == "x"
        assert figure.layout.yaxis.scaleratio == 1


def test_interactive_reuses_figure_marks_partial_and_never_integrates(monkeypatch, solution):
    graph_objects = pytest.importorskip("plotly.graph_objects")
    from relatipy import _core
    from relatipy.plotting import plot_solution_interactive

    failed = Solution(state=solution._state, integration=solution.integration,
                      status=-1, message="test numerical failure")
    before = failed.xyz.copy()
    figure = graph_objects.Figure()
    figure.add_scatter(x=[0, 1], y=[2, 3], name="caller data")
    caller_trace = figure.data[0]

    def unexpected_integration(*args, **kwargs):
        raise AssertionError("interactive rendering invoked integration")

    monkeypatch.setattr(_core, "integrate_kerr", unexpected_integration)
    monkeypatch.setattr(Orbit, "solve", unexpected_integration)
    monkeypatch.setattr(Orbit, "integrate", unexpected_integration)
    returned = plot_solution_interactive(failed, projection="xz", fig=figure)
    assert returned is figure and figure.data[0] is caller_trace
    notes = " ".join(note.text for note in figure.layout.annotations).lower()
    assert "partial" in notes and "fail" in notes
    np.testing.assert_array_equal(failed.xyz, before)
    assert not failed.xyz.value.flags.writeable


def test_plotly_is_optional_and_imported_only_when_rendering():
    script = r"""
import importlib.abc
import sys
class BlockPlotly(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == 'plotly' or fullname.startswith('plotly.'):
            raise ModuleNotFoundError('Plotly deliberately unavailable')
sys.meta_path.insert(0, BlockPlotly())
import relatipy
from relatipy.plotting import plot_solution_interactive
assert not any(name == 'plotly' or name.startswith('plotly.') for name in sys.modules)
from astropy import units as u
metric = relatipy.Kerr(mass=1*u.Msun, spin=0)
orbit = metric.orbit(a=30*metric.r_g, e=.2, inc=.7*u.rad, Omega=.4*u.rad, omega=.3*u.rad)
solution = orbit.solve(tau_eval=[0, .001]*metric._time_scale, method='dop853')
try:
    plot_solution_interactive(solution)
except ImportError as error:
    assert 'plotly' in str(error).lower()
else:
    raise AssertionError('interactive rendering succeeded with Plotly blocked')
"""
    result = subprocess.run([sys.executable, "-c", script], capture_output=True,
                            text=True, check=False, timeout=30)
    assert result.returncode == 0, result.stderr


def test_plotly_three_dimensional_view_has_native_surfaces(
    monkeypatch, solution, orbit,
):
    graph_objects = pytest.importorskip("plotly.graph_objects")
    calls = []
    original = Kerr._surface_profile

    def recorded(self, surface, theta):
        calls.append(surface)
        return original(self, surface, theta)

    monkeypatch.setattr(Kerr, "_surface_profile", recorded)
    for figure in (solution.plot(projection="3d", interactive=True),
                   orbit.preview(projection="3d", interactive=True)):
        assert isinstance(figure, graph_objects.Figure)
        surfaces = [trace for trace in figure.data if trace.type == "surface"]
        assert {trace.name for trace in surfaces} == {"Outer horizon", "Ergosphere"}
        assert figure.layout.scene.aspectmode == "data"
        # The horizon surface replaces the equatorial horizon circle.
        assert sum(trace.name == "Outer horizon" for trace in figure.data) == 1
    assert set(calls) == {"outer_horizon", "ergosurface"}


@pytest.mark.parametrize("style", (None, GRAVITATIONAL_STYLE), ids=("unit", "r_g"))
def test_horizon_surface_uses_native_cartesian_profile(solution, style):
    pytest.importorskip("plotly")
    figure = solution.plot(projection="3d", interactive=True, style=style)
    horizon = next(trace for trace in figure.data if trace.name == "Outer horizon")
    theta = np.linspace(0.0, np.pi, 72) * u.rad
    rho, z = solution._metric._surface_profile("outer_horizon", theta)
    unit = solution.xyz.unit if style is None else solution._metric.r_g
    np.testing.assert_allclose(np.hypot(horizon.x, horizon.y)[:, 0],
                               (rho / unit).decompose().value, rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(np.asarray(horizon.z)[:, 0], (z / unit).decompose().value,
                               rtol=1e-10, atol=1e-12)


def test_views_share_one_physical_scale(solution, orbit, pyplot):
    from relatipy.plotting import plot_views

    for source in (solution, orbit):
        fig, (main, top, right) = plot_views(source)
        fig.canvas.draw()
        scales = []
        for ax in (main, top, right):
            box = ax.get_window_extent()
            (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
            scales += [box.width / (x1 - x0), box.height / (y1 - y0)]
        np.testing.assert_allclose(scales, scales[0], rtol=1e-9)
        # Shared edges: the top panel has the main width, the right its height.
        assert main.get_position().width == pytest.approx(top.get_position().width)
        assert main.get_position().height == pytest.approx(right.get_position().height)
    with pytest.raises(ValueError):
        solution.plot(projection="views", interactive=True)


def test_default_projection_is_static_views(solution, orbit, pyplot):
    """Omitted projection selects three static views for both public APIs."""
    from matplotlib.figure import Figure
    from relatipy.plotting import plot_solution, preview_orbit

    for fig, axes in (solution.plot(), orbit.preview(),
                      plot_solution(solution), preview_orbit(orbit)):
        assert isinstance(fig, Figure) and len(axes) == 3
        assert all(ax.figure is fig and ax.name == "rectilinear" for ax in axes)
    with pytest.raises(ValueError, match="static figure only"):
        solution.plot(interactive=True)


def test_explicit_3d_is_interactive_and_planes_keep_their_types(solution, orbit, pyplot):
    graph_objects = pytest.importorskip("plotly.graph_objects")
    from matplotlib.figure import Figure

    for figure in (solution.plot(projection="3d"), orbit.preview(projection="3d")):
        assert isinstance(figure, graph_objects.Figure)
    for (fig, ax), name in ((solution.plot(projection="xy"), "rectilinear"),
                            (orbit.preview(projection="xz"), "rectilinear"),
                            (solution.plot(projection="3d", interactive=False), "3d")):
        assert isinstance(fig, Figure) and ax.name == name and ax.figure is fig
    assert isinstance(solution.plot(projection="yz", interactive=True), graph_objects.Figure)
    with pytest.raises(ValueError):
        solution.plot(projection="views", interactive=True)
