"""Constants-of-motion plots use the stored E, Lz and Carter Q series."""

from __future__ import annotations

import numpy as np
import pytest
from astropy import units as u
from astropy.constants import c

from relatipy.metrics.kerr import Kerr
from relatipy.plotting import plot_constants


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


def _solution(vz):
    metric = Kerr(mass=1 * u.Msun, spin=0.5)
    orbit = metric.orbit(x=20 * metric.r_g, vy=0.2 * c, vz=vz)
    return orbit.solve(tau_span=(0 * u.s, 300 * metric._time_scale))


def _ydata(ax):
    return np.asarray(ax.lines[0].get_ydata())


def test_default_plots_three_constants_against_coordinate_time(pyplot):
    solution = _solution(0.05 * c)
    fig, axes = solution.plot_constants()

    assert axes.shape == (3,)
    assert all(ax.figure is fig for ax in axes)
    assert r"$t$" in axes[-1].get_xlabel()
    np.testing.assert_array_equal(axes[0].lines[0].get_xdata(), solution.t.value)
    for ax, values, label in zip(
        axes, (solution.get_E(), solution.get_Lz(), solution.get_Q()),
        ("$E - E_0$", "$L_z - L_{z,0}$", "$Q - Q_0$"),
    ):
        # Conserved values are unresolvable by ticks: shift by the first sample.
        assert ax.get_ylabel().startswith(label)
        np.testing.assert_array_equal(_ydata(ax), values - values[0])
        assert f"{values[0]:.12g}" in ax.texts[0].get_text()


def test_resolved_series_is_plotted_unshifted(pyplot, monkeypatch):
    solution = _solution(0.05 * c)
    varying = np.linspace(1.0, 2.0, len(solution))
    monkeypatch.setattr(type(solution), "get_Q", lambda self: varying)
    _, axes = solution.plot_constants(E=False, Lz=False)
    assert axes[0].get_ylabel().startswith("$Q$")
    np.testing.assert_array_equal(_ydata(axes[0]), varying)


def test_proper_time_axis_and_switches(pyplot):
    solution = _solution(0.05 * c)
    _, axes = plot_constants(solution, time="tau", E=False)

    assert axes.shape == (2,)
    assert r"$\tau$" in axes[-1].get_xlabel()
    np.testing.assert_array_equal(axes[0].lines[0].get_xdata(), solution.tau.value)
    assert axes[0].get_ylabel().startswith("$L_z")


def test_drift_is_relative_or_absolute_for_zero_reference(pyplot):
    inclined = _solution(0.05 * c)
    _, axes = inclined.plot_constants(drift=True)
    energy = inclined.get_E()
    np.testing.assert_allclose(_ydata(axes[0]), (energy - energy[0]) / abs(energy[0]),
                               rtol=0, atol=0)
    assert _ydata(axes[0])[0] == 0.0
    assert "|E_0|" in axes[0].get_ylabel()

    equatorial = _solution(0 * c)
    _, axes = equatorial.plot_constants(drift=True, E=False, Lz=False)
    carter = equatorial.get_Q()
    assert abs(carter[0]) < 1e-12
    np.testing.assert_array_equal(_ydata(axes[0]), carter - carter[0])
    assert axes[0].get_ylabel().startswith(r"$\Delta Q$")


@pytest.mark.parametrize(
    ("kwargs", "error"),
    [
        ({"time": "lambda"}, ValueError),
        ({"energy": True}, ValueError),
        ({"E": False, "Lz": False, "Q": False}, ValueError),
        ({"Q": 1}, TypeError),
        ({"drift": "yes"}, TypeError),
    ],
)
def test_invalid_arguments(pyplot, kwargs, error):
    with pytest.raises(error):
        _solution(0.05 * c).plot_constants(**kwargs)
