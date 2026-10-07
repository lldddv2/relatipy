"""Publication format for custom figures, observables reports and corner plots.

Assertions check sizes, scales, validation and immutability rather than
colors or pixels.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest
from astropy import units as u

from relatipy.plotting import (
    DEFAULT_STYLE, MM, TARGETS, Style, figure_width, format_axes, new_figure,
    plot_corner, plot_observables, publication_style, rc_params, save_figure,
)


@pytest.fixture
def pyplot(monkeypatch):
    monkeypatch.setenv("MPLBACKEND", "Agg")
    matplotlib = pytest.importorskip("matplotlib")
    with matplotlib.rc_context():
        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt

        yield plt
        plt.close("all")


def test_style_is_immutable_and_replace_builds_variants():
    style = DEFAULT_STYLE.replace(target="aanda")
    assert style.target is TARGETS["aanda"] and DEFAULT_STYLE.target is TARGETS["thesis"]
    with pytest.raises(dataclasses.FrozenInstanceError):
        style.font = "serif"
    with pytest.raises(TypeError):
        style.lines["trajectory"] = {"color": "red"}
    with pytest.raises(ValueError):
        Style(target="nature")
    with pytest.raises(ValueError):
        Style(length_unit="km")
    custom = Style(lines={**DEFAULT_STYLE.lines, "trajectory": {"color": "#0072B2",
                                                                "linewidth": 1.0}})
    assert custom.lines["trajectory"]["color"] == "#0072B2"
    assert DEFAULT_STYLE.lines["trajectory"]["color"] == "black"


def test_new_figure_has_printed_size_and_formatted_axes(pyplot):
    fig, ax = new_figure()
    width, height = fig.get_size_inches()
    assert width == pytest.approx(TARGETS["thesis"].width_mm * MM)
    assert height == pytest.approx(0.75 * width)
    assert ax.xaxis.get_tick_params(which="major")["direction"] == "in"
    ax.plot([-2.0, 2.0], [-1.0, 1.0])
    fig.canvas.draw()
    labels = [label.get_text() for label in ax.get_xticklabels() if label.get_text()]
    assert any("\N{MINUS SIGN}" in text for text in labels)
    assert not any("-" in text for text in labels)
    fig, axes = new_figure("full", 50, nrows=2, sharex=True)
    assert fig.get_size_inches()[0] == pytest.approx(TARGETS["thesis"].full_width_mm * MM)
    assert fig.get_size_inches()[1] == pytest.approx(50 * MM)
    assert len(axes) == 2
    assert figure_width(88) == pytest.approx(88 * MM)
    with pytest.raises(ValueError):
        figure_width("half")


def test_publication_style_applies_to_plain_matplotlib(pyplot):
    import matplotlib

    before = matplotlib.rcParams["xtick.direction"]
    with publication_style() as style:
        assert matplotlib.rcParams["xtick.direction"] == style.frame.tick_direction
        assert matplotlib.rcParams["pdf.fonttype"] == 42
        fig, ax = pyplot.subplots()
        format_axes(ax)
    assert matplotlib.rcParams["xtick.direction"] == before
    assert rc_params()["axes.labelsize"] == TARGETS["thesis"].label_size


def test_publication_style_keeps_interactive_mode_set_inside(pyplot):
    # Inline backends enable interactive mode when the first figure selects
    # them; leaving the style context must not switch it back off.
    import matplotlib

    before = matplotlib.is_interactive()
    try:
        with publication_style():
            matplotlib.interactive(not before)
        assert matplotlib.is_interactive() is (not before)
    finally:
        matplotlib.interactive(before)


def test_save_figure_writes_each_format(pyplot, tmp_path):
    fig, ax = new_figure()
    ax.plot([0, 1], [0, 1])
    written = save_figure(fig, tmp_path / "figure")
    assert [path.suffix for path in written] == [".pdf", ".png"]
    assert all(path.stat().st_size > 0 for path in written)
    assert save_figure(fig, tmp_path / "only.png") == [tmp_path / "only.png"]
    assert b"FontFile2" in (tmp_path / "figure.pdf").read_bytes()  # TrueType


def _report_inputs(count=40):
    epochs = np.linspace(2000.0, 2016.0, count)
    model_epochs = np.linspace(1999.0, 2017.0, 400)
    phase = 2 * np.pi * (model_epochs - 2000.0) / 16.0
    model = {"model_ra": 40 * np.cos(phase), "model_dec": 90 + 90 * np.sin(phase),
             "model_velocity": 1000 * np.sin(phase)}
    sample = 2 * np.pi * (epochs - 2000.0) / 16.0
    return {
        "astrometry_epochs": epochs, "ra": 40 * np.cos(sample),
        "ra_err": np.full(count, 1.0), "dec": 90 + 90 * np.sin(sample),
        "dec_err": np.full(count, 1.0), "spectroscopy_epochs": epochs[::2],
        "velocity": 1000 * np.sin(sample[::2]), "velocity_err": np.full(count // 2, 30.0),
        "model_epochs": model_epochs, **model,
    }


def test_observables_report_layout_and_equal_sky_scale(pyplot):
    inputs = _report_inputs()
    fig, (sky, ra, dec, velocity) = plot_observables(
        **inputs, residual_ra=np.zeros(40), residual_dec=np.zeros(40),
        names={"data": "S2 data", "center": "Sgr A*"})
    fig.canvas.draw()
    assert fig.get_size_inches()[0] == pytest.approx(TARGETS["thesis"].full_width_mm * MM)
    box = sky.get_window_extent()
    (x0, x1), (y0, y1) = sky.get_xlim(), sky.get_ylim()
    assert x0 > x1  # RA grows to the east, to the left
    assert box.width / abs(x1 - x0) == pytest.approx(box.height / (y1 - y0), rel=1e-6)
    # Stacked time panels touch and share the epoch axis.
    assert ra.get_position().y0 == pytest.approx(dec.get_position().y1)
    assert dec.get_position().y0 == pytest.approx(velocity.get_position().y1)
    names = [text.get_text() for text in sky.get_legend().texts]
    assert names == ["S2 data", "Model", "Data − model", "Sgr A*"]
    # The legend clears every drawn point.
    legend_bottom = sky.transData.inverted().transform(
        (0, sky.get_legend().get_window_extent().y0))[1]
    assert legend_bottom > (inputs["dec"] + inputs["dec_err"]).max()


def test_observables_accept_quantities_and_reject_mismatches(pyplot):
    inputs = _report_inputs()
    quantities = {**inputs, "ra": (inputs["ra"] * u.mas).to(u.arcsec),
                  "model_velocity": inputs["model_velocity"] * u.km / u.s}
    fig, axes = plot_observables(**quantities)
    np.testing.assert_allclose(axes[1].lines[0].get_ydata(), inputs["model_ra"])
    with pytest.raises(ValueError):
        plot_observables(**{**inputs, "ra": inputs["ra"][:-1]})
    with pytest.raises(ValueError):
        plot_observables(**inputs, residual_ra=np.zeros(40))
    with pytest.raises(ValueError):
        plot_observables(**{**inputs, "dec": np.full(40, np.nan)})


def test_corner_plot_structure_offsets_and_validation(pyplot):
    rng = np.random.default_rng(1)
    samples = rng.normal([2002.32, 0.88, -3.0], [0.005, 0.0012, 4.0], size=(4000, 3))
    fig, axes = plot_corner(samples, [r"$t_\mathrm{p}$ [yr]", "$e$",
                                      r"$v_\mathrm{LSR}$ [km/s]"])
    assert fig.get_size_inches() == pytest.approx(DEFAULT_STYLE.corner.size_in)
    assert axes[0, 1] is None and axes[2, 0] is not None
    # A large constant is subtracted where it shortens the tick labels.
    assert axes[2, 0].get_xlabel().endswith("− 2000")
    assert axes[2, 1].get_xlabel() == "$e$"
    fig.canvas.draw()
    sizes = {axes[2, i].xaxis.label.get_fontsize() for i in range(3)}
    assert len(sizes) == 1 and sizes.pop() <= DEFAULT_STYLE.corner.label_size
    with pytest.raises(ValueError):
        plot_corner(samples, ["a", "b"])
    with pytest.raises(ValueError):
        plot_corner(np.full((10, 2), np.inf), ["a", "b"])
