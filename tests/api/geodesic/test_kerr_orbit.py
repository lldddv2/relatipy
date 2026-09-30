"""Public Kerr and scalar Orbit integration contracts."""

import numpy as np
import pytest
from astropy import units as u
from astropy.constants import c

from relatipy.metrics.kerr import Kerr


def _orbit():
    bh = Kerr(mass=1 * u.Msun, spin=0.5)
    return bh.orbit(x=12 * bh.r_g, vy=0.1 * c)


def test_kerr_rejects_invalid_mass_and_spin() -> None:
    with pytest.raises(ValueError):
        Kerr(mass=0 * u.kg, spin=0.5)
    with pytest.raises(ValueError):
        Kerr(mass=1 * u.Msun, spin=-0.1)
    with pytest.raises(ValueError):
        Kerr(mass=1 * u.Msun, spin=1.1)


def test_kerr_radii_and_ergosurface_are_quantities() -> None:
    bh = Kerr(mass=1 * u.Msun, spin=0.5)
    assert bh.horizons.event > bh.horizons.cauchy
    assert bh.horizons.event.unit.is_equivalent(u.m)
    polar = bh.r_ergosurface(theta=[0, np.pi / 2] * u.rad)
    assert polar.shape == (2,)
    assert polar[0] < polar[1]
    with pytest.raises(ValueError):
        bh.r_ergosurface(theta=[np.nan] * u.rad)


def test_input_families_and_physical_units() -> None:
    bh = Kerr(mass=1 * u.Msun, spin=0.5)
    with pytest.raises(ValueError, match="families"):
        bh.orbit(x=12 * bh.r_g, R=12 * bh.r_g)
    with pytest.raises(ValueError, match="requires"):
        bh.orbit(r=12 * bh.r_g, phi=0 * u.rad)
    with pytest.raises(u.UnitConversionError):
        bh.orbit(x=10 * u.Msun)
    with pytest.raises(ValueError):
        bh.orbit(x=0 * bh.r_g)


def test_initial_reset_copy_and_absolute_integrate() -> None:
    orb = _orbit()
    first = orb.initial.xyz.copy()
    assert orb.initial.tau == 0 * u.s
    assert orb.tau == orb.initial.tau
    target = orb.tau + 0.01 * orb._metric._time_scale
    assert orb.integrate(target, method="dp45") is None
    assert orb.tau == target
    copied = orb.copy()
    copied.reset()
    assert copied.xyz.unit.is_equivalent(first.unit)
    np.testing.assert_allclose(copied.xyz.to_value(first.unit), first.value)
    assert orb.tau == target
    orb.reset()
    np.testing.assert_allclose(orb.xyz.to_value(first.unit), first.value)


def test_solve_is_independent_and_samples_requested_times() -> None:
    orb = _orbit()
    times = np.array([0.0, 0.005, 0.01]) * orb._metric._time_scale
    sol = orb.solve(tau_eval=times, method="dp45")
    assert sol.status == 0
    assert sol.integration.method == "dp45"
    assert len(sol) == len(times)
    np.testing.assert_allclose(sol.tau.to_value(u.s), times.to_value(u.s))
    assert orb.tau == orb.initial.tau
    with pytest.raises(ValueError, match="strictly increasing"):
        orb.solve(tau_eval=np.array([0, 0]) * u.s)


def test_all_initial_condition_families_share_the_time_contract() -> None:
    bh = Kerr(mass=1 * u.Msun, spin=0)
    radius = 30 * bh.r_g
    angular_velocity = (0.1 * c / radius) * u.rad
    orbits = (
        bh.orbit(x=radius, vy=0.1 * c, tau=2 * u.s, t=3 * u.s),
        bh.orbit(r=radius, theta=90 * u.deg, phi=0 * u.rad,
                 vphi=angular_velocity, tau=2 * u.s, t=3 * u.s),
        bh.orbit(R=radius, Theta=90 * u.deg, Phi=0 * u.rad,
                 vPhi=angular_velocity, tau=2 * u.s, t=3 * u.s),
        bh.orbit(a=100 * bh.r_g, e=0.2, inc=30 * u.deg,
                 tau=2 * u.s, t=3 * u.s),
    )
    for orb in orbits:
        assert orb.tau == 2 * u.s
        assert orb.t == 3 * u.s
        assert orb.initial.tau == orb.tau
        assert orb.initial.t == orb.t
    np.testing.assert_allclose(
        orbits[0].xyz.to_value(u.m), orbits[1].xyz.to_value(u.m), rtol=1e-12
    )
    np.testing.assert_allclose(
        orbits[0].xyz.to_value(u.m), orbits[2].xyz.to_value(u.m), rtol=1e-12
    )


def test_requested_endpoint_is_retained_at_small_physical_times() -> None:
    orb = _orbit()
    times = np.array([0.0, 1e-8, 2e-8]) * u.s
    sol = orb.solve(tau_eval=times, method="dp45")
    assert sol.status == 0
    assert len(sol) == 3
    np.testing.assert_array_equal(sol.tau.to_value(u.s), times.to_value(u.s))


@pytest.mark.parametrize("method", ["radau", "dop853", "dp45", "projection_radau"])
def test_public_method_dispatch_reaches_requested_time(method: str) -> None:
    orb = _orbit()
    end = 0.01 * orb._metric._time_scale
    sol = orb.solve(tau_span=(0 * u.s, end), method=method)
    assert sol.status == 0
    assert sol.integration.method == method
    assert sol.tau[-1] == end


def test_default_method_is_recorded_in_lowercase() -> None:
    orb = _orbit()
    end = 0.01 * orb._metric._time_scale
    sol = orb.solve(tau_span=(0 * u.s, end))
    assert sol.integration.method == "radau"


@pytest.mark.parametrize("method", ["Radau", "DOP853", "DP45", "Projection_Radau"])
def test_legacy_method_names_are_rejected(method: str) -> None:
    orb = _orbit()
    end = 0.01 * orb._metric._time_scale
    with pytest.raises(ValueError, match="method must be one of"):
        orb.solve(tau_span=(0 * u.s, end), method=method)
    with pytest.raises(ValueError, match="method must be one of"):
        orb.integrate(end, method=method)


def test_osculating_preview_is_labeled_and_returns_axes() -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    orb = _orbit()
    fig, ax = orb.preview(projection="xy", show_horizon=True, show_isco=True)
    try:
        assert ax.figure is fig
        names = " ".join(text.get_text() for text in ax.get_legend().texts)
        assert "Kepler" in names and "horizon" in names.lower()
        # Face-on, the ISCO names are written along their circles.
        roles = {line.get_gid() for line in ax.lines}
        assert {"preview", "horizon", "isco_prograde", "isco_retrograde"} <= roles
        assert any(text.get_gid() == "isco_name" for text in ax.texts)
    finally:
        plt.close(fig)
