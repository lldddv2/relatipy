"""PerturbedKerr is a Kerr value plus one force with the Kerr orbit API."""

import copy
import pickle

import numpy as np
import pytest
from astropy import units as u
from astropy.constants import c

from relatipy import IntegrationError, Kerr, Orbit, PerturbedKerr, Solution
from relatipy.forces import Plummer

MASS = 4e6 * u.M_sun
SPIN = 0.6
KERR = Kerr(mass=MASS, spin=SPIN)
ZERO = Plummer(mass=0 * u.M_sun, radius=0.01 * u.pc)
SPHERE = Plummer(mass=4e5 * u.M_sun, radius=100 * KERR.r_g)
ORBIT = {
    "a": 200 * KERR.r_g, "e": 0.5, "inc": 0.7 * u.rad,
    "Omega": 0.3 * u.rad, "omega": 0.2 * u.rad,
}
SPAN = (0 * u.s, (3e4 * KERR.r_g / c).to(u.s))


def _solve(metric, method, **kwargs):
    return metric.orbit(**ORBIT).solve(tau_span=SPAN, method=method, **kwargs)


def test_value_semantics_repr_and_background_metric() -> None:
    metric = PerturbedKerr(mass=MASS, spin=SPIN, force=SPHERE)
    assert metric.kerr == KERR
    assert metric.force is SPHERE
    assert metric.mass == MASS and metric.spin == SPIN
    assert metric.r_g == KERR.r_g
    assert metric.horizons == KERR.horizons
    assert metric == PerturbedKerr(mass=MASS, spin=SPIN, force=SPHERE)
    assert hash(metric) == hash(PerturbedKerr(mass=MASS, spin=SPIN, force=SPHERE))
    assert metric != PerturbedKerr(mass=MASS, spin=SPIN, force=ZERO)
    assert metric != KERR
    assert repr(metric).startswith(
        "PerturbedKerr(mass=<Quantity 4000000. solMass>, spin=0.6, force=Plummer("
    )
    with pytest.raises(AttributeError, match="immutable"):
        metric._force = ZERO
    with pytest.raises(AttributeError, match="immutable"):
        del metric._kerr


def test_constructor_validates_like_kerr_and_requires_a_force() -> None:
    with pytest.raises(TypeError, match="ExternalForce"):
        PerturbedKerr(mass=MASS, spin=SPIN, force=None)
    with pytest.raises(ValueError, match="spin"):
        PerturbedKerr(mass=MASS, spin=1.5, force=ZERO)
    with pytest.raises(ValueError, match="mass"):
        PerturbedKerr(mass=-MASS, spin=SPIN, force=ZERO)


@pytest.mark.parametrize("method", ["dop853", "dp45", "radau"])
def test_zero_force_reproduces_kerr_orbit(method) -> None:
    reference = _solve(KERR, method)
    perturbed = _solve(PerturbedKerr(mass=MASS, spin=SPIN, force=ZERO), method)
    assert type(perturbed) is Solution
    assert perturbed.status == reference.status == 0
    assert perturbed.integration.method == method
    assert perturbed.tau[-1] == reference.tau[-1]
    scale = KERR.r_g.to_value(u.m)
    for name in ("x", "y", "z"):
        np.testing.assert_allclose(
            getattr(perturbed, name)[-1].to_value(u.m),
            getattr(reference, name)[-1].to_value(u.m), rtol=1e-12, atol=1e-12 * scale,
        )
    samples = np.linspace(0, 1, 7) * SPAN[1]
    np.testing.assert_allclose(
        perturbed.at(tau=samples).R.to_value(u.m),
        reference.at(tau=samples).R.to_value(u.m), rtol=1e-10,
    )


def test_force_changes_the_orbit_and_integrate_matches_solve() -> None:
    metric = PerturbedKerr(mass=MASS, spin=SPIN, force=SPHERE)
    perturbed = _solve(metric, "dop853")
    reference = _solve(KERR, "dop853")
    shift = abs(perturbed.R[-1] - reference.R[-1]) / KERR.r_g
    assert shift > 1
    orbit = metric.orbit(**ORBIT)
    orbit.integrate(SPAN[1], method="dop853")
    assert isinstance(orbit, Orbit)
    np.testing.assert_allclose(
        orbit.R.to_value(u.m), perturbed.R[-1].to_value(u.m), rtol=1e-12
    )


def test_copy_keeps_the_force() -> None:
    orbit = PerturbedKerr(mass=MASS, spin=SPIN, force=SPHERE).orbit(**ORBIT)
    clone = orbit.copy()
    first = orbit.solve(tau_span=SPAN, method="dop853")
    second = clone.solve(tau_span=SPAN, method="dop853")
    np.testing.assert_array_equal(first.x.value, second.x.value)


def test_projection_radau_is_rejected() -> None:
    orbit = PerturbedKerr(mass=MASS, spin=SPIN, force=ZERO).orbit(**ORBIT)
    with pytest.raises(ValueError, match="projection_radau.*Carter"):
        orbit.solve(tau_span=SPAN, method="projection_radau")
    with pytest.raises(ValueError, match="projection_radau"):
        orbit.integrate(SPAN[1], method="projection_radau")
    with pytest.raises(ValueError, match="projection_radau"):
        orbit.integrate(0 * u.s, method="projection_radau")


def test_numerical_failure_keeps_partial_solution_semantics() -> None:
    orbit = PerturbedKerr(mass=MASS, spin=SPIN, force=SPHERE).orbit(**ORBIT)
    step = SPAN[1] / 2e6
    solution = orbit.solve(tau_span=SPAN, method="dp45", max_step=step)
    assert solution.status == -1 and not solution.success
    assert "maximum number of steps" in solution.message
    assert solution.tau[-1] < SPAN[1]
    with pytest.raises(IntegrationError):
        orbit.integrate(SPAN[1], method="dp45", max_step=step)


@pytest.mark.parametrize(
    "roundtrip",
    [copy.copy, copy.deepcopy, lambda value: pickle.loads(pickle.dumps(value))],
)
def test_metric_orbit_and_solution_survive_copy_and_pickle(roundtrip) -> None:
    metric = PerturbedKerr(mass=MASS, spin=SPIN, force=SPHERE)
    assert roundtrip(metric) == metric
    orbit = metric.orbit(**ORBIT)
    solution = orbit.solve(tau_span=SPAN, method="dop853")
    again = roundtrip(orbit).solve(tau_span=SPAN, method="dop853")
    np.testing.assert_array_equal(again.x.value, solution.x.value)
    assert len(roundtrip(solution)) == len(solution)
