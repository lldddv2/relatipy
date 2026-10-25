"""Static axisymmetric forces conserve E' and L' through the public API.

For a potential that is static and axisymmetric about the spin axis,
``E' = -u_t exp(Phi)`` and ``L' = u_phi exp(Phi)`` are exact invariants of
the projected forced geodesic equation. The orbit is solved with
:class:`relatipy.PerturbedKerr` and a Plummer sphere compiled at run time as
a :class:`relatipy.Potential`, whose acceleration is the exact gradient of
``Phi``; ``Phi`` is evaluated on the stored samples in the black-hole
Cartesian frame.
"""

import numpy as np
import pytest
from astropy import units as u
from astropy.constants import c

from relatipy import Kerr, PerturbedKerr, Potential, _force_core, _jit
from relatipy.geodesic.solution import constants_of_motion

pytest.importorskip("sympy")


def _has_compiler() -> bool:
    try:
        _jit.compiler_command()
    except _jit.CompilerNotFoundError:
        return False
    return True


pytestmark = pytest.mark.skipif(not _has_compiler(), reason="no C compiler")

MASS = 1e6 * u.M_sun
SPIN = 0.7
KERR = Kerr(mass=MASS, spin=SPIN)
SPHERE_MASS = 0.2
SCALE = 60.0
PARAMETERS = {"m": SPHERE_MASS * MASS, "b": SCALE * KERR.r_g}
SPHERES = {
    "spherical": Potential("-m / sqrt(r**2 + b**2)", coords="spherical",
                           parameters=PARAMETERS),
    "cylindrical": Potential("-m / sqrt(R**2 + z**2 + b**2)", coords="cylindrical",
                             parameters=PARAMETERS),
}
ORBIT = {
    "a": 120 * KERR.r_g, "e": 0.4, "inc": 0.9 * u.rad,
    "Omega": 0.4 * u.rad, "omega": 1.1 * u.rad,
}
SPAN = (0 * u.s, (2e5 * KERR.r_g / c).to(u.s))


@pytest.mark.parametrize("coords", sorted(SPHERES))
@pytest.mark.parametrize("method", ["dop853", "radau", "dp45"])
def test_aligned_plummer_potential_conserves_e_prime_and_l_prime(method, coords) -> None:
    metric = PerturbedKerr(mass=MASS, spin=SPIN, force=SPHERES[coords])
    np.testing.assert_allclose(
        metric._binding().model_data, [SPHERE_MASS, SCALE], rtol=1e-12
    )
    solution = metric.orbit(**ORBIT).solve(
        tau_span=SPAN, method=method, rtol=1e-12
    )
    assert solution.status == 0 and len(solution) > 100
    state = solution[:]
    constants = constants_of_motion(state)
    energy, angular_momentum = constants[:, 0], constants[:, 1]

    rows = np.column_stack((
        np.zeros(len(solution)),
        (state.R / KERR.r_g).to_value(u.one),
        state.Theta.to_value(u.rad),
        state.Phi.to_value(u.rad),
    ))
    positions = _force_core.bl_to_cartesian(SPIN, np.ascontiguousarray(rows))
    potential = -SPHERE_MASS / np.sqrt(np.sum(positions**2, axis=1) + SCALE**2)
    weight = np.exp(potential)
    e_prime = energy * weight
    l_prime = angular_momentum * weight

    # The force is strong: the Kerr constants change by far more than the
    # drift of the invariants allowed below.
    assert np.ptp(energy) / energy[0] > 1e-4
    assert np.ptp(angular_momentum) / abs(angular_momentum[0]) > 1e-4
    tolerance = {"dop853": 1e-10, "dp45": 1e-9, "radau": 1e-8}[method]
    assert np.max(np.abs(e_prime / e_prime[0] - 1)) < tolerance
    assert np.max(np.abs(l_prime / l_prime[0] - 1)) < tolerance
