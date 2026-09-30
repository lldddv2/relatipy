"""Integrate a stable bound Kerr orbit from its parameters and Mino phases."""

from astropy import units as u

from relatipy import Kerr
from relatipy.coordinates import KerrOrbitalElements


def main() -> None:
    """Construct physical initial conditions and integrate in proper time."""
    black_hole = Kerr(mass=4e6 * u.Msun, spin=0.7)
    elements = KerrOrbitalElements(
        p=12 * black_hole.r_g,
        e=0.3,
        x=0.7,
        q_r0=0.4 * u.rad,
        q_theta0=1.2 * u.rad,
        q_phi0=0.2 * u.rad,
    )
    orbit = black_hole.orbit(elements=elements)
    solution = orbit.solve(
        tau_span=(0 * u.s, 1000 * u.s),
        method="dop853",
        rtol=1e-10,
        atol=1e-12,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    print(f"Stored {len(solution)} states; final proper time: {solution.tau[-1]}")


if __name__ == "__main__":
    main()
