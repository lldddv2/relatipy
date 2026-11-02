"""SurfaceDensityDisk and get_force: any axisymmetric Sigma(R), tabulated once."""

import pickle
import shutil

import numpy as np
import pytest
from astropy import units as u

import relatipy as rp
from relatipy import PerturbedKerrMcmcModel, _force_core
from relatipy.forces import Force, Potential, SurfaceDensityDisk, ThinDisk, get_force

pytest.importorskip("sympy")
pytestmark = pytest.mark.skipif(
    shutil.which("cc") is None and shutil.which("gcc") is None, reason="needs a C compiler")

M_BH = 4.29e6 * u.M_sun
PARSEC = (u.pc / rp.Kerr(mass=M_BH, spin=0.0).r_g).to_value(u.one)
RADII = np.geomspace(1e-4, 5, 21) * PARSEC
POINTS = np.concatenate([
    np.stack([RADII, 0 * RADII, 0 * RADII], axis=1),                                  # plane
    np.stack([0 * RADII, 0 * RADII, RADII], axis=1),                                  # axis
    np.stack([RADII * np.cos(0.3), 0 * RADII, RADII * np.sin(0.3)], axis=1),         # 17 deg
])


def _acceleration(force, points=POINTS):
    bound = force._bind(M_BH)
    return _force_core.acceleration(
        bound.callback, bound.model, np.ascontiguousarray(points), np.zeros_like(points))


def _relative(first, second):
    return np.max(np.linalg.norm(first - second, axis=1) / np.linalg.norm(second, axis=1))


@pytest.mark.parametrize(("slope", "r_in"), [(-1.0, 0.0), (-1.7, 0.0), (-3.0, 0.05), (-1.0, 0.2)])
def test_power_law_density_is_the_thin_disk(slope, r_in) -> None:
    density = get_force(
        surface_density="R**gamma", params={"gamma": slope},
        mass=1e3 * u.M_sun, r_out=1 * u.pc, r_in=r_in * u.pc)
    assert isinstance(density, SurfaceDensityDisk)
    reference = ThinDisk(mass=1e3 * u.M_sun, r_out=1 * u.pc, r_in=r_in * u.pc, slope=slope)
    assert _relative(_acceleration(density), _acceleration(reference)) < 1e-12


def _exponential_reference(points, h, eps, m_d, r_out):
    """Softened ring sum of Sigma = exp(-R / h) by scipy quad (independent of the native code)."""
    from scipy.integrate import quad
    from scipy.special import ellipe, ellipkm1

    norm = quad(lambda R: np.exp(-R / h) * R, 0, 1, epsabs=0, epsrel=1e-13)[0]
    rows = []
    for x, y, z in points / r_out:
        d = np.hypot(x, y)

        def parts(R):
            zeta2 = z * z + (eps * R) ** 2
            p2 = (d + R) ** 2 + zeta2
            q2 = (d - R) ** 2 + zeta2
            p = np.sqrt(p2)
            K, E = ellipkm1(q2 / p2), ellipe(4 * R * d / p2)
            dm = np.exp(-R / h) * R / norm
            a_d = -dm / (np.pi * p * d) * (K + E * (d * d - R * R - zeta2) / q2) if d > 1e-8 else 0.0
            return a_d, -2 * dm * E * z / (np.pi * p * q2)

        options = {"limit": 2000, "epsabs": 0, "epsrel": 1e-12, "points": [d] if 0 < d < 1 else None}
        a_d = quad(lambda R: parts(R)[0], 0, 1, **options)[0] if d > 1e-8 else 0.0
        a_z = quad(lambda R: parts(R)[1], 0, 1, **options)[0] if z else 0.0
        rows.append((a_d * x / d if d else 0.0, a_d * y / d if d else 0.0, a_z))
    return np.array(rows) * m_d / r_out**2


def test_exponential_disk_matches_an_independent_quadrature() -> None:
    disk = SurfaceDensityDisk(
        "exp(-R/h)", parameters={"h": 0.25}, mass=1e3 * u.M_sun, r_out=1 * u.pc)
    points = POINTS[[3, 8, 12, 17, 24, 30, 45, 50, 56]]
    m_d = (1e3 * u.M_sun / M_BH).to_value(u.one)
    reference = _exponential_reference(points, 0.25, 0.01, m_d, PARSEC)
    assert _relative(_acceleration(disk, points), reference) < 1e-6


def test_density_in_plane_coordinates_is_accepted() -> None:
    by_radius = SurfaceDensityDisk("exp(-R/h)", parameters={"h": 0.3}, mass=1 * u.M_sun, r_out=1 * u.au)
    by_xy = SurfaceDensityDisk(
        "exp(-sqrt(x**2 + y**2)/h)", parameters={"h": 0.3}, mass=1 * u.M_sun, r_out=1 * u.au)
    assert by_radius._shape == by_xy._shape


@pytest.mark.parametrize(
    ("text", "params", "r_in", "message"),
    [
        ("R*(1 + cos(phi))", {}, None, "axisymmetric"),
        ("x", {}, None, "axisymmetric"),
        ("R**gamma", {}, None, "unknown name"),
        ("R - 0.5", {}, None, "non-negative"),
        ("R**-2.5", {}, None, "diverges"),
    ],
)
def test_invalid_densities_are_rejected(text, params, r_in, message) -> None:
    with pytest.raises(ValueError, match=message):
        SurfaceDensityDisk(text, parameters=params, mass=1 * u.M_sun, r_out=1 * u.au, r_in=r_in)


def test_inner_edge_allows_a_steep_density() -> None:
    disk = SurfaceDensityDisk("R**-2.5", mass=1 * u.M_sun, r_out=1 * u.au, r_in=0.1 * u.au)
    assert np.all(np.isfinite(_acceleration(disk, POINTS[:5] * 1e-3)))


def test_value_semantics_pickle_and_run_time_parameters() -> None:
    disk = SurfaceDensityDisk("R**gamma", parameters={"gamma": -1.0},
                              mass=1e3 * u.M_sun, r_out=1362.5 * u.au)
    assert disk.parameter_names == ("disk_mass", "normal_theta", "normal_phi")
    copy = pickle.loads(pickle.dumps(disk))
    np.testing.assert_array_equal(_acceleration(copy), _acceleration(disk))
    heavier = disk.with_parameters(disk_mass=2e3 * u.M_sun)
    np.testing.assert_allclose(_acceleration(heavier), 2 * _acceleration(disk), rtol=1e-15)


def test_perturbed_model_accepts_a_surface_density() -> None:
    density = SurfaceDensityDisk("R**gamma", parameters={"gamma": -1.0},
                                 mass=4e3 * u.M_sun, r_out=1362.5 * u.au)
    thin = ThinDisk(mass=4e3 * u.M_sun, r_out=1362.5 * u.au)
    kerr = {"mass": 4.3e6 * u.M_sun, "spin": 0.5, "vec": (0.2, -0.3, 0.9)}
    orbit = {"a": 1000 * u.au, "e": 0.88, "inc": 2.3 * u.rad,
             "Omega": 4.0 * u.rad, "omega": 1.2 * u.rad, "t": 2018.4 * u.yr}
    sol = {"t_obs": np.linspace(2024.0, 2004.0, 15) * u.yr, "method": "dop853",
           "rtol": 1e-12, "atol": 1e-14}
    first = PerturbedKerrMcmcModel(force=density).get_ra_dec_vr(
        kerr=kerr, orbit=orbit, sol=sol, distance=8.2 * u.kpc)
    second = PerturbedKerrMcmcModel(force=thin).get_ra_dec_vr(
        kerr=kerr, orbit=orbit, sol=sol, distance=8.2 * u.kpc)
    for one, other in zip(first, second):
        np.testing.assert_allclose(one, other, rtol=0, atol=1e-10 * np.max(np.abs(other)))


def test_get_force_infers_coordinates() -> None:
    cylindrical = get_force(potential="-M/sqrt(R**2 + (a + sqrt(z**2 + b**2))**2)",
                            params={"M": 1e-3, "a": 1e4, "b": 1e3})
    assert isinstance(cylindrical, Potential) and cylindrical.coords == "cylindrical"
    spherical = get_force(potential="-M/sqrt(r**2 + b**2)", params={"M": 1e-3, "b": 1e3})
    assert spherical.coords == "spherical"
    drag = get_force(force="(-k*vx, -k*vy, -k*vz)", params={"k": 1e-3})
    assert isinstance(drag, Force) and drag.coords == "cartesian"


def test_get_force_mixes_coordinates_in_a_potential() -> None:
    mixed = get_force(potential="-M/sqrt(r**2 + b**2) + 0*R", params={"M": 1e-3, "b": 1e3})
    plain = get_force(potential="-M/sqrt(x**2 + y**2 + z**2 + b**2)", params={"M": 1e-3, "b": 1e3})
    point = np.array([[3e3, 4e3, 1e3]])
    np.testing.assert_allclose(_acceleration(mixed, point), _acceleration(plain, point), rtol=1e-14)


def test_get_force_rejects_ambiguous_calls() -> None:
    with pytest.raises(TypeError, match="exactly one"):
        get_force()
    with pytest.raises(TypeError, match="exactly one"):
        get_force(potential="-1/r", surface_density="R")
    with pytest.raises(TypeError, match="unexpected options"):
        get_force(potential="-1/r", mass=1 * u.M_sun)
    with pytest.raises(ValueError, match="mix names"):
        get_force(force="(-x, -R, 0)")
    with pytest.raises(TypeError, match="coords"):
        get_force(surface_density="R", coords="cylindrical", mass=1 * u.M_sun, r_out=1 * u.au)


def test_package_exports() -> None:
    assert rp.get_force is get_force and rp.SurfaceDensityDisk is SurfaceDensityDisk
