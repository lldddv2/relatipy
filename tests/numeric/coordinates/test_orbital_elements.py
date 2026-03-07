import numpy as np
import astropy.units as u_astropy
from relatipy.numeric.coordinates import OrbitalElements, Cartesian

# CI 1
a = 1.0
e = 0.1
inc = 0.2
Omega = 0.3
omega = 0.4
f = 0.5
mass = 1.0


# CI 2
ts = [0.0, 1.0, 2.0, 3.0]
as_ = [1.0, 2.0, 3.0, 4.0]
es = [0.1, 0.2, 0.3, 0.4]
incs = [0.2, 0.3, 0.4, 0.5]
Omegas = [0.3, 0.4, 0.5, 0.6]
omegas = [0.4, 0.5, 0.6, 0.7]
fs = [0.5, 0.6, 0.7, 0.8]
mass_2 = 1.0

# CI 3: estado cartesiano de una órbita acotada (e < 1) para que el round-trip sea válido
_mass_3 = 1e30 * u_astropy.kg  # masa central en kg explícito
_oe_ref = OrbitalElements(0.0, 1e7, 0.1, 0.2, 0.3, 0.4, 0.5, _mass_3)
_cart_ref = _oe_ref.convert_to("Cartesian")
xs = _cart_ref.xs
vs = _cart_ref.vs
mass_3 = _mass_3

class TestOrbitalElements:
    def test_orbital_elements_conservation(self):
        # CI 1
        orbital_elements = OrbitalElements(0, a, e, inc, Omega, omega, f, mass_2)
        orbital_elements_cartesian = orbital_elements.convert_to("Cartesian").convert_to("OrbitalElements", mass=mass_2)

        assert np.isclose(orbital_elements.t, orbital_elements_cartesian.t).all()
        assert np.isclose(orbital_elements.a, orbital_elements_cartesian.a).all()
        assert np.isclose(orbital_elements.e, orbital_elements_cartesian.e).all()
        assert np.isclose(orbital_elements.inc, orbital_elements_cartesian.inc).all()
        assert np.isclose(orbital_elements.Omega, orbital_elements_cartesian.Omega).all()
        assert np.isclose(orbital_elements.omega, orbital_elements_cartesian.omega).all()
        assert np.isclose(orbital_elements.f, orbital_elements_cartesian.f).all()
        assert np.isclose(orbital_elements.mass, orbital_elements_cartesian.mass).all()
        
        # CI 2
        orbital_elements = OrbitalElements(ts, as_, es, incs, Omegas, omegas, fs, mass)
        orbital_elements_cartesian = orbital_elements.convert_to("Cartesian").convert_to("OrbitalElements", mass=mass_2)

        assert np.isclose(orbital_elements.t, orbital_elements_cartesian.t).all()
        assert np.isclose(orbital_elements.a, orbital_elements_cartesian.a).all()
        assert np.isclose(orbital_elements.e, orbital_elements_cartesian.e).all()
        assert np.isclose(orbital_elements.inc, orbital_elements_cartesian.inc).all()
        assert np.isclose(orbital_elements.Omega, orbital_elements_cartesian.Omega).all()
        assert np.isclose(orbital_elements.omega, orbital_elements_cartesian.omega).all()
        assert np.isclose(orbital_elements.f, orbital_elements_cartesian.f).all()
        assert np.isclose(orbital_elements.mass, orbital_elements_cartesian.mass)

        # CI 3: round-trip Cartesian -> OrbitalElements -> Cartesian (solo índices 0-6: t,x,y,z,vx,vy,vz)
        orbital_elements = Cartesian(xs, vels=vs)
        orbital_elements_orbital_elements = orbital_elements.convert_to("OrbitalElements", mass=mass_3).convert_to("Cartesian")

        assert np.isclose(orbital_elements[0], orbital_elements_orbital_elements[0]).all()
        assert np.isclose(orbital_elements[1], orbital_elements_orbital_elements[1]).all()
        assert np.isclose(orbital_elements[2], orbital_elements_orbital_elements[2]).all()
        assert np.isclose(orbital_elements[3], orbital_elements_orbital_elements[3]).all()
        assert np.isclose(orbital_elements[4], orbital_elements_orbital_elements[4]).all()
        assert np.isclose(orbital_elements[5], orbital_elements_orbital_elements[5]).all()
        assert np.isclose(orbital_elements[6], orbital_elements_orbital_elements[6]).all()
