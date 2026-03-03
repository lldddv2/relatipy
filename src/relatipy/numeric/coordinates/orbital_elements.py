import numpy
import spiceypy as spice
from ..constants import _G, _G_SI, _c_SI
from ..utils.dimensions import validator

try:
    import astropy.units as u
except ImportError:
    u = None

def _mass_to_kg(mass):
    """Return mass as float in kg. Accepts float or astropy Quantity with mass unit."""
    if u is not None and hasattr(mass, "to_value"):
        return float(mass.to_value(u.kg))
    return float(mass)


def _state_to_spice_units(xs, vs):
    """
    Convert state to SPICE units (km, km/s, s).
    Expects xs=(t, x, y, z), vs=(vx, vy, vz) with position in [m], velocity in [v/c]
    and time in geometric units [s/c_SI] (as from the coordinate validator).
    """
    xs = numpy.asarray(xs, dtype=float)
    vs = numpy.asarray(vs, dtype=float)
    pos_m = xs[1:4]
    vel_dimless = vs  # v/c from validator
    vel_m_s = vel_dimless * _c_SI
    pos_km = pos_m / 1000.0
    vel_km_s = vel_m_s / 1000.0
    t_geom = xs[0]
    t_sec = t_geom * _c_SI
    return pos_km, vel_km_s, t_sec


class OrbitalElements:
    def __init__(self, t, a, e, inc, Omega, omega, f, mass=1):
        """
        Elementos orbitales clásicos (Keplerianos).

        Parameters
        ----------
        t : float
            Tiempo coordenado [s], por defecto 0.0
        a : float
            Semi-eje mayor [m]
        e : float
            Excentricidad (0 <= e < 1 para órbitas cerradas)
        inc : float
            Inclinación [rad]
        Omega : float
            Longitud del nodo ascendente (RAAN) [rad]
        omega : float
            Argumento del periastro [rad]
        f : float
            Anomalía verdadera [rad]
        mass : float
            Masa del cuerpo central [kg]
        """

        def _as_array_or_scalar(value):
            """Return numpy array for sequences, scalar float otherwise."""
            if isinstance(value, (list, tuple)):
                return numpy.asarray(value, dtype=float)
            arr = numpy.asarray(value)
            if arr.shape == () or arr.size == 1:
                return float(arr)
            return arr.astype(float)

        # Permitir tanto escalares como listas / arrays para los elementos orbitales
        self.t = _as_array_or_scalar(t)
        self.a = _as_array_or_scalar(a)
        self.e = _as_array_or_scalar(e)
        self.inc = _as_array_or_scalar(inc)
        self.Omega = _as_array_or_scalar(Omega)
        self.omega = _as_array_or_scalar(omega)
        self.f = _as_array_or_scalar(f)

        # La masa puede ser escalar o lista; para listas usamos el validador existente
        self.mass = _mass_to_kg(mass)
        self.mu = (_G_SI * self.mass) / 1e9  # km³/s² para SPICE
        self.name_metric = "OrbitalElements"
        self.kwargs = {"mass": mass}

        self.state_vector = numpy.array((self.t, self.a, self.e, self.inc, self.Omega, self.omega, self.f))


    # si se hace []
    def __getitem__(self, index):
        return self.state_vector[index]

    def __setitem__(self, index, value):
        self.state_vector[index] = value

    def _get_elts(self, index=None):
        """
        Construye el vector de elementos orbitales en formato SPICE (osculating elements).
        elts = [rp, ecc, inc, lnode, argp, m0, t0, mu] en km, rad, km³/s².
        """
        # Seleccionar escalar o elemento i-ésimo para el caso vectorial
        if index is not None:
            a = self.a[index]
            e = self.e[index]
            inc = self.inc[index]
            Omega = self.Omega[index]
            omega = self.omega[index]
            f = self.f[index]
            t = self.t[index] if isinstance(self.t, numpy.ndarray) else self.t
            mu = self.mu[index] if isinstance(self.mu, numpy.ndarray) else self.mu
        else:
            a = self.a
            e = self.e
            inc = self.inc
            Omega = self.Omega
            omega = self.omega
            f = self.f
            t = self.t
            mu = self.mu

        rp_km = a * (1 - e) / 1000.0  # a en m -> rp en km
        # Convertir anomalía verdadera f a anomalía media M
        # E = 2*arctan(sqrt((1-e)/(1+e)) * tan(f/2))  -> anomalía excéntrica
        # M = E - e*sin(E). Para e>=1 (órbita abierta) se evita NaN en sqrt(1-e).
        one_minus_e = numpy.maximum(1 - e, 1e-15)
        E = 2 * numpy.arctan2(
            numpy.sqrt(one_minus_e) * numpy.sin(f / 2),
            numpy.sqrt(1 + e) * numpy.cos(f / 2)
        )
        M0 = E - e * numpy.sin(E)
        # SPICE conics espera la época T0 en segundos, no en unidades geométricas
        t0_sec = float(t) * _c_SI
        return [rp_km, e, inc, Omega, omega, M0, t0_sec, mu]

    def _to_cartesian_arrays(self):
        """
        Convierte elementos orbitales a posición y velocidad cartesianas 3D
        usando spiceypy.conics.

        Devuelve:
        - r_m: array de forma (3,) y v_dimless de forma (3,) para el caso escalar
        - r_m: array de forma (N, 3) y v_dimless de forma (N, 3) para el caso vectorial
        """
        # Detectar si estamos en modo "batch" (varias órbitas)
        is_batch = isinstance(self.a, numpy.ndarray) and self.a.ndim > 0 and self.a.size > 1

        if not is_batch:
            # Comportamiento escalar original
            elts = self._get_elts()
            t_sec = float(self.t) * _c_SI  # SPICE usa segundos
            state_km = spice.conics(elts, t_sec)  # km, km/s
            r_km = numpy.array(state_km[:3])
            v_km_s = numpy.array(state_km[3:])
            r_m = r_km * 1000.0
            v_m_s = v_km_s * 1000.0
            v_dimless = v_m_s / _c_SI  # v/c
            return r_m, v_dimless

        # Modo vectorial: iterar sobre cada conjunto de elementos orbitales
        n_orbits = len(self[0])
        r_list = []
        v_list = []
        for i in range(n_orbits):
            elts_i = self._get_elts(index=i)
            t_i = self.t[i] if isinstance(self.t, numpy.ndarray) else self.t
            t_sec_i = float(t_i) * _c_SI
            state_km_i = spice.conics(elts_i, t_sec_i)
            r_km_i = numpy.array(state_km_i[:3])
            v_km_s_i = numpy.array(state_km_i[3:])
            r_list.append(r_km_i * 1000.0)
            v_list.append((v_km_s_i * 1000.0) / _c_SI)

        r_m = numpy.stack(r_list, axis=0)       # (N, 3)
        v_dimless = numpy.stack(v_list, axis=0) # (N, 3)
        return r_m, v_dimless

    def convert_to_cartesian(self):
        """Convierte a coordenadas cartesianas (4-posición + 3-velocidad)."""
        from .cartesian import Cartesian

        r_vec, v_vec = self._to_cartesian_arrays()
        # Caso escalar
        if r_vec.ndim == 1:
            xs = numpy.array([self.t, r_vec[0], r_vec[1], r_vec[2]])
            vs = numpy.array([v_vec[0], v_vec[1], v_vec[2]])
        else:
            # Caso vectorial: r_vec y v_vec tienen forma (N, 3)
            t_array = (
                self.t if isinstance(self.t, numpy.ndarray) else numpy.asarray([self.t] * r_vec.shape[0])
            )
            xs = numpy.array(
                [t_array, r_vec[:, 0], r_vec[:, 1], r_vec[:, 2]]
            )  # (4, N)
            vs = numpy.array(
                [v_vec[:, 0], v_vec[:, 1], v_vec[:, 2]]
            )  # (3, N)

        cartesian = Cartesian(xs, vels=vs, from_dxs_dt=False)

        # Propagar los parámetros orbitales para permitir comprobaciones de conservación
        cartesian.t = self.t
        cartesian.a = self.a
        cartesian.e = self.e
        cartesian.inc = self.inc
        cartesian.Omega = self.Omega
        cartesian.omega = self.omega
        cartesian.f = self.f
        cartesian.mass = self.mass

        return cartesian

    def convert_to(self, target_system, **kwargs):
        """
        Convierte a otro sistema de coordenadas pasando por cartesiana.

        Parameters
        ----------
        target_system : str
            Nombre del sistema destino (debe estar en coordinate_systems).
        """
        from . import coordinate_systems

        if target_system not in coordinate_systems:
            raise ValueError(
                f"Unsupported target coordinate system: {target_system}. "
                f"Supported systems are: {list(coordinate_systems.keys())}"
            )

        if target_system == "OrbitalElements" and "mass" not in kwargs:
            raise ValueError("mass is required to convert to OrbitalElements.")

        cartesian = self.convert_to_cartesian()
        if target_system == "Cartesian":
            return cartesian

        target_class = coordinate_systems[target_system]
        return target_class._convert_from_cartesian(cartesian.xs, cartesian.vs, **kwargs)

    @staticmethod
    def from_cartesian(xs, vs, mass, t=None):
        """
        Construye OrbitalElements desde posición y velocidad cartesianas
        usando spiceypy.oscelt.
        """
        if isinstance(xs[1], numpy.ndarray):
            return OrbitalElements.from_cartesian_vector(xs, vs, mass, t)
        return OrbitalElements.from_cartesian_scalar(xs, vs, mass, t)

    @staticmethod
    def from_cartesian_vector(xs, vs, mass, t=None):
        """
        Construye OrbitalElements desde posición y velocidad cartesianas
        usando spiceypy.oscelt.
        """
        N = len(xs[0])
        ts = numpy.zeros(N)
        as_ = numpy.zeros(N)
        es = numpy.zeros(N)
        incs = numpy.zeros(N)
        Omegas = numpy.zeros(N)
        omegas = numpy.zeros(N)
        fs = numpy.zeros(N)
        for i in range(N):
            ts[i], as_[i], es[i], incs[i], Omegas[i], omegas[i], fs[i] = OrbitalElements.from_cartesian_scalar(xs[:, i], vs[:, i], mass)
        return OrbitalElements(ts, as_, es, incs, Omegas, omegas, fs, mass)

    @staticmethod
    def from_cartesian_scalar(xs, vs, mass, t=None):
        """
        Construye OrbitalElements desde posición y velocidad cartesianas
        usando spiceypy.oscelt.

        Parameters
        ----------
        xs : array-like, shape (4,)
            4-posición [t, x, y, z] (t en unidades geométricas, posición en m)
        vs : array-like, shape (3,)
            Velocidad cartesiana [vx, vy, vz] en v/c (sin unidades)
        mass : float or Quantity
            Masa del cuerpo central [kg]
        """
        mass_kg = _mass_to_kg(mass)
        pos_km, vel_km_s, t_sec = _state_to_spice_units(xs, vs)
        if t is not None:
            t_sec = float(t) if (u is None or not hasattr(t, "to_value")) else float(t.to_value(u.s))
        state_km = numpy.concatenate([pos_km, vel_km_s]).astype(float)
        mu_km = (_G_SI * mass_kg) / 1e9  # m³/s² -> km³/s²

        elts = spice.oscelt(state_km, t_sec, mu_km)
        # elts = [rp, ecc, inc, lnode, argp, m0, t0, mu] in km, rad, km³/s²
        rp_km, e, inc, Omega, omega, M0, t0, _ = elts

        # Normalizar M0 a [0, 2π) para consistencia en el ciclo OE -> Cartesian -> OE
        M0 = numpy.mod(float(M0), 2 * numpy.pi)

        rp_m = rp_km * 1000.0
        a_m = rp_m / (1 - e)

        # Anomalía media M0 (en época t_sec) → excéntrica E → verdadera f
        E = M0
        for _ in range(50):
            E = E - (E - e * numpy.sin(E) - M0) / (1 - e * numpy.cos(E))

        one_minus_e = numpy.maximum(1 - e, 1e-15)
        f = 2 * numpy.arctan2(
            numpy.sqrt(1 + e) * numpy.sin(E / 2),
            numpy.sqrt(one_minus_e) * numpy.cos(E / 2)
        )

        t_geom = t_sec / _c_SI  # almacenar en unidades geométricas
        return OrbitalElements(t_geom, a_m, e, inc, Omega, omega, f, mass_kg)

    @staticmethod
    def _convert_from_cartesian(xs_p, vs_p, **kwargs):
        """
        Alias de from_cartesian para compatibilidad con CoordinateBase.convert_to.
        """
        mass = kwargs.get("mass")
        if mass is None:
            raise ValueError("mass is required to convert to OrbitalElements.")
        return OrbitalElements.from_cartesian(xs_p, vs_p, mass)