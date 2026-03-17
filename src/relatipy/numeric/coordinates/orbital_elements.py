import numpy
import spiceypy as spice
from ..constants import _G, _G_SI, _c_SI
from ..utils.dimensions import validator

try:
    import astropy.units as u
except ImportError:
    u = None

_M_ref_kg = 1.98892e30  # kg — 1 unidad geométrica de masa = 1 masa solar


def _mass_to_kg(mass):
    """
    Return mass as float in kg.
    - astropy Quantity con unidad de masa: convierte directamente a kg.
    - float/int: interpreta como unidades geométricas (M=1 → 1 masa solar = _M_ref_kg kg).
    """
    if u is not None and hasattr(mass, "to_value"):
        return float(mass.to_value(u.kg))
    return float(mass) * _M_ref_kg


def _state_to_spice_units(xs, vs, mass_kg):
    """
    Convert state to SPICE units (km, km/s, s).
    Expects xs=(t, x, y, z), vs=(vx, vy, vz) with position in geometric units [GM/c²],
    velocity in [v/c] and time in geometric units [GM/c³].
    """
    xs = numpy.asarray(xs, dtype=float)
    vs = numpy.asarray(vs, dtype=float)
    length_unit = _G_SI * mass_kg / _c_SI**2  # GM/c² en metros
    time_unit = _G_SI * mass_kg / _c_SI**3    # GM/c³ en segundos
    pos_m = xs[1:4] * length_unit
    vel_dimless = vs  # v/c from validator
    vel_m_s = vel_dimless * _c_SI
    pos_km = pos_m / 1000.0
    vel_km_s = vel_m_s / 1000.0
    t_geom = xs[0]
    t_sec = t_geom * time_unit
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

        length_unit = _G_SI * self.mass / _c_SI**2  # GM/c² en metros
        time_unit = _G_SI * self.mass / _c_SI**3    # GM/c³ en segundos
        rp_km = a * length_unit * (1 - e) / 1000.0  # unidades geométricas → km
        one_minus_e = numpy.maximum(1 - e, 1e-15)
        E = 2 * numpy.arctan2(
            numpy.sqrt(one_minus_e) * numpy.sin(f / 2),
            numpy.sqrt(1 + e) * numpy.cos(f / 2)
        )
        M0 = E - e * numpy.sin(E)
        t0_sec = float(t) * time_unit  # unidades geométricas → segundos
        return [rp_km, e, inc, Omega, omega, M0, t0_sec, mu]

    def _to_cartesian_arrays(self):
        """
        Convierte elementos orbitales a posición y velocidad cartesianas 3D
        usando spiceypy.conics.

        Devuelve:
        - r_geom: array de forma (3,) en unidades geométricas [GM/c²]
        - v_dimless: array de forma (3,) en v/c (adimensional)
        Para el caso vectorial:
        - r_geom: array de forma (N, 3) en unidades geométricas [GM/c²]
        - v_dimless: array de forma (N, 3) en v/c (adimensional)
        """
        # Detectar si estamos en modo "batch" (varias órbitas)
        is_batch = isinstance(self.a, numpy.ndarray) and self.a.ndim > 0 and self.a.size > 1

        length_unit = _G_SI * self.mass / _c_SI**2  # GM/c² en metros
        time_unit = _G_SI * self.mass / _c_SI**3    # GM/c³ en segundos

        if not is_batch:
            elts = self._get_elts()
            t_sec = float(self.t) * time_unit  # unidades geométricas → segundos
            state_km = spice.conics(elts, t_sec)  # km, km/s
            r_km = numpy.array(state_km[:3])
            v_km_s = numpy.array(state_km[3:])
            r_m = r_km * 1000.0
            r_geom = r_m / length_unit  # metros → unidades geométricas
            v_m_s = v_km_s * 1000.0
            v_dimless = v_m_s / _c_SI  # m/s → v/c
            return r_geom, v_dimless

        # Modo vectorial: iterar sobre cada conjunto de elementos orbitales
        n_orbits = len(self[0])
        r_list = []
        v_list = []
        for i in range(n_orbits):
            elts_i = self._get_elts(index=i)
            t_i = self.t[i] if isinstance(self.t, numpy.ndarray) else self.t
            t_sec_i = float(t_i) * time_unit  # unidades geométricas → segundos
            state_km_i = spice.conics(elts_i, t_sec_i)
            r_km_i = numpy.array(state_km_i[:3])
            v_km_s_i = numpy.array(state_km_i[3:])
            r_geom_i = (r_km_i * 1000.0) / length_unit  # km → m → unidades geométricas
            v_dimless_i = (v_km_s_i * 1000.0) / _c_SI   # km/s → m/s → v/c
            r_list.append(r_geom_i)
            v_list.append(v_dimless_i)

        r_geom = numpy.stack(r_list, axis=0)    # (N, 3)
        v_dimless = numpy.stack(v_list, axis=0) # (N, 3)
        return r_geom, v_dimless

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
            4-posición [t, x, y, z] (t en unidades geométricas [GM/c³], posición en unidades geométricas [GM/c²])
        vs : array-like, shape (3,)
            Velocidad cartesiana [vx, vy, vz] en v/c (sin unidades)
        mass : float or Quantity
            Masa del cuerpo central [kg]
        """
        mass_kg = _mass_to_kg(mass)
        pos_km, vel_km_s, t_sec = _state_to_spice_units(xs, vs, mass_kg)
        if t is not None:
            t_sec = float(t) if (u is None or not hasattr(t, "to_value")) else float(t.to_value(u.s))
        state_km = numpy.concatenate([pos_km, vel_km_s]).astype(float)
        mu_km = (_G_SI * mass_kg) / 1e9

        elts = spice.oscelt(state_km, t_sec, mu_km)
        rp_km, e, inc, Omega, omega, M0, t0, _ = elts

        M0 = numpy.mod(float(M0), 2 * numpy.pi)

        rp_m = rp_km * 1000.0
        a_m = rp_m / (1 - e)

        length_unit = _G_SI * mass_kg / _c_SI**2  # GM/c² en metros
        time_unit = _G_SI * mass_kg / _c_SI**3    # GM/c³ en segundos
        a_geom = a_m / length_unit

        E = M0
        for _ in range(50):
            E = E - (E - e * numpy.sin(E) - M0) / (1 - e * numpy.cos(E))

        one_minus_e = numpy.maximum(1 - e, 1e-15)
        f = 2 * numpy.arctan2(
            numpy.sqrt(1 + e) * numpy.sin(E / 2),
            numpy.sqrt(one_minus_e) * numpy.cos(E / 2)
        )

        t_geom = t_sec / time_unit
        return OrbitalElements(t_geom, a_geom, e, inc, Omega, omega, f, mass)

    @staticmethod
    def _convert_from_cartesian(xs_p, vs_p, **kwargs):
        """
        Alias de from_cartesian para compatibilidad con CoordinateBase.convert_to.
        """
        mass = kwargs.get("mass")
        if mass is None:
            raise ValueError("mass is required to convert to OrbitalElements.")
        return OrbitalElements.from_cartesian(xs_p, vs_p, mass)

    def _get_period(self):
        """
        Calcula el periodo orbital en unidades geométricas [GM/c³].
        """
        L = _G_SI * self.mass / _c_SI**2  # m por unidad geométrica (para esta masa)
        T = _G_SI * self.mass / _c_SI**3  # s por unidad geométrica (para esta masa)

        a_m  = self.a * L
        mu_si = _G_SI * self.mass

        T_sec = 2 * numpy.pi * numpy.sqrt(a_m**3 / mu_si)
        return T_sec / T