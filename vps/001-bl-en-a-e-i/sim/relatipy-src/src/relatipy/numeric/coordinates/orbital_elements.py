"""
Classical Keplerian orbital elements and conversion to/from Cartesian states.

This module provides :class:`OrbitalElements`, which stores the standard set
of osculating elements and uses `spiceypy` (SPICE ``conics`` / ``oscelt``) to
convert between elements and Cartesian position and velocity. Time
:math:`t` and semi-major axis :math:`a` are expressed in **geometric units**
(:math:`GM/c^3` and :math:`GM/c^2`) consistent with the rest of the numeric
stack; angular elements are given in **degrees** in the public API but are
converted internally to radians for SPICE.

Notes
-----
When ``astropy.units`` is installed, a central mass may be passed as a
:class:`astropy.units.Quantity` with dimensions of mass; otherwise a numeric
``mass`` is interpreted as geometric units where :math:`M = 1` corresponds to
one solar mass (reference mass :math:`M_\\odot \\approx 1.98892\\times 10^{30}` kg).

See Also
--------
relatipy.numeric.coordinates.cartesian.Cartesian :
    Target type for :meth:`OrbitalElements.convert_to_cartesian`.

Examples
--------
>>> import numpy as np
>>> from relatipy.numeric.coordinates.orbital_elements import OrbitalElements
>>> oe = OrbitalElements(t=0.0, a=100.0, e=0.1, inc=30.0, Omega=10.0, omega=20.0, f=45.0, mass=1.0)
>>> oe.state_vector.shape
(7,)
"""

import numpy
import spiceypy as spice
from ..constants import _G_SI, _c_SI

try:
    import astropy.units as u
except ImportError:
    u = None

# kg — one geometric mass unit = one solar mass (reference)
_M_ref_kg = 1.98892e30


def _mass_to_kg(mass):
    """
    Convert central-body mass to kilograms.

    If ``mass`` is an ``astropy.units.Quantity`` with mass dimension, it is
    converted to kg. Otherwise it is treated as a dimensionless geometric mass
    and multiplied by the reference solar mass in kg.

    Parameters
    ----------
    mass : float or astropy.units.Quantity
        Central mass. Scalar numeric values are interpreted as geometric units
        (:math:`M=1` → one solar mass in kg).

    Returns
    -------
    float
        Mass in kilograms.

    Examples
    --------
    >>> from relatipy.numeric.coordinates.orbital_elements import _mass_to_kg
    >>> bool(numpy.isclose(_mass_to_kg(1.0), 1.98892e30))
    True
    """
    if u is not None and hasattr(mass, "to_value"):
        return float(mass.to_value(u.kg))
    return float(mass) * _M_ref_kg


def _state_to_spice_units(xs, vs, mass_kg):
    """
    Map a geometric Cartesian state to SPICE physical units (km, km/s, s).

    Parameters
    ----------
    xs : array_like, shape (4,)
        Event :math:`(t, x, y, z)` with :math:`t` in geometric time
        :math:`GM/c^3` and spatial components in geometric length :math:`GM/c^2`.
    vs : array_like, shape (3,)
        Cartesian velocity components :math:`v^i/c` (dimensionless).
    mass_kg : float
        Central mass in kg, used to define the geometric–SI scale.

    Returns
    -------
    pos_km : ndarray, shape (3,)
        Position in km.
    vel_km_s : ndarray, shape (3,)
        Velocity in km/s.
    t_sec : float
        Time in seconds.

    Notes
    -----
    The conversion uses :math:`L = GM/c^2` (meters per geometric length) and
    :math:`T = GM/c^3` (seconds per geometric time) for the given ``mass_kg``.

    Examples
    --------
    >>> import numpy as np
    >>> from relatipy.numeric.coordinates.orbital_elements import _state_to_spice_units
    >>> xs = np.array([0.0, 1.0, 0.0, 0.0])
    >>> vs = np.zeros(3)
    >>> pk, vk, ts = _state_to_spice_units(xs, vs, 1.98892e30)
    >>> pk.shape
    (3,)
    """
    xs = numpy.asarray(xs, dtype=float)
    vs = numpy.asarray(vs, dtype=float)
    length_unit = _G_SI * mass_kg / _c_SI**2  # GM/c² in meters
    time_unit = _G_SI * mass_kg / _c_SI**3    # GM/c³ in seconds
    pos_m = xs[1:4] * length_unit
    vel_dimless = vs  # v/c (dimensionless)
    vel_m_s = vel_dimless * _c_SI
    pos_km = pos_m / 1000.0
    vel_km_s = vel_m_s / 1000.0
    t_geom = xs[0]
    t_sec = t_geom * time_unit
    return pos_km, vel_km_s, t_sec


class OrbitalElements:
    """
    Keplerian (osculating) orbital elements with SPICE-backed Cartesian maps.

    Stores time, semi-major axis, eccentricity, inclination, longitude of the
    ascending node, argument of periapsis, and true anomaly. Public angles are
    in degrees; time and semi-major axis use geometric units consistent with
    Cartesian coordinates in this package.

    Parameters
    ----------
    t : float or array_like, optional
        Coordinate time in geometric units :math:`GM/c^3`. Scalar or sequence
        for batch orbits. Default ``0.0``.
    a : float or array_like, optional
        Semi-major axis in geometric length units :math:`GM/c^2`. Default ``1.0``
        (avoids a degenerate :math:`a=0` conic with SPICE).
    e : float or array_like, optional
        Eccentricity; for bound orbits typically :math:`0 \\leq e < 1`.
        Default ``0.0``.
    inc : float or array_like, optional
        Inclination in degrees. Default ``0.0``.
    Omega : float or array_like, optional
        Longitude of the ascending node (RAAN) in degrees. Default ``0.0``.
    omega : float or array_like, optional
        Argument of periapsis in degrees. Default ``0.0``.
    f : float or array_like, optional
        True anomaly in degrees. Default ``0.0``.
    mass : float or astropy.units.Quantity, optional
        Central mass: geometric scalar (default ``1`` = one solar mass in kg) or
        an ``astropy`` mass quantity converted to kg.

    Attributes
    ----------
    t, a, e, inc, Omega, omega, f : float or ndarray
        Orbital elements as provided (angles in degrees).
    mass : float
        Central mass in kg.
    mu : float or ndarray
        Gravitational parameter :math:`\\mu = GM` in km³/s² for SPICE.
    name_metric : str
        Registry label for this coordinate system (``OrbitalElements``).
    kwargs : dict
        Original keyword arguments (e.g. ``{\"mass\": mass}``).
    state_vector : ndarray
        Concatenation ``(t, a, e, inc, Omega, omega, f)``.

    Examples
    --------
    >>> import numpy as np
    >>> from relatipy.numeric.coordinates.orbital_elements import OrbitalElements
    >>> oe0 = OrbitalElements()
    >>> oe0.a, oe0.e
    (1.0, 0.0)
    >>> oe = OrbitalElements(0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, mass=1.0)
    >>> len(oe.state_vector)
    7
    """

    def __init__(self, t=0.0, a=1.0, e=0.0, inc=0.0, Omega=0.0, omega=0.0, f=0.0, mass=1):
        def _as_array_or_scalar(value):
            """
            Normalize an element to a float, or to a float ndarray if sequence-like.

            Parameters
            ----------
            value : scalar or sequence
                Orbital element value.

            Returns
            -------
            float or ndarray
                Scalar float for a single value; otherwise a float array.
            """
            if isinstance(value, (list, tuple)):
                return numpy.asarray(value, dtype=float)
            arr = numpy.asarray(value)
            if arr.shape == () or arr.size == 1:
                return float(arr)
            return arr.astype(float)

        self.t = _as_array_or_scalar(t)
        self.a = _as_array_or_scalar(a)
        self.e = _as_array_or_scalar(e)
        self.inc = _as_array_or_scalar(inc)
        self.Omega = _as_array_or_scalar(Omega)
        self.omega = _as_array_or_scalar(omega)
        self.f = _as_array_or_scalar(f)
        deg2rad = numpy.pi / 180.0
        self._inc_rad = numpy.asarray(self.inc, dtype=float) * deg2rad
        self._Omega_rad = numpy.asarray(self.Omega, dtype=float) * deg2rad
        self._omega_rad = numpy.asarray(self.omega, dtype=float) * deg2rad
        self._f_rad = numpy.asarray(self.f, dtype=float) * deg2rad
        if not isinstance(self.inc, numpy.ndarray):
            self._inc_rad = float(self._inc_rad)
            self._Omega_rad = float(self._Omega_rad)
            self._omega_rad = float(self._omega_rad)
            self._f_rad = float(self._f_rad)

        self.mass = _mass_to_kg(mass)
        self.mu = (_G_SI * self.mass) / 1e9  # km³/s² for SPICE
        self.name_metric = "OrbitalElements"
        self.kwargs = {"mass": mass}

        self.state_vector = numpy.array(
            (self.t, self.a, self.e, self.inc, self.Omega, self.omega, self.f)
        )

    def __getitem__(self, index):
        """
        Return a component of the public state vector.

        Parameters
        ----------
        index : int or slice
            Index into ``state_vector``.

        Returns
        -------
        scalar or ndarray
            Selected element(s).

        Examples
        --------
        >>> from relatipy.numeric.coordinates.orbital_elements import OrbitalElements
        >>> oe = OrbitalElements(0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        >>> float(oe[1])
        10.0
        """
        return self.state_vector[index]

    def __setitem__(self, index, value):
        """
        Update a component of ``state_vector`` in place.

        Parameters
        ----------
        index : int or slice
            Index into ``state_vector``.
        value : scalar or array_like
            New value(s).

        Notes
        -----
        This does not resynchronize internal radian caches or derived fields;
        prefer constructing a new instance if elements change meaningfully.

        Examples
        --------
        >>> import numpy as np
        >>> from relatipy.numeric.coordinates.orbital_elements import OrbitalElements
        >>> oe = OrbitalElements(0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        >>> oe[1] = 20.0
        >>> float(oe[1])
        20.0
        """
        self.state_vector[index] = value

    def _get_elts(self, index=None):
        """
        Build the SPICE osculating-element vector for ``spice.conics``.

        Parameters
        ----------
        index : int, optional
            If batch elements are stored as arrays, index of the orbit to use.
            If omitted, scalar elements are used.

        Returns
        -------
        list
            ``[rp, ecc, inc, lnode, argp, m0, t0, mu]`` with ``rp`` in km,
            angles in rad, ``t0`` in s, and ``mu`` in km³/s².

        Notes
        -----
        Mean anomaly at epoch ``m0`` is obtained from the true anomaly via the
        eccentric anomaly.

        Examples
        --------
        >>> from relatipy.numeric.coordinates.orbital_elements import OrbitalElements
        >>> oe = OrbitalElements(0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, mass=1.0)
        >>> elts = oe._get_elts()
        >>> len(elts)
        8
        """
        if index is not None:
            a = self.a[index]
            e = self.e[index]
            inc = self._inc_rad[index]
            Omega = self._Omega_rad[index]
            omega = self._omega_rad[index]
            f = self._f_rad[index]
            t = self.t[index] if isinstance(self.t, numpy.ndarray) else self.t
            mu = self.mu[index] if isinstance(self.mu, numpy.ndarray) else self.mu
        else:
            a = self.a
            e = self.e
            inc = self._inc_rad
            Omega = self._Omega_rad
            omega = self._omega_rad
            f = self._f_rad
            t = self.t
            mu = self.mu

        length_unit = _G_SI * self.mass / _c_SI**2  # GM/c² en metros
        time_unit = _G_SI * self.mass / _c_SI**3    # GM/c³ en segundos
        rp_km = a * length_unit * (1 - e) / 1000.0  # geometric → km
        one_minus_e = numpy.maximum(1 - e, 1e-15)
        E = 2 * numpy.arctan2(
            numpy.sqrt(one_minus_e) * numpy.sin(f / 2),
            numpy.sqrt(1 + e) * numpy.cos(f / 2)
        )
        M0 = E - e * numpy.sin(E)
        t0_sec = float(t) * time_unit  # geometric → seconds
        return [rp_km, e, inc, Omega, omega, M0, t0_sec, mu]

    def _to_cartesian_arrays(self):
        """
        Convert orbital elements to Cartesian position and velocity arrays.

        Uses :func:`spiceypy.spiceypy.conics` in km, km/s and maps back to
        geometric length and :math:`v/c`.

        Returns
        -------
        r_geom : ndarray, shape (3,) or (N, 3)
            Position in geometric units :math:`GM/c^2`.
        v_dimless : ndarray, shape (3,) or (N, 3)
            Velocity components :math:`v^i/c`.

        Notes
        -----
        If ``a`` is a one-dimensional array with more than one entry, each
        orbit is processed in a loop (batch mode).

        Examples
        --------
        >>> import numpy as np
        >>> from relatipy.numeric.coordinates.orbital_elements import OrbitalElements
        >>> oe = OrbitalElements(0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, mass=1.0)
        >>> r, v = oe._to_cartesian_arrays()
        >>> r.shape
        (3,)
        """
        is_batch = isinstance(self.a, numpy.ndarray) and self.a.ndim > 0 and self.a.size > 1

        length_unit = _G_SI * self.mass / _c_SI**2  # GM/c² en metros
        time_unit = _G_SI * self.mass / _c_SI**3    # GM/c³ en segundos

        if not is_batch:
            elts = self._get_elts()
            t_sec = float(self.t) * time_unit  # geometric → seconds
            state_km = spice.conics(elts, t_sec)  # km, km/s
            r_km = numpy.array(state_km[:3])
            v_km_s = numpy.array(state_km[3:])
            r_m = r_km * 1000.0
            r_geom = r_m / length_unit
            v_m_s = v_km_s * 1000.0
            v_dimless = v_m_s / _c_SI
            return r_geom, v_dimless

        n_orbits = len(self[0])
        r_list = []
        v_list = []
        for i in range(n_orbits):
            elts_i = self._get_elts(index=i)
            t_i = self.t[i] if isinstance(self.t, numpy.ndarray) else self.t
            t_sec_i = float(t_i) * time_unit
            state_km_i = spice.conics(elts_i, t_sec_i)
            r_km_i = numpy.array(state_km_i[:3])
            v_km_s_i = numpy.array(state_km_i[3:])
            r_geom_i = (r_km_i * 1000.0) / length_unit
            v_dimless_i = (v_km_s_i * 1000.0) / _c_SI
            r_list.append(r_geom_i)
            v_list.append(v_dimless_i)

        r_geom = numpy.stack(r_list, axis=0)
        v_dimless = numpy.stack(v_list, axis=0)
        return r_geom, v_dimless

    def convert_to_cartesian(self):
        """
        Convert to Cartesian coordinates (four-position and three-velocity).

        Returns
        -------
        Cartesian
            Instance of :class:`~relatipy.numeric.coordinates.cartesian.Cartesian`
            with ``from_dxs_dt=False``. Orbital attributes are copied for
            downstream checks.

        Examples
        --------
        >>> from relatipy.numeric.coordinates.orbital_elements import OrbitalElements
        >>> oe = OrbitalElements(0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, mass=1.0)
        >>> cart = oe.convert_to_cartesian()
        >>> cart.name_metric
        'Cartesian'
        """
        from .cartesian import Cartesian

        r_vec, v_vec = self._to_cartesian_arrays()
        if r_vec.ndim == 1:
            xs = numpy.array([self.t, r_vec[0], r_vec[1], r_vec[2]])
            vs = numpy.array([v_vec[0], v_vec[1], v_vec[2]])
        else:
            t_array = (
                self.t if isinstance(self.t, numpy.ndarray) else numpy.asarray([self.t] * r_vec.shape[0])
            )
            xs = numpy.array(
                [t_array, r_vec[:, 0], r_vec[:, 1], r_vec[:, 2]]
            )
            vs = numpy.array(
                [v_vec[:, 0], v_vec[:, 1], v_vec[:, 2]]
            )

        cartesian = Cartesian(xs, vels=vs, from_dxs_dt=False)

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
        Convert to another registered coordinate system via Cartesian.

        Parameters
        ----------
        target_system : str
            Key in ``relatipy.numeric.coordinates.coordinate_systems``.
        **kwargs
            Passed to the target class factory (e.g. ``mass`` required when
            converting to ``\"OrbitalElements\"``).

        Returns
        -------
        CoordinateBase subclass
            Instance of the requested system.

        Raises
        ------
        ValueError
            If ``target_system`` is not registered, or if converting to
            ``\"OrbitalElements\"`` without ``mass`` in ``kwargs``.

        Examples
        --------
        >>> from relatipy.numeric.coordinates.orbital_elements import OrbitalElements
        >>> oe = OrbitalElements(0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, mass=1.0)
        >>> same = oe.convert_to("OrbitalElements", mass=1.0)
        >>> same.name_metric
        'OrbitalElements'
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
        Build :class:`OrbitalElements` from a Cartesian state using ``oscelt``.

        Dispatches to vector or scalar helpers depending on the shape of
        ``xs``.

        Parameters
        ----------
        xs : array_like
            Four-position; either shape ``(4,)`` or batch ``(4, N)``.
        vs : array_like
            Velocity ``v/c``, shape ``(3,)`` or ``(3, N)`` matching ``xs``.
        mass : float or astropy.units.Quantity
            Central mass (geometric or quantity); see :class:`OrbitalElements`.
        t : float or astropy.units.Quantity, optional
            Override time in seconds if given; otherwise taken from the geometric
            time in ``xs`` after unit conversion.

        Returns
        -------
        OrbitalElements
            Fitted osculating elements.

        Examples
        --------
        >>> import numpy as np
        >>> from relatipy.numeric.coordinates.orbital_elements import OrbitalElements
        >>> xs = np.array([0.0, 1e6, 0.0, 0.0], dtype=float)
        >>> vs = np.array([0.0, 0.0, 1e-4], dtype=float)
        >>> oe = OrbitalElements.from_cartesian(xs, vs, mass=1.0)
        >>> oe.name_metric
        'OrbitalElements'
        """
        if isinstance(xs[1], numpy.ndarray):
            return OrbitalElements.from_cartesian_vector(xs, vs, mass, t)
        return OrbitalElements.from_cartesian_scalar(xs, vs, mass, t)

    @staticmethod
    def from_cartesian_vector(xs, vs, mass, t=None):
        """
        Build :class:`OrbitalElements` for ``N`` Cartesian states (columns).

        Parameters
        ----------
        xs : array_like, shape (4, N)
            Batch of four-positions in geometric units.
        vs : array_like, shape (3, N)
            Batch of velocities :math:`v/c`.
        mass : float or astropy.units.Quantity
            Central mass, shared by all columns.
        t : float or astropy.units.Quantity, optional
            Optional per-call time override in seconds (see
            :meth:`from_cartesian_scalar`).

        Returns
        -------
        OrbitalElements
            Elements with array fields of length ``N``.

        Examples
        --------
        >>> import numpy as np
        >>> from relatipy.numeric.coordinates.orbital_elements import OrbitalElements
        >>> xs = np.array([[0.0, 0.0], [1e6, 2e6], [0.0, 0.0], [0.0, 0.0]], dtype=float)
        >>> vs = np.array([[0.0, 0.0], [0.0, 0.0], [1e-4, 1e-4]], dtype=float)
        >>> oe = OrbitalElements.from_cartesian_vector(xs, vs, mass=1.0)
        >>> np.asarray(oe.a).size
        2
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
        Build :class:`OrbitalElements` from one Cartesian state using ``oscelt``.

        Parameters
        ----------
        xs : array_like, shape (4,)
            Four-position ``(t, x, y, z)`` with ``t`` in geometric time and
            spatial part in geometric length.
        vs : array_like, shape (3,)
            Cartesian velocity :math:`v^i/c`.
        mass : float or astropy.units.Quantity
            Central mass.
        t : float or astropy.units.Quantity, optional
            If set, evaluation epoch in seconds (bypassing geometric time from
            ``xs`` after internal conversion path).

        Returns
        -------
        OrbitalElements
            Osculating elements; angles in degrees, ``t`` and ``a`` geometric.

        Notes
        -----
        True anomaly is recovered from mean anomaly with a fixed-point iteration
        on the eccentric anomaly.

        Examples
        --------
        >>> import numpy as np
        >>> from relatipy.numeric.coordinates.orbital_elements import OrbitalElements
        >>> xs = np.array([0.0, 1e6, 0.0, 0.0], dtype=float)
        >>> vs = np.array([0.0, 0.0, 1e-4], dtype=float)
        >>> oe = OrbitalElements.from_cartesian_scalar(xs, vs, mass=1.0)
        >>> hasattr(oe, "e")
        True
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

        length_unit = _G_SI * mass_kg / _c_SI**2
        time_unit = _G_SI * mass_kg / _c_SI**3
        a_geom = a_m / length_unit

        E = M0
        for _ in range(50):
            E = E - (E - e * numpy.sin(E) - M0) / (1 - e * numpy.cos(E))

        one_minus_e = numpy.maximum(1 - e, 1e-15)
        f = 2 * numpy.arctan2(
            numpy.sqrt(1 + e) * numpy.sin(E / 2),
            numpy.sqrt(one_minus_e) * numpy.cos(E / 2)
        )

        rad2deg = 180.0 / numpy.pi
        t_geom = t_sec / time_unit
        return OrbitalElements(
            t_geom,
            a_geom,
            e,
            float(inc) * rad2deg,
            float(Omega) * rad2deg,
            float(omega) * rad2deg,
            float(f) * rad2deg,
            mass,
        )

    @staticmethod
    def _convert_from_cartesian(xs_p, vs_p, **kwargs):
        """
        Bridge for ``CoordinateBase.convert_to`` (same as :meth:`from_cartesian`).

        Parameters
        ----------
        xs_p : array_like
            Four-position Cartesian state.
        vs_p : array_like
            Three-velocity :math:`v/c`.
        **kwargs
            Must include ``mass``.

        Returns
        -------
        OrbitalElements
            Elements from :meth:`from_cartesian`.

        Raises
        ------
        ValueError
            If ``mass`` is missing from ``kwargs``.

        Examples
        --------
        >>> import numpy as np
        >>> from relatipy.numeric.coordinates.orbital_elements import OrbitalElements
        >>> xs = np.array([0.0, 1e6, 0.0, 0.0], dtype=float)
        >>> vs = np.array([0.0, 0.0, 1e-4], dtype=float)
        >>> oe = OrbitalElements._convert_from_cartesian(xs, vs, mass=1.0)
        >>> oe.name_metric
        'OrbitalElements'
        """
        mass = kwargs.get("mass")
        if mass is None:
            raise ValueError("mass is required to convert to OrbitalElements.")
        return OrbitalElements.from_cartesian(xs_p, vs_p, mass)

    def _get_period(self):
        """
        Orbital period in geometric time units :math:`GM/c^3`.

        Uses Kepler's third law with semi-major axis in meters derived from
        geometric ``a`` and the central mass.

        Returns
        -------
        float or ndarray
            Period in geometric time; shape follows broadcast of ``self.a`` if
            array-like.

        Examples
        --------
        >>> from relatipy.numeric.coordinates.orbital_elements import OrbitalElements
        >>> oe = OrbitalElements(0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0, mass=1.0)
        >>> T = oe._get_period()
        >>> bool(T > 0)
        True
        """
        L = _G_SI * self.mass / _c_SI**2
        T = _G_SI * self.mass / _c_SI**3

        a_m = self.a * L
        mu_si = _G_SI * self.mass

        T_sec = 2 * numpy.pi * numpy.sqrt(a_m**3 / mu_si)
        return T_sec / T
