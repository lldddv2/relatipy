"""
Sky-frame Keplerian orbital elements with explicit BH spin orientation.

This module provides :class:`ApparentOrbitalElements`, a subclass of
:class:`~relatipy.numeric.coordinates.orbital_elements.OrbitalElements` that
accepts orbital elements expressed in the **sky frame** — the plane perpendicular
to the observer's line of sight — as reported by astrometric catalogs (e.g. for
S2 around Sgr A*).

Sky-frame convention
--------------------
:math:`\\hat{z}_{\\mathrm{sky}}` points toward the observer (line of sight),
:math:`\\hat{x}_{\\mathrm{sky}}` points toward celestial north (Dec axis), and
:math:`\\hat{y}_{\\mathrm{sky}}` points east (RA axis). The orbital elements
:math:`(i, \\Omega, \\omega)` are measured with respect to this plane.

BH spin orientation
-------------------
Two angles orient the black hole spin axis in the sky frame:

- :math:`\\zeta \\in [0^\\circ, 180^\\circ]` — polar angle of the spin axis measured
  from :math:`\\hat{z}_{\\mathrm{sky}}` (line of sight).  :math:`\\zeta = 0^\\circ`
  means the spin points directly at the observer.
- :math:`\\eta \\in [0^\\circ, 360^\\circ]` — position angle of the spin projection
  onto the sky plane, measured counter-clockwise from celestial north
  (:math:`\\hat{x}_{\\mathrm{sky}}`) toward east (:math:`\\hat{y}_{\\mathrm{sky}}`).

Rotation geometry
-----------------
The unit spin vector expressed in the sky basis is

.. math::

   \\hat{s}_{\\mathrm{sky}} = (\\sin\\zeta\\cos\\eta,\\; \\sin\\zeta\\sin\\eta,\\; \\cos\\zeta).

The rotation matrix that maps sky-frame vectors to the BH frame
(sending :math:`\\hat{s}_{\\mathrm{sky}}` to :math:`\\hat{z}_{\\mathrm{BH}}`) is

.. math::

   R = R_y(-\\zeta)\\, R_z(-\\eta) =
   \\begin{pmatrix}
     \\cos\\zeta\\cos\\eta & \\cos\\zeta\\sin\\eta & -\\sin\\zeta \\\\
     -\\sin\\eta           & \\cos\\eta             & 0            \\\\
     \\sin\\zeta\\cos\\eta & \\sin\\zeta\\sin\\eta & \\cos\\zeta
   \\end{pmatrix}.

Kerr spacetime is axially symmetric so there is no unique orientation for the
BH-frame :math:`x`-axis; the choice above (no extra twist around the spin) is the
natural default.

See Also
--------
relatipy.numeric.coordinates.orbital_elements.OrbitalElements :
    Parent class; assumes elements already expressed in the BH frame.

Examples
--------
>>> import numpy as np
>>> from relatipy.numeric.coordinates.apparent_orbital_elements import ApparentOrbitalElements
>>> aoe = ApparentOrbitalElements(0.0, 1e7, 0.1, 30.0, 40.0, 50.0, 60.0,
...                               zeta=45.0, eta=90.0, mass=1.0)
>>> aoe.state_vector.shape
(7,)
>>> aoe.zeta, aoe.eta
(45.0, 90.0)
"""

import numpy
from .orbital_elements import OrbitalElements


def _build_rotation(zeta_rad, eta_rad):
    """
    Build the 3×3 rotation matrix mapping the sky frame to the BH frame.

    Parameters
    ----------
    zeta_rad : float
        Polar angle of the spin axis from the line of sight, in radians.
    eta_rad : float
        Position angle of the spin projection on the sky plane, in radians.

    Returns
    -------
    ndarray, shape (3, 3)
        Orthogonal matrix :math:`R = R_y(-\\zeta)\\,R_z(-\\eta)`.

    Examples
    --------
    >>> import numpy as np
    >>> from relatipy.numeric.coordinates.apparent_orbital_elements import _build_rotation
    >>> R = _build_rotation(0.0, 0.0)
    >>> np.allclose(R, np.eye(3))
    True
    """
    cz, sz = numpy.cos(zeta_rad), numpy.sin(zeta_rad)
    cn, sn = numpy.cos(eta_rad), numpy.sin(eta_rad)
    return numpy.array([
        [ cz * cn,  cz * sn, -sz],
        [-sn,       cn,       0.],
        [ sz * cn,  sz * sn,  cz],
    ])


class ApparentOrbitalElements(OrbitalElements):
    """
    Sky-frame Keplerian orbital elements with explicit BH spin orientation.

    Extends :class:`OrbitalElements` so that the angles :math:`(i, \\Omega, \\omega)`
    are interpreted with respect to the **sky plane** (the plane perpendicular to
    the observer's line of sight) rather than the black hole equatorial plane.

    Two additional scalar parameters — :math:`\\zeta` and :math:`\\eta` — specify
    how the BH spin axis is oriented in the sky frame (see module docstring).
    Internally, Cartesian positions and velocities are converted to the BH frame
    via the rotation :math:`R = R_y(-\\zeta)\\,R_z(-\\eta)` before any general-relativistic
    computation.

    Parameters
    ----------
    t : float or array_like, optional
        Coordinate time in geometric units :math:`GM/c^3`. Default ``0.0``.
    a : float or array_like, optional
        Semi-major axis in geometric length units :math:`GM/c^2`. Default ``1.0``.
    e : float or array_like, optional
        Eccentricity. Default ``0.0``.
    inc : float or array_like, optional
        Inclination with respect to the sky plane, in degrees. Default ``0.0``.
    Omega : float or array_like, optional
        Longitude of the ascending node measured in the sky plane, in degrees.
        Default ``0.0``.
    omega : float or array_like, optional
        Argument of periapsis measured from the sky-plane node, in degrees.
        Default ``0.0``.
    f : float or array_like, optional
        True anomaly in degrees. Default ``0.0``.
    zeta : float, optional
        Polar angle of the BH spin axis from the line of sight, in degrees.
        Range :math:`[0^\\circ, 180^\\circ]`. Default ``0.0`` (identity sky/BH map).
    eta : float, optional
        Position angle of the spin projection on the sky plane, CCW from
        celestial north toward east, in degrees. Range :math:`[0^\\circ, 360^\\circ]`.
        Default ``0.0``.
    mass : float or astropy.units.Quantity, optional
        Central mass (geometric scalar or astropy quantity). Default ``1``.

    Attributes
    ----------
    zeta : float
        Spin polar angle from line of sight, in degrees (scalar).
    eta : float
        Spin position angle on sky plane, in degrees (scalar).
    name_metric : str
        ``"ApparentOrbitalElements"``.

    Notes
    -----
    The inherited ``state_vector`` remains a 7-element array
    ``(t, a, e, inc, Omega, omega, f)`` describing the orbit in the sky frame,
    consistent with the public API of :class:`OrbitalElements`.  The orientation
    parameters :math:`(\\zeta, \\eta)` are stored as ordinary attributes (like
    ``mass``) and are included in ``kwargs``.

    When :math:`\\zeta = \\eta = 0`, the rotation :math:`R` is the identity and
    results are identical to :class:`OrbitalElements`.

    Examples
    --------
    >>> import numpy as np
    >>> from relatipy.numeric.coordinates.orbital_elements import OrbitalElements
    >>> from relatipy.numeric.coordinates.apparent_orbital_elements import ApparentOrbitalElements
    >>> oe  = OrbitalElements(0.0, 1e7, 0.1, 30.0, 40.0, 50.0, 60.0, mass=1.0)
    >>> aoe = ApparentOrbitalElements(0.0, 1e7, 0.1, 30.0, 40.0, 50.0, 60.0,
    ...                               zeta=0.0, eta=0.0, mass=1.0)
    >>> np.allclose(oe.convert_to_cartesian().xs, aoe.convert_to_cartesian().xs)
    True
    """

    def __init__(self, t=0.0, a=1.0, e=0.0, inc=0.0, Omega=0.0, omega=0.0, f=0.0, zeta=0.0, eta=0.0, mass=1):
        super().__init__(t, a, e, inc, Omega, omega, f, mass=mass)
        self.zeta = float(zeta)
        self.eta = float(eta)
        self._zeta_rad = self.zeta * numpy.pi / 180.0
        self._eta_rad = self.eta * numpy.pi / 180.0
        self.name_metric = "ApparentOrbitalElements"
        self.kwargs = {"mass": mass, "zeta": self.zeta, "eta": self.eta}

    def _rotation_sky_to_bh(self):
        """
        Rotation matrix that maps sky-frame vectors to the BH frame.

        Returns
        -------
        ndarray, shape (3, 3)
            :math:`R = R_y(-\\zeta)\\,R_z(-\\eta)`.

        Examples
        --------
        >>> from relatipy.numeric.coordinates.apparent_orbital_elements import ApparentOrbitalElements
        >>> aoe = ApparentOrbitalElements(0.0, 1e7, 0.0, 0.0, 0.0, 0.0, 0.0,
        ...                              zeta=0.0, eta=0.0, mass=1.0)
        >>> import numpy as np
        >>> np.allclose(aoe._rotation_sky_to_bh(), np.eye(3))
        True
        """
        return _build_rotation(self._zeta_rad, self._eta_rad)

    def _to_cartesian_arrays(self):
        """
        Convert sky-frame elements to BH-frame Cartesian position and velocity.

        Calls the parent implementation to obtain sky-frame :math:`(r, v)`, then
        applies the rotation :math:`R` to map them to the BH frame.

        Returns
        -------
        r_bh : ndarray, shape (3,) or (N, 3)
            Position in the BH frame in geometric units :math:`GM/c^2`.
        v_bh : ndarray, shape (3,) or (N, 3)
            Velocity in the BH frame, dimensionless :math:`v/c`.

        Examples
        --------
        >>> import numpy as np
        >>> from relatipy.numeric.coordinates.apparent_orbital_elements import ApparentOrbitalElements
        >>> aoe = ApparentOrbitalElements(0.0, 1e7, 0.0, 0.0, 0.0, 0.0, 0.0,
        ...                              zeta=0.0, eta=0.0, mass=1.0)
        >>> r, v = aoe._to_cartesian_arrays()
        >>> r.shape
        (3,)
        """
        r_sky, v_sky = super()._to_cartesian_arrays()
        R = self._rotation_sky_to_bh()
        if r_sky.ndim == 1:
            return R @ r_sky, R @ v_sky
        return r_sky @ R.T, v_sky @ R.T

    def convert_to_cartesian(self):
        """
        Convert to BH-frame Cartesian coordinates.

        Returns
        -------
        Cartesian
            Instance of :class:`~relatipy.numeric.coordinates.cartesian.Cartesian`
            with positions and velocities expressed in the BH frame. The
            ``zeta`` and ``eta`` attributes are copied to the returned object.

        Examples
        --------
        >>> from relatipy.numeric.coordinates.apparent_orbital_elements import ApparentOrbitalElements
        >>> aoe = ApparentOrbitalElements(0.0, 1e7, 0.0, 0.0, 0.0, 0.0, 0.0,
        ...                              zeta=45.0, eta=30.0, mass=1.0)
        >>> cart = aoe.convert_to_cartesian()
        >>> cart.name_metric
        'Cartesian'
        >>> cart.zeta
        45.0
        """
        cart = super().convert_to_cartesian()
        cart.zeta = self.zeta
        cart.eta = self.eta
        return cart

    @staticmethod
    def from_cartesian(xs, vs, mass, zeta, eta, t=None):
        """
        Build :class:`ApparentOrbitalElements` from a BH-frame Cartesian state.

        Applies the inverse rotation :math:`R^\\top` to map state vectors from
        the BH frame to the sky frame, then uses SPICE ``oscelt`` to fit the
        osculating elements in the sky frame.

        Parameters
        ----------
        xs : array_like
            Four-position; shape ``(4,)`` (scalar) or ``(4, N)`` (batch).
        vs : array_like
            Velocity :math:`v/c`; shape ``(3,)`` or ``(3, N)`` matching ``xs``.
        mass : float or astropy.units.Quantity
            Central mass.
        zeta : float
            Polar angle of the BH spin from line of sight, in degrees.
        eta : float
            Position angle of the spin projection on the sky plane, in degrees.
        t : float, optional
            Time override in seconds passed to ``oscelt``.

        Returns
        -------
        ApparentOrbitalElements
            Sky-frame osculating elements with the given orientation angles.

        Examples
        --------
        >>> import numpy as np
        >>> from relatipy.numeric.coordinates.apparent_orbital_elements import ApparentOrbitalElements
        >>> aoe = ApparentOrbitalElements(0.0, 1e7, 0.1, 30.0, 40.0, 50.0, 60.0,
        ...                              zeta=45.0, eta=30.0, mass=1.0)
        >>> cart = aoe.convert_to_cartesian()
        >>> aoe2 = ApparentOrbitalElements.from_cartesian(
        ...     cart.xs, cart.vs, mass=1.0, zeta=45.0, eta=30.0)
        >>> np.allclose(aoe.state_vector, aoe2.state_vector, atol=1e-6)
        True
        """
        if isinstance(xs[1], numpy.ndarray):
            return ApparentOrbitalElements.from_cartesian_vector(xs, vs, mass, zeta, eta, t)
        return ApparentOrbitalElements.from_cartesian_scalar(xs, vs, mass, zeta, eta, t)

    @staticmethod
    def from_cartesian_scalar(xs, vs, mass, zeta, eta, t=None):
        """
        Build :class:`ApparentOrbitalElements` from one BH-frame Cartesian state.

        Parameters
        ----------
        xs : array_like, shape (4,)
            BH-frame four-position ``(t, x, y, z)`` in geometric units.
        vs : array_like, shape (3,)
            BH-frame velocity :math:`v^i/c`.
        mass : float or astropy.units.Quantity
            Central mass.
        zeta : float
            Polar angle of BH spin from line of sight, in degrees.
        eta : float
            Position angle of spin projection on sky plane, in degrees.
        t : float, optional
            Time override in seconds.

        Returns
        -------
        ApparentOrbitalElements
            Sky-frame osculating elements.

        Examples
        --------
        >>> import numpy as np
        >>> from relatipy.numeric.coordinates.apparent_orbital_elements import ApparentOrbitalElements
        >>> xs = np.array([0.0, 1e7, 0.0, 0.0])
        >>> vs = np.array([0.0, 0.0, 1e-5])
        >>> aoe = ApparentOrbitalElements.from_cartesian_scalar(xs, vs, 1.0, 0.0, 0.0)
        >>> hasattr(aoe, "zeta")
        True
        """
        R = _build_rotation(
            float(zeta) * numpy.pi / 180.0,
            float(eta) * numpy.pi / 180.0,
        )
        xs = numpy.asarray(xs, dtype=float)
        vs = numpy.asarray(vs, dtype=float)
        r_sky = R.T @ xs[1:4]
        v_sky = R.T @ vs
        xs_sky = numpy.array([xs[0], r_sky[0], r_sky[1], r_sky[2]])
        oe = OrbitalElements.from_cartesian_scalar(xs_sky, v_sky, mass, t)
        return ApparentOrbitalElements(
            oe.t, oe.a, oe.e, oe.inc, oe.Omega, oe.omega, oe.f,
            float(zeta), float(eta), mass,
        )

    @staticmethod
    def from_cartesian_vector(xs, vs, mass, zeta, eta, t=None):
        """
        Build :class:`ApparentOrbitalElements` for ``N`` BH-frame Cartesian states.

        Parameters
        ----------
        xs : array_like, shape (4, N)
            Batch of BH-frame four-positions in geometric units.
        vs : array_like, shape (3, N)
            Batch of BH-frame velocities :math:`v/c`.
        mass : float or astropy.units.Quantity
            Central mass, shared by all columns.
        zeta : float
            Polar angle of BH spin from line of sight, in degrees.
        eta : float
            Position angle of spin projection on sky plane, in degrees.
        t : float, optional
            Optional time override in seconds.

        Returns
        -------
        ApparentOrbitalElements
            Sky-frame elements with array fields of length ``N``.

        Examples
        --------
        >>> import numpy as np
        >>> from relatipy.numeric.coordinates.apparent_orbital_elements import ApparentOrbitalElements
        >>> xs = np.array([[0.0, 0.0], [1e7, 2e7], [0.0, 0.0], [0.0, 0.0]])
        >>> vs = np.array([[0.0, 0.0], [0.0, 0.0], [1e-5, 1e-5]])
        >>> aoe = ApparentOrbitalElements.from_cartesian_vector(xs, vs, 1.0, 0.0, 0.0)
        >>> np.asarray(aoe.a).size
        2
        """
        R = _build_rotation(
            float(zeta) * numpy.pi / 180.0,
            float(eta) * numpy.pi / 180.0,
        )
        xs = numpy.asarray(xs, dtype=float)
        vs = numpy.asarray(vs, dtype=float)
        r_sky = R.T @ xs[1:4]
        v_sky = R.T @ vs
        xs_sky = numpy.vstack([xs[0:1, :], r_sky])
        oe = OrbitalElements.from_cartesian_vector(xs_sky, v_sky, mass, t)
        return ApparentOrbitalElements(
            oe.t, oe.a, oe.e, oe.inc, oe.Omega, oe.omega, oe.f,
            float(zeta), float(eta), mass,
        )

    @staticmethod
    def _convert_from_cartesian(xs_p, vs_p, **kwargs):
        """
        Bridge for ``convert_to("ApparentOrbitalElements", ...)``.

        Parameters
        ----------
        xs_p : array_like
            BH-frame four-position.
        vs_p : array_like
            BH-frame three-velocity :math:`v/c`.
        **kwargs
            Must include ``mass``, ``zeta``, and ``eta``.

        Returns
        -------
        ApparentOrbitalElements
            Sky-frame elements from :meth:`from_cartesian`.

        Raises
        ------
        ValueError
            If ``mass``, ``zeta``, or ``eta`` are missing from ``kwargs``.

        Examples
        --------
        >>> import numpy as np
        >>> from relatipy.numeric.coordinates.apparent_orbital_elements import ApparentOrbitalElements
        >>> xs = np.array([0.0, 1e7, 0.0, 0.0])
        >>> vs = np.array([0.0, 0.0, 1e-5])
        >>> aoe = ApparentOrbitalElements._convert_from_cartesian(
        ...     xs, vs, mass=1.0, zeta=0.0, eta=0.0)
        >>> aoe.name_metric
        'ApparentOrbitalElements'
        """
        mass = kwargs.get("mass")
        zeta = kwargs.get("zeta")
        eta = kwargs.get("eta")
        if mass is None:
            raise ValueError("mass is required to convert to ApparentOrbitalElements.")
        if zeta is None or eta is None:
            raise ValueError("zeta and eta are required to convert to ApparentOrbitalElements.")
        return ApparentOrbitalElements.from_cartesian(xs_p, vs_p, mass, zeta, eta)
