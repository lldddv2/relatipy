"""
Equatorial sky coordinates (RA/Dec) from BH-frame position and distance.

Given a position in the **BH frame** together with the BH spin orientation
angles :math:`(\\zeta, \\eta)` and the physical distance to the observer,
this module rotates the position to the sky frame and computes the observed
right ascension and declination offsets.

Sky-frame convention
--------------------
:math:`\\hat{x}_{\\mathrm{sky}}` points north (Dec axis),
:math:`\\hat{y}_{\\mathrm{sky}}` points east (RA axis), and
:math:`\\hat{z}_{\\mathrm{sky}}` points toward the observer (line of sight).
This matches the convention used by
:class:`~relatipy.numeric.coordinates.apparent_orbital_elements.ApparentOrbitalElements`.

BH-to-sky rotation
-------------------
The rotation :math:`R = R_y(-\\zeta)\\,R_z(-\\eta)` maps sky-frame vectors to
the BH frame.  Its transpose :math:`R^\\top` reverses the mapping:

.. math::

    \\mathbf{r}_{\\mathrm{sky}} = R^\\top \\mathbf{r}_{\\mathrm{BH}}.

The north and east sky-plane components are then

.. math::

    x_N = (R^\\top \\mathbf{r}_{\\mathrm{BH}})_0, \\qquad
    y_E = (R^\\top \\mathbf{r}_{\\mathrm{BH}})_1.

Coordinate computation
----------------------
For sky-plane offsets :math:`(x_N, y_E)` and observer distance :math:`d`,
the angular offsets follow from exact spherical projection:

.. math::

    \\alpha = \\arctan\\!\\left(\\frac{y_E}{d}\\right), \\qquad
    \\delta = \\arctan\\!\\left(\\frac{x_N}{\\sqrt{d^2 + y_E^2}}\\right).

Both are returned in **decimal degrees**.

Construction paths
------------------
1. **Manual** — provide the BH-frame position ``(x, y, z)`` together with
   the spin angles ``zeta`` and ``eta`` (degrees) and ``distance``.
2. **From** :class:`~relatipy.numeric.coordinates.apparent_orbital_elements.ApparentOrbitalElements`
   — call :meth:`EquatorialCoordinate.from_apparent_orbital_elements`;
   only ``distance`` is required because ``zeta`` and ``eta`` are already
   stored in the object.

Examples
--------
>>> import astropy.units as u
>>> from relatipy.numeric.coordinates.apparent_equatorial_coordinate import EquatorialCoordinate
>>> ec = EquatorialCoordinate(0.0, 0.0, 0.0, zeta=0.0, eta=0.0,
...                           distance=8.5 * u.kpc, mass=4e6)
>>> ec.ra, ec.dec
(0.0, 0.0)
"""

import numpy
from astropy import units as u
from astropy.coordinates import Angle

from ..constants import _L_ref
from .apparent_orbital_elements import _build_rotation


def _to_physical_meters(value, mass, label):
    """
    Convert a BH-frame position component to meters.

    Parameters
    ----------
    value : float or astropy.units.Quantity
        Position component in geometric units :math:`GM_{\\odot}/c^2` (float)
        or in explicit length units (Quantity).
    mass : float or astropy.units.Quantity or None
        Central mass in solar masses (float) or with mass units (Quantity).
        Required when ``value`` is a plain float.
    label : str
        Parameter name used in error messages.

    Returns
    -------
    float
        Physical length in meters.

    Raises
    ------
    TypeError
        If ``value`` is a Quantity with incompatible units.
    ValueError
        If ``value`` is a float and ``mass`` is None, or if ``mass`` has
        incompatible units.

    Examples
    --------
    >>> import astropy.units as u
    >>> from relatipy.numeric.coordinates.apparent_equatorial_coordinate import _to_physical_meters
    >>> _to_physical_meters(1.0 * u.au, None, "x") > 0
    True
    """
    if isinstance(value, u.Quantity):
        try:
            return float(value.to(u.m).value)
        except u.UnitConversionError:
            raise TypeError(
                f"'{label}' Quantity must have length units; got '{value.unit}'."
            )

    if mass is None:
        raise ValueError(
            f"'{label}' is a plain float (geometric units). "
            "Provide 'mass' so it can be converted to physical meters."
        )

    if isinstance(mass, u.Quantity):
        try:
            mass_solar = float(mass.to(u.M_sun).value)
        except u.UnitConversionError:
            raise ValueError(
                f"'mass' Quantity must have mass units; got '{mass.unit}'."
            )
    else:
        mass_solar = float(mass)

    return float(value) * mass_solar * _L_ref


class EquatorialCoordinate:
    """
    Equatorial sky coordinates (RA/Dec) from BH-frame position and distance.

    Rotates a BH-frame Cartesian position to the sky frame using the spin
    orientation angles :math:`(\\zeta, \\eta)`, then projects the sky-plane
    offsets onto angular coordinates via exact spherical trigonometry.

    .. math::

        \\mathbf{r}_{\\mathrm{sky}} = R^\\top(\\zeta, \\eta)\\,\\mathbf{r}_{\\mathrm{BH}},

    .. math::

        \\alpha = \\arctan\\!\\left(\\frac{y_E}{d}\\right), \\qquad
        \\delta = \\arctan\\!\\left(\\frac{x_N}{\\sqrt{d^2 + y_E^2}}\\right).

    Parameters
    ----------
    x : float or astropy.units.Quantity
        BH-frame :math:`x`-component of the position. A plain float is
        interpreted in geometric units :math:`GM_{\\odot}/c^2`; a Quantity
        must carry length units.
    y : float or astropy.units.Quantity
        BH-frame :math:`y`-component. Same unit convention as ``x``.
    z : float or astropy.units.Quantity
        BH-frame :math:`z`-component. Same unit convention as ``x``.
    zeta : float
        Polar angle of the BH spin axis from the line of sight, in degrees.
        Range :math:`[0^\\circ, 180^\\circ]`.
    eta : float
        Position angle of the spin projection on the sky plane, CCW from
        celestial north toward east, in degrees.
        Range :math:`[0^\\circ, 360^\\circ]`.
    distance : astropy.units.Quantity
        Physical distance to the central object. **Mandatory** and must carry
        length units (e.g. ``u.kpc``, ``u.pc``, ``u.lyr``).
    mass : float or astropy.units.Quantity, optional
        Central mass used to convert geometric-unit floats to physical length.
        A plain float is taken as a multiple of :math:`M_{\\odot}`. Required
        when any of ``x``, ``y``, ``z`` is a plain float. Default ``None``.

    Attributes
    ----------
    ra : float
        Right ascension offset in decimal degrees. Positive east.
    dec : float
        Declination offset in decimal degrees. Positive north.
    zeta : float
        Spin polar angle from the line of sight, in degrees.
    eta : float
        Spin position angle on the sky plane, in degrees.
    distance : astropy.units.Quantity
        Physical distance passed at construction.
    name_metric : str
        ``"EquatorialCoordinate"``.

    Raises
    ------
    TypeError
        If ``distance`` is not an :class:`astropy.units.Quantity`.
    ValueError
        If ``distance`` does not have length units, or if a plain float
        position component is given without ``mass``.

    Examples
    --------
    On-axis spin (:math:`\\zeta = 0`) — rotation is the identity, so BH frame
    and sky frame coincide:

    >>> import astropy.units as u
    >>> from relatipy.numeric.coordinates.apparent_equatorial_coordinate import EquatorialCoordinate
    >>> ec = EquatorialCoordinate(0.0, 1000.0, 0.0, zeta=0.0, eta=0.0,
    ...                           distance=8.5 * u.kpc, mass=4.0e6)
    >>> ec.ra != 0.0   # y_BH → y_sky (east) → non-zero RA offset
    True
    >>> ec.dec == 0.0  # x_BH = 0 → zero Dec offset
    True
    >>> ec.name_metric
    'EquatorialCoordinate'

    From an :class:`~relatipy.numeric.coordinates.apparent_orbital_elements.ApparentOrbitalElements`
    instance (preferred path — ``zeta`` and ``eta`` are inferred automatically):

    >>> from relatipy.numeric.coordinates.apparent_orbital_elements import ApparentOrbitalElements
    >>> aoe = ApparentOrbitalElements(0.0, 1e7, 0.1, 30.0, 40.0, 50.0, 60.0,
    ...                               zeta=45.0, eta=90.0, mass=4.0e6)
    >>> ec2 = EquatorialCoordinate.from_apparent_orbital_elements(aoe, 8.5 * u.kpc)
    >>> isinstance(ec2, EquatorialCoordinate)
    True
    """

    def __init__(self, x, y, z, zeta, eta, distance, mass=None):
        if not isinstance(distance, u.Quantity):
            raise TypeError(
                "distance must be an astropy Quantity with length units "
                "(e.g. 8.5 * u.kpc)."
            )
        try:
            d_m = float(distance.to(u.m).value)
        except u.UnitConversionError:
            raise ValueError(
                f"distance must have length units (e.g. kpc, pc, lyr); "
                f"got '{distance.unit}'."
            )

        x_m = _to_physical_meters(x, mass, "x")
        y_m = _to_physical_meters(y, mass, "y")
        z_m = _to_physical_meters(z, mass, "z")

        zeta_rad = float(zeta) * numpy.pi / 180.0
        eta_rad = float(eta) * numpy.pi / 180.0
        R = _build_rotation(zeta_rad, eta_rad)

        r_bh = numpy.array([x_m, y_m, z_m])
        r_sky = R.T @ r_bh

        x_north_m = r_sky[0]
        y_east_m = r_sky[1]

        self._ra_rad = numpy.arctan2(y_east_m, d_m)
        self._dec_rad = numpy.arctan2(x_north_m, numpy.sqrt(d_m**2 + y_east_m**2))

        self.ra = float(numpy.degrees(self._ra_rad))
        self.dec = float(numpy.degrees(self._dec_rad))
        self.zeta = float(zeta)
        self.eta = float(eta)
        self.distance = distance
        self.name_metric = "EquatorialCoordinate"

    # ------------------------------------------------------------------
    # Factory
    # ------------------------------------------------------------------

    @classmethod
    def from_apparent_orbital_elements(cls, aoe, distance):
        """
        Build an :class:`EquatorialCoordinate` from an
        :class:`~relatipy.numeric.coordinates.apparent_orbital_elements.ApparentOrbitalElements`
        instance.

        Extracts the BH-frame Cartesian position and the spin orientation
        angles :math:`(\\zeta, \\eta)` directly from ``aoe``, so only the
        observer distance needs to be supplied.

        Parameters
        ----------
        aoe : ApparentOrbitalElements
            Sky-frame orbital elements with stored ``zeta``, ``eta``, and
            ``mass`` attributes.
        distance : astropy.units.Quantity
            Physical distance to the central object. Must carry length units.

        Returns
        -------
        EquatorialCoordinate
            Instance with ``ra`` and ``dec`` in decimal degrees.

        Raises
        ------
        TypeError
            If ``distance`` is not an :class:`astropy.units.Quantity`.
        ValueError
            If ``distance`` does not have length units.

        Examples
        --------
        >>> import astropy.units as u
        >>> from relatipy.numeric.coordinates.apparent_orbital_elements import ApparentOrbitalElements
        >>> from relatipy.numeric.coordinates.apparent_equatorial_coordinate import EquatorialCoordinate
        >>> aoe = ApparentOrbitalElements(0.0, 1e7, 0.1, 30.0, 40.0, 50.0, 60.0,
        ...                               zeta=45.0, eta=90.0, mass=4.0e6)
        >>> ec = EquatorialCoordinate.from_apparent_orbital_elements(aoe, 8.5 * u.kpc)
        >>> ec.zeta, ec.eta
        (45.0, 90.0)
        """
        cart = aoe.convert_to_cartesian()
        x_bh, y_bh, z_bh = cart.xs[1], cart.xs[2], cart.xs[3]
        return cls(
            x_bh, y_bh, z_bh,
            zeta=aoe.zeta,
            eta=aoe.eta,
            distance=distance,
            mass=aoe.mass,
        )

    # ------------------------------------------------------------------
    # Display
    # ------------------------------------------------------------------

    def to_sexagesimal(self):
        """
        Return RA and Dec in standard sexagesimal notation.

        Right ascension is expressed in hours (HH:MM:SS.ss) and declination
        in degrees with sign (±DD:MM:SS.s), following the IAU convention.

        Returns
        -------
        ra_str : str
            RA offset formatted as ``HH:MM:SS.ss`` (hours, minutes, seconds).
        dec_str : str
            Dec offset formatted as ``±DD:MM:SS.s`` (degrees, arcminutes,
            arcseconds).

        Notes
        -----
        The conversion uses :math:`1\\,\\mathrm{h} = 15^\\circ`, so the sign
        of ``ra_str`` reflects whether the offset is east (positive) or west
        (negative).

        Examples
        --------
        >>> import astropy.units as u
        >>> from relatipy.numeric.coordinates.apparent_equatorial_coordinate import EquatorialCoordinate
        >>> ec = EquatorialCoordinate(0.0, 1000.0, 0.0, zeta=0.0, eta=0.0,
        ...                           distance=8.5 * u.kpc, mass=4.0e6)
        >>> ra_hms, dec_dms = ec.to_sexagesimal()
        >>> isinstance(ra_hms, str) and isinstance(dec_dms, str)
        True
        """
        ra_angle = Angle(self.ra, unit=u.deg)
        dec_angle = Angle(self.dec, unit=u.deg)
        ra_str = ra_angle.to_string(unit=u.hourangle, sep=":", precision=2, pad=True)
        dec_str = dec_angle.to_string(unit=u.deg, sep=":", precision=1, alwayssign=True)
        return ra_str, dec_str

    # ------------------------------------------------------------------
    # Dunder helpers
    # ------------------------------------------------------------------

    def __repr__(self):
        return (
            f"EquatorialCoordinate("
            f"ra={self.ra:.6f} deg, "
            f"dec={self.dec:.6f} deg, "
            f"distance={self.distance})"
        )

    def __eq__(self, other):
        if not isinstance(other, EquatorialCoordinate):
            return NotImplemented
        return (
            numpy.isclose(self.ra, other.ra)
            and numpy.isclose(self.dec, other.dec)
            and self.distance == other.distance
        )
