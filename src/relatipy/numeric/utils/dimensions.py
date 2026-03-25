"""
Convert physical quantities to dimensionless geometric units.

This module provides :class:`Validator`, which normalizes
:class:`astropy.units.Quantity` inputs (and plain floats taken as already
normalized geometric values) for numerical relativity code. Reference scales
are defined in ``relatipy.numeric.constants`` (one solar mass,
:math:`L_{\\mathrm{ref}} = G M_{\\mathrm{ref}} / c^2`,
:math:`T_{\\mathrm{ref}} = G M_{\\mathrm{ref}} / c^3`).

Normalizations used:

.. math::

    \\frac{M}{M_{\\mathrm{ref}}}, \\quad
    \\frac{r}{L_{\\mathrm{ref}}}, \\quad
    \\frac{t}{T_{\\mathrm{ref}}}, \\quad
    \\frac{v}{c},

and angles in radians. Angular velocity in rad/s is multiplied by
:math:`T_{\\mathrm{ref}}` so that geometric products with radius and angle
remain consistent.

Notes
-----
For arguments that are not ``Quantity`` instances, values are assumed to
already be in the corresponding geometric normalization (no SI conversion).

The module-level :obj:`validator` is a default :class:`Validator` instance
for convenience.

Examples
--------
>>> from relatipy.numeric.utils.dimensions import validator
>>> validator.validate_mass(1.0)
1.0
"""

from astropy import units as u
from ..constants import _c_SI, _M_ref_kg, _L_ref, _T_ref
from numpy import asarray


class Validator:
    """
    Validate inputs and map physical units to dimensionless geometric form.

    Dispatch methods convert :class:`astropy.units.Quantity` values where
    needed; plain floats are returned after basic checks or unchanged,
    depending on the quantity.

    Attributes
    ----------
    associated_units_validation : dict
        Maps Astropy units to handler methods: ``u.kg``, ``u.s``, ``u.m``,
        ``u.rad``, ``u.m / u.s``, and ``u.rad / u.s``.

    Examples
    --------
    >>> from relatipy.numeric.utils.dimensions import Validator
    >>> v = Validator()
    >>> v.validate_length(2.0)
    2.0
    """

    def __init__(self):
        """
        Initialize the unit-to-method dispatch mapping.

        Examples
        --------
        >>> from relatipy.numeric.utils.dimensions import Validator
        >>> isinstance(Validator(), Validator)
        True
        """
        self.associated_units_validation = {
            u.kg: self.validate_mass,
            u.s: self.validate_time,
            u.m: self.validate_length,
            u.rad: self.validate_angle,
            u.m / u.s: self.validate_velocity,
            u.rad / u.s: self.validate_angular_velocity,
        }

    def validate_mass(self, mass):
        """
        Normalize mass to :math:`M / M_{\\mathrm{ref}}`.

        Converts ``kg`` to dimensionless mass using one solar mass as
        :math:`M_{\\mathrm{ref}}` (see ``_M_ref_kg`` in constants).

        Parameters
        ----------
        mass : float or astropy.units.Quantity
            Mass in kilograms if a ``Quantity``, or already-dimensionless
            mass if a float.

        Returns
        -------
        float
            Dimensionless mass :math:`M / M_{\\mathrm{ref}}`, or ``mass``
            unchanged when a non-negative float is passed.

        Raises
        ------
        ValueError
            If ``mass`` is a ``Quantity`` in incompatible units, or if a
            float ``mass`` is negative.

        Examples
        --------
        >>> from relatipy.numeric.utils.dimensions import validator
        >>> validator.validate_mass(1.0)
        1.0
        """
        if isinstance(mass, u.Quantity):
            try:
                return mass.to(u.kg).value / _M_ref_kg
            except Exception:
                raise ValueError(f"Invalid mass units: {mass.unit}")
        if mass < 0:
            raise ValueError("Mass must be positive")
        return mass

    def validate_length(self, length):
        """
        Normalize length to :math:`r / L_{\\mathrm{ref}}`.

        Converts meters to geometric length using
        :math:`L_{\\mathrm{ref}} = G M_{\\mathrm{ref}} / c^2`.

        Parameters
        ----------
        length : float or astropy.units.Quantity
            Length in meters if a ``Quantity``, or geometric length if a float.

        Returns
        -------
        float
            Dimensionless :math:`r / L_{\\mathrm{ref}}`, or ``length`` if a
            float was given.

        Raises
        ------
        ValueError
            If ``length`` is a ``Quantity`` in incompatible units.

        Examples
        --------
        >>> from relatipy.numeric.utils.dimensions import validator
        >>> validator.validate_length(1.0)
        1.0
        """
        if isinstance(length, u.Quantity):
            try:
                return length.to(u.m).value / _L_ref
            except Exception:
                raise ValueError(f"Invalid length units: {length.unit}")
        return length

    def validate_time(self, time):
        """
        Normalize time to :math:`t / T_{\\mathrm{ref}}`.

        Converts seconds using :math:`T_{\\mathrm{ref}} = G M_{\\mathrm{ref}} / c^3`.

        Parameters
        ----------
        time : float or astropy.units.Quantity
            Time in seconds if a ``Quantity``, or geometric time if a float.

        Returns
        -------
        float
            Dimensionless :math:`t / T_{\\mathrm{ref}}`, or ``time`` if a float.

        Raises
        ------
        ValueError
            If ``time`` is a ``Quantity`` in incompatible units.

        Examples
        --------
        >>> from relatipy.numeric.utils.dimensions import validator
        >>> validator.validate_time(0.5)
        0.5
        """
        if isinstance(time, u.Quantity):
            try:
                return time.to(u.s).value / _T_ref
            except Exception:
                raise ValueError(f"Invalid time units: {time.unit}")
        return time

    def validate_angle(self, angle):
        """
        Normalize angle to radians as a dimensionless number.

        Parameters
        ----------
        angle : float or astropy.units.Quantity
            Angle; if a ``Quantity``, it is converted to radians.

        Returns
        -------
        float
            Numeric value in radians.

        Raises
        ------
        ValueError
            If ``angle`` is a ``Quantity`` in incompatible units.

        Examples
        --------
        >>> from relatipy.numeric.utils.dimensions import validator
        >>> validator.validate_angle(0.1)
        0.1
        """
        if isinstance(angle, u.Quantity):
            try:
                return angle.to(u.rad).value
            except Exception:
                raise ValueError(f"Invalid angle units: {angle.unit}")
        return angle

    def validate_velocity(self, velocity):
        """
        Normalize velocity to :math:`v / c`.

        Parameters
        ----------
        velocity : float or astropy.units.Quantity
            Speed in m/s if a ``Quantity``, or :math:`v/c` if a float.

        Returns
        -------
        float
            Dimensionless :math:`v/c`, or ``velocity`` if a float.

        Raises
        ------
        ValueError
            If ``velocity`` is a ``Quantity`` in incompatible units, or if a
            float ``velocity`` is greater than 1 (superluminal in geometric
            units).

        Examples
        --------
        >>> from relatipy.numeric.utils.dimensions import validator
        >>> validator.validate_velocity(0.1)
        0.1
        """
        if isinstance(velocity, u.Quantity):
            try:
                return velocity.to(u.m / u.s).value / _c_SI
            except Exception:
                raise ValueError(f"Invalid velocity units: {velocity.unit}")
        if velocity > 1:
            raise ValueError("Velocity must be less than the speed of light")
        return velocity

    def validate_angular_velocity(self, angular_velocity):
        """
        Scale angular velocity to geometric units.

        SI angular velocity :math:`\\mathrm{d}\\theta/\\mathrm{d}t` in rad/s is
        mapped to a dimensionless rate by multiplying by
        :math:`T_{\\mathrm{ref}}`, so that

        .. math::

            r_{\\mathrm{geom}} \\frac{\\mathrm{d}\\theta}{\\mathrm{d}t_{\\mathrm{geom}}}
            = \\frac{r}{L_{\\mathrm{ref}}} \\,
              \\frac{\\mathrm{d}\\theta}{\\mathrm{d}t_{\\mathrm{SI}}} \\,
              T_{\\mathrm{ref}}

        stays consistent with geometric velocity conventions.

        Parameters
        ----------
        angular_velocity : float or astropy.units.Quantity
            Angular velocity in rad/s if a ``Quantity``, or the geometric
            equivalent if a float.

        Returns
        -------
        float
            ``(dθ/dt)_SI * T_ref`` for quantities, or ``angular_velocity``
            unchanged for floats.

        Raises
        ------
        ValueError
            If ``angular_velocity`` is a ``Quantity`` in incompatible units.

        Examples
        --------
        >>> from relatipy.numeric.utils.dimensions import validator
        >>> validator.validate_angular_velocity(1e-5)
        1e-05
        """
        if isinstance(angular_velocity, u.Quantity):
            try:
                return angular_velocity.to(u.rad / u.s).value * _T_ref
            except Exception:
                raise ValueError(
                    f"Invalid angular velocity units: {angular_velocity.unit}"
                )
        return angular_velocity

    def validate_scalar(self, scalar):
        """
        Dispatch a scalar ``Quantity`` to the appropriate validator by unit.

        Parameters
        ----------
        scalar : float or astropy.units.Quantity
            Value to normalize. If a ``Quantity``, its unit must appear in
            :attr:`associated_units_validation`.

        Returns
        -------
        float
            Result of the matching ``validate_*`` method, or ``scalar`` if not
            a ``Quantity``.

        Raises
        ------
        ValueError
            If ``scalar`` is a ``Quantity`` with an unsupported unit or
            conversion fails.

        Examples
        --------
        >>> import astropy.units as u
        >>> from relatipy.numeric.utils.dimensions import validator
        >>> validator.validate_scalar(3.0)
        3.0
        """
        if isinstance(scalar, u.Quantity):
            try:
                return self.associated_units_validation[scalar.unit](scalar)
            except Exception:
                raise ValueError(f"Invalid scalar units: {scalar.unit}")
        return scalar

    def validate_vector(self, vector):
        """
        Apply :meth:`validate_scalar` element-wise and return a NumPy array.

        Each entry is converted to SI if it is a ``Quantity`` (via ``.si``),
        then validated with :meth:`validate_scalar`.

        Parameters
        ----------
        vector : sequence
            Iterable of scalars or ``Quantity`` instances.

        Returns
        -------
        numpy.ndarray
            One-dimensional array of normalized floats.

        Raises
        ------
        ValueError
            Propagated from :meth:`validate_scalar` when an element cannot be
            converted.

        Examples
        --------
        >>> from relatipy.numeric.utils.dimensions import validator
        >>> validator.validate_vector([1.0, 2.0])
        array([1., 2.])
        """
        vector_si = [x.si if isinstance(x, u.Quantity) else x for x in vector]
        new_vector = []
        for i in vector_si:
            new_vector.append(self.validate_scalar(i))
        return asarray(new_vector)


validator = Validator()
