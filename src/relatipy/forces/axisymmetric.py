"""Built-in tabulated thin disk and spheroid (Piran et al. 2026).

Both models evaluate an axisymmetric, plane-symmetric field from a table of
the unit-mass distribution with unit outer radius, built once in native code
and interpolated with C2 bicubic splines. The table depends only on the
shape (slope, inner-to-outer radius ratio, softening or axis ratio); the
mass, the outer radius and the symmetry axis are written into a small
native model at bind time and can change between MCMC proposals without a
rebuild. Built tables are cached per shape in each process.
"""

from __future__ import annotations

import math
import threading
from collections import OrderedDict

import numpy as np
from astropy import units as u

from ._bound import geometric_conversion
from .core import ExternalForce

# Tables cost about 1 s and 11 MB (thin disk) to several seconds and 22 MB
# (flat spheroid); a few shapes are kept per process. Placed models keep
# their own table alive, so eviction never invalidates a bound force.
_TABLE_CACHE_SIZE = 4
_TABLES: OrderedDict[tuple, object] = OrderedDict()
_TABLES_LOCK = threading.Lock()


def _table(shape: tuple) -> object:
    """Return the built table of ``shape``, building it on first use."""
    with _TABLES_LOCK:
        table = _TABLES.get(shape)
        if table is not None:
            _TABLES.move_to_end(shape)
            return table
    from .. import _force_core

    kind, *arguments = shape
    if kind == "thin_disk":
        table = _force_core.AxisymmetricTable.thin_disk(*arguments)
    else:
        table = _force_core.AxisymmetricTable.spheroid(*arguments)
    with _TABLES_LOCK:
        table = _TABLES.setdefault(shape, table)
        _TABLES.move_to_end(shape)
        while len(_TABLES) > _TABLE_CACHE_SIZE:
            _TABLES.popitem(last=False)
    return table


def _normal_angles(normal: object) -> tuple[float, float]:
    """Return polar and azimuthal angles of a nonzero finite 3-vector."""
    vector = np.asarray(normal, dtype=np.float64)
    if vector.shape != (3,) or not np.all(np.isfinite(vector)):
        raise ValueError("normal must contain three finite components")
    norm = float(np.linalg.norm(vector))
    if not norm > 0 or not math.isfinite(norm):
        raise ValueError("normal must be nonzero")
    x, y, z = vector / norm
    return math.acos(min(1.0, max(-1.0, z))), math.atan2(y, x)


def _unit_normal(theta: float, phi: float) -> np.ndarray:
    sine = math.sin(theta)
    return np.array([sine * math.cos(phi), sine * math.sin(phi), math.cos(theta)])


def _length(value: object, name: str) -> u.Quantity:
    if not isinstance(value, u.Quantity) or value.ndim != 0:
        raise TypeError(f"{name} must be a scalar length quantity")
    length = u.Quantity(value, copy=True)
    length.to(u.m)
    if not np.isfinite(length.value):
        raise ValueError(f"{name} must be finite")
    length.flags.writeable = False
    return length


def _real(value: object, name: str) -> float:
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be a real number")
    try:
        number = float(u.Quantity(value, u.one).value)
    except (TypeError, ValueError, u.UnitsError):
        raise TypeError(f"{name} must be a dimensionless real number") from None
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite")
    return number


class _TabulatedForce(ExternalForce):
    """Shared implementation of :class:`ThinDisk` and :class:`Spheroid`."""

    __slots__ = ("_shape", "_outer", "_outer_conversion", "_options")

    _mass_name = "mass"
    _outer_name = "r_out"

    def __init__(
        self,
        mass: u.Quantity,
        outer: u.Quantity,
        shape: tuple,
        options: dict,
        normal: object,
        angles: tuple[float, float] | None,
    ) -> None:
        if not isinstance(mass, u.Quantity):
            raise TypeError("mass must be an astropy quantity")
        mass.to(u.kg)
        theta, phi = _normal_angles(normal) if angles is None else angles
        object.__setattr__(self, "_shape", shape)
        object.__setattr__(self, "_outer", outer)
        object.__setattr__(
            self, "_outer_conversion",
            geometric_conversion(outer.unit, self._outer_name),
        )
        object.__setattr__(self, "_options", options)
        super().__init__({
            self._mass_name: mass,
            "normal_theta": u.Quantity(theta, u.rad),
            "normal_phi": u.Quantity(phi, u.rad),
        })

    @property
    def normal(self) -> np.ndarray:
        """Unit symmetry axis ``(sin t cos p, sin t sin p, cos t)``.

        Returns
        -------
        numpy.ndarray
            New array of the axis built from ``normal_theta`` ``t`` and
            ``normal_phi`` ``p``, in the frame the force is used in.
        """
        theta = float(self._parameters["normal_theta"].to_value(u.rad))
        phi = float(self._parameters["normal_phi"].to_value(u.rad))
        return _unit_normal(theta, phi)

    def _check_vector(self, values: np.ndarray) -> np.ndarray:
        values = super()._check_vector(values)
        if values[0] < 0:
            raise ValueError(f"{type(self).__name__} mass must not be negative")
        return values

    def _key(self) -> tuple:
        return (
            self._shape, float(self._outer.value), self._outer.unit.to_string(),
            self._parameter_key(),
        )

    def _native_callback(self) -> tuple[int, object]:
        if self._native is None:
            table = _table(self._shape)
            object.__setattr__(self, "_native", (int(table.callback), table))
        return self._native

    def _new_storage(self) -> object:
        from .. import _force_core

        return _force_core.AxisymmetricModel(self._native_callback()[1])

    def _storage_address(self, storage: object) -> int:
        return int(storage.address)

    def _storage_data(self, storage: object) -> np.ndarray:
        return storage.parameters

    def _write(
        self,
        storage: object,
        mass: float,
        values: np.ndarray,
        rotation: np.ndarray | None,
    ) -> None:
        factors, powers = self._conversions
        model_mass = float(values[0] * factors[0] * mass ** (-powers[0]))
        factor, power = self._outer_conversion
        outer = float(self._outer.value) * factor * mass ** (-power)
        normal = _unit_normal(float(values[1]), float(values[2]))
        if rotation is not None:
            # Rotation columns are the black-hole axes in observer axes.
            normal = np.asarray(rotation).T @ normal
        storage.set(model_mass, outer, normal)

    def _replace(self, parameters: dict[str, object]) -> _TabulatedForce:
        angles = tuple(
            float(u.Quantity(parameters[name], u.rad).value)
            for name in ("normal_theta", "normal_phi")
        )
        return type(self)(
            mass=parameters[self._mass_name],
            **{self._outer_name: self._outer}, **self._options, _angles=angles,
        )

    def __reduce__(self) -> tuple:
        """Pickle the constructor arguments; the table is rebuilt lazily."""
        angles = tuple(
            float(self._parameters[name].to_value(u.rad))
            for name in ("normal_theta", "normal_phi")
        )
        return _rebuild_tabulated, (
            type(self), self._parameters[self._mass_name],
            {self._outer_name: self._outer, **self._options}, angles,
        )

    def __repr__(self) -> str:
        """Return an evaluable-style representation."""
        options = ", ".join(f"{key}={value!r}" for key, value in self._options.items())
        normal = tuple(round(float(item), 15) for item in self.normal)
        return (
            f"{type(self).__name__}(mass={self._parameters[self._mass_name]!r}, "
            f"{self._outer_name}={self._outer!r}, {options}, normal={normal!r})"
        )


class ThinDisk(_TabulatedForce):
    """Razor-thin power-law disk (Piran et al. 2026).

    Parameters
    ----------
    mass : astropy.units.Quantity
        Finite disk mass, zero or positive.
    r_out : astropy.units.Quantity
        Finite, positive outer radius ``R_out``.
    r_in : astropy.units.Quantity, optional
        Inner radius ``R_in`` in ``[0, R_out)``, default zero.
    slope : float, optional
        Power ``alpha`` of ``Sigma(R) proportional to R^alpha`` between
        ``R_in`` and ``R_out``, default ``-1``; ``R_in = 0`` needs
        ``alpha > -2``.
    softening : float, optional
        Ring softening ``eps`` in ``(0, 1)``: each ring of radius ``R'`` is
        softened with ``H = eps R'``. Default ``0.01``.
    normal : array-like, optional
        Disk normal, default ``(0, 0, 1)``. With
        :class:`~relatipy.metrics.PerturbedKerr` it is in the black-hole
        frame (``z`` along the spin); with the perturbed MCMC models it is
        in observer axes ``(Dec, RA, away)``, so the default is a disk seen
        face-on.

    Raises
    ------
    TypeError
        If a quantity or number has the wrong type.
    astropy.units.UnitConversionError
        If ``mass`` is not a mass or a radius not a length.
    ValueError
        If a value is not finite or outside its range, or ``normal`` is not a
        nonzero three-vector.

    Notes
    -----
    The run-time parameters, and the force block of the perturbed MCMC
    models, are ``disk_mass`` (in the unit of ``mass``) and the normal angles
    ``normal_theta`` and ``normal_phi`` (radians): the normal is
    ``(sin t cos p, sin t sin p, cos t)`` in the frame above, so ``t`` is
    measured from the third axis (the spin, or the direction away from the
    observer) and ``p`` from the first axis (``x``, or Dec) towards the
    second (``y``, or RA). ``r_out`` and the shape stay fixed; a different
    ``r_out`` with ``r_in = 0`` reuses the same table.

    The acceleration is the exact gradient of the softened potential, read
    from a table of the unit disk built on first use (about 1 s and 11 MB,
    then cached per shape in each process); the potential and the
    acceleration are separate interpolants that agree to about ``1e-5`` in
    the softened plane. The callback fails at the disk centre ``r = 0``.

    Examples
    --------
    The fiducial disk of Piran et al. reaches the apocentre of S301
    (``a (1 + e)`` with ``a = 83.0 mas``, ``e = 0.9832`` at ``8.277 kpc``):

    >>> from astropy import units as u
    >>> from relatipy.forces import ThinDisk
    >>> disk = ThinDisk(mass=1e3 * u.M_sun, r_out=1362.5 * u.au)
    >>> disk.parameter_names
    ('disk_mass', 'normal_theta', 'normal_phi')
    """

    __slots__ = ()
    _mass_name = "disk_mass"
    _outer_name = "r_out"

    def __init__(
        self,
        *,
        mass: u.Quantity,
        r_out: u.Quantity,
        r_in: u.Quantity | None = None,
        slope: float = -1.0,
        softening: float = 0.01,
        normal: object = (0.0, 0.0, 1.0),
        _angles: tuple[float, float] | None = None,
    ) -> None:
        """Validate the disk shape, mass, outer radius and normal."""
        outer = _length(r_out, "r_out")
        inner = _length(0 * outer.unit if r_in is None else r_in, "r_in")
        if not outer.value > 0:
            raise ValueError("r_out must be positive")
        ratio = float((inner / outer).to_value(u.one))
        if not 0 <= ratio < 1:
            raise ValueError("r_in must be in [0, r_out)")
        alpha = _real(slope, "slope")
        if ratio == 0 and not alpha > -2:
            raise ValueError("slope must exceed -2 when r_in is zero")
        eps = _real(softening, "softening")
        if not 0 < eps < 1:
            raise ValueError("softening must be in (0, 1)")
        super().__init__(
            mass, outer, ("thin_disk", alpha, ratio, eps),
            {"r_in": inner, "slope": alpha, "softening": eps}, normal, _angles,
        )


class Spheroid(_TabulatedForce):
    """Oblate power-law spheroid (Piran et al. 2026, Sect. 4).

    Parameters
    ----------
    mass : astropy.units.Quantity
        Finite total mass, zero or positive.
    s_out : astropy.units.Quantity
        Finite, positive outer label ``s_out`` (the equatorial outer radius).
    axis_ratio : float, optional
        Axis ratio ``q`` in ``(0, 1]``, default ``1`` (a sphere).
    s_in : astropy.units.Quantity, optional
        Inner label in ``[0, s_out)``, default zero.
    slope : float, optional
        Power ``alpha`` of ``rho proportional to s^alpha`` on the similar
        spheroids ``s^2 = R^2 + z^2 / q^2``, default ``-2``; ``s_in = 0``
        needs ``alpha > -3``.
    normal : array-like, optional
        Symmetry axis, in the frames described for :class:`ThinDisk`.

    Raises
    ------
    TypeError, ValueError, astropy.units.UnitConversionError
        As :class:`ThinDisk`.

    Notes
    -----
    The run-time parameters are ``mass``, ``normal_theta`` and
    ``normal_phi``, with the angle convention of :class:`ThinDisk`. With
    ``q = 1`` the field is that of a sphere, ``G M(<r) / r^2`` with
    ``M(<r) = M (r / s_out)^(alpha + 3)`` inside ``s_out``. As ``q -> 0`` at
    fixed mass it tends to a razor-thin disk with ``Sigma`` proportional to
    ``R^(alpha + 1)`` well inside the rim and without softening. The force
    is only continuous across the density edges. Tables of flat spheroids
    take longer to build (about 6 s and 22 MB at ``q = 0.05``).

    Examples
    --------
    >>> from astropy import units as u
    >>> from relatipy.forces import Spheroid
    >>> Spheroid(mass=1e3 * u.M_sun, s_out=1362.5 * u.au, axis_ratio=0.1).parameter_names
    ('mass', 'normal_theta', 'normal_phi')
    """

    __slots__ = ()
    _mass_name = "mass"
    _outer_name = "s_out"

    def __init__(
        self,
        *,
        mass: u.Quantity,
        s_out: u.Quantity,
        axis_ratio: float = 1.0,
        s_in: u.Quantity | None = None,
        slope: float = -2.0,
        normal: object = (0.0, 0.0, 1.0),
        _angles: tuple[float, float] | None = None,
    ) -> None:
        """Validate the spheroid shape, mass, outer label and axis."""
        outer = _length(s_out, "s_out")
        inner = _length(0 * outer.unit if s_in is None else s_in, "s_in")
        if not outer.value > 0:
            raise ValueError("s_out must be positive")
        ratio = float((inner / outer).to_value(u.one))
        if not 0 <= ratio < 1:
            raise ValueError("s_in must be in [0, s_out)")
        q = _real(axis_ratio, "axis_ratio")
        if not 0 < q <= 1:
            raise ValueError("axis_ratio must be in (0, 1]")
        alpha = _real(slope, "slope")
        if ratio == 0 and not alpha > -3:
            raise ValueError("slope must exceed -3 when s_in is zero")
        super().__init__(
            mass, outer, ("spheroid", alpha, q, ratio),
            {"axis_ratio": q, "s_in": inner, "slope": alpha}, normal, _angles,
        )


def _rebuild_tabulated(
    cls: type, mass: u.Quantity, options: dict, angles: tuple[float, float]
) -> _TabulatedForce:
    """Recreate a pickled tabulated force; its table is built on first use."""
    return cls(mass=mass, **options, _angles=tuple(angles))
