"""Tabulated thin disks: the power law of Piran et al. (2026) and any ``Sigma(R)``.

The model evaluates an axisymmetric, plane-symmetric field from a table of
the unit-mass distribution with unit outer radius, built once in native code
and interpolated with C2 bicubic splines. The table depends only on the
shape (surface density, inner-to-outer radius ratio and softening); the
mass, the outer radius and the symmetry axis are written into a small
native model at bind time and can change between MCMC proposals without a
rebuild. Built tables are cached per shape in each process.
"""

from __future__ import annotations

import functools
import importlib
import math
import threading
from collections import OrderedDict
from collections.abc import Mapping

import numpy as np
from astropy import units as u

from ._bound import geometric_conversion
from .core import ExternalForce

# Tables cost about 1 s and 11 MB each; a few shapes are kept per process. Placed models keep
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
    if kind == "thin_mass":
        source, name, values, ratio, eps, beta, inner_exponent = arguments
        from .compiled import _compile

        address, function = _compile(source, name)
        # The density is called only while the table is built.
        table = _force_core.AxisymmetricTable.thin_mass(
            address, np.array(values, dtype=np.float64), ratio, eps, beta, inner_exponent)
        del function
    else:
        table = _force_core.AxisymmetricTable.thin_disk(*arguments)
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
    """Tabulated axisymmetric force; base of :class:`ThinDisk`."""

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
        Disk normal, default ``(0, 0, 1)``. The perturbed MCMC models read
        it in observer axes ``(Dec, RA, away)``, so the default is a disk
        seen face-on; bound without a spin rotation it is in the black-hole
        frame (``z`` along the spin).

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


_SYMBOLIC_HINT = (
    "surface densities need SymPy; install it with `pip install relatipy[symbolic]`"
)


@functools.lru_cache(maxsize=32)
def _density_source(text: str, names: tuple[str, ...], values: tuple[float, ...]) -> tuple:
    """Generated C source, function name and SymPy expression of ``Sigma``."""
    try:
        importlib.import_module("sympy")
    except ImportError as exc:
        raise ImportError(_SYMBOLIC_HINT) from exc
    density = importlib.import_module("relatipy._codegen.density")
    generated = density.build_surface_density_source(
        text, parameters=names, values=dict(zip(names, values)))
    slope = density.central_slope(generated.expression, dict(zip(names, values)))
    return generated.source, generated.function_name, generated.expression, slope


def _check_density(expression, names: tuple[str, ...], values: tuple[float, ...], ratio: float) -> None:
    """Reject a density that is negative or not finite on ``[r_in, r_out]``."""
    import sympy as sp

    radius = sp.Symbol("R", positive=True)
    numeric = sp.lambdify(
        radius, expression.xreplace({sp.Symbol(k, real=True): v for k, v in zip(names, values)}),
        "numpy")
    lower = ratio if ratio > 0 else 1e-6
    samples = np.concatenate([np.geomspace(lower, 1.0, 200), [1.0]])
    with np.errstate(all="ignore"):
        sigma = np.broadcast_to(np.asarray(numeric(samples), dtype=np.float64), samples.shape)
    if not np.all(np.isfinite(sigma)) or np.any(sigma < 0):
        raise ValueError("surface density must be finite and non-negative on [r_in, r_out]")
    if not np.any(sigma > 0):
        raise ValueError("surface density must not vanish on the whole disk")


class SurfaceDensityDisk(_TabulatedForce):
    """Razor-thin axisymmetric disk of any surface density ``Sigma(R)``.

    Parameters
    ----------
    surface_density : str
        ``Sigma`` as an expression of ``R`` (cylindrical radius in units of
        ``r_out``; ``r`` is the same radius in the plane) or of ``x, y``,
        ``phi``, provided it does not depend on the azimuth; any
        normalization. It may use the names in ``parameters``, numbers,
        ``pi`` and the functions of :class:`~relatipy.forces.Potential`.
    mass : astropy.units.Quantity
        Finite disk mass, zero or positive.
    r_out : astropy.units.Quantity
        Finite, positive outer radius ``R_out``.
    r_in : astropy.units.Quantity, optional
        Inner radius ``R_in`` in ``[0, R_out)``, default zero.
    softening : float, optional
        Ring softening ``eps`` in ``(0, 1)``, ``H = eps R'``. Default ``0.01``.
    parameters : mapping, optional
        Dimensionless shape parameters of ``surface_density``. They fix the
        table, so they are not run-time parameters.
    normal : array-like, optional
        Disk normal, as in :class:`ThinDisk`.

    Raises
    ------
    ImportError
        If SymPy is not installed.
    TypeError
        If a quantity or number has the wrong type.
    ValueError
        If ``surface_density`` depends on the azimuth or on undeclared names,
        is negative or not finite on the disk, or has a mass that diverges at
        the centre (``Sigma ~ R^s`` with ``s <= -2`` and ``r_in = 0``).

    Notes
    -----
    ``Sigma`` is compiled to C (SymPy and a C11 compiler, as for
    :class:`~relatipy.forces.Potential`) and summed over softened rings once,
    in native code, into the same kind of table as :class:`ThinDisk`; the
    run-time parameters are ``disk_mass``, ``normal_theta`` and
    ``normal_phi``. ``surface_density="R**gamma"`` with
    ``parameters={"gamma": a}`` is the disk ``ThinDisk(slope=a)``.
    """

    __slots__ = ()
    _mass_name = "disk_mass"
    _outer_name = "r_out"

    def __init__(
        self,
        surface_density: str,
        *,
        mass: u.Quantity,
        r_out: u.Quantity,
        r_in: u.Quantity | None = None,
        softening: float = 0.01,
        parameters: Mapping[str, object] | None = None,
        normal: object = (0.0, 0.0, 1.0),
        _angles: tuple[float, float] | None = None,
    ) -> None:
        """Validate the density, shape, mass, outer radius and normal."""
        if not isinstance(surface_density, str):
            raise TypeError("surface_density must be a string expression")
        outer = _length(r_out, "r_out")
        inner = _length(0 * outer.unit if r_in is None else r_in, "r_in")
        if not outer.value > 0:
            raise ValueError("r_out must be positive")
        ratio = float((inner / outer).to_value(u.one))
        if not 0 <= ratio < 1:
            raise ValueError("r_in must be in [0, r_out)")
        eps = _real(softening, "softening")
        if not 0 < eps < 1:
            raise ValueError("softening must be in (0, 1)")
        if parameters is None:
            parameters = {}
        if not isinstance(parameters, Mapping):
            raise TypeError("parameters must be a mapping from name to value")
        names = tuple(parameters)
        values = tuple(_real(parameters[name], name) for name in names)
        source, function_name, expression, slope = _density_source(surface_density, names, values)
        _check_density(expression, names, values, ratio)
        if ratio > 0:
            beta, inner_exponent = 1.0, 3.0
        else:
            if not slope > -2:
                raise ValueError(
                    "the disk mass diverges at the centre (Sigma ~ R^s with s <= -2); "
                    "give r_in > 0")
            # A smooth density wants u = R; a cusp R^s wants u = R^(s + 2).
            smooth = math.isinf(slope) or (slope >= 0 and abs(slope - round(slope)) < 1e-9)
            beta = 1.0 if smooth else slope + 2.0
            # Steeper than R^1 at the centre, the outer disk dominates: a ~ r.
            inner_exponent = 3.0 if math.isinf(slope) else min(slope + 2.0, 3.0)
        super().__init__(
            mass, outer,
            ("thin_mass", source, function_name, values, ratio, eps, beta, inner_exponent),
            {
                "surface_density": surface_density, "r_in": inner, "softening": eps,
                "parameters": dict(zip(names, values)),
            },
            normal, _angles,
        )


def _rebuild_tabulated(
    cls: type, mass: u.Quantity, options: dict, angles: tuple[float, float]
) -> _TabulatedForce:
    """Recreate a pickled tabulated force; its table is built on first use."""
    return cls(mass=mass, **options, _angles=tuple(angles))
