"""Evaluate several stars around one Kerr black hole from one MCMC vector.

:class:`KerrMcmcMultiModel` composes one :class:`~relatipy.observables.KerrMcmcModel`
per star. The stars share distance, mass, astrometric frame, velocity offset,
spin and solver settings; each star has its own osculating Kepler elements.
:class:`KerrMcmcGaussianLikelihood` turns its predictions into an independent
Gaussian log-likelihood for samplers such as emcee.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from types import MappingProxyType

import numpy as np
from astropy import units as u

from .._validation import immutable_array
from ..geodesic import IntegrationError
from .kerr_mcmc import KerrMcmcModel, _epoch, _epochs


_MILLION_SOLAR_MASS = u.Unit(1e6 * u.M_sun)
_SHARED_PARAMETERS = (
    ("D", u.kpc),
    ("M", _MILLION_SOLAR_MASS),
    ("xS0", u.arcsec),
    ("yS0", u.arcsec),
    ("vxS0", u.arcsec / u.yr),
    ("vyS0", u.arcsec / u.yr),
    ("v_LSR", u.km / u.s),
)
_SPIN_PARAMETERS = (
    ("spin_dec", u.one),
    ("spin_ra", u.one),
    ("spin_away", u.one),
)
_STAR_PARAMETERS = (
    ("t_p", u.yr),
    ("a", u.arcsec),
    ("e", u.one),
    ("inc", u.rad),
    ("Omega", u.rad),
    ("omega", u.rad),
)
_RESERVED_BLOCKS = frozenset(("shared", "spin", "force"))
_LOG_TWO_PI = math.log(2.0 * math.pi)
_DATA_UNITS = {
    "alpha": u.arcsec,
    "delta": u.arcsec,
    "v_los": u.km / u.s,
}


def _number(value: object, unit: u.UnitBase, name: str) -> float:
    """Convert one finite parameter value, plain or a quantity, to ``unit``."""
    if isinstance(value, u.Quantity):
        if value.ndim != 0:
            raise ValueError(f"{name} must be a scalar")
        number = float(value.to_value(unit, equivalencies=u.dimensionless_angles()))
    elif isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be a real number")
    else:
        try:
            array = np.asarray(value, dtype=np.float64)
        except (TypeError, ValueError):
            raise TypeError(f"{name} must be a real number") from None
        if array.ndim != 0:
            raise ValueError(f"{name} must be a scalar")
        number = float(array)
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite")
    return number


class _ParameterLayout:
    """Immutable, picklable layout of named parameter blocks in one vector.

    ``blocks`` is a sequence of ``(block_name, ((parameter, unit), ...))``
    pairs laid out contiguously in that order.
    """

    __slots__ = ("_names", "_units", "_slices", "_index")

    def __init__(
        self, blocks: Sequence[tuple[str, Sequence[tuple[str, u.UnitBase]]]]
    ) -> None:
        names: list[str] = []
        units: list[u.UnitBase] = []
        slices: dict[str, slice] = {}
        for block, parameters in blocks:
            if block in slices:
                raise ValueError(f"duplicate parameter block {block!r}")
            start = len(names)
            for name, unit in parameters:
                names.append(name)
                units.append(unit)
            slices[block] = slice(start, len(names))
        if len(set(names)) != len(names):
            raise ValueError("parameter names must be unique")
        self._names = tuple(names)
        self._units = tuple(units)
        self._slices = slices
        self._index = {name: position for position, name in enumerate(names)}

    @property
    def size(self) -> int:
        return len(self._names)

    def pack(self, values: Mapping[str, object]) -> np.ndarray:
        """Return the parameter vector for a mapping that names every entry."""
        if not isinstance(values, Mapping):
            raise TypeError("values must be a mapping from parameter name to value")
        missing = [name for name in self._names if name not in values]
        if missing:
            raise ValueError(f"missing parameters: {', '.join(missing)}")
        unknown = sorted(str(key) for key in values.keys() - self._index.keys())
        if unknown:
            raise ValueError(f"unknown parameters: {', '.join(unknown)}")
        return np.array(
            [_number(values[name], unit, name)
             for name, unit in zip(self._names, self._units)],
            dtype=np.float64,
        )

    def unpack(self, vector: object) -> dict[str, float]:
        """Return a name-to-value mapping for one parameter vector."""
        array = np.asarray(vector, dtype=np.float64)
        if array.shape != (self.size,):
            raise ValueError(f"params must contain {self.size} numbers")
        return dict(zip(self._names, array.tolist()))


class _MultiStarBase:
    """Shared stars, epochs and parameter layout for multi-star models.

    The vector holds the shared block, the spin block, one block per star,
    then any ``extra_blocks`` in order. Subclasses implement ``_predict`` for
    a vector accepted by ``_checked``, which they may extend; ``_predict``
    then raises :class:`ValueError` only for a native initial-state
    rejection and :class:`~relatipy.geodesic.IntegrationError` for an
    integration failure.
    """

    __slots__ = (
        "_stars",
        "_astrometry_times",
        "_spectroscopy_times",
        "_reference_epoch",
        "_layout",
        "_gather",
        "_spin_slice",
    )

    def __init__(
        self,
        astrometry_times: Mapping[str, object],
        spectroscopy_times: Mapping[str, object],
        reference_epoch: object,
        extra_blocks: Sequence[tuple[str, Sequence[tuple[str, u.UnitBase]]]] = (),
    ) -> None:
        if not isinstance(astrometry_times, Mapping) or not isinstance(
            spectroscopy_times, Mapping
        ):
            raise TypeError(
                "astrometry_times and spectroscopy_times must be mappings from star name"
            )
        stars = tuple(astrometry_times)
        if not stars:
            raise ValueError("at least one star is required")
        if set(stars) != set(spectroscopy_times):
            raise ValueError(
                "astrometry_times and spectroscopy_times must name the same stars"
            )
        for star in stars:
            if not isinstance(star, str):
                raise TypeError("star names must be strings")
            if not star or "." in star or star.strip() != star:
                raise ValueError(
                    f"star name {star!r} must be nonempty, without dots or "
                    "surrounding whitespace"
                )
            if star in _RESERVED_BLOCKS:
                raise ValueError(f"star name {star!r} is reserved")
        self._stars = stars
        self._astrometry_times = {
            star: _epochs(astrometry_times[star], f"astrometry_times[{star!r}]")
            for star in stars
        }
        self._spectroscopy_times = {
            star: _epochs(spectroscopy_times[star], f"spectroscopy_times[{star!r}]")
            for star in stars
        }
        self._reference_epoch = _epoch(reference_epoch, "reference_epoch")

        blocks = [("shared", _SHARED_PARAMETERS), ("spin", _SPIN_PARAMETERS)]
        blocks.extend(
            (star, tuple((f"{star}.{name}", unit) for name, unit in _STAR_PARAMETERS))
            for star in stars
        )
        for block, _ in extra_blocks:
            if block in stars:
                raise ValueError(f"parameter block {block!r} collides with a star name")
        blocks.extend(extra_blocks)
        self._layout = _ParameterLayout(blocks)

        # One fancy index gathers every star's legacy 13-vector
        # (D, M, t_p, a, e, inc, Omega, omega, xS0, yS0, vxS0, vyS0, v_LSR).
        slices = self._layout._slices
        shared = slices["shared"]
        rows = []
        for star in stars:
            own = slices[star]
            rows.append([
                shared.start, shared.start + 1,
                *range(own.start, own.stop),
                *range(shared.start + 2, shared.stop),
            ])
        self._gather = immutable_array(np.array(rows, dtype=np.intp))
        self._spin_slice = slices["spin"]

    @property
    def stars(self) -> tuple[str, ...]:
        """Star names in parameter-block order."""
        return self._stars

    @property
    def reference_epoch(self) -> float:
        """Zero point of the linear astrometric frame drift in Julian years."""
        return self._reference_epoch

    @property
    def astrometry_times(self) -> Mapping[str, np.ndarray]:
        """Read-only astrometric epochs in Julian years for each star."""
        return MappingProxyType(self._astrometry_times)

    @property
    def spectroscopy_times(self) -> Mapping[str, np.ndarray]:
        """Read-only spectroscopic epochs in Julian years for each star."""
        return MappingProxyType(self._spectroscopy_times)

    @property
    def n_parameters(self) -> int:
        """Length of the parameter vector."""
        return self._layout.size

    @property
    def parameter_names(self) -> tuple[str, ...]:
        """Parameter names in vector order."""
        return self._layout._names

    @property
    def parameter_units(self) -> Mapping[str, u.UnitBase]:
        """Unit of each parameter, keyed by name."""
        return MappingProxyType(dict(zip(self._layout._names, self._layout._units)))

    @property
    def parameter_slices(self) -> Mapping[str, slice]:
        """Slice of each parameter block, keyed by block name."""
        return MappingProxyType(self._layout._slices)

    def pack(self, values: Mapping[str, object]) -> np.ndarray:
        """Build a parameter vector from named values.

        Parameters
        ----------
        values : mapping
            Exactly the names in :attr:`parameter_names`. Each value is a
            real number in the unit of :attr:`parameter_units` or a scalar
            :class:`~astropy.units.Quantity` convertible to it. Angles accept
            any angular unit.

        Returns
        -------
        numpy.ndarray
            A new float array ordered as :attr:`parameter_names`.

        Raises
        ------
        TypeError
            If ``values`` is not a mapping or a value is not a real number.
        ValueError
            If a name is missing or unknown, or a value is not a finite
            scalar.
        astropy.units.UnitConversionError
            If a quantity has an incompatible unit.
        """
        return self._layout.pack(values)

    def unpack(self, params: object) -> dict[str, float]:
        """Return a new mapping from parameter name to value.

        Parameters
        ----------
        params : array-like
            One vector of length :attr:`n_parameters`.

        Returns
        -------
        dict
            Plain floats in the units of :attr:`parameter_units`, in vector
            order.

        Raises
        ------
        ValueError
            If ``params`` does not have the layout length.
        """
        return self._layout.unpack(params)

    def _checked(self, params: object) -> np.ndarray:
        """Validate one vector against the shared Kerr parameter domain."""
        values = np.asarray(params, dtype=np.float64)
        if values.shape != (self._layout.size,) or not np.all(np.isfinite(values)):
            raise ValueError(
                f"params must contain {self._layout.size} finite numbers"
            )
        if values[0] <= 0 or values[1] <= 0:
            raise ValueError("D and M must be positive")
        rows = values[self._gather]
        if np.any(rows[:, 3] <= 0):
            raise ValueError("a must be positive for every star")
        eccentricity = rows[:, 4]
        if np.any(eccentricity < 0) or np.any(eccentricity >= 1):
            raise ValueError("e must be in [0, 1) for every star")
        spin = values[self._spin_slice]
        if np.linalg.norm(spin) > 1.0:
            raise ValueError("spin_vector magnitude must be at most one")
        return values

    def _predict(
        self, values: np.ndarray
    ) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
        raise NotImplementedError

    def __call__(
        self, params: object
    ) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
        return self._predict(self._checked(params))


class KerrMcmcMultiModel(_MultiStarBase):
    """Predict astrometry and redshift of several stars around one Kerr hole.

    Each star is evaluated by its own :class:`KerrMcmcModel` in the legacy
    fixed-epoch mode. All stars share the reference epoch, the solver
    settings, and the distance, mass, frame, velocity-offset and spin
    parameters of one proposal.

    Parameters
    ----------
    astrometry_times, spectroscopy_times : mapping
        Mappings from star name to observation epochs in Julian years, as
        arrays or time quantities. Both must name the same stars; the order
        of ``astrometry_times`` sets the star order. A star needs at least
        one epoch in total. Star names are nonempty strings without dots or
        surrounding whitespace and other than ``'shared'``, ``'spin'`` and
        ``'force'``.
    reference_epoch : float or astropy.units.Quantity
        Shared zero point of the linear astrometric frame drift.

    Raises
    ------
    TypeError
        If an epoch argument is not a mapping, a star name is not a string,
        or ``reference_epoch`` is not a real number or time quantity.
    ValueError
        If no star is given, the two mappings name different stars, a star
        name is invalid, an epoch series is not one-dimensional and finite,
        a star has no epoch, or ``reference_epoch`` is not finite.

    Notes
    -----
    The parameter vector is laid out as:

    1. shared block ``(D, M, xS0, yS0, vxS0, vyS0, v_LSR)`` in kpc, million
       solar masses, arcsec (two offsets), arcsec/year (two drifts) and km/s;
    2. spin block ``(spin_dec, spin_ra, spin_away)``, the dimensionless Kerr
       spin vector in observer axes ``(Dec, RA, away)`` with magnitude at most
       one;
    3. for each star ``S``, ``(S.t_p, S.a, S.e, S.inc, S.Omega, S.omega)`` in
       Julian years, arcsec, dimensionless and radians (three angles).

    :attr:`parameter_names`, :attr:`parameter_slices`, :meth:`pack` and
    :meth:`unpack` expose this layout. Each star's prediction equals a
    :class:`KerrMcmcModel` call with the vector ``(D, M, t_p, a, e, inc,
    Omega, omega, xS0, yS0, vxS0, vyS0, v_LSR)`` and the proposal spin.
    Instances are picklable and can be sent to multiprocessing pools.
    """

    __slots__ = ("_models", "_solver_set")

    def __init__(
        self,
        astrometry_times: Mapping[str, object],
        spectroscopy_times: Mapping[str, object],
        *,
        reference_epoch: object,
    ) -> None:
        """Validate per-star epochs and build one Kerr model per star."""
        super().__init__(astrometry_times, spectroscopy_times, reference_epoch)
        self._models = tuple(
            KerrMcmcModel(
                self._astrometry_times[star],
                self._spectroscopy_times[star],
                reference_epoch=self._reference_epoch,
            )
            for star in self._stars
        )
        self._solver_set = False

    def set_solver(
        self,
        *,
        method: str,
        rtol: float,
        atol: float,
        max_steps: int = 100000,
    ) -> None:
        """Set the solver of every star with :meth:`KerrMcmcModel.set_solver`.

        Parameters
        ----------
        method : {'dop853', 'radau'}
            Native integration method.
        rtol, atol : float
            Positive scalar relative and absolute tolerances.
        max_steps : int, optional
            Maximum steps per endpoint call, default ``100000``.

        Raises
        ------
        TypeError, ValueError
            As :meth:`KerrMcmcModel.set_solver`. Invalid settings leave every
            star unchanged.
        """
        # The first call validates; the others cannot fail with the same values.
        for model in self._models:
            model.set_solver(method=method, rtol=rtol, atol=atol, max_steps=max_steps)
        self._solver_set = True

    def _checked(self, params: object) -> np.ndarray:
        if not self._solver_set:
            raise ValueError("call set_solver before evaluating the model")
        return super()._checked(params)

    def _predict(
        self, values: np.ndarray
    ) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
        spin = values[self._spin_slice]
        rows = values[self._gather]
        return {
            star: model(row, spin_vector=spin)
            for star, model, row in zip(self._stars, self._models, rows)
        }

    def __call__(
        self, params: object
    ) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """Return each star's ``(alpha, delta, v_los)`` for one proposal.

        Parameters
        ----------
        params : array-like
            :attr:`n_parameters` finite numbers in the layout of the class
            notes.

        Returns
        -------
        dict
            Maps each star, in :attr:`stars` order, to read-only arrays:
            ``alpha`` and ``delta`` in arcsec at its astrometric epochs and
            ``v_los`` in km/s at its spectroscopic epochs.

        Raises
        ------
        ValueError
            If :meth:`set_solver` was not called, ``params`` has the wrong
            length or non-finite entries, ``D``, ``M`` or a star's ``a`` is
            not positive, a star's ``e`` is outside ``[0, 1)``, the spin
            magnitude exceeds one, or the native evaluator rejects a star's
            initial state.
        relatipy.geodesic.IntegrationError
            If geodesic integration or observed-time matching fails.
        ImportError
            If the compiled Kerr MCMC extension is unavailable.
        """
        return self._predict(self._checked(params))


class KerrMcmcGaussianLikelihood:
    """Independent Gaussian log-likelihood of multi-star observations.

    Parameters
    ----------
    model : KerrMcmcMultiModel
        Model whose stars and epochs the data follow.
    data : mapping
        Maps every star of ``model`` to a mapping with the keys ``alpha``,
        ``delta``, ``v_los`` and ``alpha_err``, ``delta_err``, ``v_los_err``.
        Astrometric values and errors are in arcsec and velocities in km/s,
        as plain arrays or convertible quantities. Each series has the length
        of that star's astrometric or spectroscopic epochs, values are finite
        and errors are positive and finite.

    Raises
    ------
    TypeError
        If ``model`` is not a multi-star model or ``data`` or an entry is not
        a mapping.
    ValueError
        If stars or keys are missing or unknown, or a series has the wrong
        shape, non-finite values or non-positive errors.
    astropy.units.UnitConversionError
        If a quantity has an incompatible unit.

    Notes
    -----
    Calling the instance returns

    .. math::

       \\ln L = -\\frac{1}{2} \\sum_i \\left(\\frac{y_i - m_i}{\\sigma_i}\\right)^2
       - \\sum_i \\ln \\sigma_i - \\frac{N}{2} \\ln 2\\pi

    over all stars and series. It returns ``-inf`` when the native evaluator
    rejects a star's initial state, integration or observed-time matching
    fails, or a prediction is not finite. Parameters outside the model
    domain still raise :class:`ValueError`, so a prior must exclude them.
    Instances are picklable.
    """

    __slots__ = ("_model", "_observed", "_inverse_errors", "_constant")

    def __init__(self, model: _MultiStarBase, data: Mapping[str, object]) -> None:
        """Validate and flatten per-star observations and errors."""
        if not isinstance(model, _MultiStarBase):
            raise TypeError("model must be a multi-star Kerr MCMC model")
        if not isinstance(data, Mapping):
            raise TypeError("data must be a mapping from star name")
        missing = [star for star in model.stars if star not in data]
        if missing:
            raise ValueError(f"data is missing stars: {', '.join(missing)}")
        unknown = sorted(str(star) for star in data.keys() - set(model.stars))
        if unknown:
            raise ValueError(f"data has unknown stars: {', '.join(unknown)}")
        keys = (*_DATA_UNITS, *(f"{key}_err" for key in _DATA_UNITS))
        observed = []
        inverse_errors = []
        log_errors = 0.0
        count = 0
        for star in model.stars:
            entry = data[star]
            if not isinstance(entry, Mapping):
                raise TypeError(f"data[{star!r}] must be a mapping")
            absent = [key for key in keys if key not in entry]
            if absent:
                raise ValueError(f"data[{star!r}] is missing {', '.join(absent)}")
            extra = sorted(str(key) for key in entry.keys() - set(keys))
            if extra:
                raise ValueError(f"data[{star!r}] has unknown keys: {', '.join(extra)}")
            sizes = {
                "alpha": model.astrometry_times[star].size,
                "delta": model.astrometry_times[star].size,
                "v_los": model.spectroscopy_times[star].size,
            }
            values = []
            errors = []
            for key, unit in _DATA_UNITS.items():
                value = self._series(entry[key], unit, f"data[{star!r}][{key!r}]")
                error = self._series(
                    entry[f"{key}_err"], unit, f"data[{star!r}][{key + '_err'!r}]"
                )
                for name, series in ((key, value), (f"{key}_err", error)):
                    if series.shape != (sizes[key],):
                        raise ValueError(
                            f"data[{star!r}][{name!r}] must have {sizes[key]} values"
                        )
                if np.any(error <= 0):
                    raise ValueError(f"data[{star!r}][{key + '_err'!r}] must be positive")
                values.append(value)
                errors.append(error)
            error_vector = np.concatenate(errors)
            observed.append(immutable_array(np.concatenate(values)))
            inverse_errors.append(immutable_array(1.0 / error_vector))
            log_errors += float(np.sum(np.log(error_vector)))
            count += error_vector.size
        self._model = model
        self._observed = tuple(observed)
        self._inverse_errors = tuple(inverse_errors)
        self._constant = -log_errors - 0.5 * count * _LOG_TWO_PI

    @staticmethod
    def _series(value: object, unit: u.UnitBase, name: str) -> np.ndarray:
        if isinstance(value, u.Quantity):
            array = np.asarray(value.to_value(unit), dtype=np.float64)
        else:
            array = np.asarray(value, dtype=np.float64)
        if array.ndim != 1 or not np.all(np.isfinite(array)):
            raise ValueError(f"{name} must be a finite one-dimensional series")
        return array

    @property
    def model(self) -> _MultiStarBase:
        """The model evaluated for each proposal."""
        return self._model

    def __call__(self, params: object) -> float:
        """Return the Gaussian log-likelihood of one parameter vector.

        Parameters
        ----------
        params : array-like
            One vector in the layout of :attr:`model`.

        Returns
        -------
        float
            The log-likelihood, or ``-inf`` on a native evaluation failure.

        Raises
        ------
        ValueError
            If ``params`` is invalid for the model, or the model has no
            solver.
        """
        model = self._model
        values = model._checked(params)
        try:
            predictions = model._predict(values)
        except IntegrationError:
            return -math.inf
        except ValueError:
            # ``_checked`` accepted the inputs, so this is the native
            # rejection of an initial state.
            return -math.inf
        chi_squared = 0.0
        for prediction, observed, inverse in zip(
            predictions.values(), self._observed, self._inverse_errors
        ):
            scaled = (np.concatenate(prediction) - observed) * inverse
            chi_squared += float(np.dot(scaled, scaled))
        if not math.isfinite(chi_squared):
            return -math.inf
        return self._constant - 0.5 * chi_squared
