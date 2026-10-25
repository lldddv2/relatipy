"""User-defined forces and potentials compiled to native callbacks at run time.

A definition (a Python callable or a string) is turned into the C source of
an ``rp_external_acceleration_fn`` by the private symbolic optimizer
(``relatipy._codegen``, which needs SymPy) and compiled with the system C
compiler into a cached shared library (``relatipy._jit``). The Python
definition never runs inside an integrator. Parameters are read by the
callback at run time, so changing them needs no recompilation.
"""

from __future__ import annotations

import importlib
from collections.abc import Callable, Mapping, Sequence

from astropy import units as u

from .core import ExternalForce

_COORDINATES = ("cartesian", "cylindrical", "spherical")
_SYMBOLIC_HINT = (
    "run-time forces need SymPy; install it with `pip install relatipy[symbolic]`"
)


def _generate(
    definition: Callable | str | Sequence[str],
    kind: str,
    coords: str,
    parameters: tuple[str, ...],
) -> tuple[str, str]:
    """Return the generated C source and its function name."""
    try:
        importlib.import_module("sympy")
    except ImportError as exc:
        raise ImportError(_SYMBOLIC_HINT) from exc
    try:
        codegen = importlib.import_module("relatipy._codegen.forces")
    except ModuleNotFoundError as exc:
        if exc.name is not None and exc.name.startswith("relatipy._codegen"):
            raise ImportError(
                "run-time forces need relatipy._codegen, which this build does not include"
            ) from exc
        raise
    generated = codegen.build_force_source(
        definition, kind=kind, coords=coords, parameters=parameters
    )
    return str(generated.source), str(generated.function_name)


def _compile(source: str, function_name: str) -> tuple[int, object]:
    """Compile (or load from the cache) and return the address and handle."""
    try:
        jit = importlib.import_module("relatipy._jit")
    except ModuleNotFoundError as exc:
        if exc.name == "relatipy._jit":
            raise ImportError(
                "run-time forces need relatipy._jit, which this build does not include"
            ) from exc
        raise
    function = jit.compile_function(source, function_name)
    return int(function.address), function


class _CompiledForce(ExternalForce):
    """Shared implementation of :class:`Force` and :class:`Potential`."""

    __slots__ = ("_definition", "_coords", "_source")

    _kind = ""

    def __init__(
        self,
        definition: Callable | str | Sequence[str] | None,
        *,
        coords: str = "cartesian",
        parameters: Mapping[str, object] | None = None,
        _source: tuple[str, str] | None = None,
    ) -> None:
        if coords not in _COORDINATES:
            raise ValueError("coords must be 'cartesian', 'cylindrical' or 'spherical'")
        if _source is None:
            if callable(definition) or isinstance(definition, str):
                pass
            elif (
                isinstance(definition, Sequence)
                and self._kind == "force"
                and len(definition) == 3
                and all(isinstance(item, str) for item in definition)
            ):
                definition = tuple(definition)
            else:
                expected = (
                    "a callable, a string or three strings"
                    if self._kind == "force" else "a callable or a string"
                )
                raise TypeError(f"definition must be {expected}")
        if parameters is None:
            parameters = {}
        if not isinstance(parameters, Mapping):
            raise TypeError("parameters must be a mapping from name to value")
        object.__setattr__(self, "_definition", definition)
        object.__setattr__(self, "_coords", coords)
        object.__setattr__(self, "_source", _source)
        super().__init__(parameters)

    @property
    def coords(self) -> str:
        """Coordinates the definition is written in.

        Returns
        -------
        str
            ``"cartesian"``, ``"cylindrical"`` or ``"spherical"``.
        """
        return self._coords

    @property
    def source(self) -> str:
        """Generated C source of the native callback.

        Returns
        -------
        str
            One C translation unit, generated on first access.

        Raises
        ------
        ImportError
            If SymPy or the code generator is unavailable.
        ValueError
            If the definition is invalid.
        """
        return self._generated()[0]

    def _generated(self) -> tuple[str, str]:
        if self._source is None:
            generated = _generate(
                self._definition, self._kind, self._coords, self.parameter_names
            )
            object.__setattr__(self, "_source", generated)
        return self._source

    def _native_callback(self) -> tuple[int, object]:
        if self._native is None:
            source, function_name = self._generated()
            object.__setattr__(self, "_native", _compile(source, function_name))
        return self._native

    def _key(self) -> tuple:
        return (self._coords, self._generated(), self._parameter_key())

    def _replace(self, parameters: dict[str, object]) -> _CompiledForce:
        return type(self)(
            self._definition, coords=self._coords, parameters=parameters,
            _source=self._source,
        )

    def __reduce__(self) -> tuple:
        """Pickle the generated source and parameters, not the definition."""
        return _rebuild_compiled, (
            type(self), self._coords, self._generated(), dict(self._parameters)
        )

    def __repr__(self) -> str:
        """Return a summary of the definition, coordinates and parameters."""
        definition = self._definition
        if definition is None:
            text = "<compiled source>"
        elif callable(definition):
            text = f"<callable {getattr(definition, '__qualname__', type(definition).__name__)}>"
        else:
            text = repr(definition)
        parameters = ", ".join(
            f"{name!r}: {value!r}" for name, value in self._parameters.items()
        )
        return (
            f"{type(self).__name__}({text}, coords={self._coords!r}, "
            f"parameters={{{parameters}}})"
        )


class Force(_CompiledForce):
    """Newtonian acceleration defined by an expression and compiled to C.

    Parameters
    ----------
    definition : callable, str or sequence of three str
        The three acceleration components in the orthonormal basis of
        ``coords``. A callable receives, positionally, the three coordinates,
        the three velocity components in the same basis, then the declared
        parameters, all as SymPy symbols, and returns three expressions. A
        string is a tuple expression such as ``"(-k*x, -k*y, -k*z)"``; three
        strings give one component each. Strings use the coordinate names
        ``x, y, z`` and ``vx, vy, vz`` (Cartesian), ``R, phi, z`` and
        ``v_R, v_phi, v_z`` (cylindrical), or ``r, theta, phi`` and
        ``v_r, v_theta, v_phi`` (spherical), the declared parameters, numbers,
        ``pi`` and the functions ``sqrt exp log sin cos tan asin acos atan
        atan2 sinh cosh tanh abs erf``.
    coords : {"cartesian", "cylindrical", "spherical"}, optional
        Coordinates of the black-hole Cartesian frame (``z`` along the spin)
        in which the definition is written, default ``"cartesian"``.
    parameters : mapping, optional
        Run-time parameters in the order the callable receives them. Each
        value is a real number (used as given) or a scalar quantity
        (converted to geometric units when the force is bound).

    Raises
    ------
    TypeError
        If ``definition`` or ``parameters`` has the wrong type, or a
        parameter is not a number or scalar quantity.
    ValueError
        If ``coords`` is unknown, a parameter name is not an identifier or a
        parameter is not finite.
    astropy.units.UnitConversionError
        If a parameter unit is not built from length, time, mass and angle.

    Notes
    -----
    The definition is evaluated in geometric units of the black hole with
    ``G = c = M = 1``: lengths in ``G M / c^2``, velocities in ``c`` and the
    acceleration in ``c^4 / (G M)``. Code is generated and compiled when the
    force is first bound, for example by the first
    :meth:`~relatipy.metrics.PerturbedKerr.orbit` integration; SymPy
    (``relatipy[symbolic]``) and a C11 compiler (``$CC`` or ``cc``) are
    required then, and a missing one raises :class:`ImportError` or
    :class:`RuntimeError`. Compiled libraries are cached on disk under a hash
    of the source, so a later process reuses them.

    Pickling stores the generated C source and the parameters, not the
    Python definition, so forces defined with lambdas can be sent to
    multiprocessing pools; a worker loads the cached library without SymPy.
    """

    __slots__ = ()
    _kind = "force"


class Potential(_CompiledForce):
    """Newtonian potential defined by an expression and compiled to C.

    Parameters
    ----------
    definition : callable or str
        Scalar potential ``Phi``. A callable receives, positionally, the
        three coordinates and then the declared parameters as SymPy symbols.
        A string uses the coordinate names of :class:`Force` (without
        velocities), the declared parameters, numbers, ``pi`` and the same
        functions.
    coords : {"cartesian", "cylindrical", "spherical"}, optional
        Coordinates of the definition, default ``"cartesian"``.
    parameters : mapping, optional
        Run-time parameters, as in :class:`Force`.

    Raises
    ------
    TypeError, ValueError, astropy.units.UnitConversionError
        As :class:`Force`.

    Notes
    -----
    The acceleration is ``-grad Phi``, differentiated symbolically. ``Phi``
    is in units of ``c^2`` and lengths in ``G M / c^2``. A static potential
    axisymmetric about the spin conserves ``E' = -u_t exp(Phi)`` and
    ``L' = u_phi exp(Phi)``. Compilation, caching and pickling follow
    :class:`Force`.

    Examples
    --------
    >>> from relatipy.forces import Potential
    >>> Potential("-k / sqrt(R**2 + z**2 + b**2)", coords="cylindrical",
    ...           parameters={"k": 1e-3, "b": 50.0})
    Potential('-k / sqrt(R**2 + z**2 + b**2)', coords='cylindrical', parameters={'k': <Quantity 0.001>, 'b': <Quantity 50.>})
    """

    __slots__ = ()
    _kind = "potential"


def _rebuild_compiled(
    cls: type, coords: str, generated: tuple[str, str], parameters: dict
) -> _CompiledForce:
    """Recreate a pickled compiled force from its generated source."""
    return cls(None, coords=coords, parameters=parameters, _source=tuple(generated))
