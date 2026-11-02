"""One entry point for user-defined forces: :func:`get_force`."""

from __future__ import annotations

import ast
from collections.abc import Mapping, Sequence

from .axisymmetric import SurfaceDensityDisk
from .compiled import Force, Potential
from .core import ExternalForce

__all__ = ["get_force"]

_SYSTEMS = (
    ("cartesian", frozenset({"x", "y", "z"})),
    ("cylindrical", frozenset({"R", "phi", "z"})),
    ("spherical", frozenset({"r", "theta", "phi"})),
)
_COORDINATES = frozenset().union(*(names for _, names in _SYSTEMS))
_VELOCITIES = {
    "cartesian": frozenset({"vx", "vy", "vz"}),
    "cylindrical": frozenset({"v_R", "v_phi", "v_z"}),
    "spherical": frozenset({"v_r", "v_theta", "v_phi"}),
}


def _names(texts: Sequence[str]) -> set[str]:
    names = set()
    for text in texts:
        try:
            tree = ast.parse(text.strip(), mode="eval")
        except SyntaxError as exc:
            raise ValueError(f"invalid expression {text!r}: {exc.msg}") from exc
        names |= {node.id for node in ast.walk(tree) if isinstance(node, ast.Name)}
    return names


def _infer_coords(texts: Sequence[str], velocity: bool) -> str | None:
    """Coordinate system that contains every coordinate name, or None if mixed."""
    used = _names(texts)
    for system, names in _SYSTEMS:
        allowed = names | (_VELOCITIES[system] if velocity else frozenset())
        if used & (_COORDINATES | frozenset().union(*_VELOCITIES.values())) <= allowed:
            return system
    return None


def _mixed_potential(text: str, names: tuple[str, ...]):
    """Potential written with names of several systems, as a Cartesian callable."""
    import sympy as sp

    from .._codegen.forces import parse_expression

    def definition(x, y, z, *parameters):
        cylindrical = sp.sqrt(x**2 + y**2)
        namespace = {
            "x": x, "y": y, "z": z,
            "R": cylindrical, "r": sp.sqrt(x**2 + y**2 + z**2),
            "phi": sp.atan2(y, x), "theta": sp.atan2(cylindrical, z),
        }
        namespace |= dict(zip(names, parameters))
        return parse_expression(text, namespace)

    return definition


def get_force(
    *,
    potential: str | None = None,
    force: str | Sequence[str] | None = None,
    surface_density: str | None = None,
    params: Mapping[str, object] | None = None,
    coords: str | None = None,
    **options,
) -> ExternalForce:
    """Build an external force from a potential, a force or a surface density.

    Exactly one of ``potential``, ``force`` and ``surface_density`` is given,
    as text. Coordinates are written with the fixed names ``x, y, z``
    (Cartesian), ``R, phi, z`` (cylindrical) and ``r, theta, phi``
    (spherical) of the black-hole frame; every other name must be a key of
    ``params``.

    Parameters
    ----------
    potential : str, optional
        Newtonian potential ``Phi``; compiled as :class:`Potential`. Names of
        several systems may be mixed.
    force : str or sequence of str, optional
        Acceleration components in the orthonormal basis of one coordinate
        system; compiled as :class:`Force`.
    surface_density : str, optional
        Surface density ``Sigma(R)`` of a razor-thin disk, with ``R`` in
        units of ``r_out``; it must be axisymmetric. Built as
        :class:`SurfaceDensityDisk` (summed over rings once, then
        tabulated), which needs the keyword ``options`` ``mass`` and
        ``r_out`` and accepts ``r_in``, ``softening`` and ``normal``.
    params : mapping, optional
        Parameter values. Run-time parameters for ``potential`` and
        ``force``; fixed, dimensionless shape parameters for
        ``surface_density``.
    coords : {"cartesian", "cylindrical", "spherical"}, optional
        Coordinates of ``potential`` or ``force``; inferred from the names
        used when omitted.
    **options
        Keyword arguments of :class:`SurfaceDensityDisk`.

    Returns
    -------
    ExternalForce

    Raises
    ------
    TypeError
        If not exactly one definition is given, or ``options`` are given for
        a potential or a force.
    ValueError
        If the coordinates of a force are mixed, or a definition is invalid
        (see the classes above).
    """
    given = {
        name: value for name, value in
        (("potential", potential), ("force", force), ("surface_density", surface_density))
        if value is not None
    }
    if len(given) != 1:
        raise TypeError("give exactly one of potential, force and surface_density")
    kind, definition = next(iter(given.items()))
    params = {} if params is None else params
    if not isinstance(params, Mapping):
        raise TypeError("params must be a mapping from name to value")
    if kind == "surface_density":
        if coords is not None:
            raise TypeError("coords does not apply to a surface density")
        return SurfaceDensityDisk(definition, parameters=params, **options)
    if options:
        raise TypeError(f"unexpected options for a {kind}: {', '.join(sorted(options))}")
    texts = [definition] if isinstance(definition, str) else list(definition)
    if not all(isinstance(text, str) for text in texts):
        raise TypeError(f"{kind} must be text")
    if coords is None:
        coords = _infer_coords(texts, velocity=kind == "force")
    if kind == "potential":
        if coords is None:
            return Potential(_mixed_potential(definition, tuple(params)),
                             coords="cartesian", parameters=params)
        return Potential(definition, coords=coords, parameters=params)
    if coords is None:
        raise ValueError(
            "force components mix names of several coordinate systems; "
            "use one system or give coords")
    return Force(definition, coords=coords, parameters=params)
