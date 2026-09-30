"""Internal delegation of :class:`State` properties to wrapper classes."""

from __future__ import annotations

from operator import attrgetter

from .state import State

STATE_PROPERTIES = (
    "tau",
    "txyz",
    "trqp",
    "tRQP",
    "t",
    "x",
    "y",
    "z",
    "xyz",
    "r",
    "theta",
    "phi",
    "R",
    "Theta",
    "Phi",
    "vxyz",
    "vrqp",
    "vRQP",
    "vx",
    "vy",
    "vz",
    "vr",
    "vtheta",
    "vphi",
    "vR",
    "vTheta",
    "vPhi",
    "ut",
    "uxyz",
    "urqp",
    "uRQP",
    "ux",
    "uy",
    "uz",
    "ur",
    "utheta",
    "uphi",
    "uR",
    "uTheta",
    "uPhi",
    "state_vector",
)


_DIRECT_STATE_PROPERTY_DOCS = {
    "tau": """Return stored proper time.

    Returns
    -------
    astropy.units.Quantity
        Read-only scalar or shape ``(n,)`` proper time.
    """,
    "txyz": """Return Cartesian coordinates, computing and caching them if needed.

    Returns
    -------
    CartesianCoordinates
        Immutable scalar or vectorized Cartesian coordinate view.
    """,
    "trqp": """Return spherical coordinates, computing and caching them if needed.

    Returns
    -------
    SphericalCoordinates
        Immutable scalar or vectorized spherical coordinate view.
    """,
    "tRQP": """Return stored Boyer--Lindquist coordinates.

    Returns
    -------
    BoyerLindquistCoordinates
        Immutable scalar or vectorized Boyer--Lindquist coordinate view.
    """,
    "vxyz": """Return Cartesian coordinate velocity, computing it if needed.

    Returns
    -------
    astropy.units.Quantity
        Read-only derivative ``dx^i/dt`` with respect to coordinate time,
        with shape ``(3,)`` or ``(n, 3)``.
    """,
    "vrqp": """Return spherical coordinate velocity, computing it if needed.

    Returns
    -------
    SphericalVelocity
        Immutable scalar or vectorized spherical velocity view.
    """,
    "vRQP": """Return stored Boyer--Lindquist coordinate velocity.

    Returns
    -------
    BoyerLindquistVelocity
        Immutable scalar or vectorized Boyer--Lindquist velocity view.
    """,
    "ut": """Return the stored temporal four-velocity component.

    Returns
    -------
    astropy.units.Quantity
        Read-only dimensionless scalar or shape ``(n,)`` component.
    """,
    "uxyz": """Return Cartesian spatial four-velocity, computing it if needed.

    Returns
    -------
    astropy.units.Quantity
        Read-only derivative ``dx^i/dtau`` with respect to proper time, in
        length per time, with shape ``(3,)`` or ``(n, 3)``.
    """,
    "urqp": """Return spherical spatial four-velocity, computing it if needed.

    Returns
    -------
    SphericalFourVelocity
        Immutable scalar or vectorized spherical four-velocity view.
    """,
    "uRQP": """Return stored Boyer--Lindquist spatial four-velocity.

    Returns
    -------
    BoyerLindquistFourVelocity
        Immutable scalar or vectorized Boyer--Lindquist four-velocity view.
    """,
}


def _state_property(attribute: str, name: str) -> property:
    """Create a property that delegates to the wrapped state.

    Parameters
    ----------
    attribute : str
        Name of the owner attribute holding the :class:`State`.
    name : str
        Name of the :class:`State` attribute to expose.

    Returns
    -------
    property
        Read-only, documented descriptor returning the ``name`` attribute
        of the state stored in ``attribute``.

    Notes
    -----
    Properties backed by a :class:`State` property reuse that property's
    documentation without its example block.  Direct stored fields use the
    descriptions in ``_DIRECT_STATE_PROPERTY_DOCS``.  This affects
    introspection only; the getter remains a single attribute delegation.
    """
    state_descriptor = getattr(State, name, None)
    doc = getattr(state_descriptor, "__doc__", None)
    if doc is not None:
        doc = doc.split("\n\n        Examples\n", maxsplit=1)[0]
    else:
        doc = _DIRECT_STATE_PROPERTY_DOCS.get(name)
    return property(attrgetter(f"{attribute}.{name}"), doc=doc)


def delegate_state_properties(owner: type, attribute: str) -> None:
    """Expose every :data:`STATE_PROPERTIES` name on ``owner``.

    Parameters
    ----------
    owner : type
        Class receiving read-only properties.
    attribute : str
        Name of the instance attribute holding the delegated :class:`State`.
    """
    for name in STATE_PROPERTIES:
        setattr(owner, name, _state_property(attribute, name))
