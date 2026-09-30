"""Internal layout, status decoding, and assembly of native Kerr results.

Column indices mirror ``rp_solution_reconstructed_column`` and status codes
mirror ``rp_kerr_status`` and ``rp_integrator_status`` in the C sources.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from astropy import units as u

from ..coordinates import (
    BoyerLindquistCoordinates,
    BoyerLindquistFourVelocity,
    BoyerLindquistVelocity,
    CartesianCoordinates,
    OrbitalElements,
    SphericalCoordinates,
    SphericalFourVelocity,
    SphericalVelocity,
)
from .state import State

# Reconstructed-row columns.
BL_R, BL_THETA, BL_PHI = 1, 2, 3
BL_VR, BL_VTHETA, BL_VPHI = 4, 5, 6
UT = 7
BL_UR, BL_UTHETA, BL_UPHI = 8, 9, 10
SPH_R, SPH_THETA, SPH_PHI = 11, 12, 13
SPH_VR, SPH_VTHETA, SPH_VPHI = 14, 15, 16
UX, UY, UZ = 17, 18, 19
SPH_UR, SPH_UTHETA, SPH_UPHI = 20, 21, 22
SEMIMAJOR, ECCENTRICITY = 23, 24
INCLINATION, ASCENDING_NODE, PERIAPSIS_ARGUMENT, TRUE_ANOMALY = 25, 26, 27, 28
ROW_WIDTH = 29

# Physical kind of every column used to assemble a state.
COLUMN_KIND = {
    BL_R: "length", BL_THETA: "angle", BL_PHI: "angle",
    BL_VR: "speed", BL_VTHETA: "rate", BL_VPHI: "rate",
    UT: "one",
    BL_UR: "speed", BL_UTHETA: "rate", BL_UPHI: "rate",
    SPH_R: "length", SPH_THETA: "angle", SPH_PHI: "angle",
    SPH_VR: "speed", SPH_VTHETA: "rate", SPH_VPHI: "rate",
    UX: "speed", UY: "speed", UZ: "speed",
    SPH_UR: "speed", SPH_UTHETA: "rate", SPH_UPHI: "rate",
    SEMIMAJOR: "length", ECCENTRICITY: "one",
    INCLINATION: "angle", ASCENDING_NODE: "angle",
    PERIAPSIS_ARGUMENT: "angle", TRUE_ANOMALY: "angle",
}

_KERR_STATUS = {
    1: "internal null pointer",
    2: "non-finite numerical value",
    3: "invalid parameter or outside the exterior Boyer-Lindquist chart",
    4: "Boyer-Lindquist coordinate singularity",
    5: "Kerr ring singularity",
    6: "numerical range exceeded",
    7: "non-timelike coordinate velocity",
    8: "no real null tangent",
}

_INTEGRATOR_STATUS = {
    1: "internal null pointer",
    2: "invalid argument",
    3: "non-finite numerical value",
    4: "right-hand side evaluation failed",
    5: "step size underflow",
    6: "maximum number of steps reached",
    7: "implicit solver did not converge",
    8: "stopped by an internal event",
    9: "internal event evaluation failed",
    10: "invariant projection failed",
}


def describe_reconstruction_status(code: int) -> str:
    """Return a readable reason for a native Kerr reconstruction status.

    Examples
    --------
    >>> describe_reconstruction_status(4)
    'Boyer-Lindquist coordinate singularity; native status 4'
    """
    reason = _KERR_STATUS.get(int(code), "unknown reconstruction failure")
    return f"{reason}; native status {int(code)}"


def describe_integrator_status(code: object) -> str:
    """Return a readable reason for a native integrator status.

    Examples
    --------
    >>> describe_integrator_status(5)
    'step size underflow; native status 5'
    """
    if not isinstance(code, (int, np.integer)):
        return "unknown integrator failure"
    reason = _INTEGRATOR_STATUS.get(int(code), "unknown integrator failure")
    return f"{reason}; native status {int(code)}"


def validate_reconstructed_rows(rows: np.ndarray) -> None:
    """Reject non-finite native values except approved element sentinels.

    Parameters
    ----------
    rows : numpy.ndarray
        Reconstructed native rows of shape ``(n, 29)``.

    Raises
    ------
    ValueError
        If a value is non-finite outside the approved sentinels: ``+inf``
        semi-major axis and ``nan`` undefined angles.
    """
    valid_axis = np.isfinite(rows[:, SEMIMAJOR]) | np.isposinf(rows[:, SEMIMAJOR])
    angles = rows[:, INCLINATION:ROW_WIDTH]
    if (
        not np.all(np.isfinite(rows[:, :SEMIMAJOR]))
        or not np.all(valid_axis)
        or not np.all(np.isfinite(rows[:, ECCENTRICITY]))
        or not np.all(np.isfinite(angles) | np.isnan(angles))
    ):
        raise ValueError("native state reconstruction returned non-finite values")


def assemble_state(
    column: Callable[[int], u.Quantity],
    *,
    tau: u.Quantity,
    t: u.Quantity,
    x: u.Quantity,
    y: u.Quantity,
    z: u.Quantity,
    vxyz: u.Quantity,
    state_type: type[State] = State,
) -> State:
    """Build a state from Cartesian values and reconstructed columns.

    Parameters
    ----------
    column : callable
        Returns the physical quantity for one reconstructed column index,
        already in its presentation unit.
    tau, t : astropy.units.Quantity
        Proper and coordinate times.
    x, y, z : astropy.units.Quantity
        Cartesian position components.
    vxyz : astropy.units.Quantity
        Cartesian coordinate velocity, shape ``(3,)`` or ``(n, 3)``.
    state_type : type, optional
        :class:`State` or a subclass to instantiate.

    Returns
    -------
    State
        Validated state of ``state_type``.
    """
    ux, uy, uz = column(UX), column(UY), column(UZ)
    return state_type(
        tau=tau,
        txyz=CartesianCoordinates(t, x, y, z),
        trqp=SphericalCoordinates(
            t, column(SPH_R), column(SPH_THETA), column(SPH_PHI)
        ),
        tRQP=BoyerLindquistCoordinates(
            t, column(BL_R), column(BL_THETA), column(BL_PHI)
        ),
        vxyz=vxyz,
        vrqp=SphericalVelocity(
            column(SPH_VR), column(SPH_VTHETA), column(SPH_VPHI)
        ),
        vRQP=BoyerLindquistVelocity(
            column(BL_VR), column(BL_VTHETA), column(BL_VPHI)
        ),
        ut=column(UT),
        uxyz=np.stack(
            [ux.value, uy.to_value(ux.unit), uz.to_value(ux.unit)], axis=-1
        ) * ux.unit,
        urqp=SphericalFourVelocity(
            column(SPH_UR), column(SPH_UTHETA), column(SPH_UPHI)
        ),
        uRQP=BoyerLindquistFourVelocity(
            column(BL_UR), column(BL_UTHETA), column(BL_UPHI)
        ),
        orbital_elements=OrbitalElements(
            a=column(SEMIMAJOR),
            e=column(ECCENTRICITY).to_value(u.one),
            inc=column(INCLINATION),
            Omega=column(ASCENDING_NODE),
            omega=column(PERIAPSIS_ARGUMENT),
            f=column(TRUE_ANOMALY),
        ),
    )
