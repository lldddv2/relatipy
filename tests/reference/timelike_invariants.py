"""Independent, test-only Kerr invariants in ``G = c = M = mu = 1``.

The metric is obtained by numerically inverting the inverse metric printed in
Schmidt (2002), Eq. (3), pp. 2--3, rather than calling the production kernel.
Energy and axial angular momentum follow the paragraph before Eq. (6), p. 3.
Carter Q follows Eqs. (8), (17), pp. 3--4 with fixed rest mass squared one;
it is not the alternative constant ``K = Q + (Lz - a E)**2``.  The norm target
is minus one, from Eqs. (4), (18), pp. 3, 5.  Source: arXiv:gr-qc/0202090.
Brain provenance: INF-022/C001--C004/S001, scientific manual review pending.
No formula in this module is used by the RelatiPy production backend.
"""

from __future__ import annotations

import numpy as np
from astropy import units as u
from astropy.constants import G, c


INVARIANT_NAMES = ("norm", "energy", "Lz", "CarterQ")


def kerr_metric(spin: float, coordinates: np.ndarray) -> np.ndarray:
    """Return independent covariant BL metric, with shape ``(..., 4, 4)``.

    Coordinates have shape ``(..., 4)`` and order ``(t, r, theta, phi)``.
    The caller must stay away from the BL horizon and polar singularities.
    """
    coordinates = np.asarray(coordinates, dtype=float)
    if coordinates.shape[-1] != 4:
        raise ValueError("coordinates must have final dimension four")
    radius, theta = coordinates[..., 1], coordinates[..., 2]
    sine2 = np.sin(theta) ** 2
    sigma = radius**2 + spin**2 * np.cos(theta)**2
    delta = radius**2 - 2 * radius + spin**2
    inverse = np.zeros((*coordinates.shape[:-1], 4, 4), dtype=float)
    inverse[..., 0, 0] = -((radius**2 + spin**2)**2 - delta * spin**2 * sine2) / (delta * sigma)
    inverse[..., 0, 3] = inverse[..., 3, 0] = -2 * spin * radius / (delta * sigma)
    inverse[..., 1, 1] = delta / sigma
    inverse[..., 2, 2] = 1 / sigma
    inverse[..., 3, 3] = (delta - spin**2 * sine2) / (delta * sigma * sine2)
    return np.linalg.inv(inverse)


def timelike_invariants(spin: float, canonical: np.ndarray) -> dict[str, np.ndarray]:
    """Return norm, specific E, specific Lz, and fixed-mass Carter Q.

    ``canonical`` has shape ``(..., 8)``, ordered ``(x, u)``.  Q always uses
    ``mu**2 = 1``; replacing it with the drifting measured norm would conceal
    part of the normalization error.  Values are dimensionless internally.
    """
    canonical = np.asarray(canonical, dtype=float)
    if canonical.shape[-1] != 8:
        raise ValueError("canonical state must have final dimension eight")
    coordinates, velocity = canonical[..., :4], canonical[..., 4:]
    metric = kerr_metric(spin, coordinates)
    momentum = np.einsum("...ij,...j->...i", metric, velocity)
    energy, angular_momentum = -momentum[..., 0], momentum[..., 3]
    theta = coordinates[..., 2]
    carter_q = momentum[..., 2]**2 + np.cos(theta)**2 * (
        spin**2 * (1 - energy**2) + angular_momentum**2 / np.sin(theta)**2
    )
    return {
        "norm": np.einsum("...i,...i->...", momentum, velocity),
        "energy": energy,
        "Lz": angular_momentum,
        "CarterQ": carter_q,
    }


def public_canonical_state(state, mass: u.Quantity) -> np.ndarray:
    """Convert public State/Solution BL quantities into normalized ``(x,u)``.

    This boundary conversion uses only public fields and Astropy units.
    Uppercase ``R, Theta, Phi`` identify BL coordinates in the public API.
    """
    length = (G * mass / c**2).to(u.m)
    time = (G * mass / c**3).to(u.s)
    components = (
        (state.t / time).to_value(u.one),
        (state.R / length).to_value(u.one),
        state.Theta.to_value(u.rad),
        state.Phi.to_value(u.rad),
        state.ut.to_value(u.one),
        (state.uR / c).to_value(u.one),
        (state.uTheta * time).to_value(u.rad),
        (state.uPhi * time).to_value(u.rad),
    )
    return np.stack(components, axis=-1)


def public_orbit_from_canonical(metric, canonical: np.ndarray):
    """Construct a public Orbit from supplied independent BL ``x,u``.

    Convert spatial four-velocity to coordinate velocity ``v = u / ut``;
    the production constructor independently reconstructs its future ut.
    """
    canonical = np.asarray(canonical, dtype=float)
    coordinate_velocity = canonical[5:] / canonical[4]
    time = (G * metric.mass / c**3).to(u.s)
    return metric.orbit(
        t=canonical[0] * time,
        R=canonical[1] * metric.r_g,
        Theta=canonical[2] * u.rad,
        Phi=canonical[3] * u.rad,
        vR=coordinate_velocity[0] * c,
        vTheta=coordinate_velocity[1] * u.rad / time,
        vPhi=coordinate_velocity[2] * u.rad / time,
    )
