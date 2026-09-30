"""Dimensionally valid builders shared by public contract tests.

The builders assemble records from arbitrary values; they never evaluate
physical formulae, so the resulting states are not geodesic solutions.
"""

from __future__ import annotations

import numpy as np
from astropy import units as u

from relatipy.coordinates import (
    BoyerLindquistCoordinates,
    BoyerLindquistFourVelocity,
    BoyerLindquistVelocity,
    CartesianCoordinates,
    OrbitalElements,
    SphericalCoordinates,
    SphericalFourVelocity,
    SphericalVelocity,
)
from relatipy.geodesic import IntegrationInfo, Solution, State, Termination


def make_state(n: int | None = None, cls: type[State] = State) -> State:
    """Build a dimensionally valid state without evaluating physical formulae."""
    if n is None:
        scalar_or_series = 0.0
        spatial = np.array([1.0, 2.0, 3.0])
        eccentricity = 0.2
    else:
        scalar_or_series = np.arange(n, dtype=float)
        spatial = np.arange(n * 3, dtype=float).reshape(n, 3)
        eccentricity = np.full(n, 0.2)

    return cls(
        tau=scalar_or_series * u.s,
        txyz=CartesianCoordinates(
            scalar_or_series * u.s,
            (scalar_or_series + 1) * u.km,
            (scalar_or_series + 2) * u.km,
            (scalar_or_series + 3) * u.km,
        ),
        trqp=SphericalCoordinates(
            scalar_or_series * u.s,
            (scalar_or_series + 4) * u.km,
            (scalar_or_series + 0.1) * u.rad,
            (scalar_or_series + 0.2) * u.rad,
        ),
        tRQP=BoyerLindquistCoordinates(
            scalar_or_series * u.s,
            (scalar_or_series + 5) * u.km,
            (scalar_or_series + 0.3) * u.rad,
            (scalar_or_series + 0.4) * u.rad,
        ),
        vxyz=spatial * u.km / u.s,
        vrqp=SphericalVelocity(
            (scalar_or_series + 1) * u.km / u.s,
            (scalar_or_series + 2) * u.rad / u.s,
            (scalar_or_series + 3) * u.rad / u.s,
        ),
        vRQP=BoyerLindquistVelocity(
            (scalar_or_series + 1) * u.km / u.s,
            (scalar_or_series + 2) * u.rad / u.s,
            (scalar_or_series + 3) * u.rad / u.s,
        ),
        ut=(scalar_or_series + 1) * u.one,
        uxyz=(spatial + 1) * u.km / u.s,
        urqp=SphericalFourVelocity(
            (scalar_or_series + 1) * u.km / u.s,
            (scalar_or_series + 2) * u.rad / u.s,
            (scalar_or_series + 3) * u.rad / u.s,
        ),
        uRQP=BoyerLindquistFourVelocity(
            (scalar_or_series + 1) * u.km / u.s,
            (scalar_or_series + 2) * u.rad / u.s,
            (scalar_or_series + 3) * u.rad / u.s,
        ),
        orbital_elements=OrbitalElements(
            a=(scalar_or_series + 10) * u.km,
            e=eccentricity,
            inc=(scalar_or_series + 0.1) * u.rad,
            Omega=(scalar_or_series + 0.2) * u.rad,
            omega=(scalar_or_series + 0.3) * u.rad,
            f=(scalar_or_series + 0.4) * u.rad,
        ),
    )


def make_series_state(n: int = 4) -> State:
    """Build an ``n``-sample state series with distinct values per view."""
    values = np.arange(n, dtype=float)
    t = values * u.s
    xyz = np.column_stack((values + 1, values + 2, values + 3)) * u.km
    vector = np.column_stack((values + 4, values + 5, values + 6)) * u.km / u.s
    return State(
        tau=t,
        txyz=CartesianCoordinates(t=t, x=xyz[:, 0], y=xyz[:, 1], z=xyz[:, 2]),
        trqp=SphericalCoordinates(
            t=t,
            r=(values + 10) * u.km,
            theta=(values + 0.1) * u.rad,
            phi=(values + 0.2) * u.rad,
        ),
        tRQP=BoyerLindquistCoordinates(
            t=t,
            R=(values + 11) * u.km,
            Theta=(values + 0.3) * u.rad,
            Phi=(values + 0.4) * u.rad,
        ),
        vxyz=vector,
        vrqp=SphericalVelocity(
            vr=(values + 1) * u.km / u.s,
            vtheta=(values + 1) * u.rad / u.s,
            vphi=(values + 2) * u.rad / u.s,
        ),
        vRQP=BoyerLindquistVelocity(
            vR=(values + 2) * u.km / u.s,
            vTheta=(values + 3) * u.rad / u.s,
            vPhi=(values + 4) * u.rad / u.s,
        ),
        ut=np.ones(n) * u.one,
        uxyz=vector,
        urqp=SphericalFourVelocity(
            ur=(values + 1) * u.km / u.s,
            utheta=(values + 2) * u.rad / u.s,
            uphi=(values + 3) * u.rad / u.s,
        ),
        uRQP=BoyerLindquistFourVelocity(
            uR=(values + 4) * u.km / u.s,
            uTheta=(values + 5) * u.rad / u.s,
            uPhi=(values + 6) * u.rad / u.s,
        ),
        orbital_elements=OrbitalElements(
            a=(values + 20) * u.km,
            e=values / 10,
            inc=(values + 0.1) * u.rad,
            Omega=(values + 0.2) * u.rad,
            omega=(values + 0.3) * u.rad,
            f=(values + 0.4) * u.rad,
        ),
    )


def make_solution(*, status: int = 0, termination: Termination | None = None) -> Solution:
    """Wrap :func:`make_series_state` in a solution with fixed diagnostics."""
    return Solution(
        state=make_series_state(),
        integration=IntegrationInfo(
            method="TEST",
            rtol=1e-9,
            atol=1e-12,
            max_step=None,
            first_step=None,
            n_steps=3,
            nfev=10,
        ),
        status=status,
        message="test result",
        termination=termination,
    )
