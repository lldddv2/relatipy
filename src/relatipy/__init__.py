"""RelatiPy public package.

The implementation is organized into ``coordinates``, ``geodesic``,
``metrics``, ``observables``, and ``plotting`` packages.  The most common
classes, and :func:`plot_sols` for several trajectories in one frame, are
also re-exported here.
"""

from __future__ import annotations

from . import coordinates, geodesic, metrics, observables, plotting
from .geodesic import (
    InitialState, IntegrationError, IntegrationInfo, IntegrationTerminated,
    IntegrationWarning,
    Orbit, Solution, State, Termination,
)
from .metrics import Kerr
from .observables import KerrMcmcModel
from .plotting import plot_sols

__all__ = [
    "coordinates",
    "geodesic",
    "metrics",
    "observables",
    "plotting",
    "Kerr",
    "KerrMcmcModel",
    "Orbit",
    "State",
    "InitialState",
    "Solution",
    "IntegrationInfo",
    "Termination",
    "IntegrationError",
    "IntegrationTerminated",
    "IntegrationWarning",
    "plot_sols",
]
