"""RelatiPy public package.

The implementation is organized into ``coordinates``, ``forces``,
``geodesic``, ``metrics``, ``observables``, and ``plotting`` packages.  The
most common classes, and :func:`plot_sols` for several trajectories in one
frame, are also re-exported here.
"""

from __future__ import annotations

from . import coordinates, forces, geodesic, metrics, observables, plotting
from .forces import Force, Plummer, Potential
from .geodesic import (
    InitialState, IntegrationError, IntegrationInfo, IntegrationTerminated,
    IntegrationWarning,
    Orbit, Solution, State, Termination,
)
from .metrics import Kerr, PerturbedKerr
from .observables import (
    KerrMcmcGaussianLikelihood, KerrMcmcModel, KerrMcmcMultiModel,
    PerturbedKerrMcmcModel, PerturbedKerrMcmcMultiModel,
)
from .plotting import plot_sols

__all__ = [
    "coordinates",
    "forces",
    "geodesic",
    "metrics",
    "observables",
    "plotting",
    "Kerr",
    "PerturbedKerr",
    "Force",
    "Potential",
    "Plummer",
    "KerrMcmcModel",
    "KerrMcmcMultiModel",
    "KerrMcmcGaussianLikelihood",
    "PerturbedKerrMcmcModel",
    "PerturbedKerrMcmcMultiModel",
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
