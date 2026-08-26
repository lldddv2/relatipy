"""Bank subsystem: persistent caches for expensive numeric computations.

Currently exposes :class:`PathsBank` for geodesic trajectories.
"""

from .paths_bank import Group, PathsBank, init_bank

__all__ = ["PathsBank", "init_bank", "Group"]
