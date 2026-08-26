"""
:class:`PathsBank` — content-addressed cache for geodesic trajectories.

Given a metric, an :class:`~relatipy.numeric.coordinates.base.CoordinateBase`
initial condition, a proper-time grid ``taus``, and integrator settings, a
:class:`PathsBank` instance returns the same trajectory that
``metric.geodesic.get_path(...)`` would produce — but only integrates it the
first time. Subsequent identical calls load the cached trajectory from disk.

The bank also supports named **groups** that bundle related cached paths.

Typical use:

>>> from relatipy.numeric.utils.banks import PathsBank, init_bank
>>> bank = PathsBank("notebooks/paths_bank"); init_bank(bank)         # doctest: +SKIP
>>> path = bank.get_path(metric, ic, taus, integrator="Yoshida6")     # doctest: +SKIP
"""

from ._version import BANK_FORMAT_VERSION
from .groups import Group
from .paths_bank import PathsBank, init_bank
from .reconstruct import BankError, BankIntegrityError

__all__ = [
    "PathsBank",
    "init_bank",
    "Group",
    "BankError",
    "BankIntegrityError",
    "BANK_FORMAT_VERSION",
]
