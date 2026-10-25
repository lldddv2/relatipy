"""Private symbolic optimizer that turns SymPy expressions into C99.

SymPy is an optional dependency (``relatipy[symbolic]``).  Importing this
package does not import SymPy; the first access to one of the names below
does, and raises an :class:`ImportError` naming the extra when SymPy is
missing.

Names
-----
optimize
    Optimize expressions, matrices or tensors into C99 statements.
OptimizedCode, OptimizationReport, OperationCount
    Result, report and operation-count types of :func:`optimize`.
symbolic_residuals, numeric_max_relative_error, compile_and_compare
    Verification utilities.
count_c_operations
    Arithmetic-operator count of a C source fragment.
build_force_source, ForceSource
    C source of an ``rp_external_acceleration_fn`` from a potential or a
    force given as a Python callable or a string.
"""

from __future__ import annotations

from importlib import import_module
from typing import Any

_EXPORTS = {
    "optimize": ".optimizer",
    "OptimizedCode": ".optimizer",
    "OptimizationReport": ".optimizer",
    "OperationCount": ".optimizer",
    "Definition": ".optimizer",
    "count_operations": ".optimizer",
    "reduce_angle": ".optimizer",
    "symbolic_residuals": ".verify",
    "numeric_max_relative_error": ".verify",
    "compile_and_compare": ".verify",
    "count_c_operations": ".verify",
    "build_force_source": ".forces",
    "ForceSource": ".forces",
}

__all__ = sorted(_EXPORTS)


def __getattr__(name: str) -> Any:
    module_name = _EXPORTS.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(import_module(module_name, __name__), name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
