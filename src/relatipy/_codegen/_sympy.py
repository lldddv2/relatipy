"""Lazy access to the optional SymPy dependency."""

from __future__ import annotations

from types import ModuleType

EXTRA_HINT = (
    "relatipy._codegen requires SymPy; install it with `pip install sympy`"
)


def require_sympy() -> ModuleType:
    """Import SymPy or raise an ImportError that explains how to install it.

    Returns
    -------
    module
        The imported ``sympy`` module.

    Raises
    ------
    ImportError
        If SymPy is not installed.
    """
    try:
        import sympy
    except ImportError as exc:  # pragma: no cover - exercised with monkeypatch
        raise ImportError(EXTRA_HINT) from exc
    return sympy
