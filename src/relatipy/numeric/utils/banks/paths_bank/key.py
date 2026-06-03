"""
Canonical-key construction and hashing for :class:`PathsBank`.

A cache entry is uniquely identified by:

- The metric type and its numeric parameters (``mass``, ``a``).
- The initial-conditions class, its ``xs`` and ``vs`` arrays, and any extra
  ``kwargs`` (e.g. spin ``a`` for Boyer-Lindquist, ``mass`` for orbital
  elements).
- The proper-time grid ``taus`` (full content, dtype and shape).
- The integrator name and its tolerance/step settings, with project defaults
  merged in so that omitted kwargs collapse to the same key as explicit
  defaults.

Floats are canonicalized via :func:`repr` (round-trip exact in CPython >= 3.1).
Arrays are summarized by SHA-256 of their raw bytes plus dtype and shape, so
two arrays differing only in dtype hash distinctly.

The final hash is SHA-256 of the JSON dump of the canonical dict (excluding the
``relatipy_version`` field, which is stored for diagnostics but does not
invalidate the bank when bumped).
"""

from __future__ import annotations

import hashlib
import json
from typing import Any, Mapping

import numpy as np

from ._version import BANK_FORMAT_VERSION

try:
    from relatipy import __version__ as _RELATIPY_VERSION
except Exception:  # pragma: no cover - extremely defensive
    _RELATIPY_VERSION = "unknown"


_INTEGRATOR_DEFAULTS = {
    "integrator": "Radau",
    "adaptative": True,
    "steps_per_period": 100,
    "rtol": 1e-12,
    "atol": 1e-13,
}


def _to_float(x: Any) -> float:
    """Coerce ``x`` to a plain ``float``, unwrapping ``astropy.Quantity``."""
    to_value = getattr(x, "to_value", None)
    if callable(to_value):
        try:
            return float(to_value())
        except TypeError:
            return float(to_value(getattr(x, "unit", None)))
    return float(x)


def _canonical_array(arr: Any) -> dict:
    """Summarize an array as ``{dtype, shape, sha256}`` of its raw bytes."""
    a = np.ascontiguousarray(arr)
    return {
        "dtype": str(a.dtype),
        "shape": list(a.shape),
        "sha256": hashlib.sha256(a.tobytes()).hexdigest(),
    }


def _canonical_metric(metric: Any) -> dict:
    """Reduce a metric instance to its identifying numeric parameters."""
    cls = type(metric).__name__
    if cls == "Schwarzschild":
        return {"type": cls, "mass": repr(_to_float(metric.mass))}
    if cls == "Kerr":
        return {
            "type": cls,
            "mass": repr(_to_float(metric.mass)),
            "a": repr(_to_float(metric.a)),
        }
    # Generic fallback: stable for any BaseMetric subclass with mass + kwargs.
    payload: dict = {"type": cls, "mass": repr(_to_float(metric.mass))}
    extra = getattr(metric, "kwargs", {}) or {}
    payload["kwargs"] = {k: repr(_to_float(v)) for k, v in sorted(extra.items())}
    return payload


def _canonical_ic(ic: Any) -> dict:
    """Reduce a coordinate / orbital-elements instance to its identifying state.

    Uses ``ic.state_vector`` (present on both :class:`CoordinateBase` and
    :class:`OrbitalElements`) plus ``ic.kwargs`` so the canonical key works
    uniformly across initial-condition types.
    """
    state = np.asarray(ic.state_vector, dtype=np.float64)
    kwargs = getattr(ic, "kwargs", {}) or {}
    return {
        "type": type(ic).__name__,
        "system_name": getattr(ic, "name_metric", type(ic).__name__),
        "state_vector": _canonical_array(state),
        "kwargs": {k: repr(_to_float(v)) for k, v in sorted(kwargs.items())},
    }


def _canonical_integrator(kw: Mapping[str, Any]) -> dict:
    """Merge integrator kwargs with defaults so calls collapse canonically."""
    merged = {**_INTEGRATOR_DEFAULTS, **dict(kw)}
    return {
        "integrator": str(merged["integrator"]),
        "adaptative": bool(merged["adaptative"]),
        "steps_per_period": int(merged["steps_per_period"]),
        "rtol": repr(float(merged["rtol"])),
        "atol": repr(float(merged["atol"])),
    }


def canonical_key(
    metric: Any,
    initial_conditions: Any,
    taus: Any,
    integrator_kwargs: Mapping[str, Any] | None = None,
) -> dict:
    """
    Build the canonical, JSON-friendly key dict for a cache lookup.

    Parameters
    ----------
    metric : object
        Metric instance (e.g. ``Schwarzschild`` or ``Kerr``).
    initial_conditions : CoordinateBase
        Coordinate-system instance carrying ``xs``, ``vs``, and ``kwargs``.
    taus : array_like
        Proper-time grid to be passed to ``Geodesic.get_path``.
    integrator_kwargs : mapping, optional
        Any of ``integrator``, ``adaptative``, ``steps_per_period``, ``rtol``,
        ``atol``. Missing entries are filled with project defaults.

    Returns
    -------
    dict
        Canonical key. Pass to :func:`hash_key` to obtain the SHA-256 digest.
    """
    return {
        "bank_format_version": BANK_FORMAT_VERSION,
        "relatipy_version": _RELATIPY_VERSION,
        "metric": _canonical_metric(metric),
        "initial_conditions": _canonical_ic(initial_conditions),
        "taus": _canonical_array(np.asarray(taus)),
        "integrator": _canonical_integrator(integrator_kwargs or {}),
    }


def hash_key(key: dict) -> str:
    """SHA-256 hex digest of ``key`` excluding ``relatipy_version``."""
    hashable = {k: v for k, v in key.items() if k != "relatipy_version"}
    blob = json.dumps(hashable, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()
