"""
NPZ-based storage for individual cached paths.

Each cached path lives at ``paths/<hash>.npz`` and is self-describing:
the ``meta`` field contains the integrator name, coordinate class, and
``relatipy``/bank versions, so a path can be rebuilt even if the central
manifest is lost.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from ._version import BANK_FORMAT_VERSION

try:
    from relatipy import __version__ as _RELATIPY_VERSION
except Exception:  # pragma: no cover
    _RELATIPY_VERSION = "unknown"


def save(
    npz_path: Path,
    coord: Any,
    taus: Any,
    integrator: str,
) -> dict:
    """
    Persist a CoordinateBase trajectory to ``npz_path``.

    Stores ``xs``, ``vs``, ``taus`` arrays plus the coordinate class name and
    its ``kwargs`` (floats canonicalized via :func:`repr`) so that the
    trajectory can be reconstructed without the central manifest.

    Returns the small metadata dict so the caller can persist a copy in the
    manifest entry.
    """
    npz_path = Path(npz_path)
    npz_path.parent.mkdir(parents=True, exist_ok=True)

    xs = np.asarray(coord.xs, dtype=np.float64)
    vs = np.asarray(coord.vs, dtype=np.float64)
    taus_arr = np.asarray(taus, dtype=np.float64)
    coord_kwargs = getattr(coord, "kwargs", {}) or {}
    coord_kwargs_repr = {k: repr(float(v)) for k, v in coord_kwargs.items()}

    n_points = int(xs.shape[-1]) if xs.ndim > 1 else 1
    tau_span = (
        [float(taus_arr.flat[0]), float(taus_arr.flat[-1])]
        if taus_arr.size > 0
        else [0.0, 0.0]
    )
    meta = {
        "bank_format_version": BANK_FORMAT_VERSION,
        "relatipy_version": _RELATIPY_VERSION,
        "integrator": integrator,
        "coordinate_system": getattr(coord, "name_metric", type(coord).__name__),
        "coord_class": type(coord).__name__,
        "n_points": n_points,
        "tau_span": tau_span,
    }

    np.savez_compressed(
        npz_path,
        xs=xs,
        vs=vs,
        taus=taus_arr,
        coord_class=np.asarray(type(coord).__name__),
        coord_kwargs_json=np.asarray(json.dumps(coord_kwargs_repr, sort_keys=True)),
        meta=np.asarray(json.dumps(meta, sort_keys=True)),
    )
    return meta


def load_raw(npz_path: Path) -> dict:
    """Read raw arrays + meta from an NPZ; rebuild belongs to ``reconstruct``."""
    npz_path = Path(npz_path)
    with np.load(npz_path, allow_pickle=False) as z:
        out = {
            "xs": np.array(z["xs"]),
            "vs": np.array(z["vs"]),
            "taus": np.array(z["taus"]),
            "coord_class": str(z["coord_class"]),
            "coord_kwargs": json.loads(str(z["coord_kwargs_json"])),
            "meta": json.loads(str(z["meta"])),
        }
    return out
