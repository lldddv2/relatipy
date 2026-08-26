"""
Rebuild a :class:`CoordinateBase` instance from a stored NPZ.

Reconstruction is mechanical: look up the class name in the
``coordinate_systems`` registry, decode the saved ``kwargs`` (floats stored as
``repr``-strings), and reinstantiate via ``cls(xs=..., vels=..., from_dxs_dt=False,
**kwargs)``. The proper-time grid ``taus`` is attached as an attribute so the
returned object carries everything needed for downstream analysis.
"""

from __future__ import annotations

from pathlib import Path

from . import storage


class BankError(Exception):
    """Base class for :class:`PathsBank` errors."""


class BankIntegrityError(BankError):
    """Raised when on-disk state does not match the manifest."""


def rebuild(npz_path: Path):
    """Rebuild the CoordinateBase trajectory stored in ``npz_path``."""
    # Lazy import to avoid utils -> numeric circular imports at module load.
    from relatipy.numeric.coordinates import coordinate_systems

    npz_path = Path(npz_path)
    if not npz_path.exists():
        raise BankIntegrityError(f"NPZ not found: {npz_path}")

    raw = storage.load_raw(npz_path)
    cls_name = raw["coord_class"]
    if cls_name not in coordinate_systems:
        raise BankIntegrityError(
            f"Unknown coordinate class {cls_name!r} in {npz_path}. "
            f"Known: {sorted(coordinate_systems)}"
        )
    cls = coordinate_systems[cls_name]
    kwargs = {k: float(v) for k, v in raw["coord_kwargs"].items()}

    inst = cls(xs=raw["xs"], vels=raw["vs"], from_dxs_dt=False, **kwargs)
    inst.taus = raw["taus"]
    return inst
