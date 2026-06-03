"""
:class:`PathsBank` — content-addressed disk cache for geodesic trajectories.

A bank lives in a directory containing ``manifest.json`` and ``paths/<hash>.npz``
files. The bank does NOT monkey-patch :class:`Geodesic`: the user calls
``bank.get_path(metric, ic, taus, ...)`` instead of ``metric.geodesic.get_path(
ic, taus, ...)``. On a cache miss the bank delegates to the real
``Geodesic.get_path``, persists the result, and returns it; on a hit it loads
the trajectory from disk.

Returned :class:`CoordinateBase` instances carry their cache hash as the
``_bank_hash`` attribute, so they can be passed directly to group-management
methods.
"""

from __future__ import annotations

import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

from . import key as _key
from . import storage
from .groups import Group
from .manifest import Manifest
from .reconstruct import BankError, BankIntegrityError, rebuild


def _utcnow_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _resolve_hashes(items: Iterable[Any]) -> list[str]:
    """Turn a mixed iterable of hash strings / cached objects into hashes."""
    out: list[str] = []
    for it in items:
        if isinstance(it, str):
            out.append(it)
        elif hasattr(it, "_bank_hash"):
            out.append(getattr(it, "_bank_hash"))
        else:
            raise ValueError(
                "Item is neither a hash string nor a bank-produced object "
                "(missing _bank_hash attribute)."
            )
    return out


class PathsBank:
    """
    Disk-backed cache for geodesic trajectories.

    Parameters
    ----------
    path : str or Path
        Directory holding the manifest and per-path NPZ files. Created on
        :meth:`init` if missing.
    """

    MANIFEST_FILENAME = "manifest.json"
    PATHS_SUBDIR = "paths"

    def __init__(self, path: str | os.PathLike) -> None:
        self.path = Path(path)
        self._manifest: Manifest | None = None

    # ------------------------------------------------------------------ init

    @property
    def manifest_path(self) -> Path:
        return self.path / self.MANIFEST_FILENAME

    @property
    def paths_dir(self) -> Path:
        return self.path / self.PATHS_SUBDIR

    @property
    def manifest(self) -> Manifest:
        if self._manifest is None:
            raise BankError(
                f"Bank at {self.path} not initialized; call init_bank(bank) first."
            )
        return self._manifest

    def init(self) -> "PathsBank":
        """Create or load the bank directory. Idempotent."""
        self.path.mkdir(parents=True, exist_ok=True)
        self.paths_dir.mkdir(parents=True, exist_ok=True)
        if self.manifest_path.exists():
            self._manifest = Manifest.load(self.manifest_path)
        else:
            self._manifest = Manifest.empty()
            self._manifest.save_atomic(self.manifest_path)
        return self

    def _save_manifest(self) -> None:
        self.manifest.save_atomic(self.manifest_path)

    # ------------------------------------------------------------ get_path

    def hash_of(
        self,
        metric,
        initial_conditions,
        taus,
        integrator: str = "Radau",
        adaptative: bool = True,
        steps_per_period: int = 100,
        rtol: float = 1e-12,
        atol: float = 1e-13,
    ) -> str:
        """SHA-256 hex digest identifying this configuration."""
        k = _key.canonical_key(
            metric,
            initial_conditions,
            taus,
            {
                "integrator": integrator,
                "adaptative": adaptative,
                "steps_per_period": steps_per_period,
                "rtol": rtol,
                "atol": atol,
            },
        )
        return _key.hash_key(k)

    def has(
        self,
        metric,
        initial_conditions,
        taus,
        **kwargs,
    ) -> bool:
        """Return True if a path for this configuration is cached."""
        h = self.hash_of(metric, initial_conditions, taus, **kwargs)
        return h in self.manifest.entries

    def get_path(
        self,
        metric,
        initial_conditions,
        taus,
        integrator: str = "Radau",
        adaptative: bool = True,
        steps_per_period: int = 100,
        rtol: float = 1e-12,
        atol: float = 1e-13,
    ):
        """
        Return the geodesic for this configuration, loading from cache on hit.

        On a miss, delegates to ``metric.geodesic.get_path(initial_conditions,
        taus, **integrator_kwargs)``, persists the trajectory, and returns it.
        The returned object carries its cache hash on ``._bank_hash``.
        """
        integrator_kwargs = {
            "integrator": integrator,
            "adaptative": adaptative,
            "steps_per_period": steps_per_period,
            "rtol": rtol,
            "atol": atol,
        }
        key = _key.canonical_key(metric, initial_conditions, taus, integrator_kwargs)
        h = _key.hash_key(key)

        npz_path = self.paths_dir / f"{h}.npz"

        if h in self.manifest.entries:
            if not npz_path.exists():
                raise BankIntegrityError(
                    f"Manifest references {h} but NPZ {npz_path} is missing. "
                    "Run bank.verify(prune_missing=True) to repair."
                )
            coord = rebuild(npz_path)
            coord._bank_hash = h
            return coord

        # Miss: integrate via the real Geodesic.
        coord = metric.geodesic.get_path(
            initial_conditions,
            taus,
            **integrator_kwargs,
        )

        meta = storage.save(npz_path, coord, taus, integrator=integrator)
        entry = {
            "key": key,
            "npz_filename": f"{self.PATHS_SUBDIR}/{h}.npz",
            "created_at": _utcnow_iso(),
            "integrator": integrator,
            "n_points": meta["n_points"],
            "tau_span": meta["tau_span"],
            "coordinate_system": meta["coordinate_system"],
            "relatipy_version": meta["relatipy_version"],
        }
        self.manifest.add_entry(h, entry)
        self._save_manifest()

        coord._bank_hash = h
        return coord

    # --------------------------------------------------------- introspection

    def list_paths(self) -> list[dict]:
        """Lightweight summary of every cached path."""
        out: list[dict] = []
        for h, e in self.manifest.entries.items():
            out.append(
                {
                    "hash": h,
                    "integrator": e.get("integrator"),
                    "n_points": e.get("n_points"),
                    "tau_span": e.get("tau_span"),
                    "coordinate_system": e.get("coordinate_system"),
                    "created_at": e.get("created_at"),
                }
            )
        return out

    def get_path_by_hash(self, h: str):
        if h not in self.manifest.entries:
            raise KeyError(f"No entry with hash {h!r} in bank.")
        npz_path = self.paths_dir / f"{h}.npz"
        coord = rebuild(npz_path)
        coord._bank_hash = h
        return coord

    def delete_path(self, h: str) -> None:
        """Remove a cached path; also detaches it from every group."""
        if h not in self.manifest.entries:
            raise KeyError(f"No entry with hash {h!r} in bank.")
        npz_path = self.paths_dir / f"{h}.npz"
        try:
            npz_path.unlink()
        except FileNotFoundError:
            pass
        self.manifest.remove_entry(h)
        for gname, gdict in self.manifest.groups.items():
            hashes = gdict.get("path_hashes", [])
            if h in hashes:
                gdict["path_hashes"] = [x for x in hashes if x != h]
                gdict["updated_at"] = _utcnow_iso()
        self._save_manifest()

    def verify(self, prune_missing: bool = False) -> dict:
        """
        Scan disk for missing NPZ files and orphan files; optionally prune.

        Returns a report ``{"missing": [...], "orphan_files": [...]}``.
        """
        missing: list[str] = []
        for h, e in list(self.manifest.entries.items()):
            npz_path = self.paths_dir / f"{h}.npz"
            if not npz_path.exists():
                missing.append(h)
                if prune_missing:
                    self.delete_path_safe(h)

        known = {f"{h}.npz" for h in self.manifest.entries}
        orphans: list[str] = []
        if self.paths_dir.exists():
            for f in self.paths_dir.iterdir():
                if f.is_file() and f.name not in known:
                    orphans.append(f.name)

        if prune_missing:
            self._save_manifest()

        return {"missing": missing, "orphan_files": orphans}

    def delete_path_safe(self, h: str) -> None:
        """Like :meth:`delete_path` but tolerates missing entries silently."""
        if h not in self.manifest.entries:
            return
        npz_path = self.paths_dir / f"{h}.npz"
        try:
            npz_path.unlink()
        except FileNotFoundError:
            pass
        self.manifest.remove_entry(h)
        for gdict in self.manifest.groups.values():
            hashes = gdict.get("path_hashes", [])
            if h in hashes:
                gdict["path_hashes"] = [x for x in hashes if x != h]
                gdict["updated_at"] = _utcnow_iso()

    # ----------------------------------------------------------------- groups

    def _require_known_hash(self, h: str) -> None:
        if h not in self.manifest.entries:
            raise KeyError(
                f"Hash {h!r} is not in this bank; cannot add to a group."
            )

    def save_group(
        self,
        name: str,
        *paths_or_hashes,
        notes: str | None = None,
    ) -> Group:
        """Create (or overwrite) a named group from the given hashes/objects."""
        if name in self.manifest.groups:
            raise ValueError(
                f"Group {name!r} already exists. Use add_to_group / "
                "remove_from_group, or delete it first."
            )
        hashes = _resolve_hashes(paths_or_hashes)
        for h in hashes:
            self._require_known_hash(h)
        g = Group(name=name, path_hashes=hashes, notes=notes)
        self.manifest.groups[name] = {
            "created_at": g.created_at,
            "updated_at": g.updated_at,
            "path_hashes": list(g.path_hashes),
            "notes": g.notes,
        }
        self._save_manifest()
        return g

    def list_groups(self) -> list[str]:
        return list(self.manifest.groups)

    def get_group(self, name: str) -> Group:
        if name not in self.manifest.groups:
            raise KeyError(f"No group named {name!r}.")
        return Group.from_dict(name, self.manifest.groups[name])

    def show_groups(self) -> None:
        """Pretty-print groups to stdout."""
        if not self.manifest.groups:
            print("(no groups)")
            return
        for name, g in self.manifest.groups.items():
            count = len(g.get("path_hashes", []))
            notes = g.get("notes") or ""
            print(f"- {name}  ({count} paths)  {notes}")
            for h in g.get("path_hashes", []):
                print(f"    * {h}")

    def delete_group(self, name: str) -> None:
        if name not in self.manifest.groups:
            raise KeyError(f"No group named {name!r}.")
        del self.manifest.groups[name]
        self._save_manifest()

    def rename_group(self, old: str, new: str) -> None:
        if old not in self.manifest.groups:
            raise KeyError(f"No group named {old!r}.")
        if new in self.manifest.groups:
            raise ValueError(f"Group {new!r} already exists; refusing to overwrite.")
        self.manifest.groups[new] = self.manifest.groups.pop(old)
        self.manifest.groups[new]["updated_at"] = _utcnow_iso()
        self._save_manifest()

    def add_to_group(self, name: str, *paths_or_hashes) -> Group:
        if name not in self.manifest.groups:
            raise KeyError(f"No group named {name!r}.")
        hashes = _resolve_hashes(paths_or_hashes)
        for h in hashes:
            self._require_known_hash(h)
        gdict = self.manifest.groups[name]
        existing = list(gdict.get("path_hashes", []))
        for h in hashes:
            if h not in existing:
                existing.append(h)
        gdict["path_hashes"] = existing
        gdict["updated_at"] = _utcnow_iso()
        self._save_manifest()
        return Group.from_dict(name, gdict)

    def remove_from_group(self, name: str, *paths_or_hashes) -> Group:
        if name not in self.manifest.groups:
            raise KeyError(f"No group named {name!r}.")
        targets = set(_resolve_hashes(paths_or_hashes))
        gdict = self.manifest.groups[name]
        gdict["path_hashes"] = [
            h for h in gdict.get("path_hashes", []) if h not in targets
        ]
        gdict["updated_at"] = _utcnow_iso()
        self._save_manifest()
        return Group.from_dict(name, gdict)

    def __repr__(self) -> str:
        n_paths = len(self._manifest.entries) if self._manifest is not None else "?"
        n_groups = len(self._manifest.groups) if self._manifest is not None else "?"
        return f"PathsBank(path={str(self.path)!r}, paths={n_paths}, groups={n_groups})"


def init_bank(bank: PathsBank) -> PathsBank:
    """Module-level convenience: ``init_bank(bank)`` == ``bank.init()``."""
    return bank.init()
