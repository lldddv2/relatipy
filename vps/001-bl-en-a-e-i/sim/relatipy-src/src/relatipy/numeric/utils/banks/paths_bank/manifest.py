"""
JSON manifest for :class:`PathsBank` — schema, atomic I/O, in-memory model.
"""

from __future__ import annotations

import json
import os
import tempfile
import warnings
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from ._version import BANK_FORMAT_VERSION

try:
    from relatipy import __version__ as _RELATIPY_VERSION
except Exception:  # pragma: no cover
    _RELATIPY_VERSION = "unknown"


def _utcnow_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


@dataclass
class Manifest:
    """In-memory representation of ``manifest.json``."""

    bank_format_version: str = BANK_FORMAT_VERSION
    relatipy_version: str = _RELATIPY_VERSION
    created_at: str = field(default_factory=_utcnow_iso)
    updated_at: str = field(default_factory=_utcnow_iso)
    entries: dict[str, dict] = field(default_factory=dict)
    groups: dict[str, dict] = field(default_factory=dict)

    @classmethod
    def empty(cls) -> "Manifest":
        """Create a fresh manifest with current timestamps."""
        return cls()

    @classmethod
    def load(cls, path: Path) -> "Manifest":
        """Load a manifest JSON file, warning on format/version drift."""
        path = Path(path)
        # Detect orphan tmp from a crashed write
        tmp = path.with_suffix(path.suffix + ".tmp")
        if tmp.exists() and not path.exists():
            warnings.warn(
                f"Found orphan manifest tmp {tmp}; ignoring. "
                "Bank may have crashed mid-write.",
                stacklevel=2,
            )

        with path.open("r", encoding="utf-8") as f:
            data = json.load(f)

        m = cls(
            bank_format_version=data.get("bank_format_version", BANK_FORMAT_VERSION),
            relatipy_version=data.get("relatipy_version", "unknown"),
            created_at=data.get("created_at", _utcnow_iso()),
            updated_at=data.get("updated_at", _utcnow_iso()),
            entries=data.get("entries", {}),
            groups=data.get("groups", {}),
        )
        if m.bank_format_version != BANK_FORMAT_VERSION:
            warnings.warn(
                f"Bank format version mismatch: file={m.bank_format_version!r}, "
                f"runtime={BANK_FORMAT_VERSION!r}.",
                stacklevel=2,
            )
        if m.relatipy_version != _RELATIPY_VERSION:
            warnings.warn(
                f"relatipy version drift: bank={m.relatipy_version!r}, "
                f"runtime={_RELATIPY_VERSION!r}. Cached hashes still valid.",
                stacklevel=2,
            )
        return m

    def save_atomic(self, path: Path) -> None:
        """Serialize to ``path`` via ``<path>.tmp`` + :func:`os.replace`."""
        path = Path(path)
        self.updated_at = _utcnow_iso()
        payload = asdict(self)
        directory = path.parent
        directory.mkdir(parents=True, exist_ok=True)
        fd, tmp_str = tempfile.mkstemp(
            prefix=path.name + ".", suffix=".tmp", dir=str(directory)
        )
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as f:
                json.dump(payload, f, indent=2, sort_keys=False)
                f.flush()
                os.fsync(f.fileno())
            os.replace(tmp_str, path)
        except Exception:
            try:
                os.unlink(tmp_str)
            except OSError:
                pass
            raise

    def add_entry(self, h: str, entry: dict) -> None:
        self.entries[h] = entry

    def remove_entry(self, h: str) -> dict | None:
        return self.entries.pop(h, None)
