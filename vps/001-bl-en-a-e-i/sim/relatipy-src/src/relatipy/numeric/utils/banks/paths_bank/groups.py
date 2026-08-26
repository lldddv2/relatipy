"""
Named groups of cached paths.

A :class:`Group` is a lightweight bundle of cache hashes plus optional
notes — a way to say "these paths belong together for this analysis".
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone


def _utcnow_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


@dataclass
class Group:
    """A named bundle of cached path hashes."""

    name: str
    created_at: str = field(default_factory=_utcnow_iso)
    updated_at: str = field(default_factory=_utcnow_iso)
    path_hashes: list[str] = field(default_factory=list)
    notes: str | None = None

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, name: str, data: dict) -> "Group":
        return cls(
            name=name,
            created_at=data.get("created_at", _utcnow_iso()),
            updated_at=data.get("updated_at", _utcnow_iso()),
            path_hashes=list(data.get("path_hashes", [])),
            notes=data.get("notes"),
        )

    def touch(self) -> None:
        self.updated_at = _utcnow_iso()

    def add(self, hashes: list[str]) -> None:
        for h in hashes:
            if h not in self.path_hashes:
                self.path_hashes.append(h)
        self.touch()

    def remove(self, hashes: list[str]) -> None:
        self.path_hashes = [h for h in self.path_hashes if h not in hashes]
        self.touch()

    def __len__(self) -> int:
        return len(self.path_hashes)
