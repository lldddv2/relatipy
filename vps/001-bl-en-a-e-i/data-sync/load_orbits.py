"""
Reader for the Kerr orbit bank produced by ``vps/001-bl-en-a-e-i/sim``.

Only requires numpy. The bank is a directory (``data/`` next to this file)
holding ``index.csv``, ``run.json`` and ``shards/orbits_NNNNN.npz``; see
``000-README.md`` for the layout and the physical conventions.

Examples
--------
>>> from load_orbits import OrbitBank
>>> bank = OrbitBank()                       # defaults to ./data
>>> bank.index["e"].max()                    # doctest: +SKIP
>>> orb = bank.get(1234)                     # doctest: +SKIP
>>> orb.r.shape, orb.a, orb.e, orb.inc_deg   # doctest: +SKIP
>>> x, y, z = orb.cartesian()                # doctest: +SKIP
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path

import numpy

__all__ = ["OrbitBank", "Orbit"]

_FLOAT_COLUMNS = {
    "a", "e", "inc_rad", "inc_deg", "r_p", "r_a", "period", "tau_max",
    "E", "Lz", "Q", "r_min", "r_max", "norm_err", "wall_s",
}
_INT_COLUMNS = {"orbit_id", "shard", "row", "n_samples", "bound", "status"}


@dataclass
class Orbit:
    """
    One trajectory of the bank in Boyer-Lindquist coordinates.

    Attributes
    ----------
    orbit_id : int
        Row of the Latin-Hypercube design that generated this orbit.
    a, e, inc_rad : float
        Keplerian initial conditions (semi-major axis in ``M``, eccentricity,
        inclination in radians). These are labels, not integration outputs.
    spin : float
        Dimensionless Kerr spin the orbit was integrated in.
    tau : ndarray, shape (n_samples,)
        Proper time, uniform from ``0`` to ``tau_max``.
    t, r, theta, phi : ndarray, shape (n_samples,)
        Boyer-Lindquist coordinates. ``t`` is float64, the rest float32.
    u : ndarray, shape (4, n_samples)
        Four-velocity ``dx^mu/dtau``, rows ordered ``(t, r, theta, phi)``.
    E, Lz, Q : float
        Specific energy, axial angular momentum and Carter constant.
    period : float
        Keplerian period ``2 pi a**1.5`` used to set ``tau_max``.
    bound : bool
        ``E < 1``. Read the note in ``000-README.md`` before using the
        unbound orbits.
    status : int
        ``0`` if the integration succeeded, ``1`` if the arrays are NaN.
    """

    orbit_id: int
    a: float
    e: float
    inc_rad: float
    spin: float
    tau: numpy.ndarray
    t: numpy.ndarray
    r: numpy.ndarray
    theta: numpy.ndarray
    phi: numpy.ndarray
    u: numpy.ndarray
    E: float
    Lz: float
    Q: float
    period: float
    bound: bool
    status: int

    @property
    def inc_deg(self) -> float:
        return float(numpy.degrees(self.inc_rad))

    def cartesian(self):
        """
        Pseudo-Cartesian embedding of the Boyer-Lindquist path, for plotting.

        Returns
        -------
        tuple of ndarray
            ``(x, y, z)`` with
            ``x = sqrt(r^2 + a_spin^2) sin(theta) cos(phi)``,
            ``y = sqrt(r^2 + a_spin^2) sin(theta) sin(phi)``,
            ``z = r cos(theta)``. This is the oblate-spheroidal embedding
            the Kerr horizon is a sphere in; it is a visualisation aid, not a
            flat-space position.
        """
        rho = numpy.sqrt(self.r.astype(float) ** 2 + self.spin ** 2)
        s = numpy.sin(self.theta.astype(float))
        return (rho * s * numpy.cos(self.phi),
                rho * s * numpy.sin(self.phi),
                self.r.astype(float) * numpy.cos(self.theta.astype(float)))


class OrbitBank:
    """
    Lazy accessor over the sharded orbit bank.

    Parameters
    ----------
    root : path-like, optional
        Directory holding ``index.csv`` and ``shards/``. Defaults to ``data``
        next to this module.

    Attributes
    ----------
    index : dict of str to ndarray
        Column-oriented view of ``index.csv``; every column has one entry per
        orbit, aligned with ``index["orbit_id"]``.
    run : dict
        Contents of ``run.json`` (configuration and conventions of the run).
    """

    def __init__(self, root=None):
        self.root = Path(root) if root is not None else Path(__file__).parent / "data"
        if not (self.root / "index.csv").is_file():
            raise FileNotFoundError(
                f"{self.root}/index.csv not found — has the sync run yet? "
                "See 000-README.md."
            )
        self.index = self._read_index(self.root / "index.csv")
        run_path = self.root / "run.json"
        self.run = json.loads(run_path.read_text(encoding="utf-8")) if run_path.is_file() else {}
        self._cache = {}
        self._cache_id = None

    @staticmethod
    def _read_index(path):
        with path.open("r", newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        cols = {}
        for name in rows[0] if rows else []:
            raw = [row[name] for row in rows]
            if name in _FLOAT_COLUMNS:
                cols[name] = numpy.array(raw, dtype=float)
            elif name in _INT_COLUMNS:
                cols[name] = numpy.array(raw, dtype=numpy.int64)
            else:
                cols[name] = numpy.array(raw, dtype=object)
        return cols

    def __len__(self):
        return len(self.index["orbit_id"])

    def select(self, **bounds):
        """
        Boolean mask over the index.

        Parameters
        ----------
        **bounds
            ``column=(lo, hi)`` for a range, or ``column=value`` for equality.
            Bounds are inclusive and ``None`` means open.

        Returns
        -------
        ndarray of bool
            Mask aligned with ``index``.

        Examples
        --------
        >>> mask = bank.select(bound=1, e=(0.5, 0.9), inc_deg=(0, 30))  # doctest: +SKIP
        >>> ids = bank.index["orbit_id"][mask]                          # doctest: +SKIP
        """
        mask = numpy.ones(len(self), dtype=bool)
        for name, spec in bounds.items():
            col = self.index[name]
            if isinstance(spec, tuple):
                lo, hi = spec
                if lo is not None:
                    mask &= col >= lo
                if hi is not None:
                    mask &= col <= hi
            else:
                mask &= col == spec
        return mask

    def shard(self, shard_id):
        """Load one shard ``.npz``, keeping only the most recent one in memory."""
        if self._cache_id != shard_id:
            path = self.root / "shards" / f"orbits_{int(shard_id):05d}.npz"
            with numpy.load(path) as z:
                self._cache = {k: z[k] for k in z.files}
            self._cache_id = shard_id
        return self._cache

    def get(self, orbit_id) -> Orbit:
        """Return the :class:`Orbit` with the given ``orbit_id``."""
        where = numpy.nonzero(self.index["orbit_id"] == int(orbit_id))[0]
        if where.size == 0:
            raise KeyError(f"orbit_id {orbit_id} is not in the index")
        k = int(where[0])
        z = self.shard(self.index["shard"][k])
        row = int(self.index["row"][k])

        ns = int(z["n_samples"])
        tau = numpy.linspace(0.0, float(z["tau_max"][row]), ns)
        xsp = z["xsp"][row]
        return Orbit(
            orbit_id=int(self.index["orbit_id"][k]),
            a=float(z["ci"][row, 0]),
            e=float(z["ci"][row, 1]),
            inc_rad=float(z["ci"][row, 2]),
            spin=float(z["spin"]),
            tau=tau,
            t=z["t"][row],
            r=xsp[0], theta=xsp[1], phi=xsp[2],
            u=z["u"][row],
            E=float(z["constants"][row, 0]),
            Lz=float(z["constants"][row, 1]),
            Q=float(z["constants"][row, 2]),
            period=float(z["period"][row]),
            bound=bool(z["bound"][row]),
            status=int(z["status"][row]),
        )

    def iter_orbits(self, orbit_ids=None):
        """
        Iterate orbits shard by shard, so each ``.npz`` is opened once.

        Parameters
        ----------
        orbit_ids : array_like of int, optional
            Subset to iterate. Defaults to every orbit in the index.

        Yields
        ------
        Orbit
        """
        ids = (self.index["orbit_id"] if orbit_ids is None
               else numpy.asarray(orbit_ids, dtype=numpy.int64))
        order = numpy.argsort([self.index["shard"][
            int(numpy.nonzero(self.index["orbit_id"] == i)[0][0])] for i in ids])
        for i in numpy.asarray(ids)[order]:
            yield self.get(int(i))


if __name__ == "__main__":
    bank = OrbitBank()
    idx = bank.index
    ok = idx["status"] == 0
    print(f"bank: {len(bank)} orbits, {int(ok.sum())} integrated, "
          f"{int(idx['bound'].sum())} bound (E < 1)")
    print(f"  a   [{idx['a'].min():.1f}, {idx['a'].max():.1f}] M")
    print(f"  e   [{idx['e'].min():.4f}, {idx['e'].max():.4f}]")
    print(f"  inc [{idx['inc_deg'].min():.2f}, {idx['inc_deg'].max():.2f}] deg")
    print(f"  r_p [{idx['r_p'].min():.2f}, {idx['r_p'].max():.1f}] M")
    print(f"  max |g_uu u^u u^v - 1| = {idx['norm_err'][ok].max():.2e}")
