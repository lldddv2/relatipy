"""
End-to-end smoke test:

- Build a real Schwarzschild geodesic with a small Yoshida6 run.
- Run twice through the bank.
- Assert the second call never reaches ``Geodesic.get_path``.
- Assert trajectories are bit-for-bit identical.
- Bundle the path in a group and re-read it.
"""

import astropy.units as u
import numpy as np

from relatipy.numeric.coordinates import OrbitalElements
from relatipy.numeric.geodesic.geodesic import Geodesic
from relatipy.numeric.metrics import Schwarzschild
from relatipy.numeric.utils.banks import PathsBank, init_bank


def test_round_trip_through_bank(tmp_path, monkeypatch):
    bank = init_bank(PathsBank(tmp_path / "bank"))
    M = 1.989e30 * u.kg
    metric = Schwarzschild(M)
    ic = OrbitalElements(a=50.0, e=0.2, mass=M)
    taus = np.linspace(0.0, 100.0, 120)
    kw = dict(integrator="Yoshida6", adaptative=False, rtol=1e-9, atol=1e-12)

    calls = {"n": 0}
    real = Geodesic.get_path

    def counted(self, *a, **kw_):
        calls["n"] += 1
        return real(self, *a, **kw_)

    monkeypatch.setattr(Geodesic, "get_path", counted)

    p1 = bank.get_path(metric, ic, taus, **kw)
    assert calls["n"] == 1
    p2 = bank.get_path(metric, ic, taus, **kw)
    assert calls["n"] == 1, "second call must hit the cache, not the integrator"

    assert np.array_equal(p1.xs, p2.xs)
    assert np.array_equal(p1.vs, p2.vs)

    g = bank.save_group("smoke", p1, notes="end-to-end")
    assert g.path_hashes == [p1._bank_hash]
    g2 = bank.get_group("smoke")
    assert g2.path_hashes == [p1._bank_hash]

    # Reloading the bank from disk should also see the entry.
    bank2 = init_bank(PathsBank(tmp_path / "bank"))
    assert p1._bank_hash in bank2.manifest.entries
    assert "smoke" in bank2.list_groups()
