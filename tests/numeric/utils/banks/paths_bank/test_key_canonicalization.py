"""Tests for canonical key construction and hashing."""

import numpy as np
import pytest

from relatipy.numeric.coordinates import OrbitalElements
from relatipy.numeric.metrics import Kerr, Schwarzschild
from relatipy.numeric.utils.banks.paths_bank import key as K


def _ic():
    return OrbitalElements(a=50.0, e=0.2, mass=1.0)


def _kw(**over):
    base = {
        "integrator": "Yoshida6",
        "adaptative": True,
        "steps_per_period": 100,
        "rtol": 1e-12,
        "atol": 1e-13,
    }
    base.update(over)
    return base


def test_same_inputs_same_hash():
    m = Schwarzschild(mass=1.0)
    taus = np.linspace(0.0, 1.0, 10)
    h1 = K.hash_key(K.canonical_key(m, _ic(), taus, _kw()))
    h2 = K.hash_key(K.canonical_key(m, _ic(), taus, _kw()))
    assert h1 == h2


def test_different_mass_different_hash():
    taus = np.linspace(0.0, 1.0, 10)
    h1 = K.hash_key(K.canonical_key(Schwarzschild(1.0), _ic(), taus, _kw()))
    h2 = K.hash_key(K.canonical_key(Schwarzschild(2.0), _ic(), taus, _kw()))
    assert h1 != h2


def test_different_metric_type_different_hash():
    taus = np.linspace(0.0, 1.0, 10)
    h1 = K.hash_key(K.canonical_key(Schwarzschild(1.0), _ic(), taus, _kw()))
    h2 = K.hash_key(K.canonical_key(Kerr(1.0, 0.5), _ic(), taus, _kw()))
    assert h1 != h2


def test_taus_dtype_matters():
    m = Schwarzschild(1.0)
    t64 = np.linspace(0.0, 1.0, 10, dtype=np.float64)
    t32 = t64.astype(np.float32)
    h64 = K.hash_key(K.canonical_key(m, _ic(), t64, _kw()))
    h32 = K.hash_key(K.canonical_key(m, _ic(), t32, _kw()))
    assert h64 != h32


def test_taus_content_matters():
    m = Schwarzschild(1.0)
    h1 = K.hash_key(K.canonical_key(m, _ic(), np.linspace(0, 1, 10), _kw()))
    h2 = K.hash_key(K.canonical_key(m, _ic(), np.linspace(0, 1, 11), _kw()))
    assert h1 != h2


def test_integrator_default_collapse():
    """Omitting kwargs equals passing project defaults explicitly."""
    m = Schwarzschild(1.0)
    taus = np.linspace(0.0, 1.0, 10)
    h_empty = K.hash_key(K.canonical_key(m, _ic(), taus, {}))
    h_explicit = K.hash_key(
        K.canonical_key(
            m,
            _ic(),
            taus,
            {
                "integrator": "Radau",
                "adaptative": True,
                "steps_per_period": 100,
                "rtol": 1e-12,
                "atol": 1e-13,
            },
        )
    )
    assert h_empty == h_explicit


def test_rtol_changes_hash():
    m = Schwarzschild(1.0)
    taus = np.linspace(0.0, 1.0, 10)
    h1 = K.hash_key(K.canonical_key(m, _ic(), taus, _kw(rtol=1e-12)))
    h2 = K.hash_key(K.canonical_key(m, _ic(), taus, _kw(rtol=1e-11)))
    assert h1 != h2


def test_relatipy_version_excluded_from_hash():
    m = Schwarzschild(1.0)
    taus = np.linspace(0.0, 1.0, 10)
    k = K.canonical_key(m, _ic(), taus, _kw())
    assert "relatipy_version" in k
    h_orig = K.hash_key(k)
    k_bumped = dict(k); k_bumped["relatipy_version"] = "0.0.0-pre"
    assert K.hash_key(k_bumped) == h_orig
