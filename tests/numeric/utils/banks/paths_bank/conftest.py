"""Shared fixtures for PathsBank tests."""

import astropy.units as u
import numpy as np
import pytest

from relatipy.numeric.coordinates import OrbitalElements
from relatipy.numeric.metrics import Schwarzschild
from relatipy.numeric.utils.banks import PathsBank, init_bank


_M_SUN = 1.989e30 * u.kg


@pytest.fixture
def schw():
    """Solar-mass Schwarzschild black hole."""
    return Schwarzschild(_M_SUN)


@pytest.fixture
def ic_orbital():
    """Slightly eccentric orbit around a solar-mass BH."""
    return OrbitalElements(a=50.0, e=0.2, inc=0.0, mass=_M_SUN)


@pytest.fixture
def taus():
    return np.linspace(0.0, 100.0, 80)


@pytest.fixture
def bank(tmp_path):
    b = PathsBank(tmp_path / "bank")
    return init_bank(b)
