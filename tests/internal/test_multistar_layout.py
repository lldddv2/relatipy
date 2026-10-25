"""Checks for the private multi-star layout reused by perturbed models."""

import pickle

import numpy as np
import pytest
from astropy import units as u

from relatipy.observables.kerr_mcmc_multi import _MultiStarBase, _ParameterLayout


FORCE = ("force", (("log_M_disk", u.one), ("disk_normal_x", u.one)))


class _ForcedModel(_MultiStarBase):
    __slots__ = ()

    def __init__(self, astrometry, spectroscopy, reference_epoch):
        super().__init__(astrometry, spectroscopy, reference_epoch, (FORCE,))

    def _predict(self, values):
        return {"rows": values[self._gather], "force": values[self.parameter_slices["force"]]}


def test_extra_blocks_are_appended_without_moving_kerr_blocks():
    model = _ForcedModel({"S2": [2010.0], "S301": [2020.0]}, {"S2": [], "S301": []}, 2000.0)
    assert model.parameter_names[-2:] == ("log_M_disk", "disk_normal_x")
    assert model.parameter_slices["force"] == slice(22, 24)
    assert model.parameter_slices["S301"] == slice(16, 22)
    vector = np.arange(1.0, 25.0) * 1e-2
    vector[0:2] = (8.0, 4.0)
    vector[[12, 18]] = 0.5
    result = model(vector)
    np.testing.assert_array_equal(result["force"], vector[22:])
    np.testing.assert_array_equal(
        result["rows"][1], vector[[0, 1, 16, 17, 18, 19, 20, 21, 2, 3, 4, 5, 6]]
    )
    restored = pickle.loads(pickle.dumps(model))
    assert restored.parameter_names == model.parameter_names


def test_layout_rejects_duplicates_and_block_collisions():
    with pytest.raises(ValueError, match="duplicate parameter block"):
        _ParameterLayout((("a", (("x", u.one),)), ("a", (("y", u.one),))))
    with pytest.raises(ValueError, match="unique"):
        _ParameterLayout((("a", (("x", u.one),)), ("b", (("x", u.one),))))
    with pytest.raises(ValueError, match="collides"):
        _MultiStarBase({"disk": [2010.0]}, {"disk": []}, 2000.0,
                       (("disk", (("m", u.one),)),))
    with pytest.raises(ValueError, match="reserved"):
        _ForcedModel({"force": [2010.0]}, {"force": []}, 2000.0)
