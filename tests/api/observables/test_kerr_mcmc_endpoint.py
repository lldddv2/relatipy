"""Regression tests for integrations whose last step ends one ulp short."""

from __future__ import annotations

import astropy.units as u
import numpy as np
import pytest

from relatipy import KerrMcmcModel

# S2 osculating at 2010.35 (GRAVITY Collaboration 2024, Table B.1). With
# dop853 at rtol 1e-11 the step sequence to these epochs stops one ulp before
# the endpoint, which the integrator used to report as a step underflow.
S2_CALL = dict(
    kerr=dict(mass=4.2996e6 * u.Msun, spin=0.7, vec=(0.0, 0.0, 1.0)),
    orbit=dict(
        a=0.12502 * 8.2759e3 * u.au,
        e=0.88444,
        inc=134.67 * u.deg,
        Omega=228.21 * u.deg,
        omega=66.279 * u.deg,
        f=3.1413949029953905 * u.rad,
        t=2010.35 * u.yr,
    ),
    distance=8.2759 * u.kpc,
)


@pytest.mark.parametrize("epoch", [2010.94, 2010.9471, 2010.9470999999608])
def test_final_sliver_step_does_not_fail(epoch: float) -> None:
    model = KerrMcmcModel()
    solutions = [
        model.get_ra_dec_vr(
            sol=dict(
                t_obs=np.array([epoch]) * u.yr,
                method="dop853",
                rtol=rtol,
                atol=rtol * 1e-2,
            ),
            **S2_CALL,
        )
        for rtol in (1e-11, 1e-12)
    ]
    for observable, reference in zip(*solutions):
        np.testing.assert_allclose(
            np.asarray(observable), np.asarray(reference), rtol=1e-9
        )
