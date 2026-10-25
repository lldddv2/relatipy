"""Evaluate several stars around one Kerr black hole with an external force.

:class:`PerturbedKerrMcmcMultiModel` is the perturbed counterpart of
:class:`~relatipy.observables.KerrMcmcMultiModel`: the same shared, spin and
per-star blocks, followed by one force block. All stars share one bound
force whose parameters are rewritten in place for each proposal.
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from ..forces import ExternalForce
from ..forces._bound import BoundForce
from .kerr_mcmc import _spin_frame
from .kerr_mcmc_multi import _MultiStarBase
from .perturbed_kerr_mcmc import (
    _KG_PER_MILLION_SOLAR_MASS,
    PerturbedKerrMcmcModel,
    _require_force,
)


class PerturbedKerrMcmcMultiModel(_MultiStarBase):
    """Predict astrometry and redshift of several stars under one force.

    Parameters
    ----------
    astrometry_times, spectroscopy_times : mapping
        Per-star epochs, as for
        :class:`~relatipy.observables.KerrMcmcMultiModel`.
    reference_epoch : float or astropy.units.Quantity
        Shared zero point of the linear astrometric frame drift.
    force : relatipy.forces.ExternalForce
        Force felt by every star. Its parameters form the force block and
        its directions are read in observer axes ``(Dec, RA, away)``.

    Raises
    ------
    TypeError
        If ``force`` is not an external force, or as ``KerrMcmcMultiModel``.
    ValueError
        As ``KerrMcmcMultiModel``.

    Notes
    -----
    The parameter vector is the layout of ``KerrMcmcMultiModel`` (shared
    block, spin block, one block per star) followed by the force block
    ``force.<name>`` for each of the force's
    :attr:`~relatipy.forces.ExternalForce.parameter_names`, in its
    :attr:`~relatipy.forces.ExternalForce.parameter_units`. For example a
    :class:`~relatipy.forces.Plummer` sphere adds ``force.mass`` and
    ``force.radius``.

    Each proposal binds nothing anew: the single force handle shared by all
    stars is updated in place with the proposal mass, force block and spin
    orientation, then each star is evaluated as a
    :class:`~relatipy.observables.PerturbedKerrMcmcModel` legacy call. With
    an identically zero force the predictions equal those of
    ``KerrMcmcMultiModel`` bit for bit.
    :class:`~relatipy.observables.KerrMcmcGaussianLikelihood` accepts this
    model. Instances are picklable and can be sent to multiprocessing
    pools; each process binds its own handle on first use.
    """

    __slots__ = ("_models", "_solver_set", "_force", "_bound", "_force_slice")

    def __init__(
        self,
        astrometry_times: Mapping[str, object],
        spectroscopy_times: Mapping[str, object],
        *,
        reference_epoch: object,
        force: ExternalForce,
    ) -> None:
        """Validate per-star epochs and the force, and build one model per star."""
        force = _require_force(force)
        super().__init__(
            astrometry_times, spectroscopy_times, reference_epoch,
            extra_blocks=(("force", force._mcmc_parameters()),),
        )
        self._models = tuple(
            PerturbedKerrMcmcModel(
                self._astrometry_times[star],
                self._spectroscopy_times[star],
                reference_epoch=self._reference_epoch,
                force=force,
            )
            for star in self._stars
        )
        self._solver_set = False
        self._force = force
        self._bound: BoundForce | None = None
        self._force_slice = self._layout._slices["force"]

    def __getstate__(self) -> dict:
        """Return picklable state without the native force handle."""
        return {
            name: getattr(self, name)
            for cls in type(self).__mro__
            for name in getattr(cls, "__slots__", ())
            if name != "_bound"
        }

    def __setstate__(self, state: dict) -> None:
        """Restore a pickled model; the force is bound again on first use."""
        for name, value in state.items():
            setattr(self, name, value)
        self._bound = None

    @property
    def force(self) -> ExternalForce:
        """Force felt by every star; its parameters form the force block."""
        return self._force

    def set_solver(
        self,
        *,
        method: str,
        rtol: float,
        atol: float,
        max_steps: int = 100000,
    ) -> None:
        """Set the solver of every star, as :meth:`KerrMcmcMultiModel.set_solver`.

        Parameters
        ----------
        method : {'dop853', 'radau'}
            Native integration method.
        rtol, atol : float
            Positive scalar relative and absolute tolerances.
        max_steps : int, optional
            Maximum steps per endpoint call, default ``100000``.

        Raises
        ------
        TypeError, ValueError
            As :meth:`KerrMcmcModel.set_solver`. Invalid settings leave every
            star unchanged.
        """
        for model in self._models:
            model.set_solver(method=method, rtol=rtol, atol=atol, max_steps=max_steps)
        self._solver_set = True

    def _checked(self, params: object) -> np.ndarray:
        if not self._solver_set:
            raise ValueError("call set_solver before evaluating the model")
        values = super()._checked(params)
        self._force._check_vector(values[self._force_slice])
        if self._bound is None:
            # Bind here, outside the likelihood's error handling, so that a
            # missing compiler or an invalid definition raises.
            self._bound = self._force._handle()
        return values

    def _predict(
        self, values: np.ndarray
    ) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
        spin, rotation = _spin_frame(values[self._spin_slice])
        bound = self._bound
        if bound is None:
            bound = self._force._handle()
            self._bound = bound
        bound.update(
            values[1] * _KG_PER_MILLION_SOLAR_MASS, values[self._force_slice], rotation
        )
        rows = values[self._gather]
        return {
            star: model._evaluate_legacy(row, spin, rotation, bound)
            for star, model, row in zip(self._stars, self._models, rows)
        }

    def __call__(
        self, params: object
    ) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """Return each star's ``(alpha, delta, v_los)`` for one proposal.

        Parameters
        ----------
        params : array-like
            :attr:`n_parameters` finite numbers in the layout of the class
            notes.

        Returns
        -------
        dict
            Maps each star, in :attr:`stars` order, to read-only arrays:
            ``alpha`` and ``delta`` in arcsec at its astrometric epochs and
            ``v_los`` in km/s at its spectroscopic epochs.

        Raises
        ------
        ValueError
            As :meth:`KerrMcmcMultiModel.__call__`, or if a force parameter
            is outside the force's domain (for example a negative Plummer
            mass).
        relatipy.geodesic.IntegrationError
            If integration, the force callback or arrival matching fails.
        ImportError, RuntimeError
            If the extension is unavailable, or a run-time force needs SymPy
            or a C compiler that is missing.
        """
        return self._predict(self._checked(params))
