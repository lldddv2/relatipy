"""Define structured exceptions for orbital integration operations.

These exception types carry integration failures and structured termination data.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from .integration import Termination

if TYPE_CHECKING:
    from astropy import units as u

    from .state import State


class IntegrationError(RuntimeError):
    """Report a numerical failure during an orbit integration.

    Parameters
    ----------
    *args : object
        Positional arguments forwarded to :class:`RuntimeError`, normally a
        human-readable failure message.

    Examples
    --------
    >>> error = IntegrationError("numerical failure")
    >>> str(error)
    'numerical failure'
    """


class IntegrationWarning(UserWarning):
    """Warn about integration settings that likely suit the orbit poorly.

    Emitted with :func:`warnings.warn`; it never stops an integration and
    can be filtered with :func:`warnings.filterwarnings`.

    Examples
    --------
    >>> import warnings
    >>> warnings.simplefilter("ignore", IntegrationWarning)
    """


class IntegrationTerminated(RuntimeError):
    """Report an early, non-failing termination by an internal event.

    Parameters
    ----------
    reason : str
        Stable identifier for the internal termination condition.
    tau : astropy.units.Quantity
        Scalar proper time of the last valid exterior state.  It
        must equal ``state.tau`` after unit conversion.
    state : State
        Last valid physical state supplied by the caller.
    message : str or None, optional
        Optional human-readable description.  If omitted, ``reason`` is used.

    Attributes
    ----------
    termination : Termination
        Validated terminal-event record.
    reason : str
        Stable identifier for the terminal condition.
    tau : astropy.units.Quantity
        Read-only proper time of the last valid state.
    state : State
        Last valid scalar state.

    Raises
    ------
    TypeError, astropy.units.UnitConversionError, ValueError
        If the arguments are rejected by :class:`Termination`.
    """

    __slots__ = ("_termination",)

    def __init__(
        self,
        reason: str,
        tau: u.Quantity,
        state: State,
        message: str | None = None,
    ) -> None:
        """Validate termination data and initialize the runtime error."""
        self._termination = Termination(reason, tau, state)
        super().__init__(message if message is not None else reason)

    @property
    def termination(self) -> Termination:
        """Return the validated terminal-event record.

        Returns
        -------
        Termination
            Immutable record holding :attr:`reason`, :attr:`tau`, and
            :attr:`state`.
        """
        return self._termination

    @property
    def reason(self) -> str:
        """Return the stable identifier for the terminal event.

        Returns
        -------
        str
            Event identifier supplied at construction.
        """
        return self._termination.reason

    @property
    def tau(self) -> u.Quantity:
        """Return the proper time of the last valid physical state.

        Returns
        -------
        astropy.units.Quantity
            Read-only proper time of the last valid state in the input unit.
        """
        return self._termination.tau

    @property
    def state(self) -> State:
        """Return the last valid physical state supplied at construction.

        Returns
        -------
        State
            Immutable scalar state at :attr:`tau`.
        """
        return self._termination.state
