"""
Cartesian spacetime coordinates.

This module provides :class:`Cartesian`, a concrete implementation of
:class:`~relatipy.numeric.coordinates.base.CoordinateBase` for the standard
Cartesian chart :math:`x^\\mu = (t, x^1, x^2, x^3)`. In this chart the stored
spatial velocities :math:`v^i` coincide component-wise with the coordinate
derivatives :math:`\\mathrm{d}x^i/\\mathrm{d}\\tau`, so the hooks
:meth:`Cartesian._get_vs_from_dxs_dt` and :meth:`Cartesian._get_dxs_dt_from_vs`
implement the identity map on the three spatial components.

See Also
--------
relatipy.numeric.coordinates.base.CoordinateBase :
    Abstract base class defining the coordinate state vector and conversions.

Examples
--------
>>> import numpy as np
>>> from relatipy.numeric.coordinates.cartesian import Cartesian
>>> xs = np.array([0.0, 1.0, 2.0, 3.0])
>>> c = Cartesian(xs, vels=np.array([0.1, 0.0, 0.0]), from_dxs_dt=False)
>>> np.allclose(c.dxs_dt, c.vs)
True
"""

import numpy
from .base import CoordinateBase
from ..constants import _c


class Cartesian(CoordinateBase):
    """
    Cartesian coordinates :math:`(t, x^1, x^2, x^3)` with spatial state.

    Spatial velocities :math:`v^i` and derivatives
    :math:`\\mathrm{d}x^i/\\mathrm{d}\\tau` are represented by the same
    three-component arrays in this chart.

    Parameters
    ----------
    xs : array_like of shape (4,)
        Event :math:`(x^0, x^1, x^2, x^3)` with :math:`x^0` typically time.
    vels : array_like of shape (3,) or None, optional
        Either :math:`v^i` (``from_dxs_dt=False``) or
        :math:`\\mathrm{d}x^i/\\mathrm{d}\\tau` (``from_dxs_dt=True``).
        If None, both are zeros like the spatial part of ``xs``.
    from_dxs_dt : bool, optional
        If True, ``vels`` is interpreted as :math:`\\mathrm{d}x^i/\\mathrm{d}\\tau`
        and :math:`v^i` is filled via :meth:`_get_vs_from_dxs_dt`. If False,
        ``vels`` is :math:`v^i` and derivatives come from
        :meth:`_get_dxs_dt_from_vs`.

    Attributes
    ----------
    name_metric : str
        Set to ``"Cartesian"``.
    xs : ndarray
        Four-position, shape ``(4,)``.
    vs : ndarray
        Spatial velocities :math:`v^i`, shape ``(3,)``.
    dxs_dt : ndarray
        :math:`\\mathrm{d}x^i/\\mathrm{d}\\tau`, same shape as ``vs``; equal
        to ``vs`` component-wise in this chart.
    state_vector : ndarray
        Concatenation ``(x^\\mu, v^i)``, length 7.

    Raises
    ------
    ValueError
        If ``xs`` does not have length 4 (from the base class).

    Examples
    --------
    >>> import numpy as np
    >>> from relatipy.numeric.coordinates.cartesian import Cartesian
    >>> c = Cartesian([0.0, 0.0, 0.0, 0.0], vels=[1.0, 0.0, 0.0])
    >>> c.convert_to_cartesian() is c
    True
    """

    def __init__(self, xs, vels=None, from_dxs_dt=False):
        """
        Construct Cartesian coordinates and velocity state.

        Parameters
        ----------
        xs : array_like of shape (4,)
            Four-position :math:`(x^0, x^1, x^2, x^3)`.
        vels : array_like of shape (3,) or None, optional
            Velocities or coordinate derivatives, depending on ``from_dxs_dt``.
        from_dxs_dt : bool, optional
            Whether ``vels`` stores :math:`\\mathrm{d}x^i/\\mathrm{d}\\tau`.

        Examples
        --------
        >>> import numpy as np
        >>> from relatipy.numeric.coordinates.cartesian import Cartesian
        >>> Cartesian(np.zeros(4), vels=np.zeros(3)).name_metric
        'Cartesian'
        """
        super().__init__(
            xs, vels=vels, from_dxs_dt=from_dxs_dt, system_name="Cartesian"
        )

    def _get_vs_from_dxs_dt(self):
        """
        Set spatial velocities equal to coordinate derivatives (Cartesian chart).

        Returns
        -------
        ndarray of shape (3,)
            :math:`v^i = \\mathrm{d}x^i/\\mathrm{d}\\tau` for :math:`i \\in \\{1,2,3\\}`.

        Notes
        -----
        Requires ``self.dxs_dt`` to be set by the base constructor.

        Examples
        --------
        >>> import numpy as np
        >>> from relatipy.numeric.coordinates.cartesian import Cartesian
        >>> c = Cartesian(np.zeros(4), vels=[1.0, 2.0, 3.0], from_dxs_dt=True)
        >>> np.allclose(c.vs, [1.0, 2.0, 3.0])
        True
        """
        x1_prime_dot = self.dxs_dt[0]
        x2_prime_dot = self.dxs_dt[1]
        x3_prime_dot = self.dxs_dt[2]

        return numpy.array([x1_prime_dot, x2_prime_dot, x3_prime_dot])

    def _get_dxs_dt_from_vs(self):
        """
        Set coordinate derivatives equal to spatial velocities (Cartesian chart).

        Returns
        -------
        ndarray of shape (3,)
            :math:`\\mathrm{d}x^i/\\mathrm{d}\\tau = v^i` for :math:`i \\in \\{1,2,3\\}`.

        Notes
        -----
        Requires ``self.vs`` to be set by the base constructor.

        Examples
        --------
        >>> import numpy as np
        >>> from relatipy.numeric.coordinates.cartesian import Cartesian
        >>> c = Cartesian(np.zeros(4), vels=[0.5, -0.5, 0.0], from_dxs_dt=False)
        >>> np.allclose(c.dxs_dt, c.vs)
        True
        """
        dx1_dt = self.vs[0]
        dx2_dt = self.vs[1]
        dx3_dt = self.vs[2]
        return numpy.array([dx1_dt, dx2_dt, dx3_dt])

    def convert_to_cartesian(self):
        """
        Return this instance unchanged (already Cartesian).

        Returns
        -------
        Cartesian
            ``self``.

        Examples
        --------
        >>> from relatipy.numeric.coordinates.cartesian import Cartesian
        >>> c = Cartesian([0.0, 1.0, 0.0, 0.0], vels=[0.0, 0.0, 0.0])
        >>> c.convert_to_cartesian() is c
        True
        """
        return self
