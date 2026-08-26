"""
Cylindrical coordinates for the relatipy numeric coordinate stack.

This module implements :class:`Cylindrical`, a chart with four-position

.. math::

    x^\\mu = (t, \\rho, \\phi, z),

where :math:`\\rho \\ge 0` is the radial distance from the :math:`z` axis,
:math:`\\phi` is the azimuth about that axis, and :math:`z` is the usual
Cartesian height. Spatial velocities are stored as
:math:`(v_\\rho, v_\\phi, v_z)` with :math:`v_\\phi = \\rho\\, \\mathrm{d}\\phi/\\mathrm{d}\\tau`
(momentum-like azimuthal component), while :attr:`~relatipy.numeric.coordinates.base.CoordinateBase.dxs_dt`
holds :math:`(\\mathrm{d}\\rho/\\mathrm{d}\\tau, \\mathrm{d}\\phi/\\mathrm{d}\\tau, \\mathrm{d}z/\\mathrm{d}\\tau)`.

The map to Cartesian spatial coordinates is

.. math::

    x = \\rho \\cos\\phi, \\quad y = \\rho \\sin\\phi,

with corresponding linear velocity components as implemented in
:meth:`Cylindrical._convert_to_cartesian`.

Notes
-----
The filename ``cilindrical`` is a historical spelling; the public class name
is :class:`Cylindrical`.

See Also
--------
relatipy.numeric.coordinates.base.CoordinateBase :
    Base class defining state vectors and integration interface.

Examples
--------
>>> import numpy as np
>>> from relatipy.numeric.coordinates import Cylindrical
>>> xs = np.array([0.0, 1.0, 0.0, 0.0])
>>> c = Cylindrical(xs, vels=np.zeros(3), from_dxs_dt=False)
>>> c.name_metric
'Cylindrical'
"""

import numpy
from numpy import sin, cos, sqrt, arctan2

from .base import CoordinateBase


class Cylindrical(CoordinateBase):
    """
    Cylindrical chart :math:`(t, \\rho, \\phi, z)` with linked velocities.

    Inherits storage and validation from :class:`CoordinateBase`. Spatial
    derivatives :math:`\\mathrm{d}x^i/\\mathrm{d}\\tau` are
    :math:`(\\dot\\rho, \\dot\\phi, \\dot z)`; spatial velocities are
    :math:`(v_\\rho, v_\\phi, v_z)` with :math:`v_\\phi = \\rho \\dot\\phi`.

    Parameters
    ----------
    xs : array_like of shape (4,)
        Event coordinates :math:`(t, \\rho, \\phi, z)`.
    vels : array_like of shape (3,) or None, optional
        Either :math:`(v_\\rho, v_\\phi, v_z)` if ``from_dxs_dt`` is False,
        or :math:`(\\dot\\rho, \\dot\\phi, \\dot z)` if True. If None, zeros
        are used.
    from_dxs_dt : bool, optional
        If True, interpret ``vels`` as spatial coordinate derivatives and
        compute :math:`v^i` via :meth:`_get_vs_from_dxs_dt`. If False, the
        converse via :meth:`_get_dxs_dt_from_vs`.

    Attributes
    ----------
    name_metric : str
        Set to ``\"Cylindrical\"`` (see base class).

    Examples
    --------
    >>> import numpy as np
    >>> from relatipy.numeric.coordinates import Cylindrical
    >>> xs = np.array([0.0, 2.0, np.pi / 2, 1.0])
    >>> c = Cylindrical(xs, vels=np.array([0.0, 0.0, 0.0]))
    >>> float(c.xs[1])  # rho
    2.0
    """

    def __init__(self, xs, vels=None, from_dxs_dt=False):
        super().__init__(
            xs, vels=vels, from_dxs_dt=from_dxs_dt, system_name="Cylindrical"
        )

    def _get_vs_from_dxs_dt(self):
        """
        Map spatial derivatives to cylindrical velocity components.

        Uses :math:`v_\\rho = \\dot\\rho`, :math:`v_\\phi = \\rho\\,\\dot\\phi`,
        :math:`v_z = \\dot z`.

        Returns
        -------
        ndarray of shape (3,)
            Array ``[v_rho, v_phi, v_z]``.

        Examples
        --------
        >>> import numpy as np
        >>> from relatipy.numeric.coordinates import Cylindrical
        >>> xs = np.array([0.0, 1.0, 0.5, 0.0])
    >>> c = Cylindrical(xs, vels=np.array([0.1, 0.2, 0.3]), from_dxs_dt=True)
    >>> float(c.vs[1])  # rho * phi_dot == 1.0 * 0.2
    0.2
        """
        v_rho = self.dxs_dt[0]
        v_phi = self.dxs_dt[1] * self.xs[1]
        v_z = self.dxs_dt[2]
        return numpy.array([v_rho, v_phi, v_z])

    def _get_dxs_dt_from_vs(self):
        """
        Map cylindrical velocities to spatial coordinate derivatives.

        Uses :math:`\\dot\\rho = v_\\rho`,
        :math:`\\dot\\phi = v_\\phi / \\rho`, :math:`\\dot z = v_z`.

        Returns
        -------
        ndarray of shape (3,)
            Array ``[drho_dt, dphi_dt, dz_dt]``.

        Notes
        -----
        For :math:`\\rho = 0`, :math:`\\dot\\phi = v_\\phi/\\rho` is undefined
        in the continuum model; the implementation performs floating-point
        division as given.

        Examples
        --------
        >>> import numpy as np
        >>> from relatipy.numeric.coordinates import Cylindrical
        >>> xs = np.array([0.0, 2.0, 0.0, 0.0])
    >>> c = Cylindrical(xs, vels=np.array([0.0, 1.0, 0.0]), from_dxs_dt=False)
    >>> float(c.dxs_dt[1])  # v_phi / rho == 0.5
    0.5
        """
        drho_dt = self.vs[0]
        dphi_dt = self.vs[1] / self.xs[1]
        dz_dt = self.vs[2]
        return numpy.array([drho_dt, dphi_dt, dz_dt])

    @staticmethod
    def _convert_to_cartesian(xs, vs):
        """
        Convert cylindrical position and velocities to Cartesian spatial form.

        The time component is copied unchanged: :math:`t' = t`. Spatial
        position follows :math:`x = \\rho\\cos\\phi`, :math:`y = \\rho\\sin\\phi`,
        :math:`z' = z`. Velocities use

        .. math::

            v_x = \\cos\\phi\\, v_\\rho - \\sin\\phi\\, v_\\phi, \\quad
            v_y = \\sin\\phi\\, v_\\rho + \\cos\\phi\\, v_\\phi, \\quad
            v_z' = v_z.

        Parameters
        ----------
        xs : ndarray
            Four-position :math:`(t, \\rho, \\phi, z)`, same shape as used by
            the instance ``xs`` array (typically shape ``(4,)``).
        vs : ndarray
            Cylindrical velocities :math:`(v_\\rho, v_\\phi, v_z)` (typically
            shape ``(3,)``).

        Returns
        -------
        xs_p : ndarray
            Cartesian four-position :math:`(t, x, y, z)`.
        vs_p : ndarray
            Cartesian spatial velocities :math:`(v_x, v_y, v_z)`.

        Examples
        --------
        >>> import numpy as np
        >>> from relatipy.numeric.coordinates.cilindrical import Cylindrical
        >>> xs = np.array([0.0, 1.0, 0.0, 0.0])
        >>> vs = np.array([1.0, 0.0, 0.0])
        >>> xp, vp = Cylindrical._convert_to_cartesian(xs, vs)
        >>> float(xp[1]), float(vp[0])
        (1.0, 1.0)
        """
        xs_p = numpy.zeros_like(xs)
        vs_p = numpy.zeros_like(vs)

        rho, phi, z = xs[1], xs[2], xs[3]
        sin_phi = sin(phi)
        cos_phi = cos(phi)

        xs_p[0] = xs[0]
        xs_p[1] = rho * cos_phi
        xs_p[2] = rho * sin_phi
        xs_p[3] = z

        vs_p[0] = cos_phi * vs[0] - sin_phi * vs[1]
        vs_p[1] = sin_phi * vs[0] + cos_phi * vs[1]
        vs_p[2] = vs[2]

        return xs_p, vs_p

    @staticmethod
    def _convert_from_cartesian(xs_p, vs_p):
        """
        Build a :class:`Cylindrical` instance from Cartesian state.

        Sets :math:`\\rho = \\sqrt{x^2 + y^2}`,
        :math:`\\phi = \\operatorname{atan2}(y, x)`, :math:`z` unchanged, and
        inverts the velocity map of :meth:`_convert_to_cartesian`.

        Parameters
        ----------
        xs_p : ndarray
            Cartesian four-position :math:`(t, x, y, z)`.
        vs_p : ndarray
            Cartesian spatial velocities :math:`(v_x, v_y, v_z)`.

        Returns
        -------
        Cylindrical
            New instance with ``from_dxs_dt=False`` and velocities in
            :math:`(v_\\rho, v_\\phi, v_z)` form.

        Examples
        --------
        >>> import numpy as np
        >>> from relatipy.numeric.coordinates.cilindrical import Cylindrical
        >>> xs_p = np.array([0.0, 1.0, 0.0, 0.0])
        >>> vs_p = np.array([0.0, 1.0, 0.0])
        >>> c = Cylindrical._convert_from_cartesian(xs_p, vs_p)
        >>> float(c.xs[2])  # phi
        0.0
        >>> float(c.vs[1])  # v_phi
        1.0
        """
        xs = numpy.zeros_like(xs_p)
        vs = numpy.zeros_like(vs_p)

        xs[0] = xs_p[0]
        xs[1] = sqrt(xs_p[1] ** 2 + xs_p[2] ** 2)
        xs[2] = arctan2(xs_p[2], xs_p[1])
        xs[3] = xs_p[3]

        sin_phi = sin(xs[2])
        cos_phi = cos(xs[2])

        vs[0] = cos_phi * vs_p[0] + sin_phi * vs_p[1]
        vs[1] = -sin_phi * vs_p[0] + cos_phi * vs_p[1]
        vs[2] = vs_p[2]

        coordinate = Cylindrical(xs, vels=vs, from_dxs_dt=False)
        return coordinate
