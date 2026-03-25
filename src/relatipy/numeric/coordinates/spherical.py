"""
Spherical coordinate chart for the numeric coordinate stack.

This module defines :class:`Spherical`, a subclass of
:class:`~relatipy.numeric.coordinates.base.CoordinateBase` that stores an event
as :math:`x^\\mu = (t, r, \\theta, \\phi)` with polar angle :math:`\\theta`
and azimuth :math:`\\phi`, consistent with the spatial map

.. math::

    (x, y, z) = (r \\sin\\theta \\cos\\phi,\\;
                 r \\sin\\theta \\sin\\phi,\\;
                 r \\cos\\theta).

Spatial velocities :math:`v^i` are taken as the physical orthonormal components
:math:`(v_r,\\, r\\dot{\\theta},\\, r\\sin\\theta\\,\\dot{\\phi})` in the same
units as the coordinate derivatives :math:`\\dot{r},\\, \\dot{\\theta},\\,
\\dot{\\phi}` (up to the overall affine-parameter normalization used by the
integrator). Conversions to and from Cartesian :math:`(t, x, y, z)` and the
map between :math:`v^i` and :math:`\\mathrm{d}x^i/\\mathrm{d}\\tau` are
implemented here.

Notes
-----
At :math:`r = 0` or :math:`\\sin\\theta = 0`, the relations between
:math:`v^i` and :math:`\\mathrm{d}x^i/\\mathrm{d}\\tau` involve division by
:math:`r` or :math:`r\\sin\\theta`; callers should avoid these coordinate
singularities or handle non-finite values explicitly.

See Also
--------
relatipy.numeric.coordinates.base.CoordinateBase :
    Base class API and state vector layout.
relatipy.numeric.coordinates :
    Package overview and spherical convention summary.

Examples
--------
>>> import numpy as np
>>> from relatipy.numeric.coordinates.spherical import Spherical
>>> xs = np.array([0.0, 1.0, np.pi / 2, 0.0])
>>> s = Spherical(xs, vels=np.zeros(3), from_dxs_dt=False)
>>> s.xs.shape
(4,)
"""

import numpy
from numpy import sin, cos, sqrt, arctan2, arccos

from .base import CoordinateBase


class Spherical(CoordinateBase):
    """
    Spherical coordinates :math:`(t, r, \\theta, \\phi)` with Cartesian map.

    Subclass of :class:`~relatipy.numeric.coordinates.base.CoordinateBase`
    using the standard polar angle :math:`\\theta` and azimuth :math:`\\phi`
    as in the module docstring.

    Parameters
    ----------
    xs : array_like of shape (4,)
        Event :math:`(x^0, r, \\theta, \\phi)`; typically :math:`x^0 = t` in
        geometric units.
    vels : array_like of shape (3,) or None, optional
        Either :math:`(v_r,\\, r\\dot{\\theta},\\, r\\sin\\theta\\,\\dot{\\phi})`
        when ``from_dxs_dt`` is False, or
        :math:`(\\mathrm{d}r/\\mathrm{d}\\tau,\\,
        \\mathrm{d}\\theta/\\mathrm{d}\\tau,\\,
        \\mathrm{d}\\phi/\\mathrm{d}\\tau)` when True. If None, velocities
        default to zeros (see base class).
    from_dxs_dt : bool, optional
        If True, ``vels`` are coordinate derivatives and :math:`v^i` are filled
        via :meth:`_get_vs_from_dxs_dt`. If False, the inverse map uses
        :meth:`_get_dxs_dt_from_vs`.

    Attributes
    ----------
    All attributes are those of
    :class:`~relatipy.numeric.coordinates.base.CoordinateBase` (including
    ``xs``, ``vs``, ``dxs_dt``, ``state_vector``, and ``name_metric``).

    Examples
    --------
    >>> import numpy as np
    >>> from relatipy.numeric.coordinates.spherical import Spherical
    >>> xs = np.array([0.0, 2.0, np.pi / 2, np.pi / 2])
    >>> Spherical(xs, vels=np.zeros(3), from_dxs_dt=False).name_metric
    'Spherical'
    """

    def __init__(self, xs, vels=None, from_dxs_dt=False):
        """
        Construct a spherical coordinate state.

        Parameters
        ----------
        xs : array_like of shape (4,)
            :math:`(t, r, \\theta, \\phi)`.
        vels : array_like of shape (3,) or None, optional
            Velocities or coordinate derivatives, per ``from_dxs_dt``.
        from_dxs_dt : bool, optional
            Interpretation flag for ``vels`` (default False).

        Examples
        --------
        >>> import numpy as np
        >>> from relatipy.numeric.coordinates.spherical import Spherical
        >>> Spherical(np.zeros(4), vels=np.zeros(3)).xs.shape
        (4,)
        """
        super().__init__(
            xs, vels=vels, from_dxs_dt=from_dxs_dt, system_name="Spherical"
        )

    def _get_vs_from_dxs_dt(self):
        """
        Map coordinate derivatives to spherical velocity components.

        With :math:`\\dot{x}^i = \\mathrm{d}x^i/\\mathrm{d}\\tau`,

        .. math::

            v_r = \\dot{r}, \\quad
            v_\\theta = r \\dot{\\theta}, \\quad
            v_\\phi = r \\sin\\theta \\, \\dot{\\phi}.

        Returns
        -------
        ndarray of shape (3,)
            Stack :math:`(v_r,\\, r\\dot{\\theta},\\, r\\sin\\theta\\,\\dot{\\phi})`.

        Examples
        --------
        >>> import numpy as np
        >>> from relatipy.numeric.coordinates.spherical import Spherical
        >>> xs = np.array([0.0, 1.0, np.pi / 2, 0.0])
        >>> s = Spherical.__new__(Spherical)
        >>> s.xs, s.dxs_dt = xs, np.array([0.0, 1.0, 0.0])
        >>> Spherical._get_vs_from_dxs_dt(s)
        array([0., 1., 0.])
        """
        x1_prime_dot = self.dxs_dt[0]
        x2_prime_dot = self.dxs_dt[1] * self.xs[1]
        x3_prime_dot = self.dxs_dt[2] * self.xs[1] * sin(self.xs[2])

        return numpy.array([x1_prime_dot, x2_prime_dot, x3_prime_dot])

    def _get_dxs_dt_from_vs(self):
        """
        Map spherical velocity components to coordinate derivatives.

        Inverts :meth:`_get_vs_from_dxs_dt`:

        .. math::

            \\dot{r} = v_r, \\quad
            \\dot{\\theta} = \\frac{v_\\theta}{r}, \\quad
            \\dot{\\phi} = \\frac{v_\\phi}{r \\sin\\theta}.

        Returns
        -------
        ndarray of shape (3,)
            :math:`(\\dot{r},\\, \\dot{\\theta},\\, \\dot{\\phi})`.

        Examples
        --------
        >>> import numpy as np
        >>> from relatipy.numeric.coordinates.spherical import Spherical
        >>> xs = np.array([0.0, 1.0, np.pi / 2, 0.0])
        >>> s = Spherical.__new__(Spherical)
        >>> s.xs, s.vs = xs, np.array([0.0, 1.0, 0.0])
        >>> Spherical._get_dxs_dt_from_vs(s)
        array([0., 1., 0.])
        """
        dx1_dt = self.vs[0]
        dx2_dt = self.vs[1] / self.xs[1]
        dx3_dt = self.vs[2] / (self.xs[1] * sin(self.xs[2]))
        return numpy.array([dx1_dt, dx2_dt, dx3_dt])

    @staticmethod
    def _convert_to_cartesian(xs, vs):
        """
        Convert spherical position and velocities to Cartesian.

        The time component is unchanged: :math:`t' = t`. Spatial Cartesian
        coordinates follow the standard spherical map; velocity components
        :math:`v_x, v_y, v_z` are obtained from
        :math:`(v_r,\\, r\\dot{\\theta},\\, r\\sin\\theta\\,\\dot{\\phi})`
        by the usual orthonormal change of basis.

        Parameters
        ----------
        xs : ndarray of shape (4,) or broadcast-compatible
            :math:`(t, r, \\theta, \\phi)`.
        vs : ndarray of shape (3,) or broadcast-compatible
            :math:`(v_r,\\, r\\dot{\\theta},\\, r\\sin\\theta\\,\\dot{\\phi})`.

        Returns
        -------
        xs_p : ndarray
            Cartesian four-position :math:`(t, x, y, z)`, same shape as ``xs``.
        vs_p : ndarray
            Cartesian spatial velocities :math:`(v_x, v_y, v_z)`, same shape as
            ``vs``.

        Examples
        --------
        >>> import numpy as np
        >>> from relatipy.numeric.coordinates.spherical import Spherical
        >>> xs = np.array([0.0, 1.0, np.pi / 2, 0.0])
        >>> xp, vp = Spherical._convert_to_cartesian(xs, np.zeros(3))
        >>> np.allclose(xp, [0.0, 1.0, 0.0, 0.0])
        True
        """
        xs_p = numpy.zeros_like(xs)
        vs_p = numpy.zeros_like(vs)

        xs_p[0] = xs[0]
        xs_p[1] = xs[1] * sin(xs[2]) * cos(xs[3])
        xs_p[2] = xs[1] * sin(xs[2]) * sin(xs[3])
        xs_p[3] = xs[1] * cos(xs[2])

        sin_x2 = sin(xs[2])
        cos_x2 = cos(xs[2])
        sin_x3 = sin(xs[3])
        cos_x3 = cos(xs[3])

        vs_p[0] = sin_x2 * cos_x3 * vs[0] + cos_x2 * cos_x3 * vs[1] - sin_x3 * vs[2]

        vs_p[1] = sin_x2 * sin_x3 * vs[0] + cos_x2 * sin_x3 * vs[1] + cos_x3 * vs[2]

        vs_p[2] = cos_x2 * vs[0] - sin_x2 * vs[1]

        return xs_p, vs_p

    @staticmethod
    def _convert_from_cartesian(xs_p, vs_p):
        """
        Build a :class:`Spherical` instance from Cartesian state.

        Recovers :math:`r = \\sqrt{x^2 + y^2 + z^2}`,
        :math:`\\theta = \\arccos(z/r)`, :math:`\\phi = \\mathrm{atan2}(y, x)`
        (with the usual branch), and the inverse velocity transformation.

        Parameters
        ----------
        xs_p : ndarray of shape (4,) or broadcast-compatible
            Cartesian :math:`(t, x, y, z)`.
        vs_p : ndarray of shape (3,) or broadcast-compatible
            Cartesian spatial velocities :math:`(v_x, v_y, v_z)`.

        Returns
        -------
        Spherical
            Spherical coordinate object with ``from_dxs_dt=False``.

        Examples
        --------
        >>> import numpy as np
        >>> from relatipy.numeric.coordinates.spherical import Spherical
        >>> xs_p = np.array([0.0, 1.0, 0.0, 0.0])
        >>> sph = Spherical._convert_from_cartesian(xs_p, np.zeros(3))
        >>> np.isclose(sph.xs[1], 1.0) and np.isclose(sph.xs[3], 0.0)
        True
        """
        xs = numpy.zeros_like(xs_p)
        vs = numpy.zeros_like(vs_p)

        xs[0] = xs_p[0]
        xs[1] = sqrt(xs_p[1] ** 2 + xs_p[2] ** 2 + xs_p[3] ** 2)
        xs[2] = arccos(xs_p[3] / xs[1])
        xs[3] = arctan2(xs_p[2], xs_p[1])

        sin_x2 = sin(xs[2])
        cos_x2 = cos(xs[2])
        sin_x3 = sin(xs[3])
        cos_x3 = cos(xs[3])

        vs[0] = sin_x2 * cos_x3 * vs_p[0] + sin_x2 * sin_x3 * vs_p[1] + cos_x2 * vs_p[2]

        vs[1] = cos_x2 * cos_x3 * vs_p[0] + cos_x2 * sin_x3 * vs_p[1] - sin_x2 * vs_p[2]

        vs[2] = -sin_x3 * vs_p[0] + cos_x3 * vs_p[1]

        return Spherical(xs, vels=vs, from_dxs_dt=False)
