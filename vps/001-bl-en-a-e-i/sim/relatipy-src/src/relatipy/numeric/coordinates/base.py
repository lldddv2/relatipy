"""
Abstract base class for spacetime coordinate representations.

This module defines :class:`CoordinateBase`, which stores a four-position
:math:`x^\\mu = (t, x^1, x^2, x^3)` together with spatial coordinate velocities
:math:`v^i` and coordinate derivatives :math:`\\mathrm{d}x^i/\\mathrm{d}\\tau`
(or the analogous evolution parameter used by the subclass). Subclasses
implement conversions to and from Cartesian coordinates and the map between
:math:`v^i` and :math:`\\mathrm{d}x^i/\\mathrm{d}\\tau`.

The affine parameter contraction

.. math::

    g_{\\mu\\nu} u^\\mu u^\\nu

for the four-velocity ansatz :math:`u^0 = 1`, :math:`u^i = \\mathrm{d}x^i/\\mathrm{d}\\tau`
is exposed for consistency checks along timelike geodesics.

Notes
-----
Quantities are validated through :obj:`~relatipy.numeric.utils.dimensions.validator`
when ``astropy.units`` quantities are passed.

See Also
--------
relatipy.numeric.coordinates :
    Concrete coordinate systems and the ``coordinate_systems`` registry.

Examples
--------
>>> import numpy as np
>>> from relatipy.numeric.coordinates import Cartesian
>>> xs = np.array([0.0, 1.0, 0.0, 0.0])
>>> c = Cartesian(xs, vels=np.zeros(3), from_dxs_dt=False)
>>> c.state_vector.shape
(7,)
"""

from numpy import array, zeros_like, concatenate, ones, einsum

from ..utils.dimensions import validator


class CoordinateBase:
    """
    Base class for coordinate charts with position and velocity state.

    Stores :math:`x^\\mu`, spatial velocities :math:`v^i`, derivatives
    :math:`\\mathrm{d}x^i/\\mathrm{d}\\tau`, and the concatenated state vector
    used by integrators. Subclasses must implement static conversion helpers
    and the relations between :math:`v^i` and :math:`\\mathrm{d}x^i/\\mathrm{d}\\tau`.

    Parameters
    ----------
    xs : array_like of shape (4,)
        Event coordinates :math:`(x^0, x^1, x^2, x^3)` (typically
        :math:`x^0 = t` in geometric units).
    vels : array_like of shape (3,) or None, optional
        Either spatial velocities :math:`v^i` (when ``from_dxs_dt`` is False)
        or :math:`\\mathrm{d}x^i/\\mathrm{d}\\tau` (when True). If None, both
        are set to arrays of zeros matching the spatial part of ``xs``.
    from_dxs_dt : bool, optional
        If True, interpret ``vels`` as :math:`\\mathrm{d}x^i/\\mathrm{d}\\tau`
        and compute :math:`v^i` via :meth:`_get_vs_from_dxs_dt`. If False,
        interpret ``vels`` as :math:`v^i` and compute derivatives via
        :meth:`_get_dxs_dt_from_vs`.
    system_name : str, optional
        Label stored in :attr:`name_metric` for the coordinate system.
    **kwargs
        Additional attributes set on the instance (e.g. mass parameters for
        Boyer-Lindquist coordinates).

    Attributes
    ----------
    name_metric : str
        Name of the coordinate system.
    xs : ndarray
        Four-position, shape ``(4,)`` or broadcastable batch shape as used by
        metrics (see subclass documentation).
    vs : ndarray
        Spatial velocities :math:`v^i`, shape ``(3,)`` or batched.
    dxs_dt : ndarray
        Coordinate derivatives :math:`\\mathrm{d}x^i/\\mathrm{d}\\tau`, same
        shape as ``vs``.
    state_vector : ndarray
        Concatenation ``(x^\\mu, v^i)``, length 7 for a single particle.
    kwargs : dict
        Copy of extra keyword arguments passed to the constructor.

    Raises
    ------
    ValueError
        If ``xs`` does not have length 4.

    Examples
    --------
    >>> import numpy as np
    >>> from relatipy.numeric.coordinates import Cartesian
    >>> xs = np.array([0.0, 2.0, 0.0, 0.0])
    >>> base = Cartesian(xs, vels=np.zeros(3))
    >>> float(base[0])
    0.0
    """

    def __init__(
        self, xs, vels=None, from_dxs_dt=False, system_name="CoordinateBase", **kwargs
    ):
        self.name_metric = system_name

        xs = validator.validate_vector(xs)
        if vels is not None:
            vels = validator.validate_vector(vels)

        self.kwargs = kwargs
        for key, value in kwargs.items():
            setattr(self, key, value)

        if len(xs) != 4:
            raise ValueError(
                "xs must be a list or array of length 4 representing (t, x1, x2, x3), not {}".format(
                    len(xs)
                )
            )

        self.xs = array(xs)

        if vels is None:
            self.vs = zeros_like(xs[1:])
            self.dxs_dt = zeros_like(xs[1:])
        else:
            if not from_dxs_dt:
                self.vs = vels
                self.dxs_dt = self._get_dxs_dt_from_vs()
            else:
                self.dxs_dt = vels
                self.vs = self._get_vs_from_dxs_dt()

        self.state_vector = concatenate((self.xs, self.vs))

    def convert_to_cartesian(self):
        """
        Transform this state to Cartesian coordinates.

        Returns
        -------
        Cartesian
            Instance of :class:`~relatipy.numeric.coordinates.cartesian.Cartesian`
            with positions and velocities in Cartesian form.

        Raises
        ------
        NotImplementedError
            If the subclass does not implement :meth:`_convert_to_cartesian`.

        Examples
        --------
        >>> from relatipy.numeric.coordinates import Cartesian
        >>> c = Cartesian([0.0, 1.0, 0.0, 0.0], vels=[0.0, 0.0, 0.0])
        >>> c_cart = c.convert_to_cartesian()
        >>> float(c_cart.xs[1])
        1.0
        """
        from .cartesian import Cartesian

        xs_p, vs_p = self._convert_to_cartesian(self.xs, self.vs, **self.kwargs)
        coordinate = Cartesian(xs_p, vels=vs_p, from_dxs_dt=False)
        return coordinate

    def convert_to(self, target_system, **kwargs):
        """
        Convert to another registered coordinate system.

        The transformation is performed by mapping to Cartesian coordinates
        first, then constructing the target system via its
        ``_convert_from_cartesian`` hook.

        Parameters
        ----------
        target_system : str
            Key in :obj:`relatipy.numeric.coordinates.coordinate_systems`
            (e.g. ``"Cartesian"``, ``"Spherical"``, ``"Cylindrical"``).
        **kwargs
            Forwarded to the target class's ``_convert_from_cartesian``.

        Returns
        -------
        CoordinateBase
            Instance of the requested coordinate class.

        Raises
        ------
        ValueError
            If ``target_system`` is not a registered name.

        Examples
        --------
        >>> from relatipy.numeric.coordinates import Cartesian
        >>> c = Cartesian([0.0, 3.0, 4.0, 0.0], vels=[0.0, 0.0, 0.0])
        >>> sph = c.convert_to("Spherical")
        >>> float(sph.xs[1])  # radial coordinate r
        5.0
        """
        from . import coordinate_systems

        if target_system not in coordinate_systems:
            raise ValueError(
                f"Unsupported target coordinate system: {target_system}. Supported systems are: {list(coordinate_systems.keys())}"
            )

        cartesian = self.convert_to_cartesian()
        if target_system == "Cartesian":
            return cartesian

        target_class = coordinate_systems[target_system]
        return target_class._convert_from_cartesian(cartesian.xs, cartesian.vs, **kwargs)

    def _get_vs_from_dxs_dt(self):
        """
        Compute spatial velocities :math:`v^i` from :math:`\\mathrm{d}x^i/\\mathrm{d}\\tau`.

        Returns
        -------
        ndarray
            Components :math:`v^i` in the chart of the subclass.

        Raises
        ------
        NotImplementedError
            Always on the base class; subclasses must override.

        Examples
        --------
        >>> CoordinateBase._get_vs_from_dxs_dt(None)  # doctest: +SKIP
        """
        raise NotImplementedError(
            "Subclasses must implement _get_vs_from_dxs_dt method."
        )

    def _get_dxs_dt_from_vs(self):
        """
        Compute :math:`\\mathrm{d}x^i/\\mathrm{d}\\tau` from spatial velocities :math:`v^i`.

        Returns
        -------
        ndarray
            Derivatives :math:`\\mathrm{d}x^i/\\mathrm{d}\\tau`.

        Raises
        ------
        NotImplementedError
            Always on the base class; subclasses must override.

        Examples
        --------
        >>> CoordinateBase._get_dxs_dt_from_vs(None)  # doctest: +SKIP
        """
        raise NotImplementedError(
            "Subclasses must implement _get_dxs_dt_from_vs method."
        )

    @staticmethod
    def _convert_to_cartesian(xs, vs, **kwargs):
        """
        Map position and velocity from this chart to Cartesian coordinates.

        Parameters
        ----------
        xs : array_like
            Four-position in the subclass chart.
        vs : array_like
            Spatial velocities in the subclass chart.
        **kwargs
            Extra parameters required by the chart (subclass-specific).

        Returns
        -------
        xs_p : ndarray
            Cartesian four-position.
        vs_p : ndarray
            Cartesian spatial velocities.

        Raises
        ------
        NotImplementedError
            Always on the base class; subclasses must override.

        Examples
        --------
        >>> CoordinateBase._convert_to_cartesian(None, None, None, None)  # doctest: +SKIP
        """
        raise NotImplementedError(
            "Subclasses must implement _convert_to_cartesian static method."
        )

    @staticmethod
    def _convert_from_cartesian(xs_p, vs_p, **kwargs):
        """
        Construct chart coordinates from Cartesian position and velocity.

        Parameters
        ----------
        xs_p : array_like
            Cartesian four-position.
        vs_p : array_like
            Cartesian spatial velocities.
        **kwargs
            Extra parameters required by the chart (subclass-specific).

        Returns
        -------
        xs : ndarray
            Four-position in the subclass chart.
        vs : ndarray
            Spatial velocities in the subclass chart.

        Raises
        ------
        NotImplementedError
            Always on the base class; subclasses must override.

        Examples
        --------
        >>> CoordinateBase._convert_from_cartesian(None, None, None, None)  # doctest: +SKIP
        """
        raise NotImplementedError(
            "Subclasses must implement _convert_from_cartesian static method."
        )

    def __getitem__(self, index):
        """
        Access components of the state vector :math:`(x^\\mu, v^i)`.

        Indices ``0``–``3`` refer to :math:`x^\\mu`, and ``4``–``6`` to
        spatial velocities. Tuple indexing delegates to the underlying array.

        Parameters
        ----------
        index : int or tuple of int
            Component index or NumPy-style advanced index into
            :attr:`state_vector`.

        Returns
        -------
        float or ndarray
            Selected component(s).

        Raises
        ------
        IndexError
            If ``index`` is an integer outside ``0`` … ``6``.

        Examples
        --------
        >>> import numpy as np
        >>> from relatipy.numeric.coordinates import Cartesian
        >>> c = Cartesian(np.array([1., 2., 3., 4.]), vels=[0., 0., 0.])
        >>> float(c[0]), float(c[6])
        (1.0, 0.0)
        """
        if isinstance(index, tuple):
            return self.state_vector[index]
        if index > 6 or index < 0:
            raise IndexError(
                "Index out of range. Must be between 0 and 6. [x0,x1,x2,x3,v1,v2,v3]"
            )
        return self.state_vector[index]

    def _get_ds_dtau(self, metric):
        """
        Evaluate :math:`g_{\\mu\\nu} u^\\mu u^\\nu` for the stored four-velocity ansatz.

        The implementation uses :math:`u^0 = 1` and :math:`u^i = \\mathrm{d}x^i/\\mathrm{d}\\tau`
        (batched over trailing dimensions). For a timelike geodesic in the
        project's conventions, this scalar should agree with :math:`-c^2` (or
        the equivalent in geometric units).

        Parameters
        ----------
        metric : object
            Metric provider with a ``metric(xs)`` method returning
            :math:`g_{\\mu\\nu}` with shape compatible with ``xs`` (e.g.
            ``(..., 4, 4)`` for batched positions).

        Returns
        -------
        ndarray
            Contraction :math:`g_{\\mu\\nu} u^\\mu u^\\nu`, one value per batch
            element when applicable.

        Examples
        --------
        ``xs`` and ``dxs_dt`` must be batched with a trailing sample axis
        (e.g. shapes ``(4, N)`` and ``(3, N)``) so that ``metric.metric(xs)``
        has shape ``(4, 4, N)``.

        >>> import numpy as np
        >>> from relatipy.numeric.coordinates import Cartesian
        >>> from relatipy.numeric.metrics.Minkowski_metric import Minkowski
        >>> m = Minkowski()
        >>> x = Cartesian(
        ...     np.array([[0.0], [0.0], [0.0], [0.0]]),
        ...     vels=np.zeros((3, 1)),
        ...     from_dxs_dt=False,
        ... )
        >>> val = x._get_ds_dtau(m)
        >>> np.asarray(val).shape
        (1,)
        """
        g = metric.metric(self.xs)

        u = ones((4, len(self.dxs_dt[0])))
        u[1:, :] = self.dxs_dt

        return einsum("ijn,in,jn->n", g, u, u)

    def _get_Lz(self, mass_particle=1.0):
        """
        Azimuthal angular momentum :math:`L_z` in cylindrical coordinates.

        Uses :math:`L_z = m \\, \\rho \\, v_\\phi` with :math:`\\rho` the
        cylindrical radius and :math:`v_\\phi` the azimuthal velocity component
        returned by the cylindrical chart.

        Parameters
        ----------
        mass_particle : float, optional
            Particle mass factor :math:`m` (default 1).

        Returns
        -------
        float or ndarray
            :math:`L_z` in the same shape as broadcast by ``convert_to``.

        Examples
        --------
        >>> import numpy as np
        >>> from relatipy.numeric.coordinates import Cartesian
        >>> c = Cartesian([0.0, 3.0, 4.0, 0.0], vels=[0.0, 0.0, 1.0])
        >>> lz = c._get_Lz(mass_particle=1.0)
        >>> bool(abs(lz) < 1e9)
        True
        """
        cylindrical = self.convert_to("Cylindrical")
        rho = cylindrical.xs[1]
        v_phi = cylindrical.vs[1]
        return mass_particle * rho * v_phi

    def _get_E(self, metric):
        """
        Conserved energy-like quantity :math:`E = -g_{0\\mu} u^\\mu`.

        With :math:`u^0 = 1` and :math:`u^i = \\mathrm{d}x^i/\\mathrm{d}\\tau`,
        this reduces to a contraction of the metric's time row with the
        four-velocity.

        Parameters
        ----------
        metric : object
            Metric provider with ``metric(xs)`` returning :math:`g_{\\mu\\nu}`,
            shape ``(4, 4, N)`` or compatible for einsum with ``dxs_dt``.

        Returns
        -------
        ndarray
            Values of :math:`-g_{0\\mu} u^\\mu` per batch index.

        Examples
        --------
        Same batching requirements as :meth:`_get_ds_dtau`.

        >>> import numpy as np
        >>> from relatipy.numeric.coordinates import Cartesian
        >>> from relatipy.numeric.metrics.Minkowski_metric import Minkowski
        >>> m = Minkowski()
        >>> x = Cartesian(
        ...     np.array([[0.0], [0.0], [0.0], [0.0]]),
        ...     vels=np.zeros((3, 1)),
        ...     from_dxs_dt=False,
        ... )
        >>> e = x._get_E(m)
        >>> np.asarray(e).shape
        (1,)
        """
        from numpy import einsum, ones

        g = metric.metric(self.xs)

        u = ones((4, len(self.dxs_dt[0])))
        u[1:, :] = self.dxs_dt

        return -einsum("jn,jn->n", g[0, :, :], u)
