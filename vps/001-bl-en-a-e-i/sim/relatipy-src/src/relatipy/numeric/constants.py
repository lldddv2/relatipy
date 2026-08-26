r"""
Physical and geometric-unit constants for numeric relativity.

Quantities :math:`G` and :math:`c` are set to unity in geometric units for
internal calculations. This module also provides SI values for :math:`G` and
:math:`c`, a solar-mass reference in kilograms, and conversion factors from
geometric length and time to SI when the mass scale is one solar mass.

The reference length and time scales are

.. math::

    L_{\mathrm{ref}} = \frac{G_{\mathrm{SI}} M_{\odot}}{c_{\mathrm{SI}}^2}, \qquad
    T_{\mathrm{ref}} = \frac{G_{\mathrm{SI}} M_{\odot}}{c_{\mathrm{SI}}^3}.

Constants
---------
The following module-level attributes are floats.

_G : float
    Gravitational constant in geometric units (:math:`G=1`).
_c : float
    Speed of light in geometric units (:math:`c=1`).
_G_SI : float
    Newton's gravitational constant in SI
    (:math:`\mathrm{m}^3\,\mathrm{kg}^{-1}\,\mathrm{s}^{-2}`).
_c_SI : float
    Speed of light in SI (m/s).
_M_ref_kg : float
    Solar mass :math:`M_{\odot}` used as the reference mass (kg).
_L_ref : float
    Meters per geometric length unit at mass scale :math:`M_{\odot}`.
_T_ref : float
    Seconds per geometric time unit at mass scale :math:`M_{\odot}`.

Examples
--------
>>> from relatipy.numeric.constants import _G, _c, _G_SI
>>> (_G, _c)
(1.0, 1.0)
>>> _G_SI > 0
True
"""

_G = 1.0
_c = 1.0

_G_SI = 6.67430e-11  # m^3 kg^-1 s^-2
_c_SI = 299792458.0  # m/s

_M_ref_kg = 1.98892e30  # kg, solar mass (reference mass unit)
_L_ref = _G_SI * _M_ref_kg / _c_SI**2  # m per geometric length unit
_T_ref = _G_SI * _M_ref_kg / _c_SI**3  # s per geometric time unit
