"""
Two-dimensional plotting helpers for Relatipy.

This subpackage exposes Matplotlib-based utilities styled for publication-quality
figures (serif fonts, inward ticks, optional minor ticks, thin spines), aligned
with common astronomy/orbits figure conventions.

Notes
-----
The public surface of this module is intentionally small: it re-exports
:class:`~relatipy.visualization._2D.sci.SciSubplot`, which configures
``matplotlib`` defaults and provides a figure/axes pair (or subplot layout)
for 2D plots.

Examples
--------
>>> from relatipy.visualization._2D import SciSubplot
>>> workspace = SciSubplot(figsize=(4, 3))
>>> _ = workspace.ax.plot([0, 1], [0, 1])
"""

from .sci import SciSubplot

__all__ = [
    "SciSubplot",
]
