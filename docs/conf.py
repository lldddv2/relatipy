"""Sphinx configuration for the RelatiPy documentation."""

from importlib.metadata import PackageNotFoundError, version


project = "RelatiPy"
author = "lldddv2"
copyright = "2026, lldddv2"

try:
    release = version("relatipy")
except PackageNotFoundError:
    release = "0.0.0"

version = release

needs_sphinx = "8.2"
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_reredirects",
    "myst_nb",
]

# Tutorial notebooks are committed with their outputs; the docs build does not
# execute them. Re-run a notebook locally after changing it.
nb_execution_mode = "off"
myst_enable_extensions = ["dollarmath", "amsmath"]

# Pages moved into section folders; keep the old URLs working.
redirects = {
    "installation": "user-guide/installation.html",
    "usage": "user-guide/usage.html",
    "kerr-mcmc": "user-guide/kerr-mcmc.html",
    "integration-roadmap": "user-guide/scope-and-limits.html",
    "contributing": "development/contributing.html",
    "developer-architecture": "development/architecture.html",
    "peer-validation": "development/peer-validation.html",
    "api": "reference/index.html",
    "metrics-reference": "reference/metrics.html",
    "geodesic-reference": "reference/geodesic.html",
    "coordinates-reference": "reference/coordinates.html",
    "observables-reference": "reference/observables.html",
    "plotting-reference": "reference/plotting.html",
}

root_doc = "index"
language = "en"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**/.ipynb_checkpoints"]

# Every Python cross-reference must resolve; the strict build uses ``-n``.
nitpicky = True
nitpick_ignore = [
    # Plotly publishes no Sphinx object inventory, so its figure class cannot
    # be resolved through intersphinx.
    ("py:class", "plotly.graph_objects.Figure"),
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "astropy": ("https://docs.astropy.org/en/stable", None),
    "matplotlib": ("https://matplotlib.org/stable", None),
}

autosummary_generate = True
autodoc_member_order = "bysource"
# Types are documented in the NumPy-style docstrings; rendering annotations as
# well would duplicate them with unresolved ``np.``/``u.`` aliases.
autodoc_typehints = "none"
autodoc_default_options = {
    "members": True,
    "show-inheritance": True,
}

napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_param = True
napoleon_use_rtype = True
# Convert NumPy-style type strings ("float, optional", ``{"xy", "xz"}``) into
# cross-references and literals so nitpicky mode checks every real type name.
napoleon_preprocess_types = True
napoleon_type_aliases = {
    "Quantity": "~astropy.units.Quantity",
    "u.Quantity": "~astropy.units.Quantity",
    "u.UnitBase": "~astropy.units.UnitBase",
    "np.ndarray": "~numpy.ndarray",
    "ndarray": "~numpy.ndarray",
    "array-like": ":term:`array-like <numpy:array_like>`",
    "array_like": ":term:`array_like <numpy:array_like>`",
    "mapping": "~collections.abc.Mapping",
    "sequence": "~collections.abc.Sequence",
    "Figure": "~matplotlib.figure.Figure",
    "Axes": "~matplotlib.axes.Axes",
    # Short names used in type fields of public docstrings.
    "Orbit": "~relatipy.geodesic.Orbit",
    "State": "~relatipy.geodesic.State",
    "Solution": "~relatipy.geodesic.Solution",
    "OrbitalElements": "~relatipy.coordinates.OrbitalElements",
    "KerrOrbitalElements": "~relatipy.coordinates.KerrOrbitalElements",
    "CartesianCoordinates": "~relatipy.coordinates.CartesianCoordinates",
    "CartesianStateVector": "~relatipy.coordinates.CartesianStateVector",
    "SphericalCoordinates": "~relatipy.coordinates.SphericalCoordinates",
    "SphericalVelocity": "~relatipy.coordinates.SphericalVelocity",
    "SphericalFourVelocity": "~relatipy.coordinates.SphericalFourVelocity",
    "BoyerLindquistCoordinates": "~relatipy.coordinates.BoyerLindquistCoordinates",
    "BoyerLindquistVelocity": "~relatipy.coordinates.BoyerLindquistVelocity",
    "BoyerLindquistFourVelocity": "~relatipy.coordinates.BoyerLindquistFourVelocity",
}

html_theme = "sphinx_rtd_theme"
html_title = f"RelatiPy {release} documentation"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_logo = "_static/logo-dark.svg"
html_favicon = "_static/logo-mark.svg"
html_theme_options = {
    "logo_only": True,
    "collapse_navigation": False,
    "navigation_depth": 4,
}
