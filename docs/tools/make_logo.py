"""Generate the RelatiPy logo: the curved space around a Kerr black hole.

The mark is the exact isometric embedding of the equatorial spatial slice
(t = const, theta = pi/2) of Kerr with spin a = 0.998, in units G = c = M = 1.
Its induced metric is

    dl^2 = (r^2 / Delta) dr^2 + R(r)^2 dphi^2,
    Delta = r^2 - 2 r + a^2,   R(r)^2 = r^2 + a^2 + 2 a^2 / r,

embedded in flat 3-space as a surface of revolution with cylindrical radius
R(r) and height z(r) = int sqrt(r^2 / Delta - (dR/dr)^2) dr, integrated from
the outer horizon r_+ (where R = 2 exactly). For a = 0.998 the integrand is
positive everywhere outside r_+, so the whole exterior slice is embedded.
The substitution r = r_+ + s^2 removes the 1/sqrt(r - r_+) singularity at the
throat before integrating.

A square Cartesian grid, x = R cos(phi) and y = R sin(phi), is laid on that
surface and drawn in perspective from straight above the throat (camera on
the spin axis, 12 r_g above mid-height, close enough that depth shows as
strong foreshortening toward the hole). The z-buffer for hidden lines is kept
for oblique views; from above nothing is hidden. The grid has no border:
stroke opacity falls smoothly to zero between the superellipses
(|x|^4 + |y|^4)^(1/4) = 2 and 7.5 r_g. The mark stays nearly square, yet every
line fades along its own length, so no edge line is left.

To make the spin visible, every grid point is carried by the zero-angular-
momentum observer (ZAMO) at its radius for a coordinate time of 6 r_g / c:
it turns by omega(r) * 6, with omega the equatorial frame-dragging rate. The
twist is about 2.8 rad at the horizon (omega_H = a / (2 r_+)) and falls off
roughly as 2 a / r^3, which gives the swirl. The embedding itself is static;
the swirl is a picture of frame dragging, not of the spatial geometry.

The grey ring is the outer horizon (the throat, R = 2); in oblique views its
hidden part is dashed. The grey trail is one turn of the prograde circular
equatorial geodesic at r = 3.5 r_g, integrated by RelatiPy with ``dop853`` and
lifted onto the surface through its Boyer-Lindquist (r, phi). Its opacity
rises linearly from 0.2 at the start of the turn to 1 halfway round and falls
back to 0.2 at the end. Around it are nine synthetic data
points with 1-sigma error bars along the image axes: projected samples of the
geodesic plus Gaussian scatter (sigma between 0.15 and 0.35 r_g at the orbit's
image scale, truncated at 1.2 sigma, fixed seed), applied in the image plane
so they are undistorted. They evoke
fitting an orbit to observations; they are not real measurements.

The background is transparent. Two colour variants are written: logo.svg and
logo-mark.svg (dark strokes, for light pages) and logo-dark.svg and
logo-mark-dark.svg (light strokes, for dark pages). Grid, horizon and
wordmark are greys, the orbit a slightly darker grey and the data points crimson. The
wordmark uses JetBrains Mono Light (SIL Open Font License), converted to
outlines so the SVG does not depend on installed fonts.

Run from the repository root: ``python docs/tools/make_logo.py [font.ttf]``.
"""

import sys
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.constants import c
from fontTools.pens.svgPathPen import SVGPathPen
from fontTools.pens.transformPen import TransformPen
from fontTools.ttLib import TTFont
from scipy.ndimage import minimum_filter

from relatipy import Kerr
from relatipy.coordinates import KerrOrbitalElements

SPIN = 0.998  # Thorne's limit for an accreting hole
ELEVATION = np.radians(90.0)
AZIMUTH = np.radians(0.0)
CAMERA = 12.0  # camera distance from the throat, in r_g
R_FADE = (2.0, 7.5)  # opacity 1 inside the first superellipse radius, 0 beyond the second
SQUARENESS = 4.0  # exponent of the superellipse norm (2 = circle, large = square)
CELL = 0.9  # grid spacing, in r_g
WIDTH = 0.45  # grid stroke width, px
# Colours are placeholders filled per variant: light pages get dark strokes
# (logo*.svg) and dark pages light strokes (logo*-dark.svg).
PALETTES = {
    "": {"@grid@": "#5c5c5c", "@text@": "#3d3d3d", "@data@": "#c8102e", "@orbit@": "#6b6b6b"},
    "-dark": {"@grid@": "#b0b0b0", "@text@": "#d6d6d6", "@data@": "#e8405e", "@orbit@": "#c4c4c4"},
}
GREY = "@grid@"  # grid and horizon
ORBIT = "@orbit@"  # grey for the geodesic
DATA = "@data@"  # synthetic data points and error bars (crimson)
DATA_POINTS = 9
DRAG_TIME = 6.0  # coordinate time the grid is carried by ZAMOs, in r_g / c
ORBIT_RADIUS = 3.5  # prograde circular equatorial orbit, in r_g (ISCO is 1.24)
FONT = sys.argv[1] if len(sys.argv) > 1 else "/usr/share/fonts/TTF/JetBrainsMono-Light.ttf"
OUT = Path(__file__).resolve().parents[1] / "_static"

# Embedding profile R(r), z(r).
a = SPIN
r_plus = 1.0 + np.sqrt(1.0 - a**2)
r_minus = 1.0 - np.sqrt(1.0 - a**2)
bh = Kerr(mass=1 * u.Msun, spin=SPIN)
assert np.isclose((bh.horizons.event / bh.r_g).decompose().value, r_plus)

s = np.linspace(0.0, np.sqrt(60.0 - r_plus), 200001)
r = r_plus + s**2
R = np.sqrt(r**2 + a**2 + 2 * a**2 / r)
dR_dr = (r - a**2 / r**2) / R
# r^2 / Delta * (dr/ds)^2 = 4 r^2 / (r - r_minus), finite at the throat.
dz_ds = np.sqrt(4 * r**2 / (r - r_minus) - (dR_dr * 2 * s) ** 2)
assert np.all(np.isfinite(dz_ds)), "slice not embeddable"
z = np.r_[0.0, np.cumsum(0.5 * (dz_ds[1:] + dz_ds[:-1]) * np.diff(s))]


def height(radius: np.ndarray) -> np.ndarray:
    return np.interp(radius, R, z)


# Perspective view. Returns raw image-plane coordinates and camera depth.
Z_MID = 0.5 * height(R_FADE[1])


def project(x: np.ndarray, y: np.ndarray):
    zz = height(np.hypot(x, y)) - Z_MID
    xr = x * np.cos(AZIMUTH) - y * np.sin(AZIMUTH)
    yr = x * np.sin(AZIMUTH) + y * np.cos(AZIMUTH)
    depth = CAMERA + yr * np.cos(ELEVATION) - zz * np.sin(ELEVATION)
    up = yr * np.sin(ELEVATION) + zz * np.cos(ELEVATION)
    return CAMERA * xr / depth, -CAMERA * up / depth, depth


def alpha(half_width: np.ndarray) -> np.ndarray:
    """Opacity from the superellipse radius (|x|^p + |y|^p)^(1/p)."""
    t = np.clip((half_width - R_FADE[0]) / (R_FADE[1] - R_FADE[0]), 0.0, 1.0)
    return (1.0 - t * t * (3.0 - 2.0 * t)) ** 2  # squared: softer tail


# Z-buffer of the surface out to the fade square's corners, used to hide lines.
N = 1600
rho, phi = np.meshgrid(
    np.interp(np.linspace(0, height(R_FADE[1] * 1.5), 1500), z, R),
    np.linspace(0, 2 * np.pi, 3000),
)
X, Y, D = project((rho * np.cos(phi)).ravel(), (rho * np.sin(phi)).ravel())
assert D.min() > 0, "camera must stay above the drawn surface"
lo = np.array([X.min(), Y.min()])
span = max(X.max() - lo[0], Y.max() - lo[1]) * 1.001


def pixel(px: np.ndarray, py: np.ndarray):
    return (
        np.clip(((py - lo[1]) / span * N).astype(int), 0, N - 1),
        np.clip(((px - lo[0]) / span * N).astype(int), 0, N - 1),
    )


zbuf = np.full((N, N), np.inf)
np.minimum.at(zbuf, pixel(X, Y), D)
zbuf = minimum_filter(zbuf, size=3)

LEVELS = 40  # opacity is quantised so each polyline carries one value
strokes = []  # grid: (points, opacity)
solid = []  # horizon and orbit: (points, (colour, width))


def drag(x: np.ndarray, y: np.ndarray):
    """Carry grid points along with zero-angular-momentum observers.

    A ZAMO at Boyer-Lindquist radius r turns at the frame-dragging rate
    omega = 2 a r / ((r^2 + a^2)^2 - a^2 Delta) on the equator, so after
    DRAG_TIME its azimuth has grown by omega * DRAG_TIME. Grid points keep their
    circumferential radius R and are rotated by that angle.
    """
    radius = np.hypot(x, y)
    rr = np.interp(radius, R, r)  # R(r) is monotonic outside the horizon
    omega = 2 * a * rr / ((rr**2 + a**2) ** 2 - a**2 * (rr**2 - 2 * rr + a**2))
    angle = np.arctan2(y, x) + omega * DRAG_TIME
    return radius * np.cos(angle), radius * np.sin(angle)


def draw(x: np.ndarray, y: np.ndarray, fade: bool = True, style=None, hidden_style=None) -> None:
    """Add a curve on the surface; hidden parts are dropped or use ``hidden_style``."""
    px, py, depth = project(x, y)
    # Seen from straight above the surface z(R) cannot hide itself.
    visible = (depth <= zbuf[pixel(px, py)] + 0.25) | np.isclose(ELEVATION, np.pi / 2)
    if not fade:
        for run in np.split(np.arange(len(px)), np.where(np.diff(visible.astype(int)))[0] + 1):
            if len(run) > 1 and (visible[run[0]] or hidden_style):
                solid.append((np.c_[px[run], py[run]], style if visible[run[0]] else hidden_style))
        return
    norm = (np.abs(x) ** SQUARENESS + np.abs(y) ** SQUARENESS) ** (1 / SQUARENESS)
    level = np.round(alpha(norm) * LEVELS).astype(int)
    key = np.where(visible, level, 0)
    cuts = np.where(np.diff(key))[0] + 1
    for run in np.split(np.arange(len(px)), cuts):
        if key[run[0]] > 0:
            run = np.r_[run, run[-1] + 1] if run[-1] + 1 < len(px) else run
            if len(run) > 1:
                strokes.append((np.c_[px[run], py[run]], key[run[0]] / LEVELS))


extent = R_FADE[1] + CELL
t = np.linspace(-extent, extent, 700)
for c0 in np.arange(-extent, extent + 1e-9, CELL):
    for x, y in ((np.full_like(t, c0), t), (t, np.full_like(t, c0))):
        outside = np.hypot(x, y) >= 2.0  # the throat R = 2 is the horizon
        for run in np.split(np.arange(len(t)), np.where(np.diff(outside.astype(int)))[0] + 1):
            if outside[run[0]] and len(run) > 1:
                draw(*drag(x[run], y[run]))

# Outer horizon: the throat circle R = 2 at the bottom of the embedding.
ring = np.linspace(0, 2 * np.pi, 720)
draw(2.0 * np.cos(ring), 2.0 * np.sin(ring), fade=False,
     style=(GREY, 1.1, None), hidden_style=(GREY, 0.6, "2 2"))

# A bound prograde equatorial geodesic integrated by RelatiPy, lifted onto the
# surface: Boyer-Lindquist (r, phi) maps to (R(r) cos phi, R(r) sin phi).
elements = KerrOrbitalElements(
    p=ORBIT_RADIUS * bh.r_g, e=0.0, x=1.0,
    q_r0=0 * u.rad, q_theta0=0 * u.rad, q_phi0=0 * u.rad,
)
# Native accepted steps are used directly (dense via max_step): Solution.at()
# currently misinterpolates between widely spaced steps on this orbit.
t_g = (bh.r_g / c).to(u.s)
sol = bh.orbit(elements=elements).solve(
    tau_span=(0 * u.s, 120 * t_g),
    method="dop853", rtol=1e-11, atol=1e-13, max_step=0.05 * t_g,
)
r_orbit = (sol.R / bh.r_g).decompose().value
phi_orbit = sol.Phi.to_value(u.rad)
stop = int(np.argmax(phi_orbit - phi_orbit[0] >= 2 * np.pi)) + 1  # one turn
R_orbit = np.sqrt(r_orbit**2 + a**2 + 2 * a**2 / r_orbit)[:stop]
# From straight above nothing hides the orbit, so it is projected directly.
px, py, _ = project(R_orbit * np.cos(phi_orbit[:stop]), R_orbit * np.sin(phi_orbit[:stop]))
trail = np.c_[px, py]

# Fit the visible drawing into the 256 px box.
# Framing ignores the faintest tails (opacity < 0.04), which are practically invisible.
allpts = np.vstack([pts for pts, op in strokes if op >= 0.04] + [pts for pts, _ in solid] + [trail])
bmin, bmax = allpts.min(axis=0), allpts.max(axis=0)
fit = 244 / (bmax - bmin).max()
offset = 128 - 0.5 * (bmin + bmax) * fit
MARK = "".join(
    '<polyline points="'
    + " ".join(f"{x:.1f},{y:.1f}" for x, y in pts * fit + offset)
    + f'" fill="none" stroke="{GREY}" stroke-opacity="{op:.3f}" stroke-width="{WIDTH}" '
    'stroke-linecap="round" stroke-linejoin="round"/>'
    for pts, op in strokes
) + "".join(
    '<polyline points="'
    + " ".join(f"{x:.1f},{y:.1f}" for x, y in pts * fit + offset)
    + f'" fill="none" stroke="{colour}" stroke-width="{width}" '
    + (f'stroke-dasharray="{dash}" ' if dash else "")
    + 'stroke-linecap="round" stroke-linejoin="round"/>'
    for pts, (colour, width, dash) in solid
)

# Orbit trail: opacity rises linearly from 0.2 at the start of the turn to 1 at
# its midpoint and falls back to 0.2 at the end; butt caps keep adjacent chunks
# from overlapping.
TRAIL_CHUNKS = 64
edges = np.linspace(0, len(trail) - 1, TRAIL_CHUNKS + 1).astype(int)
fitted = trail * fit + offset
MARK += "".join(
    '<polyline points="'
    + " ".join(f"{x:.1f},{y:.1f}" for x, y in fitted[i0 : i1 + 1])
    + f'" fill="none" stroke="{ORBIT}" stroke-opacity="{0.2 + 0.8 * (1 - abs(2 * (k + 0.5) / TRAIL_CHUNKS - 1)):.3f}" '
    'stroke-width="1.5" stroke-linecap="butt" stroke-linejoin="round"/>'
    for k, (i0, i1) in enumerate(zip(edges[:-1], edges[1:]))
)
# Synthetic "observations" of the orbit, as a telescope would record them:
# positions in the image plane with Gaussian scatter and 1-sigma error bars
# along the image axes (fixed seed). Scatter and bars are applied after the
# projection, so they are not distorted by the embedding or the perspective;
# sigma in r_g is converted with the image scale at the orbit's radius. They
# are illustrative, not real data.
rng = np.random.default_rng(7)
samples = np.linspace(0, stop - 1, DATA_POINTS + 2).astype(int)[1:-1]
true = fitted[samples]
# Image scale at the orbit: pixels per unit of circumferential radius R there.
R_probe = R_orbit[0] + np.array([0.0, 0.01])
px_per_rg = abs(np.diff(project(R_probe, np.zeros(2))[0])[0]) / 0.01 * fit
sigma = rng.uniform(0.15, 0.35, size=(DATA_POINTS, 2)) * px_per_rg  # px
# Scatter is truncated at 1.2 sigma so no single point strays far from the orbit.
observed = true + sigma * np.clip(rng.standard_normal((DATA_POINTS, 2)), -1.2, 1.2)
for (cx, cy), (sx, sy) in zip(observed, sigma):
    MARK += (
        f'<path d="M{cx - sx:.1f},{cy:.1f}H{cx + sx:.1f}M{cx:.1f},{cy - sy:.1f}V{cy + sy:.1f}" '
        f'stroke="{DATA}" stroke-width="1.2" stroke-linecap="round"/>'
        f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="3.0" fill="{DATA}"/>'
    )


def wordmark(text: str, x0: float, baseline: float, size: float) -> tuple[str, float]:
    font = TTFont(FONT)
    glyphs, cmap = font.getGlyphSet(), font.getBestCmap()
    k = size / font["head"].unitsPerEm
    paths, x = [], x0
    for ch in text:
        name = cmap[ord(ch)]
        pen = SVGPathPen(glyphs)
        glyphs[name].draw(TransformPen(pen, (k, 0, 0, -k, x, baseline)))
        paths.append(pen.getCommands())
        x += glyphs[name].width * k
    return " ".join(paths), x


def write(name: str, width: float, body: str, height: int = 256) -> None:
    (OUT / name).write_text(
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width:.0f} {height}" '
        f'role="img" aria-label="RelatiPy">{body}</svg>\n'
    )


def paint(svg: str, palette: dict) -> str:
    for token, colour in palette.items():
        svg = svg.replace(token, colour)
    return svg


# In the banner the mark is drawn at 125 % (320 px tall), taller than the
# wordmark, which is centred vertically beside it.
text, end = wordmark("RelatiPy", 336, 198, 118)
BANNER = f'<g transform="scale(1.25)">{MARK}</g><path d="{text}" fill="@text@"/>'
for suffix, palette in PALETTES.items():
    write(f"logo-mark{suffix}.svg", 256, paint(MARK, palette))
    write(f"logo{suffix}.svg", end + 24, paint(BANNER, palette), height=320)
print(f"r_+ = {r_plus:.6f}, throat R = {R[0]:.6f}; wrote logo[-dark].svg, logo-mark[-dark].svg")
