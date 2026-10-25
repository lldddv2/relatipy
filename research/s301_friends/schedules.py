"""Past and future observation schedules for S301, S2, S29, S38 and S55.

A schedule is a dict of 1-D numpy arrays:
    astrometry: t [yr], sigma_ra [mas], sigma_dec [mas], instrument, origin
    rv:         t [yr], sigma [km/s], instrument, origin
origin is one of
    real_epoch real epoch and real error (S2 CSV; S301 Extended Data Table 1)
    flagged    real-looking point of doubtful provenance (S301 2023.577 digitized)
    past_synth plausible past epoch built to match the published epoch counts
    future     planned GRAVITY+ / ERIS / ELT epoch
    peri       extra epoch of a pericentre campaign
    peri_edge  pericentre-campaign epoch moved to the nearest season edge because
               the pericentre itself is outside the visibility window

All epochs (past synthetic and future) are restricted to the Galactic Centre
visibility window (default 1 Mar - 30 Sep, ScheduleConfig.season).  Grids are
deterministic: n epochs are spread evenly over the observable time of an interval.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np

from forward import load_stars, period_yr

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
S2_CSV = DATA / "s2_astrometry_gravity2020_fig1.csv"
S301_CSV = DATA / "s301_radec_digitized_nature2026_fig.csv"
YEAR_D = 365.25

# arXiv:2607.12664 Extended Data Table 1: the 19 GRAVITY epochs of S301
S301_ED_TABLE1 = np.array([
    2017.2389, 2021.5759, 2023.2604, 2023.4162, 2023.4298, 2023.5063, 2024.2325,
    2024.2380, 2024.2405, 2024.3281, 2024.4073, 2024.4756, 2024.4837, 2024.5438,
    2025.2070, 2025.3602, 2025.4393, 2025.5159, 2025.5978])
S301_SIGMA_PAST = 0.207  # mas per coordinate (pixel-grid noise floor, Methods p.19)
S301_UNMATCHED_T = 2023.577267  # digitized point with no Extended Data Table 1 epoch

SEASON_MAR_SEP = (60, 274)   # day of year [start, end): 1 Mar .. 30 Sep
SEASON_LONG = (32, 310)      # 1 Feb .. 6 Nov (arXiv:2609.06936 blackout 6 Nov - 1 Feb)

FRIENDS = ("S2", "S29", "S38", "S55")


@dataclass
class ScheduleConfig:
    s301_solution: str = "A"
    include_past: bool = True
    include_future: bool = True
    t_future: float = 2026.0           # start of the GRAVITY+ era
    t_end: float = 2041.0
    season: tuple = SEASON_MAR_SEP
    # astrometry (GRAVITY+), mas
    ast_per_yr: float = 10.0
    s301_sigma_ast: float = 0.100
    friend_sigma_ast: dict = field(default_factory=lambda: {"S2": 0.030, "S29": 0.040,
                                                            "S38": 0.040, "S55": 0.050})
    s301_n_peri_ast: int = 10
    s301_peri_half: float = 0.06       # yr (+-22 d)
    friend_n_peri_ast: int = 10
    friend_peri_half: float = 0.15     # yr
    # S301 RV (ELT): None/0 -> no RV
    s301_sigma_rv: float | None = None
    s301_rv_start: float = 2030.0
    s301_rv_per_yr: float = 4.0
    s301_n_peri_rv: int = 10
    # friends RV (ERIS-like)
    friend_sigma_rv: float = 10.0
    friend_rv_per_yr: float = 3.0
    friend_n_peri_rv: int = 5
    # past data options
    s301_include_unmatched: bool = False   # keep the 2023.577 point (flagged)
    s2_drop_duplicates: bool = True        # two GRAVITY rows of the S2 CSV are exact duplicates
    past_gravity_sigma: float = 0.050      # mas, synthetic past GRAVITY epochs
    past_naco_sigma: float = 0.50          # mas, synthetic past NACO epochs (friends)


# ----------------------------------------------------------------------------- time grids
def _season_bounds(season):
    return (season[0] - 1) / YEAR_D, (season[1] - 1) / YEAR_D


def in_season(t, season=SEASON_MAR_SEP):
    t = np.asarray(t, float)
    s0, s1 = _season_bounds(season)
    fr = t - np.floor(t)
    return (fr >= s0) & (fr < s1)


def _observable_intervals(t0, t1, season):
    s0, s1 = _season_bounds(season)
    out = []
    for y in range(int(np.floor(t0)) - 1, int(np.ceil(t1)) + 1):
        a, b = max(t0, y + s0), min(t1, y + s1)
        if b > a:
            out.append((a, b))
    return out


def season_grid(t0, t1, n, season=SEASON_MAR_SEP):
    """n epochs evenly spread over the observable time inside [t0, t1]."""
    n = int(round(n))
    iv = _observable_intervals(t0, t1, season)
    if n <= 0 or not iv:
        return np.empty(0)
    lens = np.array([b - a for a, b in iv])
    cum = np.concatenate([[0.0], np.cumsum(lens)])
    s = (np.arange(n) + 0.5) * cum[-1] / n
    k = np.clip(np.searchsorted(cum, s, side="right") - 1, 0, len(iv) - 1)
    return np.array([iv[j][0] + (x - cum[j]) for j, x in zip(k, s)])


def regular(t0, t1, per_yr, season=SEASON_MAR_SEP):
    """~per_yr epochs per calendar year, all within the season."""
    s0, s1 = _season_bounds(season)
    obs = sum(b - a for a, b in _observable_intervals(t0, t1, season))
    return season_grid(t0, t1, per_yr * obs / (s1 - s0), season)


def pericentres(el, M6, R0, t0, t1):
    P = period_yr(el["a_mas"], M6, R0)
    k = np.arange(np.floor((t0 - el["t_p"]) / P) - 1, np.ceil((t1 - el["t_p"]) / P) + 2)
    tp = el["t_p"] + k * P
    return tp[(tp >= t0) & (tp <= t1)]


def peri_campaign(tp, n, half, season, t0, t1):
    """n epochs on a quarter-offset grid in [tp - half, tp + half] (nearest epoch ~1/4 spacing
    from tp).  Epochs outside the season are moved into the last/first `half` yr of the
    nearest season edge (origin 'peri_edge')."""
    tt = tp - half + (np.arange(n) + 0.25) * 2 * half / n
    ok = in_season(tt, season) & (tt >= t0) & (tt <= t1)
    t_in, n_out = tt[ok], int((~ok & (tt >= t0) & (tt <= t1)).sum())
    t_edge = np.empty(0)
    if n_out:
        s0, s1 = _season_bounds(season)
        y = np.floor(tp)
        cands = np.array([y - 1 + s1, y + s1, y + s0, y + 1 + s0])
        edge = cands[np.argmin(np.abs(cands - tp))]
        is_end = np.isclose(edge - np.floor(edge), s1)
        lo, hi = (edge - half, edge) if is_end else (edge, edge + half)
        t_edge = lo + (np.arange(n_out) + 0.5) * (hi - lo) / n_out
        t_edge = t_edge[(t_edge >= t0) & (t_edge <= t1) & in_season(t_edge, season)]
    return t_in, t_edge


# ----------------------------------------------------------------------------- containers
def _ast(t, sig, instrument, origin, sig_dec=None):
    t = np.atleast_1d(np.asarray(t, float))
    sr = np.broadcast_to(np.asarray(sig, float), t.shape).copy()
    sd = sr.copy() if sig_dec is None else np.broadcast_to(np.asarray(sig_dec, float), t.shape).copy()
    return dict(t=t, sigma_ra=sr, sigma_dec=sd,
                instrument=np.array([instrument] * t.size if isinstance(instrument, str) else instrument, dtype=object),
                origin=np.array([origin] * t.size, dtype=object))


def _rv(t, sig, instrument, origin):
    t = np.atleast_1d(np.asarray(t, float))
    return dict(t=t, sigma=np.broadcast_to(np.asarray(sig, float), t.shape).copy(),
                instrument=np.array([instrument] * t.size, dtype=object),
                origin=np.array([origin] * t.size, dtype=object))


def _cat(parts, kind):
    keys = ["t", "sigma_ra", "sigma_dec", "instrument", "origin"] if kind == "ast" else \
           ["t", "sigma", "instrument", "origin"]
    parts = [p for p in parts if p["t"].size]
    if not parts:
        return {k: np.empty(0, dtype=object if k in ("instrument", "origin") else float) for k in keys}
    out = {k: np.concatenate([p[k] for p in parts]) for k in keys}
    o = np.argsort(out["t"], kind="stable")
    return {k: v[o] for k, v in out.items()}


# ----------------------------------------------------------------------------- real data
def load_s2_csv(drop_duplicates=True):
    rows = list(csv.DictReader(open(S2_CSV)))
    seen, out = set(), []
    for r in rows:
        key = tuple(r.values())
        if drop_duplicates and key in seen:
            continue
        seen.add(key)
        out.append(r)
    f = lambda k: np.array([float(r[k]) for r in out])
    return dict(t=f("t_yr"), ra=f("ra_mas"), dec=f("dec_mas"), sigma_ra=f("ra_err_mas"),
                sigma_dec=f("dec_err_mas"), instrument=np.array([r["instrument"] for r in out], dtype=object))


def load_s301_digitized():
    """Digitized S301 positions matched to Extended Data Table 1.  Returns arrays with t
    (table epoch; digitized year for the unmatched point), ra, dec, digitized error-bar
    half-widths, sigma = 0.207 mas, and flag ('matched' / 'unmatched_2023.577')."""
    rows = list(csv.DictReader(open(S301_CSV)))
    ty = np.array([float(r["year_digitized"]) for r in rows])
    t, flag = ty.copy(), np.empty(ty.size, dtype=object)
    used = set()
    for k, x in enumerate(ty):
        j = int(np.argmin(np.abs(S301_ED_TABLE1 - x)))
        if abs(S301_ED_TABLE1[j] - x) < 0.004 and j not in used:
            t[k], flag[k] = S301_ED_TABLE1[j], "matched"
            used.add(j)
        else:
            flag[k] = "unmatched_2023.577"
    assert len(used) == 19 and (flag != "matched").sum() == 1
    f = lambda key: np.array([float(r[key]) for r in rows])
    return dict(t=t, t_digitized=ty, ra=f("ra_mas"), dec=f("dec_mas"),
                halfwidth_ra=f("ra_errorbar_halfwidth_mas"), halfwidth_dec=f("dec_errorbar_halfwidth_mas"),
                sigma=np.full(ty.size, S301_SIGMA_PAST), flag=flag)


# ----------------------------------------------------------------------------- past schedules
# Published epoch counts (GRAVITY 2024, arXiv:2409.12261 p.2), data to 2022.7.
PAST_COUNTS = {
    "S29": dict(naco=(66, 2002.3, 2019.7), gravity=(29, 2019.6, 2022.7),
                rv=[("SINFONI", 17, 2007.6, 2019.7, 20.0), ("GNIRS", 2, 2021.2, 2021.45, 20.0)]),
    "S38": dict(naco=(110, 2004.2, 2018.6), gravity=(23, 2021.2, 2022.7),
                rv=[("SINFONI", 8, 2008.3, 2019.6, 25.0), ("Keck", 2, 2010.3, 2012.6, 30.0),
                    ("ERIS", 1, 2022.3, 2022.45, 10.0)]),
    "S55": dict(naco=(42, 2004.5, 2013.6), gravity=(27, 2021.2, 2022.6),
                rv=[("SINFONI", 2, 2014.3, 2014.7, 30.0)]),
}
EXT_START = 2022.75   # unpublished GRAVITY/ERIS monitoring 2022.75 - t_future (plausible)
EXT_GRAVITY_PER_YR = 6.0
EXT_ERIS_PER_YR = 3.0


def past_s2(cfg):
    d = load_s2_csv(cfg.s2_drop_duplicates)
    ast = [_ast(d["t"], d["sigma_ra"], d["instrument"], "real_epoch", d["sigma_dec"])]
    ng = int((d["instrument"] == "GRAVITY").sum())
    # GRAVITY 2024 has 82 GRAVITY epochs to 2022.7; the CSV stops at 2019.74
    ast.append(_ast(season_grid(2019.75, 2022.7, 82 - ng, cfg.season), cfg.past_gravity_sigma, "GRAVITY", "past_synth"))
    ast.append(_ast(regular(EXT_START, cfg.t_future, EXT_GRAVITY_PER_YR, cfg.season), cfg.past_gravity_sigma, "GRAVITY", "past_synth"))
    # RV: 92 SINFONI, 3 Keck, 2 NACO, 4 GNIRS, 3 ERIS (span 1992.2-2022.7); errors assumed
    rv = [_rv(season_grid(2003.3, 2019.7, 62, cfg.season), 15.0, "SINFONI", "past_synth"),
          _rv(season_grid(2018.15, 2018.75, 30, cfg.season), 10.0, "SINFONI", "past_synth"),
          _rv(season_grid(2000.4, 2002.6, 3, cfg.season), 30.0, "Keck", "past_synth"),
          _rv(season_grid(2002.3, 2003.3, 2, cfg.season), 50.0, "NACO", "past_synth"),
          _rv(season_grid(2018.2, 2018.7, 4, cfg.season), 20.0, "GNIRS", "past_synth"),
          _rv(season_grid(2022.2, 2022.7, 3, cfg.season), 10.0, "ERIS", "past_synth"),
          _rv(regular(EXT_START, cfg.t_future, EXT_ERIS_PER_YR, cfg.season), 10.0, "ERIS", "past_synth")]
    return _cat(ast, "ast"), _cat(rv, "rv")


def past_friend(name, cfg):
    c = PAST_COUNTS[name]
    ast = [_ast(season_grid(c["naco"][1], c["naco"][2], c["naco"][0], cfg.season), cfg.past_naco_sigma, "NACO", "past_synth"),
           _ast(season_grid(c["gravity"][1], c["gravity"][2], c["gravity"][0], cfg.season), cfg.past_gravity_sigma, "GRAVITY", "past_synth"),
           _ast(regular(EXT_START, cfg.t_future, EXT_GRAVITY_PER_YR, cfg.season), cfg.past_gravity_sigma, "GRAVITY", "past_synth")]
    rv = [_rv(season_grid(t0, t1, n, cfg.season), s, ins, "past_synth") for ins, n, t0, t1, s in c["rv"]]
    rv.append(_rv(regular(EXT_START, cfg.t_future, EXT_ERIS_PER_YR, cfg.season), 10.0, "ERIS", "past_synth"))
    return _cat(ast, "ast"), _cat(rv, "rv")


def past_s301(cfg):
    parts = [_ast(S301_ED_TABLE1, S301_SIGMA_PAST, "GRAVITY", "real_epoch")]
    if cfg.s301_include_unmatched:
        parts.append(_ast([S301_UNMATCHED_T], S301_SIGMA_PAST, "GRAVITY", "flagged"))
    return _cat(parts, "ast"), _cat([], "rv")


# ----------------------------------------------------------------------------- future schedules
def _future(el, bh, t0, t1, per_yr, n_peri, half, season):
    t_reg = regular(t0, t1, per_yr, season)
    t_in, t_edge = [], []
    for tp in pericentres(el, bh["M_1e6Msun"], bh["R0_kpc"], t0 - half, t1 + half):
        a, b = peri_campaign(tp, n_peri, half, season, t0, t1)
        t_in.append(a); t_edge.append(b)
    cat = lambda x: np.concatenate(x) if x else np.empty(0)
    return t_reg, cat(t_in), cat(t_edge)


def future_s301(el, bh, cfg):
    tr, tpk, te = _future(el, bh, cfg.t_future, cfg.t_end, cfg.ast_per_yr, cfg.s301_n_peri_ast,
                          cfg.s301_peri_half, cfg.season)
    s = cfg.s301_sigma_ast
    ast = _cat([_ast(tr, s, "GRAVITY+", "future"), _ast(tpk, s, "GRAVITY+", "peri"),
                _ast(te, s, "GRAVITY+", "peri_edge")], "ast")
    rv = _cat([], "rv")
    if cfg.s301_sigma_rv:
        ins = "ELT/HARMONI" if cfg.s301_sigma_rv >= 15 else "ELT/MICADO"
        t0 = max(cfg.s301_rv_start, cfg.t_future)
        tr, tpk, te = _future(el, bh, t0, cfg.t_end, cfg.s301_rv_per_yr, cfg.s301_n_peri_rv,
                              cfg.s301_peri_half, cfg.season)
        sr = cfg.s301_sigma_rv
        rv = _cat([_rv(tr, sr, ins, "future"), _rv(tpk, sr, ins, "peri"), _rv(te, sr, ins, "peri_edge")], "rv")
    return ast, rv


def future_friend(name, el, bh, cfg):
    tr, tpk, te = _future(el, bh, cfg.t_future, cfg.t_end, cfg.ast_per_yr, cfg.friend_n_peri_ast,
                          cfg.friend_peri_half, cfg.season)
    s = cfg.friend_sigma_ast[name]
    ast = _cat([_ast(tr, s, "GRAVITY+", "future"), _ast(tpk, s, "GRAVITY+", "peri"),
                _ast(te, s, "GRAVITY+", "peri_edge")], "ast")
    tr, tpk, te = _future(el, bh, cfg.t_future, cfg.t_end, cfg.friend_rv_per_yr, cfg.friend_n_peri_rv,
                          cfg.friend_peri_half, cfg.season)
    sr = cfg.friend_sigma_rv
    rv = _cat([_rv(tr, sr, "ERIS", "future"), _rv(tpk, sr, "ERIS", "peri"), _rv(te, sr, "ERIS", "peri_edge")], "rv")
    return ast, rv


# ----------------------------------------------------------------------------- public API
def star_names(cfg):
    return ("S301" + cfg.s301_solution,) + FRIENDS


def build(name, cfg=None, stars=None):
    """Return (astrometry, rv) schedules of one star."""
    cfg = cfg or ScheduleConfig()
    stars = stars or load_stars()
    g = stars["global"]
    bh = dict(M_1e6Msun=g["M_1e6Msun"], R0_kpc=g["R0_kpc"])
    el = stars["stars"][name]
    parts_a, parts_r = [], []
    if cfg.include_past:
        a, r = past_s301(cfg) if name.startswith("S301") else past_s2(cfg) if name == "S2" else past_friend(name, cfg)
        parts_a.append(a); parts_r.append(r)
    if cfg.include_future:
        a, r = future_s301(el, bh, cfg) if name.startswith("S301") else future_friend(name, el, bh, cfg)
        a = {k: v[a["t"] <= cfg.t_end] for k, v in a.items()}
        r = {k: v[r["t"] <= cfg.t_end] for k, v in r.items()}
        parts_a.append(a); parts_r.append(r)
    return _cat(parts_a, "ast"), _cat(parts_r, "rv")


def summary(cfg=None):
    """Text table: epoch counts per star / instrument / origin, and pericentres."""
    cfg = cfg or ScheduleConfig()
    stars = load_stars()
    g = stars["global"]
    lines = []
    for name in star_names(cfg):
        a, r = build(name, cfg, stars)
        el = stars["stars"][name]
        tps = pericentres(el, g["M_1e6Msun"], g["R0_kpc"], 2026.0, cfg.t_end)
        lines.append(f"{name}: pericentres 2026-{cfg.t_end:g}: " +
                     (", ".join(f"{x:.2f}{'' if in_season(x, cfg.season) else ' (out of season)'}" for x in tps) or "none"))
        for lab, d, sk in (("ast", a, "sigma_ra"), ("rv", r, "sigma")):
            keys = sorted(set(zip(d["instrument"], d["origin"])), key=lambda k: d["t"][(d["instrument"] == k[0]) & (d["origin"] == k[1])].min())
            for ins, org in keys:
                m = (d["instrument"] == ins) & (d["origin"] == org)
                lines.append(f"  {lab:3s} {ins:12s} {org:10s} n={m.sum():4d}  {d['t'][m].min():.2f}-{d['t'][m].max():.2f}"
                             f"  sigma {np.median(d[sk][m]):.3g}")
    return "\n".join(lines)


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description="Print schedule summaries.")
    p.add_argument("--rv", type=float, default=None, help="S301 RV error [km/s] from 2030")
    p.add_argument("--t-end", type=float, default=2041.0)
    p.add_argument("--season", choices=["mar-sep", "long"], default="mar-sep")
    a = p.parse_args()
    print(summary(ScheduleConfig(s301_sigma_rv=a.rv, t_end=a.t_end,
                                 season=SEASON_MAR_SEP if a.season == "mar-sep" else SEASON_LONG)))
