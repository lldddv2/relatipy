"""Mock GRAVITY/GRAVITY+ astrometry and ERIS/SINFONI/ELT radial velocities.

Python API
    ds = make_dataset(Scenario(...), forward=KerrForward())
    write_dataset(ds, outdir)

CLI
    python3 mock.py --list
    python3 mock.py A_rv1 B_rv1            # selected scenarios
    python3 mock.py --all                  # every predefined scenario
    python3 mock.py A_rv1 --seed 7 --chi 0.9 --spin-dir par --out /tmp/x
    python3 mock.py A_rv1 --model mymodule:make_model   # pluggable forward model

Output: data/mock/<scenario>/
    <star>_astrometry.csv   t, ra, dec, sigma_ra, sigma_dec, instrument, origin  (mas, yr)
    <star>_rv.csv           t, v, sigma, instrument, origin                      (km/s, yr)
    nonoise/<star>_*.csv    the same, noise-free
    scenario.json           configuration, truth (elements, M, R0, spin), seed, model, git commit

The truth for every star is the published element set of stars.json (S301 solution
A or B), M and R0 from GRAVITY 2024, spin chi along a named direction.  Past epochs
are mock values at the real epochs: the real measured positions are in data/ and are
not mixed into the mocks.  No reference-frame offsets/drifts, no confusion, white
Gaussian noise with the scheduled 1-sigma errors.
"""
from __future__ import annotations

import argparse
import csv
import importlib
import json
import subprocess
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

import schedules as sch
from forward import KerrForward, load_stars, spin_vector

HERE = Path(__file__).resolve().parent
MOCK_DIR = HERE / "data" / "mock"
STAR_ORDER = ("S301A", "S301B", "S2", "S29", "S38", "S55")  # fixes the per-star RNG streams


@dataclass
class Scenario:
    name: str
    description: str = ""
    schedule: dict = field(default_factory=dict)   # ScheduleConfig overrides
    chi: float = 0.7
    spin_dir: str | list = "b45l45"                # SPIN_DIRS name (S301A frame) or (Dec, RA, away)
    seed: int = 20261025


def _sc(name, desc, **kw):
    return Scenario(name=name, description=desc, schedule=kw)


SCENARIOS = {s.name: s for s in [
    _sc("past", "Data to 2026.0 only: S301 19 real epochs (207 uas, no RV); friends real/plausible past.",
        include_future=False),
    _sc("A_rvnone", "Truth S301A, GRAVITY+ 2026-2041, no S301 RV.", s301_sigma_rv=None),
    _sc("A_rv20", "Truth S301A, S301 RV 20 km/s from 2030 (HARMONI-like).", s301_sigma_rv=20.0),
    _sc("A_rv5", "Truth S301A, S301 RV 5 km/s from 2030.", s301_sigma_rv=5.0),
    _sc("A_rv1", "Truth S301A, S301 RV 1 km/s from 2030 (ELT/MICADO).", s301_sigma_rv=1.0),
    _sc("A_rv1_2035", "As A_rv1 but data end 2035.0 (one future pericentre).", s301_sigma_rv=1.0, t_end=2035.0),
    _sc("A_rv1_longseason", "As A_rv1 with the 1 Feb - 6 Nov window: the 2031.81 pericentre is observable.",
        s301_sigma_rv=1.0, season=list(sch.SEASON_LONG)),
    _sc("A_rv1_s301at30", "As A_rv1 with S301 future astrometry at 30 uas.", s301_sigma_rv=1.0, s301_sigma_ast=0.030),
    _sc("B_rv1", "Truth S301B (z-mirror of A), S301 RV 1 km/s from 2030.", s301_sigma_rv=1.0, s301_solution="B"),
    _sc("B_rvnone", "Truth S301B, no S301 RV.", s301_sigma_rv=None, s301_solution="B"),
]}


def schedule_config(sc: Scenario) -> sch.ScheduleConfig:
    kw = dict(sc.schedule)
    if "season" in kw:
        kw["season"] = tuple(kw["season"])
    return sch.ScheduleConfig(**kw)


def truth_bh(sc: Scenario, stars: dict) -> dict:
    g = stars["global"]
    vec = spin_vector(sc.spin_dir, stars["stars"]["S301A"])
    return dict(M_1e6Msun=g["M_1e6Msun"], R0_kpc=g["R0_kpc"], chi=float(sc.chi), spin_vec=vec.tolist())


def make_dataset(sc: Scenario, forward=None, stars=None) -> dict:
    """Return {'scenario', 'bh', 'stars': {name: {'ast': {...}, 'rv': {...}, 'el': {...}}}}."""
    forward = forward or KerrForward()
    stars = stars or load_stars()
    cfg = schedule_config(sc)
    bh = truth_bh(sc, stars)
    out = {}
    for name in sch.star_names(cfg):
        el = stars["stars"][name]
        a, r = sch.build(name, cfg, stars)
        ra, de, v = forward(el, a["t"], r["t"], bh)
        rng = np.random.default_rng([sc.seed, STAR_ORDER.index(name)])
        n_ra, n_de = rng.standard_normal(a["t"].size), rng.standard_normal(a["t"].size)
        n_v = rng.standard_normal(r["t"].size)
        a.update(ra_true=ra, dec_true=de, ra=ra + a["sigma_ra"] * n_ra, dec=de + a["sigma_dec"] * n_de)
        r.update(v_true=v, v=v + r["sigma"] * n_v)
        out[name] = dict(ast=a, rv=r, el=el)
    return dict(scenario=sc, cfg=cfg, bh=bh, stars=out)


def _write_csv(path, cols, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        w.writerows(rows)


def _fmt(x, nd):
    return f"{x:.{nd}f}"


def write_dataset(ds: dict, outdir: Path, model_name: str = "forward.KerrForward") -> None:
    outdir = Path(outdir)
    (outdir / "nonoise").mkdir(parents=True, exist_ok=True)
    for name, d in ds["stars"].items():
        a, r = d["ast"], d["rv"]
        for sub, kra, kde, kv in (("", "ra", "dec", "v"), ("nonoise", "ra_true", "dec_true", "v_true")):
            base = outdir / sub
            _write_csv(base / f"{name}_astrometry.csv",
                       ["t", "ra", "dec", "sigma_ra", "sigma_dec", "instrument", "origin"],
                       [[_fmt(t, 6), _fmt(x, 6), _fmt(y, 6), _fmt(sx, 4), _fmt(sy, 4), ins, org]
                        for t, x, y, sx, sy, ins, org in zip(a["t"], a[kra], a[kde], a["sigma_ra"],
                                                            a["sigma_dec"], a["instrument"], a["origin"])])
            _write_csv(base / f"{name}_rv.csv", ["t", "v", "sigma", "instrument", "origin"],
                       [[_fmt(t, 6), _fmt(v, 4), _fmt(s, 3), ins, org]
                        for t, v, s, ins, org in zip(r["t"], r[kv], r["sigma"], r["instrument"], r["origin"])])
    try:
        commit = subprocess.run(["git", "-C", str(HERE), "rev-parse", "--short", "HEAD"],
                                capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        commit = "unknown"
    sc = ds["scenario"]
    meta = dict(
        scenario=asdict(sc), schedule_config=asdict(ds["cfg"]), bh_truth=ds["bh"],
        elements_truth={k: v["el"] for k, v in ds["stars"].items()},
        counts={k: dict(astrometry=int(v["ast"]["t"].size), rv=int(v["rv"]["t"].size)) for k, v in ds["stars"].items()},
        forward_model=model_name, relatipy_version=_relatipy_version(), git_commit_at_generation=commit,
        units=dict(t="decimal Julian yr (observed time)", ra_dec="mas, RA positive east, Dec positive north, "
                   "relative to Sgr A*", v="km/s, positive = receding (redshift-equivalent)"),
        noise="independent Gaussian, sigma as listed; RNG numpy default_rng([seed, star index in "
              + str(list(STAR_ORDER)) + "])",
    )
    with open(outdir / "scenario.json", "w") as fh:
        json.dump(meta, fh, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))


def _relatipy_version():
    try:
        from importlib.metadata import version
        return version("relatipy")
    except Exception:
        return "unknown"


def _load_model(spec: str | None):
    if not spec:
        return KerrForward(), "forward.KerrForward"
    mod, fn = spec.split(":")
    return getattr(importlib.import_module(mod), fn)(), spec


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("scenarios", nargs="*")
    p.add_argument("--all", action="store_true")
    p.add_argument("--list", action="store_true")
    p.add_argument("--seed", type=int)
    p.add_argument("--chi", type=float)
    p.add_argument("--spin-dir")
    p.add_argument("--model", help="module:factory returning a forward-model callable")
    p.add_argument("--out", type=Path, default=MOCK_DIR)
    a = p.parse_args(argv)
    if a.list:
        for s in SCENARIOS.values():
            print(f"{s.name:18s} {s.description}")
        return
    names = list(SCENARIOS) if a.all else a.scenarios
    if not names:
        p.error("give scenario names, --all or --list")
    model, model_name = _load_model(a.model)
    stars = load_stars()
    for n in names:
        sc = SCENARIOS[n]
        if a.seed is not None: sc.seed = a.seed
        if a.chi is not None: sc.chi = a.chi
        if a.spin_dir is not None: sc.spin_dir = a.spin_dir
        ds = make_dataset(sc, model, stars)
        write_dataset(ds, a.out / n, model_name)
        tot = {k: (v["ast"]["t"].size, v["rv"]["t"].size) for k, v in ds["stars"].items()}
        print(f"{n}: " + ", ".join(f"{k} {x}/{y}" for k, (x, y) in tot.items()) + f"  -> {a.out / n}")


if __name__ == "__main__":
    main()
