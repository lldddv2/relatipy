# S301 and friends: inputs and mock data

These are the inputs for a joint fit of S301 with S2, S29, S38 and S55, using exact Kerr plus
(later) an extended-mass disk, to forecast the spin of Sgr A*. The data are meant to look like
GRAVITY/GRAVITY+ astrometry, ERIS/SINFONI radial velocities and, later, ELT/MICADO radial
velocities. The background is in `../purposes_2026-10/HANDOFF_S301_DISK.md`,
`../purposes_2026-10/literature/extract_papers_s301_disk.md` and the prototype in
`../purposes_2026-10/disk/` (`stars.py`, `realism_past.py`, `REPORT.md`).

Python: `/usr/bin/python3`, which has relatipy 1.2.0, numpy, scipy and astropy. Pandas is not needed.

| file | content |
|---|---|
| `data/` | the two real CSVs (S2 astrometry, S301 digitized positions), `data/README.md` (provenance, columns, units, flags) and `data/mock/` |
| `stars.json` | elements with 1-sigma errors, element (osculation) epoch and source for S301A, S301B, S2, S29, S38 and S55; global M and R0 with errors; NACO frame parameters; the A/B mirror note |
| `forward.py` | forward-model interface, default `KerrForward` (relatipy `KerrMcmcModel.get_ra_dec_vr`), Kepler helpers, spin directions |
| `schedules.py` | past and future observation schedules (`ScheduleConfig`, `build(star, cfg)`) and the real-data loaders |
| `mock.py` | dataset generator and CLI, with predefined scenarios |
| `validate.py` → `validation.md` | conventions check and residuals against the real S2 and S301 data |

## Regenerate

```
cd research/s301_friends
/usr/bin/python3 validate.py            # writes validation.md (~5 s)
/usr/bin/python3 schedules.py --rv 1    # prints the schedule summary (counts, spans, errors)
/usr/bin/python3 mock.py --list
/usr/bin/python3 mock.py --all          # writes data/mock/<scenario>/ (~2 s)
/usr/bin/python3 mock.py A_rv1 --seed 7 --chi 0.9 --spin-dir par --out /tmp/x
```

## Conventions (checked in validation.md)

- Elements follow Gillessen+09, which is the convention of GRAVITY 2024 Table B.1 and of the S301
  paper. The observer frame is x = +Dec (north), y = +RA (east), z = away from us. Omega is counted
  from north through east, the ascending node is where the star recedes through the sky plane, and
  i > 90 deg means clockwise on the sky. relatipy's observer axes `(Dec, RA, away)` are the same
  frame, so published elements go in unchanged. Its outputs are alpha = +RA and delta = +Dec in
  arcsec, and v > 0 means receding (redshift-equivalent, including the gravitational redshift and the
  transverse Doppler shift). Foschi+ (arXiv:2607.26713) use a z-flipped frame (Omega + 180, omega +
  180); do not mix it with this one.
- Osculation epoch: `KerrForward` converts the Kepler elements to the initial state at
  `element_epoch`. It takes the true anomaly there from Kepler's equation with GM = M and the
  published t_p.
  - S2, S29, S38 and S55 use the GRAVITY 2024 epochs 2010.35, 1977, 2000 and 2012.
  - The S301 epoch is not published. We use the preceding apocentre, t_p − P/2 = 2018.785. The
    digitized points accept any epoch before the 2023.13 pericentre and reject any epoch after it
    (RMS 3 mas), because the osculating omega jumps by the ~2 deg Schwarzschild advance at
    pericentre.
  - Setting the state at t_p is badly conditioned for e = 0.98 (30 mas errors for S301, 1 to 4 mas
    for S2).
- Spin: `chi` is the magnitude and `spin_vec` the direction in `(Dec, RA, away)`. Named directions
  (`par`, `perpN`, `perpM`, `b45l45`, `b120l200`, `b70l310`) are defined relative to the S301A
  orbital frame (node N, L × N, L), as in the prototype. The orbital pole is
  L = (sin Ω sin i, −cos Ω sin i, cos i).
- Light: straight-line propagation with the Rømer delay. There is no Shapiro delay and no lensing,
  as in relatipy. Epochs are observed times in decimal Julian years.

## Stars and globals (`stars.json`)

- M = 4.2996 ± 0.0118 × 10⁶ Msun and R0 = 8.2759 ± 0.0086 kpc (GRAVITY 2024 Table B.1). The S301
  paper fixed 4.297 and 8.277, which are stored as `global.alternative`.
- S301A and S301B: Extended Data Table 2 of arXiv:2607.12664.
- S2, S29, S38 and S55: GRAVITY 2024 Table B.1. The t_p error of S29 is lost in the text layer and
  stored as null.
- B is essentially the z-mirror of A: Ω_B − Ω_A = 183.1 deg and ω_B − ω_A = −178.3 deg, with i, e, a
  and t_p equal within the errors. On the sky the two tracks differ by at most 0.57 mas over
  2017-2041 (median 0.25 mas away from pericentre), while their RVs have opposite sign. Only RVs
  break the degeneracy.

## Schedules (`schedules.py`, `ScheduleConfig`)

Visibility window: 1 March to 30 September by default. `SEASON_LONG` is 1 February to 6 November,
from the blackout quoted in arXiv:2609.06936. Every synthetic and future epoch lies inside the
window. Grids are deterministic: n epochs are spread evenly over the observable time.

**Past (to 2026.0)**

| star | astrometry | RV |
|---|---|---|
| S301 | the 19 Extended Data Table 1 epochs at 207 µas (`real_epoch`); the 2023.577 point is optional and flagged | none |
| S2 | the CSV, 128 NACO/SHARP and 55 GRAVITY epochs, with real errors (`real_epoch`); plus 27 GRAVITY 2019.75-2022.7 to reach GRAVITY 2024's 82, and 6/yr 2022.75-2026, at 50 µas | 92 SINFONI (62 over 2003.3-2019.7 at 15 km/s, 30 in the 2018 season at 10 km/s), 3 Keck (30), 2 NACO (50), 4 GNIRS (20), 3 ERIS 2022 (10), then ERIS 3/yr (10 km/s) |
| S29 | 66 NACO 2002.3-2019.7 (0.5 mas), 29 GRAVITY 2019.6-2022.7 (50 µas), then 6/yr | 17 SINFONI (20 km/s), 2 GNIRS 2021 (20), then ERIS 3/yr |
| S38 | 110 NACO 2004.2-2018.6, 23 GRAVITY 2021.2-2022.7, then 6/yr | 8 SINFONI (25), 2 Keck (30), 1 ERIS 2022 (10), then ERIS 3/yr |
| S55 | 42 NACO 2004.5-2013.6, 27 GRAVITY 2021.2-2022.6, then 6/yr | 2 SINFONI 2014 (30), then ERIS 3/yr |

The friends' counts and spans are those of GRAVITY 2024 (p.2). The epochs themselves are plausible
(`past_synth`), and so are the RV and NACO errors, which the paper does not give.

**Future (2026.0 to t_end = 2041.0)**

- GRAVITY+ astrometry, 10 epochs/yr. S301 at 100 µas, S2 at 30, S29 and S38 at 40, S55 at 50 µas.
- Pericentre campaigns add extra epochs on a quarter-offset grid. S301 gets 10 astrometric and 10
  RV epochs within ±0.06 yr (±22 d) of 2031.81 and 2040.49. The friends get 10 astrometric and 5 RV
  epochs within ±0.15 yr; this applies to S2 at 2034.43 and S55 at 2033.91. S29 and S38 have no
  pericentre before 2041.
- An out-of-season pericentre has its epochs moved into the last (or first) half-window of the
  nearest season edge (`peri_edge`).
  - S301 at 2031.81, which is 22 October, gets 10 epochs on 9-28 September 2031.
  - S55 at 2033.91 gets its epochs in September 2033.
  - With `SEASON_LONG`, the 2031.81 pericentre itself is observed.
- S301 RV: none, 20, 5 or 1 km/s (`s301_sigma_rv`), from 2030.0 (`s301_rv_start`), 4/yr plus the
  pericentre campaigns. The instrument label is `ELT/HARMONI` for ≥ 15 km/s and `ELT/MICADO`
  otherwise.
- Friends' RV: ERIS-like, 10 km/s, 3/yr plus pericentre campaigns.

## Mocks (`mock.py`)

- `make_dataset(Scenario, forward=...)` and `write_dataset(...)`. A forward model is any callable
  `model(elements, t_ast, t_rv, bh) -> (ra_mas, dec_mas, v_kms)`. The default is
  `forward.KerrForward` (exact Kerr, homogeneous, dop853 with rtol 1e-11). A Kerr + disk model plugs
  in through the same signature, from Python or with `--model module:factory`.
- Truth: the published central elements (S301A or S301B), M and R0 from GRAVITY 2024, and
  χ = 0.7 along `b45l45` by default (the same observer-frame vector in every scenario).
- Noise: independent Gaussian noise with the scheduled σ. The generator is
  `numpy.random.default_rng([seed, star index])` with seed 20261025 by default. Each star has its own
  stream, so the past epochs and the friends get identical noise across scenarios. A noise-free copy
  goes to `nonoise/`.
- The past epochs carry mock values at the real epochs. The real measurements live only in `data/*.csv`.
- No reference-frame offsets or drifts, no confusion, no correlated noise.
- Files in `data/mock/<scenario>/`:
  - `<star>_astrometry.csv`: t, ra, dec, sigma_ra, sigma_dec, instrument, origin (yr, mas)
  - `<star>_rv.csv`: t, v, sigma, instrument, origin (yr, km/s)
  - `nonoise/`: the same files without noise
  - `scenario.json`: config, truth, counts, model, relatipy version and git commit

| scenario | content | S301 ast / RV epochs |
|---|---|---|
| `past` | data to 2026.0 only | 19 / 0 |
| `A_rvnone`, `A_rv20`, `A_rv5`, `A_rv1` | truth A, 2026-2041, S301 RV none/20/5/1 km/s from 2030 | 189 / 0 or 64 |
| `A_rv1_2035` | as A_rv1, data end 2035.0 | 119 / 30 |
| `A_rv1_longseason` | as A_rv1, 1 Feb - 6 Nov window | 189 / 64 |
| `A_rv1_s301at30` | as A_rv1, S301 future astrometry at 30 µas | 189 / 64 |
| `B_rv1`, `B_rvnone` | truth B | 189 / 64 or 0 |

The friends have 388/163 (S2), 263/73 (S29), 301/65 (S38) and 247/61 (S55) astrometric/RV epochs
to 2041.

Size of the spin signal at fixed elements (χ = 0.7 minus χ = 0, noise-free, A_rv1): S301 up to
370 µas and 42 km/s; the friends 2 to 8 µas and below 0.15 km/s.

## Validation summary (details in `validation.md`)

Residual RMS per coordinate, data minus model, with χ = 0 and the published elements:

| data | RMS RA | RMS Dec | chi2/N |
|---|---|---|---|
| S2 GRAVITY 2016-2019 (55) | 0.128 mas | 0.082 mas | 1.37 |
| S2 NACO/SHARP 1992-2019 (128) | 1.02 mas | 1.61 mas | 0.91 |
| S2 NACO/SHARP from 2002.5 (118) | 0.53 mas | 0.40 mas | |
| S301A digitized (19 matched) | 0.23 mas | 0.47 mas | 3.2 (at 207 µas) |
| S301B digitized (19 matched) | 0.21 mas | 0.46 mas | 3.0 |

## Open issues

- The S301 osculation epoch is an assumption (anything before the 2023.13 pericentre fits).
- The S29 t_p error is missing from stars.json.
- The S301 digitized file still needs manual review: the Dec scatter is twice the RA scatter, and the
  2023.577 point is unexplained.
- relatipy 1.2.0 occasionally fails observed-time matching at isolated epochs (step-sequence
  dependent). One reproducer: S2, χ = 0.7, vec (0, 0, 1), t_obs = 2010.9470999999608, rtol 1e-11.
  `KerrForward` retries with other tolerances, which agree to < 0.001 µas, and counts the retries in
  `n_retries`. The library was not modified.
- The published elements come from 1PN fits (× f_SP for GRAVITY 2024); we use them as exact-Kerr
  truth. The mocks are self-consistent, but the real-data residuals include that model mismatch.
- No NACO frame parameters in the mocks. If the past NACO data are fitted, add x0, y0, vx0, vy0.
- The synthetic past epochs and the RV and NACO errors of the friends are plausible, not real.
