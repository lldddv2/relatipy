# Kerr+disk prototype: final report (Q0–Q3)

Forward model: PN (1PN + 1.5PN spin-orbit; 2PN omitted in truth and model) plus a
razor-thin disk (Σ ∝ R^-1, Piran+ fiducial) from a tabulated ring sum. S301 solution A
(B as check), truth χ = 0.7, 6 spin orientations. Spin errors along S301's apse line (P),
in-plane minor axis (Q) and orbital angular momentum (L, parallel component).

Inputs: S301 19 past epochs at 207 µas, future 100 µas, RV from 2030 (none/20/5/1 km/s);
M, R0 priors from GRAVITY 2024; friends with published element errors as priors and future
30 µas + 10 km/s. Fit: M, R0, spin vector, 6 elements per star, disk mass + normal.

## Validation
- PN nodal precession of S301 is 0.87% above exact Kerr (relatipy): 6.908e-4 vs
  6.848e-4 rad/orbit; scales as 1/p (missing 2.5PN spin-orbit, omitted on both sides).
- Ring formula vs brute force: 1e-13. Table interpolation error ≤ 8e-6 (5e-7 near plane/edge).
- Piran+ Table 1 disk/LT ratios reproduced to 1–5% (S301 i=30°, ω=90°: 1.44 vs 1.4; S2 14.0 vs 14).
- Fisher step-size changes (×0.3, ×3): results stable to 0.3%.

## Q0 calibration (Nature forecast)
χ = 1 aligned, 100 µas + 1 km/s, 2026–2035: σ(χ) = 0.10, orientation 27°/61°.
GRAVITY: σ < 0.2, ±30°. Agreement within a factor 2 (nuisance: M, R0, 6 elements).

## Q1: Fisher σ (P / Q / L), median over spin orientations

| case | 2035, RV 1 | 2041, no RV | 2041, RV 5 | 2041, RV 1 |
|---|---|---|---|---|
| a: S301, no disk | 1.58 / 0.21 / 0.13 | 1.82 / 0.20 / 0.14 | 1.17 / 0.18 / 0.08 | 0.40 / 0.10 / 0.07 |
| b: S301, disk 1e3 (CW) | 2.54 / 0.51 / 0.40 | 11 / 1.4 / 2.6 | 4.3 / 0.54 / 0.78 | 1.60 / 0.22 / 0.28 |
| c: joint, disk 1e3 (CW) | 2.01 / 0.34 / 0.16 | 1.96 / 0.23 / 0.10 | 1.47 / 0.21 / 0.09 | 1.04 / 0.14 / 0.07 |
| joint, no disk | 1.55 / 0.20 / 0.12 | 0.91 / 0.19 / 0.10 | 0.78 / 0.17 / 0.08 | 0.37 / 0.10 / 0.07 |

- σ(M_disk) case c: 234 Msun (2041), ~840 Msun (2035); case b unconstrained (> 3000).
- M_disk = 300 or Piran geometry: spin σ barely change. RV 20 km/s ≈ no RV.
- 30 µas (joint, 2041, RV 1): perpendicular RMS 0.55, parallel 0.03. Solution B: 1.53/0.15/0.07.

## Q2: full orbit vs secular
- LT: 99.5–100% of the per-orbit plane change within ±5% of P around pericentre;
  disk: 12–15% (63% within ±25%).
- No disk: secular ≈ full. S301 alone with disk: secular rank-deficient (2–3 of 6), full finite.
- Joint: full ~4× tighter on perpendicular spin than secular node+inclination
  (0.50/0.92 vs 1.97/3.91); with apsidal rate 1.7–3.5×. Disk mass ~4× better (234 vs ~900).

## Q3: bias when the disk is ignored (noise-free, 2041, RV 1)

| truth disk | fit | P | Q | L | Δχ² |
|---|---|---|---|---|---|
| 1e3 (CW / 2 random) | S301 alone | −4.0 / +5.8 / +8.2σ | ≤ 3.6σ | ≤ 1.2σ | ≤ 5 |
| 1e3 | joint | −2.2 / +5.9 / +9.0σ | ≤ 3.5σ | ≤ 0.4σ | 49 / 6.6 / 12 |
| 300 | either | ≤ 2.6σ | ≤ 1σ | ≤ 1σ | — |

Data to 2035 only: all ≤ 0.8σ. Joint without S301 RV (CW): up to 4.4σ.

## Gaussianity checks
- S301 alone, no disk: profile likelihood Gaussian. With disk: curved along disk
  orientation (Δχ² at −2σ: 4.0–5.5).
- emcee (case c, RV 1, 38 parameters, 160 walkers, 3000 steps, τ ≈ 190, not fully
  converged): P 0.54, Q 0.10, L 0.064, M_disk 142 Msun, disk angles 6°/10°
  (Fisher: 1.03, 0.14, 0.068, 234 Msun, 25°/25°). Fisher is conservative with a disk.

## Caveats
2PN omitted (must be modelled; the parallel component comes from apsidal precession and
needs f_SP fixed). No light bending/Shapiro; RV only v²/2 + M/r. Disk profile and edge
fixed; no granularity. RV epochs ≥ 10 days from pericentre. No seasonal visibility
(the 2031.81 pericentre is in October, Galactic Centre not observable).

## Verdict
Joint full-orbit fit with RV beats the secular strategy (~4× on perpendicular spin, breaks
the single-star degeneracy, disk mass ~4× better). For χ = 0.7 and ≤ 5 km/s RV, two
pericentres (to 2041): parallel ±0.07–0.09 (~10σ), Q ±0.14–0.2; apse-line component stays
σ ≈ 0.5–1 unless 30 µas. Ignoring a 1e3 Msun disk can bias the perpendicular spin by up to
9σ with no detectable misfit for S301 alone; the friends reveal it (Δχ² up to 49).
