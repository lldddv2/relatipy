# NOTES (survive compaction) — relatipy purposes hunt
## Budget: user said $3/$100 at ~04:00 UTC (Opus 5). Model now Opus 5.5. User asleep since ~04:10.
## Rules: no Sonnet. Haiku retrieval; Opus subagents allowed for judgment. Interrupt kills bg agents; nohup jobs survive.
## Lib facts: timelike Kerr only; no photons/Shapiro/lensing; no external forces; no GIL release; KerrMcmcModel gives RA/Dec(arcsec)+vlos(km/s) with Romer (t_obs inversion).
## RESULTS
R1 speed: S2 fwd model 200 ep/34yr: ~6.6ms (150/s/core); 10 ep 1.2ms (811/s); rtol 1e-4..1e-13 same cost; t_obs +30%. radau ~2x slower than dop853.
R2 PN spin-orbit sign validated vs EXACT Kerr (a=3000rg,e=.3,i=1,a*=.9,20 orbits):
   a_LT=(2/r^3)[3(n.S)(n x v)+(v x S)] -> dOmega/orbit 7.8593e-5 vs exact 7.8583e-5 (1.3e-4 rel). analytic LT 7.9288e-5 (0.9% off). other variants off x2, x0.5, x-0.5.
   apsidal: exact 6.7688e-3 vs 1.5PN 6.7287e-3 rad/orbit (0.6% short at a=3000rg)
R3 scan2 (inject a*=.9 axis 45deg, refit a*=0 all 8 params; 20yr; first peri t=2yr; 4 ep/yr + 30 ep per peri):
   S2-like a=1035 e=.88: ~0.1-0.4 uas rms, max 0.7-2 uas -> spin invisible to S2
   a=600 e=.976 (rp 340rg): 17 uas rms, 68 max; a=1035 e=.976 (rp 586): 3.7 rms 18.7 max
   a=400 e=.976 (rp 226): 33 rms 146 max ; a=250 e=.95 (rp 295): 21 rms 121 max
   -> signal set by rp & #orbits. See scan2.json
R4 refit.py (older, 40yr, sparse peri): S2 a*=.9: 1.3-2.4 uas rms, max 4.4-6.5
BUGS FIXED: (i) peri not sampled when IC at peri & baseline=P; (ii) PN fit p0 f=0 vs truth f0.
R5 PN bias (bias.py; truth exact Kerr, model numba PN, astrometry only, 20yr, all nuisance free):
   S2 a*=0: Newtonian 38.0 rms/218 max uas; 1PN 0.91/3.66; 1PN x f_SP free -> 0.20/0.40 with f_SP=1.0424 (+4.2% BIAS!)
   S4716-like a*=0: Newtonian 95/404; 1PN 1.22/3.1
   HYPOTHESIS: f_SP bias = BL vs harmonic gauge (r_BL=r_h+1 rg) since sky map = coordinate projection w/o photon propagation. Test: gauge.log
   S2 a*=0.9: 1PN(no spin) 1.2/4.5
R6 robustness (robust.py, 150 random perturbations each): t_obs 0 fails everywhere; t_eval dop853 1/150 (S4716-like), 3/150 (S29); radau 0. Fisher hit a t_obs failure at rtol 1e-13/1e-12 for S4716-like M-1e-5 -> sporadic.
R7 numba PN (pn_fast.py) validated: 1.6e-6 rg error vs scipy 1e-13 at 2000 steps/orbit; ~10-20ms/eval. Exact Kerr relatipy 6.6ms -> exact GR no costlier than PN.
R8 Fisher S2 alone 50uas/10kms a*=.9: s(a*)=6.9 (unconstrained)
R9 GAUGE TEST (bias_gauge.py): S2 1PN harmonic 0.91/3.66 uas -> with r->r+1 (BL) 0.010/0.037; f_SP harmonic 1.0424 -> BL 1.00045.
   S62-like: harmonic 0.40/2.8 -> BL 0.012/0.064; f_SP BL 1.0028. => 2PN dynamics negligible; observable mapping (gauge/light prop) dominates.
R10 spin bias (harmonic PN, astrometry only): S2 a* free -> 2.39 (true .9) garbage; S62-like 0.911; S4716-like 1.043.
R11 Fisher a*=.9 (axis 60/30), 20yr: s(a*) S2: 6.9/4.4/1.5/0.44 at 50/30/10/3 uas; S62-like: 0.17/0.10/0.035/0.011 (axis ~45/27/9/3 deg);
   S4716-like 1.7/1.6/0.73/0.26; S29 useless; all4: 0.126/0.078/0.027/0.008
R12 RV: lib uses 1+z=u^t+u_away (flat photon momentum). Omitted metric Doppler factor (static frame, no bending): S2 2.4, S4716 3.0, S62 42, rp300 70, rp100 373 km/s
R13 NETWORK: proxy blocks arxiv.org, export.arxiv.org, cdsarc/vizier, ui.adsabs, IOP, A&A, zenodo, doi.org -> no real data; WebSearch only.
   fix: env settings > Network access > Allowed domains (keep package managers ticked). docs https://code.claude.com/docs/en/cloud-environments#network-access
R14 C core already has rp_kerr_null_geodesic_rhs (used w/ arg 1.0 for timelike) -> photon extension is close.
USER: reminded 4 integrators incl projection_radau -> do invariant drift + phase benchmark.
R15 INTEGRATORS (integ_bench2.py, a*=.9 p=10 e=.5 x=.5, 1000 radial orbits, accepted steps, vs KerrGeoPy analytic; CPU contended):
   dop853 1e-12: 0.51s dE 3.5e-12 dLz 2e-9 dQ 9e-10 dphi 1.6e-5 | radau 1e-12: 5.9s dE 9e-10 dLz 2e-7 dQ 2e-7 dphi 5.8e-3
   projection_radau 1e-12: 6.0s dE/dLz/dQ ~5e-16..9e-16 dphi 1.9e-6 (best) | dp45 1e-12: 2.0s dphi 6.5e-4
   tau_eval-sampled invariants dominated by INTERPOLATION error (1e-5) -> don't use tau_eval to judge drift.
   projection_radau FAILS (status 10 'invariant projection failed') at rtol 1e-8,1e-9,1e-10, always at r-r_p=+6e-3..1.3e-2 (radial turning point). rtol 1e-11 OK.
   -> radau (DEFAULT) is dissipative: 10x slower & ~400x worse phase than dop853. Recommend default dop853; fix projection near turning points.
R16 MCMC S2+S62-like 10uas/5kms (17 par): step1000 a*=.943+-.020 cos_th=.38+-.09 ph=19.8+-5.2deg (truth .9/.5/30); tau~30-40.
R17 KerrGeoPy speed (analytic, python): S2-like 200 samples 18.6ms; +t-inversion 51ms/eval (19/s). relatipy full obs+Romer 6.6-8.6ms -> 6-8x faster.
R18 Fisher over 60 isotropic spin axes, S62-like alone 20yr: 30uas: median s(a*)=.27 (p16-84 .10-.48), only ~15% <0.1;
    10uas/5kms: median .09 (.03-.16), ~55% <0.1; axis err median 11deg(a*=.9) 16(a*=.6) 33(a*=.3). => GRAVITY+ era, not current.
IDEA P5: QPE timing (orbiter a~300rg around 4e5 Msun SMBH, P~18h, ~1000s of orbits over baseline) -> LT nodal 1.2e-3 rad/orbit, 2PN apsidal ~1e-3 rad/orbit -> PN fails; agent lit_H_qpe.md
R19 SOFTWARE (lit_B): no pip exact-Kerr timelike + astrometry/RV + Romer package exists. KerrGeoPy (arXiv:2406.01413, JOSS 10.21105/joss.06587) analytic, no obs, excludes polar/marginal/extremal.
   PN obs code public: GC-OrbitFit (github pmplewa, 1PN). orbitize! Keplerian. Grould+2017 (arXiv:1709.04492, GYOTO full GR S2) & Zhang+2015 (arXiv:1508.06293) no public code.
   PyGRO arXiv:2504.20152 (Della Monica). No tool does light-travel-time inversion. EinsteinPy uses FANTASY.
R20 MCMC1 (40 walkers) collapsed: widths 0.021/0.087/4.5deg < Fisher 0.035/0.137/9.75. MCMC2 (96 walkers, Fisher init): step1500 a*=.951+-.026 cos=.374+-.129 ph=19.6+-7.1 -> consistent w/ Fisher.
R21 numba PN prealloc: 0.9 ms/orbit; relatipy exact Kerr 0.25-0.39 ms/orbit at a=300rg -> exact faster than PN.
R22 LIT D: flares GRAVITY2018 circ Kepler+bending ~7rg i~160; Baubock+2020 arXiv:2002.08374 Kerr hotspot out-of-plane; Yfantis+2024 arXiv:2408.07120 ~12rg 1.3x Kepler (super-K), mm companion arXiv:2310.07762 a*>0.8.
   pulsars: Liu+2012 arXiv:1112.2151; Psaltis+2016 arXiv:1510.00394; Zhang&Saha2017 ApJ849,33 arXiv:1709.08341 numerical Kerr; Hu,Wang&Shao2026 ApJS arXiv:2602.19546 2PN+SO+Q; Ben-Salem&Hackmann2022 arXiv:2203.10931 PN unreliable tight orbits.
R23 LIT H QPE: Zhou+ I arXiv:2401.11190, II 2405.06429, III 2411.18046 (osculating, decay+disk precession); QPE-FIT Chakraborty+2025 arXiv:2508.20162 ApJ: 3PN Kerr geodesic exprs + GPU, ~10% M,a,e; spin NOT recovered a>=100rg 10-100 orbits; Arcodia+2026 arXiv:2604.09788 eRO-QPE2 32 eruptions ~325 cycles.
R24 LIT C: GRAVITY2024 A&A692,A242 arXiv:2409.12261 ext mass <1200Msun 1sig (S2,S29,S38,S55); scalar cloud MNRAS524,1075 arXiv:2306.17215; DellaMonica&deMartino2023 PRD108 L101303.
   Merritt+2010 arXiv:0911.4718: cluster swamps frame dragging beyond ~0.5 mpc -> CAVEAT for spin purpose; S301 exists (snippet, cluster precession paper).
R25 LIT A: S301 REAL: Nature 2026 DOI 10.1038/s41586-026-10894-w arXiv:2607.12664 (+2607.24931, 2609.06936, 2609.06977), ESO eso2612: K~19.3 P~8.7yr e~0.982-0.983 rp~272-284 rg (11.5-12.1AU), next peri 2031; spin forecast w/ GRAVITY+ 100uas(assumed) + MICADO RV 1 km/s to 2035; "spin within ~a decade".
   S300 (A&A 657 A82 2022, arXiv:2112.07477): a~171mas e~0.6 P~26yr (unconfirmed). GRAVITY2024 f_SP .911-.918+-.13 (S2), 1.135+-.110 combined; Abd El Dayem+2025 arXiv:2505.24789 spin+quadrupole to 2PN.
   Waisberg+2018 arXiv:1802.08198: need a(1-e^2)^(3/4) <~300 R_S (4yr,120ep,10uas,chi=.9). S301: 651 R_S. Grould+2017 lensing ~20uas peri, LT 10-40uas/3 orbits (snippets).
R26 LIT E: no paper fits exact-GR mock with PN and reports bias (not proof). Radau IIA for geodesics: none found. Projection = standard (Hairer-Lubich-Wanner). Wu+2021 ApJ914,63 explicit symplectic Kerr.
R27 QPE a=300 e=.2 a*=.7: 1.5PN(BL) fit to exact over 100 orbits: rms 0.088 s, max 0.43 s (M=4e5), a*_fit=0.652 (bias -0.05). 
R28 S301 Fisher (a*=.7, 60 iso axes, S2 angles placeholder, 4ep/yr+30/peri): 2026-35 (1 peri) 100uas/1kms med s(a*)=.60; 30/10 .82; 10/1 .21.
    2026-40 (2 peri): 100/1 .19 (.12-.28, axis 24deg); 30/10 .33; 30/1 .12; 10/5 .12; 10/1 .079 (axis 11deg). RV 1km/s carries most info.
R29 S301 peri RV raw errors vs static-frame ref (no bending): 1PN additive 256, 1PN+cross 100, SRxGR mult 84, relatipy 77 km/s; bending est ~80 km/s. MICADO 1 km/s -> need null geodesics.
REFEREE (Opus) launched with full brief.
