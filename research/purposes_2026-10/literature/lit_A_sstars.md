# Literature survey A: relativistic orbit fitting of S-stars around Sgr A*
Status key: [snip] = taken from WebSearch result summaries (full text NOT read); [unverified] = not confirmed in this session.
Access: arxiv.org, export.arxiv.org, doi.org, ui.adsabs.harvard.edu, www.aanda.org, www.eso.org, pubmed.ncbi.nlm.nih.gov return HTTP 403 at proxy CONNECT (session egress policy). WebFetch returns getaddrinfo ENOTFOUND for every host tried. No primary full text was read; all numbers below are search-snippet level and need checking against the papers.

## Q1. GRAVITY Collaboration S2 papers (2018, 2020, 2022, 2024, later)

Equations of motion (as far as snippets show): Keplerian orbit plus post-Newtonian-type relativistic terms (transverse Doppler "SRT", gravitational redshift "RS", Romer delay) with a single scaling factor f (2018) or f_SP for Schwarzschild precession (2020); f=0 Newtonian, f=1 GR [snip]. No snippet mentions a Shapiro delay term, Kerr/spin, or a full geodesic integration [unverified for all 2018-2024 papers].

1. GRAVITY Collaboration (Abuter et al.) 2018, A&A 615, L15. "Detection of the gravitational redshift in the orbit of the star S2 near the Galactic centre massive black hole". arXiv:1807.09409 (v1, submitted 25 Jul 2018).
- Joint detection of gravitational redshift + transverse Doppler, z ~ 200 km/s/c (data up to and including pericentre) [snip].
- f = 0.90 +/- 0.09 (stat) +/- 0.15 (sys) [snip].
- Pericentre passage May 2018 at ~2800 r_g, with r_g = GM/c^2 = 6.1e11 cm ~ 5 uas [snip].
- Free-parameter count, sampler, RV precision, per-epoch astrometric precision: [unverified].
- Extended-mass limit: none in this paper [unverified].

2. GRAVITY Collaboration (Abuter et al.) 2020, A&A 636, L5. "Detection of the Schwarzschild precession in the orbit of the star S2 near the Galactic centre massive black hole". arXiv:2004.07187 (submitted 15 Apr 2020).
- 14 parameters fitted jointly: six orbital elements, distance, f_SP, and others (full list unverified) [snip].
- f_SP = 1.1 with 1-sigma error 0.19 (fit including flare data) [snip]. Your citation value 1.10 +/- 0.19 matches this snippet.
- Significance ~5-6 sigma, stable under posterior fitting vs MCMC [snip].
- Sampler: "MCMC" named, emcee not named [snip].
- Astrometry: 54 GRAVITY interferometric separations 2016.7-2019.7, rms ~65 uas [snip]. NACO: 75 points 2003.3-2019.7, rms ~1.7 mas [snip]. A 2019 review (arXiv:1901.08343) quotes ~30 uas precision for GRAVITY 2018 [snip]. These two numbers differ (rms of residuals vs quoted precision).
- RV: SINFONI spectra each assigned ~12 km/s; 92 SINFONI spectra 2003.3-2019.45 [snip]. RV precision for NIRC2/NACO spectra not given [unverified].
- Extended-mass limit: none [unverified].

3. GRAVITY Collaboration 2022, A&A 657, L12, DOI 10.1051/0004-6361/202142465. "Mass distribution in the Galactic Center based on interferometric astrometry of multiple stellar orbits". arXiv:2112.07478 (submitted 14 Dec 2021).
- Four stars for the GR fit (S2, S29, S38, S55) plus earlier AO data, data to ~2021.6; 13 outer stars for the mass-distribution analysis [snip].
- Schwarzschild precession ~7 sigma [snip, via 2024 follow-up].
- M = 4.30e6 M_sun to ~0.25% [snip]. Review arXiv:2211.07008 quotes 4.297 +/- 0.012 x 10^6 M_sun [snip].
- R0 = 8.277 +/- 0.009 kpc [snip, review arXiv:2211.07008]. MPE page: 8275 pc, +/-9 stat, +/-33 sys, attributed to a 2021 GRAVITY paper [snip]. Fig. 4 annotations show R0 fixed at 8277 pc and 8279.5 pc, fit not identified [snip].
- Extended mass inside S2 apocentre (~0.23 arcsec): <~3000 M_sun (1 sigma), ~0.1% of M [snip]. Tentative excess mass beyond ~2.5 arcsec with the 13 outer stars [snip].
- Astrometric precision: ~50 uas for stars within 0.1 arcsec of Sgr A* [snip].
- Free-parameter count, sampler, RV precision, density profile used in fit: [unverified].

4. Abd El Dayem, K. et al. (GRAVITY Collaboration) 2024, A&A 692, A242, published online 17 Dec 2024. "Improving constraints on the extended mass distribution in the Galactic Center with stellar orbits". arXiv:2409.12261.
- Data to September 2022 [snip].
- f_SP, S2 alone: 0.918 +/- 0.128 with NACO zero-point prior; MCMC run gives 0.911 +/- 0.131 [snip]. A separate snippet quotes a combined-star f_SP of 1.135 +/- 0.110 [snip, attribution unclear]. These two snippets are not reconciled.
- Schwarzschild precession ~10 sigma [snip]. One snippet says the precession was measured from GRAVITY astrometry + RVs alone, with NACO dropped [snip], which conflicts with the NACO zero-point prior above.
- Extended mass: power-law and Plummer profiles both consistent with zero M(<S2 orbit); 1-sigma upper limit ~1200 M_sun [snip].
- Sampler: MCMC, 200,000 realizations (Fig. B.1) [snip].
- Free-parameter count, astrometric precision, RV precision: [unverified].

5. Post-2024 items (authors/titles only partly retrieved):
- arXiv:2507.01510 (July 2025 preprint, authors not retrieved): cluster-object perturbations shift in-plane precession by up to ~1.5 arcmin (~13% of ~12 arcmin per orbit) and make Sgr A* wander by up to ~6 uas; simulations put 2026 apocentre Dec residuals above 30 uas in 35-60% of runs [snip; numbers unverified].
- A&A 2025, article aa54223-25 (Sep 2025 issue; title/authors not retrieved): precession ~12 arcmin per orbit, redshift and prograde precession treated as established [snip].
- Nature 2026, "Discovery of a star sensitive to the spin of Sagittarius A*", K. Abd El Dayem et al. with GRAVITY+ Collaboration; online 19 Aug 2026; DOI 10.1038/s41586-026-10894-w [snip]; preprint arXiv:2607.12664 [snip]. One press-coverage snippet names Felix Mang as lead author; the Nature-page snippet names K. Abd El Dayem. Unresolved. Details in Q2.
- Also cited by snippet: GRAVITY 2019 (A&A 625, L10 per my memory [unverified]; arXiv:1904.05721, "A geometric distance measurement to the Galactic Center black hole with 0.3% uncertainty") with f_RS = 1.04 +/- 0.05 [snip].
- Do et al. 2019 (Science 365, 664; arXiv:1907.10731), Keck data: Upsilon = 0.88 +/- 0.17, Newtonian excluded at 5 sigma [snip]. See Q4.

Q1 bottom line: f_SP (2020) = 1.10 +/- 0.19; f_SP (2024, S2 alone) = 0.918 +/- 0.128 (NACO prior) or 0.911 +/- 0.131 (MCMC); combined-star value 1.135 +/- 0.110 (unreconciled). Extended-mass 1-sigma limit 3000 M_sun (2022) -> 1200 M_sun (2024). Equations, sampler and RV precision per paper are NOT verified.


## Q7 (coordinator extra question). Do S300 and S301 exist in the literature?
Scale used for derived numbers: 1 mas = 8.277 AU (R0 = 8.277 kpc); r_g = GM/c^2 = 0.04241 AU for M = 4.297e6 M_sun; R_S = 2GM/c^2 = 0.08483 AU. Derived values are labelled [derived] and are my arithmetic from the snippet orbital elements, not paper values.

Citation check: the requested "GRAVITY Collaboration 2021, A&A 645, A127, Deep images of the Galactic center with GRAVITY" does NOT match what the searches returned. The deep-images paper appears as arXiv:2112.07477, published A&A 657, A82 (2022), article number aa42459-21 [snip]. A&A 645 A127 is therefore unverified and probably a wrong reference; please check.

### S300: exists, first reported in the GRAVITY deep-images paper (arXiv:2112.07477, A&A 657, A82, 2022 [snip])
- Discovery: faint, fast-moving object in 2021 GRAVITY data (observations March-July 2021), found with the G^R imaging method; no known S-star or Gaia match, so named S300 [snip].
- Proper motion ~66 mas/yr westward (the paper's arXiv version); 2-sigma acceleration bound -23 mas/yr^2 [snip]. Transverse speed ~2600 km/s at 8.277 kpc [derived from 66 mas/yr].
- K-band magnitude: ~19.0-19.3 by flux comparison with S29; bracketing ranges 19.3-19.7 (March) to 20.1-20.5 (July) [snip].
- Orbit (one retrieval only): central solution a ~171 mas, e ~0.6, P ~26 yr [snip]. A second retrieval of the same paper did not show these values and cut off before the orbit section, so they are unverified.
- Pericentre [derived from a = 171 mas, e = 0.6]: a(1-e) = 68 mas = 566 AU = ~13,350 r_g = ~6,670 R_S. Kepler check P = 25.7 yr, consistent with the 26 yr.
- Bound-orbit constraint on |z| = 50-137 mas; foreground Galactic-disk object not excluded [snip].
- Whether GRAVITY confirmed or refuted S62 / S4711-S4716 in this paper: S62's expected inward path is corroborated [snip]. S4711-S4716 are not mentioned in the retrieved text [unverified]. S38, S42, S60, S63 and S29, S55 get refreshed orbits [snip].

### S301: exists, reported in 2026 (Nature, GRAVITY+ Collaboration with K. Abd El Dayem et al.; preprint arXiv:2607.12664 [snip]); not in the 2021 deep-images paper [snip]
- Magnitude: K ~ 19.3 [snip]; spectral type ~F1.5 and mass ~1.1-1.5 M_sun are snippet-level inferences [unverified].
- Period: 8.68-8.7 yr [snip]; Kepler check with a = 687 AU gives 8.69 yr [derived].
- Eccentricity: 0.982-0.983 [snip]. Semimajor axis ~83 mas [snip] = 687 AU [derived].
- Pericentre: 136-142 R_S from Sgr A* [snip] = ~11.5-12.1 AU [derived with R_S = 0.0848 AU] = ~270-285 r_g [derived]. Using a(1-e) directly gives 11.7-12.4 AU, consistent.
- Peak speed ~25,000 km/s (arXiv follow-up: 25,600 km/s) [snip]. S2 at pericentre ~7650 km/s for comparison [snip].
- Predicted GR precession ~1.9-2.0 deg per orbit [snip]. Radial velocity not yet measured, so two candidate orientations; arXiv follow-up adopts one [snip].
- Related preprints (titles and IDs from search results; contents not read): arXiv:2607.24931 "S301 and friends: Measuring the spin of Sgr A*"; arXiv:2609.06936 "Catching Up with the Fastest Star in the Galaxy: A Two-Passage GRAVITY+ Campaign to Detect S301's Schwarzschild Precession"; arXiv:2609.06977 "A New Window for Testing the PPN in the Strong-Field Regime: Stellar Rotation of S301".
- Attribution conflict: one snippet names Felix Mang as lead in press coverage; the Nature snippet names K. Abd El Dayem et al.; a sci.news snippet attributes the study to Abd El Dayem's group at Observatoire de Paris [snip, unresolved].

### Other names in the same snippets (for Q3 cross-check)
- Earlier "fastest star" claim named S4714 at 23,900 +/- 8,800 km/s [snip]. An S62 record-distance claim was later walked back after checks put S62 much farther from Sgr A* [snip; unverified].
- arXiv:2207.02142, "Observation of S4716 - a Star with a 4 yr Orbit around Sgr A*", ApJ, DOI 10.3847/1538-4357/ac752f [snip] (probably Peissker et al. 2022, see Q3).

Q7 bottom line: S300 and S301 both exist in the recent literature. S300: P ~ 26 yr, e ~ 0.6, r_p ~ 566 AU (~13,350 r_g), K ~ 19.0-19.3 (orbit values not confirmed by a second retrieval). S301: P ~ 8.7 yr, e ~ 0.982-0.983, r_p ~ 11.5-12.4 AU (136-142 R_S), K ~ 19.3. Neither is in GRAVITY 2022 Table for S2/S29/S38/S55.


## Q2. Spin-from-stars forecasts (historical + 2020-2026)
Conventions: r_g = GM/c^2 (= 0.04241 AU for M = 4.297e6 M_sun); R_S = 2 r_g. Derived numbers marked [derived].

Citation status:
- Will 2008, ApJL 674, L25-L28, arXiv:0711.1677 (submitted 11 Nov 2007). Title "Testing the general relativistic 'no-hair' theorems using the Galactic center black hole Sgr A*" [snip, verified by search]. Author: C. M. Will (Washington Univ. in St. Louis). INSPIRE drops the minus sign in the quadrupole relation; arXiv has it [snip].
- Merritt, Alexander, Mikkola & Will 2010, PRD 81, 062002 (DOI 10.1103/PhysRevD.81.062002, 17 pp, March 2010), title "Testing properties of the Galactic center black hole using stellar orbits" [snip, consistent across two records]. arXiv ID not retrieved. One reference list gives PRD 75, 044007 (2007) for a related citation; not reconciled.
- Zhang, Lu & Yuan 2015, ApJ 809, 127: NOT matched in any search. Closest title match: "Possibility of measuring spin precession of the nearest supermassive black hole by using S stars", arXiv:1404.5160, Research in Astronomy and Astrophysics 14(11), 1415-1422 (2014) per a library record [snip]. Author list and the ApJ citation are unverified; please check the citation.
- Grould, Vincent, Paumard & Perrin 2017, "General relativistic effects on the orbit of the S2 star with GRAVITY", A&A, DOI 10.1051/0004-6361/201731148 (article aa31148-17), arXiv:1709.04492 [snip]. "A&A 608, A60" per your citation, consistent with the DOI but not directly confirmed.
- Waisberg et al. 2018, arXiv:1802.08198, title "What stellar orbit is needed to measure the spin of the Galactic center black hole from astrometric data?" [snip]. Authors: I. Waisberg with R. Genzel listed in one snippet, not confirmed. Journal MNRAS 476, 3600 not confirmed in the retrieved text.
- Fragione & Loeb 2020, "An upper limit on the spin of Sgr A* based on stellar orbits in its vicinity", ApJL 901, L32, DOI 10.3847/2041-8213/abb9b4 (online 1 Oct 2020); arXiv:2008.11734 (v1, 26 Aug 2020) [snip].

Predicted LT / spin signal, required orbit, degeneracies, method:

1. Will 2008 (arXiv:0711.1677). Hypothetical stars on tight, highly elongated orbits; orbital-plane precession up to ~10 uas/yr from frame dragging + quadrupole; largest values need spin > 0.5 of maximum [snip]. Two or more stars tracked with precise positions could give spin vector and quadrupole and test Q = -a^2 M [snip]. Orbit: unspecified. Degeneracy: spin and quadrupole jointly. Method: not retrieved [unverified]; probably analytic secular rates.

2. Merritt, Alexander, Mikkola & Will 2010 (PRD 81, 062002). Astrometric monitoring at milli-pc separations could reveal spin and quadrupole; several candidate star and stellar-remnant layouts evaluated; open question whether relativistic effects can be separated from stellar perturbations [snip]. Numbers, orbits and method: [unverified].

3. Zhang et al. (2015 citation / RAA 2014 title match). Assumed astrometry 1-30 uas and redshift precision 0.1-10 km/s in the forecasts that Waisberg et al. (2018) cite [snip]. Other numbers [unverified].

4. Grould et al. 2017 (arXiv:1709.04492). Tests whether GRAVITY can detect several relativistic signals in S2 astrometry and RVs, LT included. Figures from snippets (inconsistent, all need checking against the A&A text):
- LT shift ~12 uas per orbit, detectable after one or two revolutions (a later review, not identified) [snip].
- Conference slides: LT astrometric effect ~40 uas and ~1 km/s spectroscopic, for a = 0.99 over three orbits (47 yr) [snip].
- Proceedings: astrometric ~15 uas at S2 passage; spectroscopic < 3 km/s [snip].
- Preprint judges constraint possible with GRAVITY + SINFONI; tight limits need long monitoring or a closer star [snip].
- Orbit used: S2 (pericentre ~2800 r_g = ~120 AU [derived from the 2018 letter snippet]). Degeneracies and method: [unverified].

5. Waisberg et al. 2018 (arXiv:1802.08198). Reference case: 4 yr baseline, 120 epochs, 10 uas precision, chi = 0.9: bound on a_orb(1-e^2)^(3/4) of ~300 R_S (formula garbled in snippet; check) [snip]. Bound relaxes with longer baseline, more epochs, better precision, higher spin [snip]. Stars with K < 19 in the central arcsec meeting these conditions: ~0.1 expected (30 epochs/yr over 4-10 yr, 10 uas) [snip]. Redshift precision needed is out of reach (~30 km/s current for S2) [snip]. Method: analytic bound (per snippet) [partly verified].

6. Fragione & Loeb 2020 (arXiv:2008.11734). Upper limit chi <~ 0.1 from the S-star layout in two nearly edge-on disks; a fast spin would have reoriented the orbital angular momenta through frame dragging [snip]. Not an astrometric forecast; analytic/dynamical argument.

7. Iorio 2020, arXiv:2009.01158, "The short-period S-stars S4711, S62, S4714 and the Lense-Thirring effect due to the spin of Sgr A*" [snip]. Argues earlier node-precession estimates require spin aligned with the line of sight, which it judges unlikely [snip]. Quoted LT rates for S4714 of ~7-14 arcsec/yr [snip; garbled, unverified].

8. arXiv:2207.02142 [authors per your brief Peissker et al. 2022; not confirmed by any snippet] ("Observation of S4716 - a Star with a 4 yr Orbit around Sgr A*", ApJ, DOI 10.3847/1538-4357/ac752f). LT precession rate ~1.11 +/- 0.13 arcsec/yr for spin parameter 0.5, a model estimate; LT scales as Gamma^2 while Schwarzschild precession and redshift scale as Gamma; Gamma(S4716) ~ 8.0e-4 [snip]. Per-orbit astrometric conversion: node shift ~2 uas per orbit [derived: ~4.4 arcsec per 4-yr orbit x apocentre ~95 mas; order-of-magnitude only].

9. Iorio (probable author, unconfirmed) 2023, arXiv:2306.17432, "Is it possible to measure the Lense-Thirring orbital shifts of the short-period S-star S4716 orbiting Sgr A*?" Title confirmed; verdict and numbers not retrieved [unverified]. Attribution to Iorio (2023) comes from a 2024 snippet (arXiv:2406.09916).

10. Abd El Dayem et al. 2025, arXiv:2505.24789, "The effects of the spin and quadrupole moment of Sgr A* on the orbits of S stars" [snip]. Post-Newtonian to 2PN order [snip]. Polar configurations give largest out-of-plane shift [snip]. A hypothetical "S2/10" star (PN parameter 10x S2) gives LT shifts ~10^3 times S2's over one S2 period [snip]. Degeneracy with quadrupole stated [snip].

11. Chinese Physics C 2025, DOI 10.1088/1674-1137/ad9304, "Inferring a spinning black hole in an expanding universe via the S2 star around the galactic center" (Chongqing University group per snippet; authors not confirmed) [snip]. Monte Carlo (MCMC) fits to public S2 data leave spin unconstrained; tighter-orbit simulations: LT precession axis follows spin direction; 1-sigma bound on Lambda ~7.3e-34 km^-2 [snip].

12. Discovery of S301 spin sensitivity, Nature 2026 (K. Abd El Dayem et al. with GRAVITY+ Collaboration; online 19 Aug 2026; DOI 10.1038/s41586-026-10894-w; preprint arXiv:2607.12664) [snip]. Star S301: P ~8.7 yr, e ~0.982-0.983, pericentre 136-142 R_S = 272-284 r_g = 11.5-12.1 AU [snip + derived]. Mock GRAVITY+ astrometry (assumed 100 uas) + ELT/MICADO RVs (assumed 1 km/s) for 2026-2035 at maximal spin [snip]. Spin uncertainty at 2035 not retrieved [unverified]. Detection within ~a decade [snip, press].
13. arXiv:2607.24931, "S301 and friends: Measuring the spin of Sgr A*" (A&A manuscript S301_LT, dated 29 Jul 2026 in the header) [snip]. Nodal LT shift ~9.6e-4 rad per orbit at maximal spin (~0.055 deg) [snip]. Derived: x apocentre 164.6 mas gives ~160 uas per orbit (order of magnitude). Disk-like mass can mimic the nodal signal; benchmarking against S2, S55, S38 proposed [snip]. Apsidal LT term 0.114 deg per orbit at chi = 1 quoted in a companion campaign paper [snip; attribution unclear].
14. arXiv:2609.06936, "Catching Up with the Fastest Star in the Galaxy: A Two-Passage GRAVITY+ Campaign to Detect S301's Schwarzschild Precession" [title snip; contents not read].
15. arXiv:2609.06977, "A New Window for Testing the PPN in the Strong-Field Regime: Stellar Rotation of S301" [title snip]. One summary says LT nodal precession of S301 ~30x S2 [snip; attribution unclear].
16. Earlier estimate for comparison: S2 LT per orbit ~12 uas [snip, see item 4].

Q2 summary for the table: predicted LT astrometric signal: S2 ~12-15 uas per orbit (Grould 2017 snippets); S301 ~0.16 mas per orbit at maximal spin (nodal, derived from 9.6e-4 rad, order of magnitude); S4716 ~2 uas per orbit (derived). Required orbit: Waisberg et al. 2018 reference bound a(1-e^2)^(3/4) of ~300 R_S for 4 yr, 120 epochs, 10 uas, chi = 0.9 (formula garbled in snippet, not a pericentre bound); S301 pericentre 136-142 R_S, about ten times closer than S2's ~1415 R_S [snip + derived]. Degeneracies: quadrupole, stellar disk/cluster perturbers, line-of-sight alignment. Method: most are analytic or perturbative; the 2026 forecasts use simulated data and fits (numerical; details unverified).


## Q3. Short-period stars S62, S4711-S4716 (Peissker et al. 2020, ApJ 899, 50; 2022, ApJ 933, 49). STATUS: INCOMPLETE
WebSearch budget for this turn was exhausted during Q3 (200 calls, shared by all agents). The searches for Peissker 2020 and 2022 and for disputes failed with "budget used up". No further searching was done, per the limit notice. Items below are what the earlier snippets support.

- Peissker et al. 2020 (ApJ 899, 50) and 2022 (ApJ 933, 49): not retrieved. Titles and authors are from your brief only [unverified].
- Only paper retrieved for S4716: arXiv:2207.02142, "Observation of S4716 - a Star with a 4 yr Orbit around Sgr A*", ApJ, DOI 10.3847/1538-4357/ac752f [snip; authors not in snippet]. Period ~4 yr [snip, from title].
- Orbital elements (period, a in arcsec and AU, e, inclination, pericentre in AU and r_g) for S62, S4711, S4712, S4713, S4714, S4715, S4716: UNVERIFIED for all. No table was retrieved.
- Relativistic parameter for S4716: Gamma ~ 8.0e-4 [snip]. Its definition (Gamma = r_g/(a(1-e^2))?) was not checked, so no pericentre is derived here.
- Fastest-star claim: S4714 at 23,900 +/- 8,800 km/s, an earlier press claim later superseded by S301 (~25,000 km/s) [snip; attribution unverified].

Disputed or doubtful claims (what was found):
1. Measurability of Lense-Thirring node precession for S4711, S62, S4714: Iorio 2020, arXiv:2009.01158, "The short-period S-stars S4711, S62, S4714 and the Lense-Thirring effect due to the spin of Sgr A*" [snip]. Iorio argues that earlier node-precession estimates need the spin axis to be aligned with the line of sight, which he judges unlikely [snip]. Quoted LT rates for S4714 of ~7-14 arcsec/yr [snip, garbled, unverified].
2. S62 distance: a record close-approach claim for S62 was later walked back after checks placed it much farther from Sgr A* than assumed [snip, unverified; checker not identified].
3. S4716 LT detectability: Iorio, arXiv:2306.17432, "Is it possible to measure the Lense-Thirring orbital shifts of the short-period S-star S4716 orbiting Sgr A*?" Title confirmed; the verdict was NOT retrieved [unverified].
4. Context (not a dispute): Burkert et al. 2023, arXiv:2306.02076, "The Orbital Structure and Selection Effects of the Galactic Center S-Star Cluster": a tight, roughly exponential pericentre-eccentricity relation across the cluster [snip].

Q3 bottom line: orbital elements of S62 and S4711-S4716 are NOT verified. Disputes found: Iorio (2020) on LT measurability (alignment assumption); an unverified S62 distance revision. Peissker 2020/2022 numbers need a follow-up search.


## Q4. Elements of S2, S29, S38, S55 (GRAVITY 2022 table) and Sgr A* M, R0, r_g. STATUS: PARTIAL
- GRAVITY 2022 table (a, e, i, Omega, omega, t_peri, P for S2, S29, S38, S55): NOT retrieved. All elements UNVERIFIED.
  Partial S2 values from snippets (not from the table): P ~ 16 yr [snip, arXiv:1807.09409 summary]; e ~ 0.88 [snip, arXiv:2004.07187 summary]; apocentre ~0.23 arcsec [snip, arXiv:2112.07478 summary]; pericentre passage May 2018 at ~120 AU, ~2800 r_g [snip, 2018 letter]. No values were retrieved for S29, S38, S55.
- Sgr A* mass (GRAVITY 2022): 4.297 +/- 0.012 x 10^6 M_sun [snip, from review arXiv:2211.07008]. Another snippet on arXiv:2112.07478 says "about 4.30 x 10^6 M_sun, roughly 0.25%" [snip]. Both are consistent with each other.
- R0 (GRAVITY 2022): 8.277 +/- 0.009 kpc [snip, review arXiv:2211.07008]. MPE page: 8275 pc, +/-9 (stat) +/-33 (sys), attributed to a 2021 GRAVITY paper [snip]. The 2022 value looks to omit the systematic term; this is an inference, unverified.
- Do et al. 2019, Science 365, 664 (arXiv:1907.10731; DOI 10.1126/science.aav8137): redshift parameter Upsilon = 0.88 +/- 0.17, consistent with GR and Newton excluded at 5 sigma [snip]. Mass and R0 from Do et al. are NOT retrieved [unverified].
- r_g = GM/c^2 for Sgr A* (derived, my arithmetic, G/c/M_sun constants standard): 0.04241 AU (M = 4.297e6 M_sun). Angular size at R0 = 8.277 kpc: 5.124 +/- 0.015 micro-arcsec (propagated M and R0 errors; 5.125 at 8.275 kpc). R_S = 2 r_g = 10.25 micro-arcsec. Consistent with the 2018 letter's "R_g ~ 5 uas" [snip].
- Derived: 120 AU (S2 pericentre) = 2829 r_g = 1415 R_S.

Q4 bottom line: r_g = 5.12 +/- 0.02 uas (derived). Mass 4.297 +/- 0.012 x 10^6 M_sun and R0 8.277 +/- 0.009 kpc (GRAVITY 2022, via review snippet). Table elements unverified.


## Q5. Forecast astrometric precision (GRAVITY+, ELT/MICADO) and spin-detection date. STATUS: PARTIAL
Precision assumptions and forecasts found (all snippet-level):
- GRAVITY+ astrometry: the discovery preprint arXiv:2607.12664 assumes ~100 uas GRAVITY+ astrometry in its 2026-2035 mock datasets [snip]. A 2018 GRAVITY-related paper (arXiv:1802.08198, Waisberg et al.) describes GRAVITY's target as ~10-100 uas [snip]. The official GRAVITY+ design figure was NOT retrieved [unverified].
- ELT/MICADO radial velocity: the same preprint assumes ~1 km/s accuracy [snip]. MICADO astrometric precision: NOT retrieved [unverified].
- Current measured precision: GRAVITY ~30 uas per 2018 (2019 review arXiv:1901.08343, [snip]); 54 GRAVITY separations, rms ~65 uas, 2016.7-2019.7 (arXiv:2004.07187, [snip]); ~50 uas for stars within 0.1 arcsec of Sgr A* (arXiv:2112.07478, [snip]).
- Threshold: GRAVITY's ~30 uas detection threshold for 2026 apocentre residuals; simulations give Dec residuals above 30 uas in 35-60% of runs (arXiv:2507.01510, [snip, numbers unverified]).

Spin-detection date forecasts found:
- S301 (Nature 2026, K. Abd El Dayem et al. with GRAVITY+; online 19 Aug 2026; DOI 10.1038/s41586-026-10894-w): two full revolutions needed; the next close approach is in 2031 [snip, ESO release eso2612, https://www.eso.org/public/news/eso2612/]. Press quotes a spin measurement "within roughly a decade" (attributed to Felix Mang) [snip].
- Mock data 2026-2035 (arXiv:2607.12664) [snip]. The 2035 spin uncertainty was NOT retrieved [unverified].
- arXiv:2607.24931 ("S301 and friends: Measuring the spin of Sgr A*") and arXiv:2609.06936 (two-passage GRAVITY+ campaign to detect S301's Schwarzschild precession): titles only [snip]. Forecast details NOT retrieved.
- Talk, not a paper: S2 alone needs ~47 yr to reach +/-0.1 (1 sigma) on the spin; 6- and 17-month detection times for hypothetical S2/5 and S2/10 scenarios (K. Abd El Dayem, 2023 talk slides, [snip]).

Q5 bottom line: GRAVITY+ 100 uas and ELT/MICADO 1 km/s are assumptions in the 2026 preprint, not verified official specs. Spin detection: "within about a decade" (press, from 2026); the 2031 pericentre is the first test; mock forecasts run to 2035.


## Q6. Verbatim statements that computational cost forced a PN / orbit-averaged model instead of full GR integration in an MCMC. STATUS: NOT FOUND
- No verbatim quote was retrieved. No snippet contained a statement that computational cost forced a post-Newtonian or orbit-averaged model in an MCMC. Nothing is quoted here, to avoid inventing text.
- Related but NOT qualifying: arXiv:2409.12261 (2024) MCMC with 200,000 realizations [snip, no cost statement]; arXiv:2505.24789 "carried out to second post-Newtonian order" [snip, no cost statement retrieved]; Waisberg et al. 2018 analytic bound [snip, no cost statement].
- Suggested next step after budget reset: search the full texts of arXiv:2004.07187, arXiv:2112.07478, arXiv:2409.12261, arXiv:2505.24789 and arXiv:2607.24931 for "computational", "cost", "post-Newtonian", "orbit-averaged", "MCMC", "full general relativistic integration". Access to arXiv is currently blocked at the egress proxy (403), so this needs a route the session allows.

Q6 bottom line: 0 verified verbatim quotes.

## Survey-wide notes
- Access: arxiv.org, export.arxiv.org, doi.org, ui.adsabs.harvard.edu, www.aanda.org, www.eso.org, pubmed.ncbi.nlm.nih.gov, iopscience.iop.org return 403 at proxy CONNECT (README: do not route around). WebFetch: getaddrinfo ENOTFOUND for every host tried. WebSearch snippets were the only source.
- WebSearch budget: 200 calls per turn shared by all agents; exhausted during Q3. Q3-Q6 are incomplete; a follow-up turn can continue.
- Citation discrepancies to check: (a) "GRAVITY 2021 A&A 645 A127" for deep images is not matched; the deep-images paper is arXiv:2112.07477, A&A 657, A82 (2022) [snip]. (b) "Zhang, Lu & Yuan 2015 ApJ 809 127" is not matched; the title match arXiv:1404.5160 is RAA 14, 1415 (2014) [snip]. (c) "Peissker et al. 2020/2022" not retrieved. (d) Sign of quadrupole relation in Will 2008 (INSPIRE vs arXiv) [snip]. (e) Merritt et al. 2010 vs PRD 75, 044007 (2007) reference [snip].
- Attribution conflicts unresolved: the 2024 f_SP values (0.918 +/- 0.128 S2-only vs 1.135 +/- 0.110 combined); Nature 2026 first author (K. Abd El Dayem vs Felix Mang); S300 orbit values (retrieved once, not in a second retrieval).

