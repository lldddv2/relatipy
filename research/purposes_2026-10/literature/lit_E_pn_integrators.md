# Literature survey E: PN orbit models vs numerical geodesics; PN bias; integrators
Date of survey: 2026-10-25.

## ACCESS STATUS (read first)
- Egress proxy returned 403 (CONNECT denied) for: arxiv.org, export.arxiv.org, ar5iv.labs.arxiv.org, www.aanda.org, iopscience.iop.org, ui.adsabs.harvard.edu, academic.oup.com, en.wikipedia.org.
- WebFetch returned getaddrinfo ENOTFOUND for the same hosts.
- Consequence: NO primary full text opened yet. Every number below comes from WebSearch result snippets (search-index summaries) and is tagged [snippet]. Items tagged [unverified] are garbled or single-source.
- Not retried against the proxy (policy denial). Environment network policy must be changed by the owner to proceed (see final report).

## Q1. PN orbit vs full numerical geodesic (Schwarzschild/Kerr, S-stars); astrometric difference; PN bias

Q1-a. Grould, Vincent, Paumard & Perrin 2017, A&A 608, A60, "General relativistic effects on the orbit of the S2 star with GRAVITY". arXiv:1709.04492 (id from search index URL). ADS 2017A&A...608A..60G.
  - Mock data: fully GR, ray-traced with GYOTO [snippet]. Not PN.
  - Effect-size table: NOT retrieved as a verified table. Snippet-level values, all [unverified]:
    * Lensing: ~20 uas near pericentre, ~2 uas typical, ~1 uas near apocentre; detectable in a few years at 10 uas astrometry and 10 km/s spectroscopy.
    * Lense-Thirring: ~10-40 uas over 3 orbital periods, orientation dependent (spin a=0.99, i'=45 deg, Omega'=160 deg).
    * Romer delay: detectable in ~4 months under models A and B (10 uas, 100 km/s).
    * Transverse Doppler: ~2 months at 10 uas for 10 km/s (model A).
    * Gravitational redshift: ~4 months at 10 uas for 100 km/s (model A).
    * Pericentre advance: "a few years" at 10 uas and 1-10 km/s (abstract, snippet).
    * Shapiro delay: no figure found.
    * Model table (models A-F vs full GR), per-period differences: ~3, 8, 16 mas for A-C and 7-40 uas for D-F. Snippet says the mapping of model to effect is unclear. Treat as unreliable.
  - Spin constraint (snippet): 1-sigma spin error ~0.4 after 16 yr falling to ~0.1 after 47 yr (10 uas, 10 km/s, a=0.99). Caveat: drops all photon-path effects except lensing. After three orbital periods: spin-axis inclination +-25 deg, node +-40 deg (snippet).
  - Where the numbers should be: section "Constraint on various effects" of the A&A full HTML (aanda.org aa31148-17) [not accessible].

Q1-b. Preto & Saha 2009, ApJ 703, 1743, "On post-Newtonian orbits and the Galactic-center stars". arXiv:0906.2226 (id from search index).
  - Kerr PN corrections written as Hamiltonian; corrections integrated numerically with an adaptive-step symplectic scheme [snippet].
  - S2 GR astrometric signal: tenths of a mas (hundreds of uas) in position; tens of km/s in velocity [snippet].
  - Error check: analytic pericentre-shift expressions (second-order Schwarzschild + first-order spin) [snippet]. NO comparison with full geodesic integration found in snippets.
  - Caveat in snippet: extended Newtonian mass (stars, gas, dark matter) can mimic GR for orbits wider than ~0.01 pc; separable because GR terms vary differently in time.

Q1-c. Hyman, Psaltis & Ozel 2023, ApJ 952, 35 (volume from a citing snippet) "Analytic Post-Newtonian Astrometric and Spectroscopic Models of Orbits around Black Holes". arXiv:2208.01655. DOI 10.3847/1538-4357/accb52.
  - First-order PN analytic astrometry + spectroscopy model, motivated by cost of numerical geodesic integration [snippet].
  - Argue numerical geodesic runs bundle same-order effects into one observable; disentangling needs on/off runs [snippet].
  - Shapiro effect possibly detectable for an S-star within ~a decade (30-m class telescopes) [snippet].
  - No uas accuracy vs numerical geodesic found in snippets [unverified].

Q1-d. Angelil & Saha 2010, ApJ 711, 157, "Relativistic redshift effects and the Galactic-center stars". arXiv:1001.1957 (snippet URL).
  - Orbital PN redshift and light-path (post-Minkowskian) redshift contributions of similar size, each a few km/s near S2 pericentre [snippet].

Q1-e. Gainutdinov & Baryshev 2020, Universe 6, 177, doi 10.3390/universe6100177, arXiv:2010.06386, "Relativistic effects in orbital motion of the S-stars at the Galactic Center". Constrains PPN beta, gamma from S-star data [snippet]. No geodesic comparison found.
  - Related: Gainutdinov arXiv:2002.12598 (PPN motion of S-stars). Snippet beta~0.97, gamma~0.81 with asymmetric errors ~0.4-0.65 [unverified, do not cite numbers].

Q1-f. Khlghatyan 2026, MNRAS 551(1), doi 10.1093/mnras/stag1407, arXiv:2609.05172, "PPN-spin degeneracies in mock S62-like stellar-orbit inference" (sole author per arXiv listing [snippet]).
  - Model: 1PN + spin-orbit Hamiltonian; synthetic astrometry + RV; Bayesian inference [snippet].
  - Upsilon (periapsis advance) recoverable at current precision; gamma and beta degenerate; Upsilon-chi strongly correlated [snippet].
  - Pairing a wide and a tight orbit separates spin from PPN [snippet].
  - Whether mocks were generated by exact geodesic integration: NOT FOUND.

Q1-g. de Mora Losada, Della Monica, de Martino & De Laurentis, arXiv:2410.22864, "Future prospects for measuring 1PPN parameters using observations of S2 and S62 at the Galactic Center". Published A&A 2025 (aanda aa52814-24, Feb 2025 URL) [snippet]; arXiv listing snippet says "accepted by MNRAS". VENUE CONFLICT, unresolved.
  - Mock S2 dataset + MCMC; GRAVITY + SINFONI, multi-period [snippet].
  - Model: PPN expansion of Schwarzschild [snippet]. Mock likely generated from the same PPN-type model (same-model injection/recovery would hide truncation bias) [inference].
  - S2+S62 forecast: gamma to ~1%, beta to ~2% [snippet].
  - S2 results in Section 4 of arXiv HTML [snippet]. Bias numbers not seen.

Q1-h. Liodis et al. 2026, arXiv:2605.20574, "Post-Newtonian orbital mechanics around a black hole in modified gravity" (MOG; first author Ioannis Liodis per snippet). PN S-star dynamics in MOG. Precession is the key effect [snippet]. No geodesic comparison in snippets.

DERIVED BY ME (not from a paper; assumptions stated):
  - S2 pericentre ~120 AU [snippet]. With M = 4.0e6 Msun (assumed), r_g = GM/c^2 = 5.9e9 m = 0.0395 AU, so r_p ~ 3.0e3 r_g (snippet says ~1400 R_S = 2800 r_g; consistent).
  - (v_p/c)^2 at pericentre with v_p ~ 7650 km/s [snippet]: ~6.5e-4 (matches snippet value for first-order correction).

## Q1 key answers (status)
- Bias question (fit exact-GR mock with PN model, report parameter bias): NOT FOUND in 5 targeted searches. Closest: 2410.22864 (same-model mock, no bias), 2609.05172 (mock source unknown), 2402.00333 (Yukawa; snippet claims simplified-metric fits off by up to ~1/5 vs f(R) metric [unverified]). Cannot assert nobody did it; search coverage limited and primary texts blocked.
- PN breakdown pericentre (r_g) or precision: NOT FOUND explicitly.
- Grould 2017 full table with all numbers: NOT VERIFIED (primary text blocked).

## Q2. Source/equation numbers for 1PN test-particle acceleration and 1.5PN spin-orbit acceleration; GRAVITY f_SP definition

Q2-a. 1PN test-particle acceleration (user form):
  a = -(GM/r^2) n + (GM/(c^2 r^2)) [ n (4GM/r - v^2) + 4 rdot v ]   (harmonic gauge)
  Literature status: formula NOT located in any retrieved text. Equation numbers: UNVERIFIED.
  Candidate sources (existence per search snippets; content unverified):
    * Will 2008, ApJ 674 L25, DOI 10.1086/528847 (snippet: arXiv 0711.1677 per snippet; unverified). Abstract is about orbital-plane (nodal) precession of short-period, high-e stars, 10 uas/yr order, Q2 = -J^2/M no-hair test. Snippet does NOT show PN EOM.
    * Merritt, Alexander, Mikkola & Will 2010, PRD 81, 062002 "Testing properties of the Galactic center black hole using stellar orbits" (authors/journal/volume confirmed by snippet). Abstract: stellar masking of frame dragging beyond ~0.5 mpc; quadrupole needs stars within ~0.2 mpc (snippet). EOM not seen.
    * Poisson & Will 2014, "Gravity: Newtonian, Post-Newtonian, Relativistic", CUP. Snippets: Will course outline ties PN celestial mechanics to Chapters 8-10 (unverified). Equation numbers: NOT FOUND.
    * Blanchet & Iyer 2003, CQG 20, 755 (arXiv gr-qc/0209089 per snippet). Harmonic-gauge relative EOM with a = -(Gm/r^2)[(1+A) n + B v].
  My consistency check (no source): with nu -> 0 the standard 1PN relative-acceleration coefficients A = v^2/c^2 - 4GM/(c^2 r), B = -4 rdot/c^2 give exactly the user's formula. A (1PN, general nu) = -(3/2) nu rdot^2 + v^2 + 3 nu v^2 - (Gm/r)(4+2 nu) is quoted in snippet of arXiv 2505.01278 (4th-order ISCO paper; unverified). B = -2(2-nu) rdot is from memory, unverified.

Q2-b. 1.5PN spin-orbit (test particle, central spin J):
  a_SO = (2G/(c^2 r^3)) [ 3 (n.J)(n x v) + (v x J) ]
  My own derivation (not sourced): stationary linear gravitomagnetic metric g_0i = -(2G/(c^3 r^3)) (J x x)_i. Test-particle Lagrangian L ~ m A.v + m U gives a = grad U + v x (curl A). Curl gives B = (2/r^3)[J - 3(J.n) n]. Then a = v x B = (2/r^3)[v x J + 3 (J.n)(n x v)]. Coefficient and sign confirmed.
  Independent check: equatorial prograde circular orbit, Kerr frequency Omega = sqrt(M)/(r^{3/2} + a sqrt(M)) (G=c=1; BPT form, standard). Expansion gives outward radial extra acceleration 2 a M^{3/2}/r^{7/2} = 2 G J sqrt(GM)/(c^2 r^{7/2}). The user formula gives |a_SO| = (2G/(c^2 r^3)) v J with v = sqrt(GM/r): identical, and outward (repulsive) for prograde. Consistent.
  CAUTION: my reconstruction of Kidder 1995 Eq. (2.2c) (test-mass reduction) gave -v x J + 3 rdot (n x J) and disagrees with the above; I do NOT trust my memory of Kidder (2.2c). Do not cite its numbers.
  Literature source candidates (existence only):
    * Kidder 1995, PRD 52, 821, arXiv gr-qc/9506022 (from a JHEP reference-list snippet). Equation (2.2c) NOT verified.
    * Merritt et al. 2010 and Will 2008: frame-dragging / nodal precession discussions (abstract level only).
  Equation numbers for any spin-orbit EOM: UNVERIFIED.

Q2-c. GRAVITY f_SP (S2 Schwarzschild precession, 2020):
  Paper: GRAVITY Collaboration 2020, A&A 636, L5, "Detection of the Schwarzschild precession in the orbit of the star S2 near the Galactic centre massive black hole". arXiv:2004.07187 (from search index URL). ADS 2020A&A...636L...5G.
  - Definition: dimensionless f_SP; 0 = Newtonian, 1 = GR [snippet, consistent across 2 snippets].
  - Multiplies: CONFLICTING snippet statements.
      (i) "multiplies the Schwarzschild precession term" (pericentre advance per orbit = f_SP x 3 pi R_S/[a(1-e^2)] = f_SP x 6 pi GM/[c^2 a (1-e^2)]); snippet: 12.1 arcmin per orbit for S2 at f_SP=1.
      (ii) "multiplies every Schwarzschild-metric term" (snippet 'ad hoc factor').
      Also snippet: f_SP=0 comparison orbit still carries Romer, SRT and 'RS' (gravitational redshift) terms. If accurate, f_SP does not touch redshift. The exact equation number: NOT VERIFIED.
  - Value (snippet): f_SP = 1.10 +- 0.19 (MCMC, data through 2019; 'agrees with GR'). Two snippets (ESO-hosted manuscript and summary) agree.
  - Later values (snippet only, source unclear, unverified): 0.85 +- 0.16 (S2-only, revised zero points); 0.997 +- 0.144 (S2+S29+S38+S55); 1.016 +- 0.226 (2024 S2 study). Do NOT cite without primary.
  - Mix-up warning: one snippet attributed f = 0.90 +- 0.09 (stat) +- 0.15 (sys) to the precession paper. That value belongs to the 2018 redshift paper (see below).

Q2-d. GRAVITY 2018 redshift factor (for comparison; not the same parameter):
  GRAVITY Collaboration 2018, A&A 615, L15, arXiv:1807.09409 (snippet URL), "Detection of the gravitational redshift in the orbit of the star S2 near the Galactic centre massive black hole". Snippet: factor f (0 Newtonian, 1 GR) multiplies the combined gravitational-redshift + transverse-Doppler PN contribution. f = 0.90 +- 0.09 (stat) +- 0.15 (sys) (snippet). Symbol 'f_RS' NOT confirmed in snippet.
  Do et al. 2019 (Science 365, 664; arXiv:1907.10731 per snippet): redshift parameter 0.88 +- 0.17 (snippet).

Q2-e. Poisson & Will 2014 textbook equation numbers: NOT FOUND.

## Q2 status
- 1PN test-particle EOM: formula verified by my consistency check against the standard relative-acceleration coefficients; literature equation numbers UNVERIFIED.
- Spin-orbit: formula derived + checked against Kerr prograde circular frequency; literature equation numbers UNVERIFIED; Kidder recollection unreliable.
- f_SP: value 1.10 +- 0.19 (snippet, 2 consistent). Exact term(s) multiplied: UNRESOLVED from snippets (conflict). Need A&A full text (aa37813-20) equation.

## Q3. Structure-preserving integrators for geodesics

Q3-a. Christian & Chan 2021, ApJ 909, 67, "FANTASY: User-friendly Symplectic Geodesic Integrator for Arbitrary Metrics with Automatic Differentiation". DOI 10.3847/1538-4357/abdc28. arXiv:2010.02237 (search URL). ADS 2021ApJ...909...67C.
  - Abstract (snippet) calls it symplectic; metric derivatives via automatic differentiation; "two coupled copies of phase space" (snippet, unverified).
  - CONFLICT resolved by primary GitHub README (github.com/pierrechristian/FANTASY, fetched 2026-10-25): "Is technically semi-symplectic instead of symplectic." README cites "Christian and Chan, 2021"; no accuracy numbers, no orbit-count benchmark. Qualitative only: "The smaller omega is, the smaller the integration error."
  - EinsteinPy wraps FANTASY (orders 2,4,6,8; snippet from EinsteinPy docs).

Q3-b. Explicit symplectic integrator series (Wang/Wu et al.), 'Construction of explicit symplectic integrators in general relativity':
  - Part I, Schwarzschild: Wang Y., Sun W., Liu F., Wu X., ApJ 907, 66 (2021) [from a reference-list snippet; not verified against journal]. Arxiv id not seen.
  - Part II, Reissner-Nordstrom: Wang Y., Sun W., Liu F., Wu X., ApJ 909, 22 (2021) [reference-list snippet]. arXiv:2103.02864 (URL in search index).
  - Part III, RN-(A)dS: ApJS 254, 8 (2021) [reference-list snippet]. Not ApJ.
  - Part IV, Kerr black holes: Wu X., Wang Y., Sun W., Liu F., ApJ 914, 63 (2021), DOI 10.3847/1538-4357/abfc45, arXiv:2106.12356 (IOP page; published 2021 June 16). Journal volume/page from snippet.
    * Method: time transformation function; Kerr Hamiltonian split into FIVE pieces, each with closed-form solution in new time (snippet).
    * Carter constant (snippet): EP4 most accurate; RK4 the only compared scheme failing to conserve it; S4, IE4, EP4 keep it. Unstable spherical-photon test: RK4 Hamiltonian and Carter errors peak near t ~ 116; other schemes conserve both through t = 2000 (time unit not stated in snippet). Rounding error becomes significant when truncation error ~1e-8 (snippet).
  - Review: Wu X., Wang Y., Sun W., Liu F.-Y., Han W.-B., ApJ 940, 166 (2022), "Explicit symplectic methods in black hole spacetimes", DOI 10.3847/1538-4357/ac9c5d, arXiv:2210.13185 (snippet). Rotating black ring: 13 pieces; time-transformed Kerr-Newman with disformal parameter: 8 pieces (snippet).
  - arXiv:2412.01045 (2024; 'Explicit symplectic integrators with adaptive time steps in curved spacetimes'): direct Kerr split tops out at two pieces because Sigma appears in denominators (snippet).
  - Note: the user's 'Wang, Wu et al.' - first author is Ying Wang for Parts I-III and Xin Wu for Part IV and the 2022 review.
  - 2021-2023 PRD listings: NOT found.

Q3-c. Seyrich & Lukes-Gerakopoulos 2012, PRD 86, 124013 (APS DOI 10.1103/PhysRevD.86.124013; ADS 2012PhRvD..86l4013S), title "Symmetric integrator for nonintegrable Hamiltonian relativistic systems", arXiv:1207.3175 (search index).
  - Content (snippet): time-symmetric symplectic integrator with newly designed step-size control; claims reversibility and conservation of constants of motion; targets long-run geodesic motion in nonintegrable spacetimes; claims it beat standard symplectic schemes and an adaptive RK on speed and accuracy.
  - NOT a Gauss-collocation paper for Kerr. The user's description ('Gauss collocation for Kerr') does not match the snippet.
  - Seyrich 2013: arXiv:1302.4580 'Gauss collocation methods for efficient structure preserving integration of post-Newtonian equations of motion' (PN EOM, not Kerr geodesics; acknowledges Lukes-Gerakopoulos and Lubich per snippet; journal unverified).

Q3-d. Radau IIA / Gauss-type implicit RK for geodesics:
  - No paper found using Radau IIA for geodesics (3 targeted searches).
  - Closest: arXiv:1801.02378 (OSTI 1542035) 'Generalized, energy-conserving numerical simulations of particles in general relativity. I. Time-like and null geodesics' (authors NOT retrieved). Benchmarks against analytic Schwarzschild and Kerr orbits; implicit midpoint baseline (= 1-stage Gauss-Legendre, standard fact).
  - Background (standard numerical analysis, not sourced here): 3-stage Radau IIA has order 5 (Radau IIA(5)); it is stiffly accurate and L-stable but NOT symplectic. s-stage Gauss-Legendre has order 2s and is symplectic. Cooper (1987, IMA J. Numer. Anal.) shows symplectic RK methods preserve quadratic invariants with CONSTANT coefficients. Carter's constant is quadratic in momenta but its coefficients depend on position, so that theorem does not transfer directly (background; unverified citation details).
  - Cross-check of user code: RelatiPy exposes methods radau (default), dop853, dp45, projection_radau in src/relatipy/geodesic/orbit.py (lines 33, 613-626), and kerr_mcmc.py accepts dop853 or radau (line 46). (Local repo grep, 2026-10-25.)

Q3-e. Projection onto constants of motion (norm, E, Lz, Carter Q): standard or novel?
  - No source found in retrieved text. Background: projection (post-step correction onto the invariant manifold) is standard in geometric integration (Hairer, Lubich & Wanner, 'Geometric Numerical Integration', Springer, 2nd ed. 2006, chapter on first integrals / projection; section numbers unverified).
  - Known criticisms (background, unverified): projection generally destroys symplecticity and time-reversibility; can add secular drift in other invariants if the projection is not onto all invariants consistently; needs the Jacobian of the invariants; can make the step non-smooth. Novelty of RelatiPy's projection_radau: NOT ASSESSED (no literature hit).

Q3-f. Carter-constant drift over how many orbits: NOT FOUND.
  - 2106.12356 reports t = 2000 (unit not stated), not orbit counts.
  - arXiv:2607.27129 (2026, 'Quadrupolar tidal effects destroy the integrability of black hole geodesics') reports drift vs proper time up to 2546 for a perturbed Kerr orbit (snippet; not orbits; Carter not conserved in that perturbed model).
  - Search budget exhausted before a direct query (see status).

## Q3 status
- FANTASY symplectic status: resolved by primary README: semi-symplectic.
- Seyrich/Lukes-Gerakopoulos 2012: exists, but is NOT Gauss collocation for Kerr (user premise incorrect per snippet).
- Wang/Wu series: Parts I-IV and 2022 review located with IDs/DOIs (volume/page from snippets).
- Radau IIA for geodesics: no literature hit.
- Projection: standard background; no novelty or criticism sources retrieved.
- Carter drift over orbits: not found.

## Q4. PyGRO and Della Monica / de Martino S2 papers (2021-2025)

Q4-a. PyGRO. arXiv:2504.20152 'PyGRO: a Python Integrator for General Relativistic Orbits', 2025. Primary GitHub README (github.com/rdellamonica/pygro, fetched):
  - BibTeX: author = {Riccardo Della Monica} ONLY; year 2025; eprint 2504.20152; primaryClass gr-qc; journal not in BibTeX. 'de Martino' does NOT appear on the page. (The user's 'Della Monica & de Martino' attribution is not supported by this README.)
  - Journal: A&A, accepted 2025, aa54300-25 per sibling-notes search snippet (UNVERIFIED; aanda blocked).
  - License GPL-3.0. Version: README badge 1.0.2; master pyproject.toml says version "1.0.1" (discrepancy).
  - Dependencies (master pyproject.toml, fetched): numpy>=2.2.2, scipy>=1.15.1, sympy>=1.11.1.
  - Scope (README): numerical integration of geodesics for timelike and null paths; any metric given as an analytic expression (sympy).
  - Integrator (pygro/integrators.py, fetched via WebFetch summary of code): custom explicit adaptive embedded Runge-Kutta, keys: dp45 (Dormand-Prince 5(4)), dp853 (DormandPrince853, 13 stages), rkf45, rkf78, ck45. Error norm atol + rtol*||y||; step rejection on err > 1; atol = rtol = 1e-10 default (accuracy_goal = precision_goal = 10); hmax 1e16, hmin 1e-16; max_iter 1000; step factor clamps 0.2 and 10. NOT symplectic; NO Radau; NO scipy solve_ivp call (per code). Docstring gives no order or conservation claim.
  - Metrics shipped: not stated in README (folder 'default_metrics' exists in pygro/, not inspected).
  - MCMC / fitting: NOT in README. Compute cost / timing: NOT stated. Likelihood counts: NOT stated.
  - Minor code notes (not asked): CashKarp45.integrate uses self.verbose never assigned; CashKarp45.next_step computes Qw = min(Q2, 10*QUIT2) but uses Q2 (apparent typo). Per fetch summary.
  - GitHub profile (rdellamonica): only pygro visible as pinned repo; other 8 repos not visible in fetched page.

Q4-b. Della Monica and/or de Martino S2 papers 2021-2025 (snippet-level only; authors marked where unverified):
  1. arXiv:2410.22864, 'Future prospects for measuring 1PPN parameters using observations of S2 and S62 at the Galactic Center'. Authors (snippet): de Mora Losada, Della Monica, de Martino, De Laurentis. Venue: A&A 2025 (aa52814-24 URL) vs 'accepted MNRAS' (arXiv snippet): CONFLICT. Method: mock S2 dataset; MCMC fitting; PPN expansion of Schwarzschild (static, spherically symmetric) (snippet). NOT a numerical geodesic integration per snippets. Forecast: gamma ~1%, beta ~2% with S2+S62 (snippet). Bias numbers not seen.
  2. arXiv:2306.17215, 'Using the motion of S2 to constrain scalar clouds around SgrA*' (2023). Authors UNVERIFIED.
  3. arXiv:2312.02653, 'Using the motion of S2 to constrain vector clouds around SgrA*' (2023). Authors UNVERIFIED.
  4. arXiv:2504.02908, 'Exploring the presence of a fifth force at the Galactic Center' (2025). Authors UNVERIFIED.
  5. arXiv:2602.18994, 'Testing the Spacetime Geometry of Sgr A* with the Relativistic Orbit of S2 star' (2026, v2 html). Authors UNVERIFIED.
  6. 'Unveiling the nature of SgrA* with the geodesic motion of S-stars' (ResearchGate record 356792460; ~2021). Authors/venue UNVERIFIED. An inspirehep-hosted PDF carries 'Citation: Della Monica, R.; de Martino, I.; de Laurentis, M.' (MDPI-style), so this may be a three-author paper; unverified.
  7. 'Relativistic orbits of S2 star in the presence of scalar field', EPJC (2024), DOI-suffix s10052-024-12477-3 (Springer). Authors UNVERIFIED.
  8. arXiv:2402.00333, 'Constraining the Yukawa gravity with post-Newtonian approximation using S-star orbits around the supermassive black hole in our Galactic Center' (2024). Authors UNVERIFIED. Snippet: mock-data section; claims simplified-metric fits can be off from f(R)-derived metric by up to ~1/5 even at small |kappa| (unverified).
  - Integrator, MCMC setup, likelihood evaluation counts, compute cost for items 2-8: NOT RETRIEVED (search budget exhausted).

## Q4 status
- PyGRO integrator: identified from primary code (custom explicit adaptive RK, not symplectic/Radau).
- PyGRO authorship: Della Monica only (BibTeX and README).
- Per-paper MCMC/integrator/evaluation counts: NOT FOUND.

## Session-level status (for final report)
- Primary-source access blocked by environment network policy: arxiv.org, export.arxiv.org, ar5iv.labs.arxiv.org, www.aanda.org, iopscience.iop.org, ui.adsabs.harvard.edu, academic.oup.com, en.wikipedia.org (403 CONNECT). WebFetch: getaddrinfo ENOTFOUND for these.
- GitHub (github.com web pages) and search snippets were usable.
- WebSearch budget for this turn exhausted (200 per turn, shared across agents) after ~32 of my own calls; remaining items (Carter drift vs orbits; Radau/geodesic papers; Q4 per-paper details; Grould table numbers) need a follow-up search.
