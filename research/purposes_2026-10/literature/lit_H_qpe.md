# Literature survey: timing models for QPEs and repeating pTDEs (Kerr geodesics)

Status legend: VERIFIED = seen in a source title/listing; SUMMARY = paraphrase from a search-snippet summary, primary text NOT read; UNVERIFIED = could not confirm.
Access: arXiv, ADS, CDS, IOP, A&A, doi.org, zenodo blocked for fetch. Semantic Scholar, INSPIRE, Astrobites, Wikipedia, RG, emergentmind, par.nsf.gov, indico: DNS failure in WebFetch. WebSearch budget (200/turn) was exhausted during Q1-Q2 work; Q3-Q5 are therefore partial.

## Q1. Zhou et al. paper series (QPE orbit timing)

Citations
- Paper I: Zhou, C., Huang, L., Guo, K., Li, Y.-P., Pan, Z., "Probing orbits of stellar mass objects deep in galactic nuclei with quasi-periodic eruptions", arXiv:2401.11190 (VERIFIED title/ID/authors via arXiv listing snippet; v1 Jan 2024, revision May 2024 per snippet). Journal: "accepted by Physical Review D" (SUMMARY of Astrobites post). Volume/article id: UNVERIFIED.
- Paper II: "... II: population analysis", arXiv:2405.06429 (VERIFIED title/ID). Authors per snippet: Zhou, C., Zhong, B., Zeng, Y., Huang, L., Pan, Z. Journal: UNVERIFIED.
- Paper III: "... III: Long term evolution", arXiv:2411.18046 (VERIFIED title/ID; v1 28 Nov 2024; authors Zhou, C., Zeng, Y., Pan, Z). Later version titled "Secular evolution of quasi-periodic eruptions" (arXiv v3). An IOP listing titled "Secular Evolution of Quasiperiodic Eruptions" with DOI 10.3847/1538-4357/adcee2 surfaced in search; match to Paper III is probable but not opened. The DOI prefix 1538-4357 is the ApJ ISSN.

Method (all SUMMARY, primary text not read)
- Paper I: flare = stellar-mass object striking a disk left by an earlier TDE. Orbit constrained from flare onset times. For GSN 069: e < 0.1 (3 sigma), a of order 10^2 r_g, wet-EMRI channel favoured, loss-cone and Hills origins disfavoured. Whether a Kerr metric is used: UNVERIFIED (snippets mention EMRI equations of motion; a Oct 2024 talk by Pan reportedly says "geodesic").
- Paper II: a later paper's summary says Paper II "solved the geodesic equation head-on". Sampler: UNVERIFIED (a summary infers nested sampling from "live points"; inference only).
- Paper III: adds orbital decay from disk collisions and disk precession. Uses an osculating-trajectory scheme; abstract claims speed-up "by orders of magnitude". One summary gives "about 52 times faster" for GSN 069 sampling; another summary of the same paper did NOT find that number, so treat 52x as UNVERIFIED. Listed tools: bilby, nessai, pymultinest. Bayesian support for orbital decay in GSN 069 and disk precession in eRO-QPE2. Narrows M_SMBH, QPE radiative efficiency, orbiter type, disk surface density and viscosity. Orbits nearly circular (e of order 10^-2, a of order 10^2 r_g). Full-sample analysis: eccentricity splits sources into two groups; period-mass exponent about 0.8 in the low-e group.

Not found (UNVERIFIED / not stated in what I could read)
- Exact free-parameter count and list; sampler settings; likelihood-evaluation count; CPU time per evaluation; treatment of light-travel time, Shapiro delay, spin; spin constraints; which of eRO-QPE1/2/3/4, RX J1301.9, AT2019qiz, Ansky were fitted.
- Cost as a stated limitation: Paper III's osculating scheme is motivated by the cost of long-term orbit integration (qualitative, SUMMARY only).

Related 2026 papers surfaced (SUMMARY only)
- arXiv:2610.10110, "Quasi-periodic eruption spectral-timing: EMRIs crossing warped accretion disks" (dated 7 Oct 2026; authors not confirmed). Kerr-geodesic frequency model + quiescent disk spectrum; joint fit of GSN 069 and eRO-QPE2; 17 free parameters, emcee; posteriors depend on starting positions; eRO-QPE2 alone admits six comparably good modes.
- arXiv:2610.21844, "Secular Evolution of Quasiperiodic Eruptions: eRO-QPE2": eccentricity upper limits 0.026-0.033 (95%); uses a forced orbit, not a Kerr geodesic.
- arXiv:2608.05534, "No evidence for a supermassive black hole binary in GSN 069": light-travel-time companion search; fitted amplitude consistent with zero.

## Q2. Other QPE timing / orbital analyses

Found (with status)
- Xian et al. 2021: arXiv:2110.10855 (VERIFIED ID via search link). SUMMARY: argues eruption timing could gauge BH spin and test the no-hair theorem. No numerical spin forecast seen. Journal/volume: UNVERIFIED.
- Linial & Metzger 2023: arXiv:2303.16231, "EMRI + TDE = QPE: Periodic X-ray Flares from Star-Disk Collisions in Galactic Nuclei" (VERIFIED title/ID; author and ApJ 957, 34 are from the brief, UNVERIFIED).
- Miniutti et al. 2019: arXiv:1909.04693, "Nine-hour X-ray quasi-periodic eruptions from a low-mass black hole galactic nucleus" (VERIFIED title/ID). Nature 573, 381 (2019), DOI 10.1038/s41586-019-1556-x: a RePEc record with this DOI surfaced in search (not opened).
- Miniutti et al. 2025 (GSN 069 timing): SUMMARY only. Reports correlated delays that conflict with the anti-correlation expected from EMRI apsidal precession. Citation/journal UNVERIFIED.
- Chakraborty et al. 2025: arXiv:2508.20162, "Prospects for EMRI/MBH parameter estimation using Quasi-Periodic Eruption timings: short-timescale analysis" (VERIFIED title/ID; ApJ per listing; DOI 10.3847/1538-4357/ae003b per IOP URL in search). Author lead Joheen Chakraborty (SUMMARY). Method (SUMMARY): analytic Kerr-geodesic expressions expanded to 3PN, GPU-accelerated Bayesian inference ("QPE-FIT"). Results (SUMMARY): about 10% precision on M_BH, EMRI semi-major axis and eccentricity within tens of orbits for e ~ 0.1-0.3. Spin NOT recoverable for a >= 100 r_g with 10-100 orbit spans. Disk precession constrained to 10-50% when disk misaligned and precessing. Chakraborty et al. 2021: not found.
- Arcodia et al. 2026: arXiv:2604.09788, "Even a precessing clock is right twice per orbit - The super-periods of eRO-QPE2 and challenges for quasi-periodic eruption orbital models" (VERIFIED title/ID; first author R. Arcodia per SUMMARY; ApJ May 2026 per SUMMARY). Numbers (SUMMARY, check against text): 32 eruptions over about one month (about 325 cycles). Periodic clock near 2.24 h with super-periods about 4.4 d and about 95 d. Period-derivative bound about 2e-6 s/s (3 sigma). Long-modulation lower limit about 918 periods at 1 sigma in the body vs about 1000 in the abstract (inconsistent in summary; check). Disk-instability fits as a damped random walk. Colliding-EMRI models survive only if one eruption per orbit is seen.
- Arcodia et al. 2024 (eRO-QPE3/4 discovery): arXiv:2401.17275 (ID from search; authorship and journal UNVERIFIED). Summary says small-number statistics limit the population work.
- Arcodia et al. 2024, "Ticking away: the long-term X-ray timing and spectral evolution of eRO-QPE2": title from Coventry pure-portal listing; arXiv ID and journal UNVERIFIED.
- Di Wang 2024: arXiv:2402.07177, "The orbital evolution of the tidally stripped star and disk-driven stable mass transfer for QPEs in GSN069" (VERIFIED title/ID; author per SUMMARY). SUMMARY: mass transfer controls the orbit, so the stripped object is not destroyed.
- Liu, Liu, Pan, Deng, Shen, Yu 2026: arXiv:2603.00226, "Quasi-periodic Eruptions from Stellar-mass Black Holes Impacting Accretion Disks in Galactic Nuclei" (VERIFIED title, authors, arXiv submission 27 Feb 2026). Simulation paper, not a timing fit. SUMMARY: stellar impactors give one bright flare per orbit (front flare more than 10x the rear); BH impactors give balanced ejecta; stellar impacts do not reproduce alternating long-short recurrence; BH collisions favoured for the 10^44-10^48 erg range.
- Not searched / not found: Franchini et al. 2023 (A&A 675 A100?) - not identified in any result; an A&A PDF "aa46565-23" surfaced but was not identified. Arcodia et al. 2021 (Nature 592, 704). Pasham et al. 2024 (AT2019qiz, Nature). Nicholl et al. 2024 (Nature 634, 804). Kejriwal, Witzany, Zajacek, Pasham, Chua 2024 (LISA EMRI counterparts). All UNVERIFIED (no search budget left).

Spin forecast: the only numeric statement is the SUMMARY of arXiv:2508.20162 above (spin not recovered at a >= 100 r_g). No other spin precision forecast found.

## Q3. Repeating partial TDEs (ASASSN-14ko, AT2018fyk, AT2020vdq, AT2022dbl, eRASSt J045650-203750)

Status: NOT SEARCHED. The WebSearch budget for this turn was exhausted during Q1-Q2, and none of the search results touched repeating pTDEs. No citations, periods, or Kerr-timing results for these sources are reported here.
Needed follow-up queries: Payne et al. 2021 (ASASSN-14ko, ApJ 910, 125 per brief, UNVERIFIED); ASASSN-14ko period about 114 d and negative period derivative (brief, UNVERIFIED); AT2018fyk rebrightening predictions; AT2020vdq; AT2022dbl; eRASSt J045650-203750 (Liu et al. per brief, UNVERIFIED). Also: any Kerr-geodesic or precession timing fit of these sources (pericentre in r_g, spin, Lense-Thirring/relativistic precession).

## Q4. Machine-readable QPE eruption-time tables

Result: NO machine-readable eruption-time table identified in this run. Partial leads only.
- Observations review: arXiv:2607.15359, "Observations of X-ray quasi-periodic eruptions" (2026; title and ID VERIFIED from search; SUMMARY says it counts 13 known QPE sources and covers population and selection effects). Whether it tabulates eruption times: UNVERIFIED.
- Arcodia et al. 2026 (arXiv:2604.09788): 32 eRO-QPE2 eruption timings used (SUMMARY). Whether they are released as a table or data file: UNVERIFIED.
- eRO-QPE3/4 discovery paper: arXiv:2401.17275 (ID only; tables UNVERIFIED).
- Rate study: arXiv:2403.17059, "Cosmic hide and seek: the volumetric rate of X-ray quasi-periodic eruptions" (title VERIFIED; uses four eROSITA QPE sources; not a timing table).
- Not checked (no budget): Miniutti et al. 2019 (arXiv:1909.04693), Pasham et al. 2024, Nicholl et al. 2024, Chakraborty et al. 2025 (arXiv:2508.20162) data availability. Zenodo blocked. A GitHub repository search was not run in this turn.
- Sources listed in the brief (GSN 069, eRO-QPE1/2, RX J1301.9): no tables checked.

## Q5. Statements that a fast exact Kerr timelike-geodesic solver with event detection is needed or is the bottleneck

No explicit statement found that event detection (equatorial or disk-plane crossings) is required. Cost statements found (all SUMMARY unless marked):
- arXiv:2508.20162 (Chakraborty et al. 2025, ApJ per listing): the cost bottleneck is the trajectory step. GPU port gives about 2.4 s (CPU) vs 0.012 s (GPU) per trajectory step with 1000 walkers (about 200x). Single likelihood evaluation about 0.26 s -> 0.015 s at 100 walkers (about 17x), about 105x at 1000 walkers. A summarizer's arithmetic put the trajectory step at 85-90% of CPU time (UNVERIFIED). The model uses analytic Kerr-geodesic expressions expanded to 3PN, which trades exactness for speed; the accuracy against exact Kerr integration is UNVERIFIED. Statement that Bayesian fits need "millions of model calls": attribution UNVERIFIED.
- arXiv:2411.18046 (Zhou, Zeng, Pan; Paper III): the osculating-trajectory method is motivated by the cost of long-term EMRI orbit calculation, with speed-up "by orders of magnitude" (abstract, SUMMARY). Paper II reportedly solved the geodesic equation head-on (SUMMARY, UNVERIFIED). This is an implicit statement that direct geodesic integration is the bottleneck, not an explicit one.
- Public code: "KerrGeoPy: A Python Package for Computing Timelike Geodesics in Kerr Spacetime" (ResearchGate title VERIFIED; authors, date, method, event-detection support: UNVERIFIED; ResearchGate fetch failed).
- Inference, not a stated claim: the collision timing model needs disk-plane crossings, so an event finder is implied by the model structure.
