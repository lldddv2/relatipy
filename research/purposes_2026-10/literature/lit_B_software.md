# Competitive landscape: Kerr timelike geodesic software (lit_B_software)

Context: RelatiPy (pypi relatipy; github lldddv2/relatipy). Timelike Kerr geodesics in C
(Radau IIA(5), DOP853, DP45, projection integrator), Astropy-units Python API, KerrMcmcModel
(sky offsets arcsec + LOS velocity km/s, Romer delay). Claimed ~150 full S2 forward models/s/core.

Fields per tool: language | timelike/null | analytic/numerical | generic inclined eccentric bound Kerr |
pip-installable | published speed | observables (astrometry/redshift/LTT) | reference

## Status log
- [start] Loaded WebSearch/WebFetch via ToolSearch. Beginning per-tool survey.

### [1] KerrGeoPy (BHPT, Park & Nasipak, NASA Goddard)
- Lang: Python (C? unverified). Timelike: YES (bound stable + plunging). Yachim fork (github.com/Yachim/KerrGeoPy) also claims lightlike (fork relation to official unverified).
- Method: analytic (Mino-time elliptic/Weierstrass solutions; Fujita & Hikida 2009 + Dyson & van de Meent 2023 plunges, arXiv:2302.03704).
- Generic inclined eccentric bound Kerr: YES (per search snippet: four-velocities, fundamental frequencies, constants of motion, separatrix). Needs paper confirmation.
- pip: YES (`pip install kerrgeopy`, PyPI 0.9.3; conda-forge). Repo github.com/blackholeperturbationtoolkit/kerrgeopy; paper arXiv:2406.01413.
- Speed numbers: NOT YET (check arXiv:2406.01413).
- Observables (astrometry/redshift/LTT): NOT in snippets; likely none (verify).
- Refs: arXiv:2406.01413 (KerrGeoPy); Fujita & Hikida 2009; Dyson & van de Meent 2023 arXiv:2302.03704.

### [2] BHPT KerrGeodesics (Mathematica)
- github.com/BlackHolePerturbationToolkit/KerrGeodesics, "Code to compute quantities related to bound timelike Kerr geodesics". MIT license. Wolfram Language.
- Timelike: YES (bound). Plunging analytic solutions reportedly incorporated (arXiv:2302.03704 says so; unverified in repo itself).
- Separatrix algorithm: Stein & Warburton arXiv:1912.07609 says new separatrix-locating algorithm ~45x faster than prior robust approach.
- pip: NO (Mathematica). Generic inclined eccentric: likely YES (to verify). Observables: none expected (geodesic quantities; verify). Speed (absolute): NOT YET.
- Latest tag v0.9.0 (2023-07-03) per snippet; repo listing shows updates later (unverified).

### [3] FastEMRIWaveforms "few" (Katz, Chua, Speri, Warburton, Hughes 2021)
- Ref: Phys. Rev. D 104, 064047 (2021); arXiv:2104.04582; DOI 10.1103/PhysRevD.104.064047. (NOTE: 5th author = Scott A. Hughes; my query wrongly said "Stein".) MIT open-access copy exists.
- Project: github.com/BlackHolePerturbationToolkit/FastEMRIWaveforms ; bhptoolkit.org/FastEMRIWaveforms. Snippet says later reference lists cite v2.0.0.
- Lang: Python front end + C/C++ (+CUDA GPU). pip: `fastemriwaveforms` (verify name). Waveform generator = fluxes/inspirals (not astrometry).
- Timelike: YES (EMRI, massive secondary). Kerr eccentric equatorial fast fluxes; generic inclination scope needs verification.
- Speed: published numbers NOT YET (paper says "rapid" model; check specifics).
- Observables: GW strain (no astrometry, no redshift/LTT on S2-type observer model) - verify.

### [4] Fujita & Hikida 2009 (analytic)
- Title: "Analytical solutions of bound timelike geodesic orbits in Kerr spacetime", CQG 26 135002 (2009), DOI 10.1088/0264-9381/26/13/135002, arXiv:0906.1420 [gr-qc], OSTI 21375492.
- Method: analytic; elliptic integrals, Mino time; claims first analytic expressions for fundamental frequencies. Generic bound orbits: YES (inclined eccentric, fundamental freq). Paper only (no code). Not pip.
- Used by KerrGeoPy.

### [2b] KerrGeodesics README (fetched github.com/BlackHolePerturbationToolkit/KerrGeodesics)
- About: "Code to compute quantities related to bound timelike Kerr geodesics"; README: "provides a set of functions for calculating properties of bound timelike geodesic orbits about a Kerr black hole."
- Named quantities: NONE listed in README (freq/Mino/separatrix not named there). Generic inclined eccentric: not stated in README (paclet; Kernel/ folder + PacletInfo.wl => Wolfram paclet). Install: not stated in README. License MIT, "Copyright 2017-2020". Docs in Documentation/English; homepage not fetched.
- Timing numbers: NOT stated. Observables (redshift/astrometry/LTT/coord time): NOT stated.
- => Need the Documentation or the 1912.07609 / 2302.03704 papers for scope & speed.

### [1b] arXiv fetch attempt (WebFetch) FAILED: "getaddrinfo ENOTFOUND arxiv.org" for arxiv.org/abs/2406.01413 and arxiv.org/abs/2104.04582. Retrying via alternate hosts (zenodo, yachim.cz docs, github README, researchgate, hackernoon).

### [1c] KerrGeoPy README (fetched github.com/blackholeperturbationtoolkit/kerrgeopy) -- ANSWERS for KerrGeoPy
- Lang: "a python implementation of the KerrGeodesics Mathematica library" (pure Python vs compiled: not stated; uses scipy>=1.8 funcs, matplotlib>=3.7 for plots).
- Install: `pip install kerrgeopy` and `conda install -c conda-forge kerrgeopy` (YES pip).
- Timelike: YES. Lightlike: not stated in README (JOSS/arXiv title = timelike).
- Bound/plunging: StableOrbit(a,p,e,x) [bound, eccentric 0<e<1, x=cos(incl) in [-1,1]] and PlungingOrbit(a,E,L,Q). Orbit(a,x0,u0) from initial 4-velocity.
- Generic inclined eccentric bound: YES (x in (-1,1)), BUT polar (x=0) NOT supported; marginally bound NOT supported; extremal Kerr (a=1) NOT supported. Equatorial is x=+/-1 -- needs care (unverified in README).
- API outputs: trajectory() -> callables t(lambda), r(lambda), theta(lambda), phi(lambda) as functions of MINO time (t,r in units of M); constants_of_motion() -> dimensionless E,L,Q; mino_frequencies() -> Upsilon_r, Upsilon_theta, Upsilon_phi, Gamma; fundamental_frequencies() -> Omega_r,Omega_theta,Omega_phi in Boyer-Lindquist time. Separatrix / 4-velocity along orbit / observer-frame: NOT stated.
- Performance numbers: NOT stated (none in README).
- Observables (redshift / astrometry / light-travel time): NOT stated; stated purpose = EMRI trajectories. => No astrometric/spectroscopic/LTT observables.
- Citation: JOSS doi 10.21105/joss.06587 (Park & Nasipak); Zenodo 10.5281/zenodo.8417442; arXiv 2406.01413. Method refs: Dyson & van de Meent arXiv:2302.03704 (plunges); Fujita & Hikida arXiv:0906.1420 (stable). License MIT.
- Docs: kerrgeopy.readthedocs.io (not fetched yet).

### [1d] Fetch status
- arxiv.org: WebFetch DNS fail; curl via proxy -> "CONNECT tunnel failed, response 403" (proxy denies arxiv). yachim.cz: DNS fail (WebFetch). Next: check proxy README for per-tool fix, try JOSS (joss.theoj.org) and readthedocs.

### [PROXY] Access policy findings (do NOT retry or route around; README says report blocked hosts)
- Proxy (/root/.ccr/README.md + __agentproxy/status) returns 403 CONNECT policy denials for: arxiv.org:443, export.arxiv.org:443, www.aanda.org:443, iopscience.iop.org:443, b2find.eudat.eu:443 (also cdn.jsdelivr.net:443 once).
- CONSEQUENCE: arXiv full texts (KerrGeoPy 2406.01413, FEW 2104.04582, Dyson&vdM 2302.03704, Stein&Warburton 1912.07609, Fujita&Hikida 0906.1420, Gralla&Lupsasca, van de Meent) and A&A/IOP full texts are NOT fetchable. Use WebSearch snippets + GitHub/PyPI/docs pages for content; mark anything from snippets only as "snippet".
- Reachable so far: github.com (WebFetch OK), zenodo (in search results), osti.gov (search result), DNS-failing: yachim.cz.

### [3b] FEW README (fetched github.com/BlackHolePerturbationToolkit/FastEMRIWaveforms)
- Lang: Python front end ("import few"), C++ compiler needed for source build, CUDA optional (nvcc). Python >=3.12 for source; wheels 3.12/3.13/3.14.
- pip: YES: `pip install fastemriwaveforms` (CPU); `fastemriwaveforms-cuda12x`; `fastemriwaveforms-cuda13x`. conda-forge CPU (Linux/macOS/Windows), CUDA plugin Linux-only.
- Kerr generic inclined eccentric geodesic module: NOT stated in README; no geodesic/trajectory utility documented there (module list: trajectory, amplitudes, summation, utilities). Dataset dir "dataset-preprocessing/KerrEccEq" hints at Kerr eccentric equatorial. => generic Kerr geodesic: NOT CONFIRMED.
- Speed: only tagline "Blazingly fast EMRI waveforms"; no numbers (test-suite time 71.5 s for 20 tests, not a waveform benchmark).
- Observables: GW waveforms only (strain/flux/phase). NO astrometry, redshift, LTT stated.
- Citation: cite FEW papers + Zenodo 10.5281/zenodo.3969004. License MIT. SemVer.
- Paper-level speed numbers (paper PRD 104 064047): NOT retrievable (arXiv blocked). Only snippet: "rapid, accurate, fully relativistic model" (Chua et al. PRL 2021 lineage).

### [3c] JOSS fetch FAILED: "getaddrinfo ENOTFOUND joss.theoj.org" (DNS). JOSS paper 10.21105/joss.06587 NOT readable -> KerrGeoPy perf numbers still unknown (none in README).

### [5] GYOTO (Vincent, Paumard, Gourgoulhon, Perrin 2011, CQG 28 225011; arXiv:1109.4769; OSTI 22021240)
- Lang: C++ core (Python bindings: UNVERIFIED; check gyoto.obspm.fr). 
- Timelike: YES per abstract snippet ("integrates null and timelike geodesic equations in the Kerr metric and in numerically computed metrics"). Null: yes.
- Method: numerical integration (not closed form). Generic inclined eccentric Kerr: presumably yes (general metric integrator) - UNVERIFIED specifics.
- pip: NOT known (UNVERIFIED). Speed: NOT found in snippets.
- Observables: ray tracing images (null) + massive-particle paths; authors intended stellar orbit fitting around Kerr (snippet). Redshift/LTT: via ray-tracing images / photon paths; not verified for massive-particle astrometry.
- 2016 validation: Vincent et al. A&A 2016 "Testing the validity of the ray-tracing code GYOTO" (arXiv:1605.04195) - validation vs Geokerr is null-geodesic; no massive-particle check found in snippet.
- Cite: arXiv:1109.4769; CQG 28 225011 (2011) (volume/number from snippet PDF header, verify).

### [7] van de Meent 2020 (requested item) -- NOT FOUND
- Targeted search returned no 2020 van de Meent paper on Kerr timelike geodesics. Closest hits (snippets only, not verified): Fujita & Hikida 2009 (arXiv:0906.1420); van de Meent & Warburton 2018 "isofrequency pairing of geodesic orbits in Kerr geometry" (ResearchGate, venue/year unverified); van de Meent & Shah (PRD; metric perturbations eccentric equatorial); Dyson & van de Meent arXiv:2302.03704 (plunges, 2023).
- => Status: "unverified" -- do NOT cite a 2020 van de Meent paper until identified. Ask / re-search with alternative title words.
- Mino-time frequency papers 2025: arXiv:2508.00516 (coordinate-time frequency = Mino freq / common factor), arXiv:2511.06684 (analytic bound geodesics in EOB-deformed Kerr, Mino time).

### [6] geokerr (Dexter & Agol 2009, ApJ 696, 1616)
- Lang: Fortran 77 (releases 3/4/2009, 5/18/2009 on Agol UW page: faculty.washington.edu/agol/geokerr/index.html). ASCL ascl:1011.015. Python f2py wrapper exists per snippet (unverified).
- Timelike: NO -- null (photon) geodesics only (ASCL: "Fortran program for computing null geodesics in the Kerr metric").
- Method: semi-analytic, Carlson elliptic integrals (Hamilton-Jacobi); needs turning points pre-set (per YNOGK authors).
- pip: NO. Speed: none in snippets. Observables: photon orbit only (used inside grtrans for radiative transfer); no astrometric/redshift for massive particles.
- Exact ApJ title NOT confirmed by search (memory suggests "A Fast New Public Code for Computing Photon Orbits in a Kerr Spacetime" -- UNVERIFIED).

### [8] YNOGK (Yang & Wang 2013, ApJS 207, 6)
- DOI 10.1088/0067-0049/207/1/6; arXiv:1305.1250; ASCL ascl:1305.008. Fortran code, free download from authors.
- Timelike: NO (null). Method: Weierstrass + Jacobi elliptic functions; no turning points needed (improvement over geokerr).
- pip: NO. Speed: NOT in snippets. Observables: photon paths only.
- FOLLOW-UP found (not in user list, relevant!): ynogkm -- arXiv:1311.4436 "ynogkm: A New Public Code For Calculating time-like Geodesics In The Kerr-Newman Spacetime"; A&A 2014 (aa22565-13, full html at aanda.org -- BLOCKED host, so only snippet). Timelike: YES (Kerr-Newman, Kerr as charge=0). Analytic elliptic-function: per snippet "same team carries approach over". Generic bound Kerr-Newman: UNVERIFIED. Fortran. pip: NO (UNVERIFIED). Observables: none noted.

### [9] Gralla & Lupsasca 2020 (null)
- "The Null Geodesics of the Kerr Exterior": arXiv:1910.12881; Phys. Rev. D 101, 044032 (2020) (NOTE: article no. is 044032, not 044031 as user wrote).
- PRD 101 044031 is most likely "Lensing by Kerr black holes" (arXiv:1910.12873) per one reference-list snippet (UNVERIFIED pairing).
- Timelike: NO (null only). Analytic: YES (elliptic integrals in Legendre form; Jacobi-function inversion; all null geodesics outside horizon). No code pip. Observables: lensing/images (null).

### [10] kgeo (null solver)
- Repo NOT found by search. arXiv:2507.15775 (2025, "Learning Null Geodesics for Gravitational Lensing Rendering in GR") describes kgeo as a rapid, precise elliptic-integral solver for light paths in Kerr, based on Gralla & Lupsasca; benchmarked vs GravLensX. => NULL ONLY. Timelike: NO. pip/lang/speed: UNVERIFIED (not in snippet). Need direct repo link (not found).
- Other candidates seen: github.com/Yachim/KerrGeoPy (fork; timelike+lightlike; bundles shadow renderer). einsteinpy-geodesics repo (python wrapper of Julia Kerr solver per snippet -- verify).
- LEAD: github.com/openjournals/joss-reviews/issues/6587 (KerrGeoPy JOSS review thread; may contain benchmark discussion) -- fetch via WebFetch.

### [1e] KerrGeoPy JOSS review thread (github.com/openjournals/joss-reviews/issues/6587)
- No performance numbers, no accuracy comparisons to other codes, no observables, no limitations raised (review unfinished; checklist only). Version under review v0.9.3; Zenodo archive DOI 10.5281/zenodo.11386563.
- => KerrGeoPy published speed numbers: NONE found in README, JOSS thread, or snippets. (arXiv full text blocked.)

### [11] PyGRO (Della Monica & ...)
- arXiv:2504.20152 "PyGRO: a Python Integrator for General Relativistic Orbits" (A&A, accepted 2025; aa54300-25). PyPI: pygro; docs rdellamonica.github.io style (deepwiki rdellamonica/pygro).
- Timelike AND null paths; numerical integrator; general analytically-given metric. Lang: Python (C? unverified).
- Kerr support: UNCONFIRMED (only Schwarzschild demo seen in snippets). Generic inclined eccentric Kerr: UNCONFIRMED. Speed: none found. Observables: not stated in snippet (no astrometry/LTT seen).
- Author list: user wrote "Della Monica & de Martino"; snippet confirms Della Monica only, "de Martino" NOT confirmed (UNVERIFIED).
- pip: YES (PyPI pygro, per snippet "code is on PyPI"; version unverified).

### [12] kerrgeodesic_gw (BHPT; Gourgoulhon, Le Tiec, Vincent, Warburton)
- github.com/BlackHolePerturbationToolkit/kerrgeodesic_gw ; SageMath package for GW emission from bodies orbiting Kerr BH. License GPL v2+.
- Lang: Python on SageMath (SageManifolds). pip: `sage -pip install kerrgeodesic_gw` (PyPI kerrgeodesic-gw; 0.5 latest seen; manual says 0.4.1). Needs SageMath.
- Kerr geodesics: class KerrGeodesic -> SageMath's general-purpose geodesic integrator (= numerical, SLOW by design). Timelike AND null notebooks. GW from circular orbits notebook. Generic inclined eccentric bound: not confirmed in snippets.
- Observables: GW (not astrometry/redshift/LTT). Galactic Center GW application cited (2019 A&A by same authors; per slide deck; UNVERIFIED).
- Speed: none. Not MCMC-oriented.

### [13] EinsteinPy (geodesics module)
- einsteinpy-project.readthedocs.io (stable/latest api/geodesic). pip: `einsteinpy` (YES, per docs). Lang: Python (autodiff metric derivatives).
- Timelike + Null classes; metric options Schwarzschild, Kerr, KerrNewman; Kerr via metric_params=(a,). Numerical: FANTASY integrator (Christian & Chan 2021, ApJ, per docs) with orders 2,4,6,8; default steps 50 vs 100 doc conflict; rtol/atol 1e-2 default.
- Generic inclined eccentric bound Kerr: YES in principle (general initial 4-pos/4-mom numerical integration) -- UNVERIFIED for orbit classes but architecture is generic.
- Speed: NONE found. Observables: none (no astrometry/redshift/LTT; only trajectories).
- Paper: academia.edu copies only (unofficial); EinsteinPy paper citation UNVERIFIED.

### [14] REBOUNDx gr / gr_full / gr_potential (Tamayo et al. 2020; arXiv:1908.05634, App. B)
- Lang: C library + Python (pip `reboundx`; PyPI 2.2.1 seen). Numerical N-body + 1PN forces.
- Timelike: massive bodies only (no null). Analytic: NO (numerical symplectic/N-body integration with PN accelerations); gr_potential integrates analytically per snippet (position-only potential, apsidal precession).
- Kerr / spin / Lense-Thirring: NONE in snippets -> effectively Schwarzschild-PN (1PN) only. Generic inclined eccentric Kerr: NO (not Kerr). Inclination via N-body orbits possible but no Kerr metric.
- Observables: orbital elements only; no astrometry/redshift/LTT. Speed: none published in snippets.
- => Not a competitor for exact Kerr.

### [15] orbitize! (Blunt et al.; arXiv:1910.01756 per ar5iv link; journal unverified, Caltech repo dates paper Mar 2020)
- Python, pip `orbitize` (libraries.io: 3.1.0 seen). Methods: OFTI + parallel-tempered MCMC, Keplerian solver; optional n-body backend (use_rebound).
- GR: NO (no Kerr/PN). Observables: relative astrometry, RV, Hipparcos IAD, Gaia astrometry. No redshift/LTT modelling.
- Speed: "under ten lines" (usability claim; no speed numbers). Citations per docs: Blunt+2019, Blunt+2017 (OFTI), 2024 paper (unverified contents).
- => Not relevant for GR; included as usability/benchmark context for orbit-fit packages.

### [16] RAPTOR (Bronzwaer, Davelaar, Younsi et al.)
- RAPTOR I: A&A 613, A2 (2018), DOI 10.1051/0004-6361/201732149; time-dependent radiative transfer in arbitrary spacetimes. RAPTOR II: polarized, A&A 641, A126 (2020) (per snippet). Lang: C++ (UNVERIFIED). pip: NO.
- Timelike: NO (light paths via geodesic eqs in general metric -> null). Kerr: general metric supported (analytic or numerical metrics per 2018 paper snippet). Speed: NOT in snippets. Observables: images/spectra/polarization of plasma (null), not stellar astrometry.

### [17] ipole (Moscibrodzka & Gammie 2018, MNRAS 475, 43; AFD-Illinois/ipole)
- Lang: C (per snippet). pip: NO. Null only: YES-by-purpose (polarized covariant ray tracing for EHT); timelike: NO. Kerr support: not explicitly confirmed in snippets (repo check pending).
- Observables: images/polarization (null). Speed: not in snippets. Semi-analytic radiative transfer scheme (paper title).

### [18] GC/S2 orbit-fitting public codes (found)
- pmplewa/GC-OrbitFit (github.com/pmplewa/GC-OrbitFit): Python S2 orbit fit; 1PN post-Newtonian corrections + gravitational redshift + time-correlated astrometric noise; astrometry+RV (NACO/OSIRIS/SINFONI data). pip from GitHub (per snippet). Sampler: emcee -> dynesty (nested sampling). NO full Kerr geodesic (1PN approx). Speed: author says readability > performance. Kerr: NO (unverified repo).
- cPhysPlus/gr-orbit-toolkit (github.com/cPhysPlus/gr-orbit-toolkit; arXiv:2511.19442 "gr-Orbit-Toolkit: a Python-based software for simulating and visualizing relativistic orbits", Yachay Tech). Simulation (not fitting?). Kerr/timelike/speed: UNVERIFIED -> fetch repo README.
- PyGRO (A&A 2025, aa54300-25; arXiv:2504.20152): forward geodesic integration with user-supplied symbolic metric; S2 example Schwarzschild; relation to fitting: forward modeling. aanda.org BLOCKED.
- ExoSOFT (Eastman? IOP AJ 2017, doi 10.3847/1538-3881/aa5e4a): astrometry+RV Keplerian fitter (Cython core). NO GR.
- GRAVITY 2020 (via Wikipedia snippet): S2 consistent with Schwarzschild geodesics (GR precession/redshift detection). Not a code.

### [18b] fetched repo READMEs (GitHub is reachable)
- cPhysPlus/gr-orbit-toolkit: Python; `git clone` + `pip install -e .` (no PyPI page seen); MIT; authors Milagros Delgado & Wladimir Banda-Barragan (Springer chapter doi 10.1007/978-3-032-08366-1_20 per About sidebar). SIMULATES/VISUALIZES only (CSV histories, GIF). Integrator RK3 example; "relativity=False" flag; spacetime NOT named; Kerr: NOT stated. No fitting, no observables (astrometry/redshift/LTT), no speed numbers. => NOT a Kerr competitor.
- pmplewa/GC-OrbitFit: Python; sample S2 orbit fit; 1st-order PN ("takes into account first-order post-Newtonian corrections"), gravitational redshift effect, time-correlated astrometric noise; orbit integration via REBOUND/REBOUNDx; sampler dynesty (latest) / emcee (earlier). Data: astrometry_NACO/SHARP.csv, velocity_NACO/OSIRIS/SINFONI.csv (Gillessen+2017 Tab. 5 per README; doi 10.3847/1538-4357/aa5c41). Install: `pip install git+https://github.com/pmplewa/GC-OrbitFit.git`. Kerr/spin/LT: NOT stated. Speed: "readability ... not performance" (no numbers). MIT per sidebar. => Schwarzschild-PN fitter, NOT exact Kerr.
- kerrgeopy.readthedocs.io: DNS FAIL (not reachable from WebFetch). inspirehep.net: DNS FAIL (WebFetch) -> PyGRO paper unreadable via that link.

### [11b] PyGRO README (github.com/rdellamonica/pygro) -- ANSWERS
- AUTHOR: Riccardo Della Monica ONLY (README BibTeX: "author={Riccardo Della Monica}", eprint 2504.20152, year 2025). "de Martino" is NOT listed => user's attribution "Della Monica & de Martino" NOT supported by README.
- Lang: Python ("PyGRO is a Python library"). pip: `pip install pygro`; version badge decodes to 1.0.2 (UNVERIFIED vs PyPI). License GPL-3.0.
- Scope: "numerical integration of the geodesic equations describing a particle or photon orbit" -> timelike AND null; "in any metric theory of gravity, given an analytic expression of the metric tensor" (general metric). Kerr: NOT named/demonstrated in README. Generic bound orbits: not stated. Integrator: not named (RK/symplectic unknown). Speed: none. Observables: none (geodesic only -- no redshift/astrometry/LTT).
- => Generic numerical integrator framework (not Kerr-specific, no MCMC-fast claims).

### [19] Analytic-vs-numeric speed evidence (search summaries; underlying pages mostly iopscience/aanda = BLOCKED, so snippet-only)
- Dexter & Agol 2009 (ApJ 696, 1616; iopscience page not fetchable): summary says an earlier Agol (1997) analytic code ran ~100x faster than direct integration for photon paths from far away to a thin disk (problem-structure dependent). SNIPPET-ONLY, null geodesics.
- Cieslik, "Kerr Geodesics in Terms of Weierstrass Elliptic Functions", arXiv:2305.07771 (2023) -- analytic Weierstrass form; speed numbers NOT seen in snippet (candidate for part b; fetch if possible).
- Symplectic integrator for Kerr: "Construction of explicit symplectic integrators in general relativity. IV. Kerr black holes", arXiv:2106.12356 (authors not in snippet) -- numerical alternative; speed unverified.
- YNOGK paper has a speed-test section (arXiv:1305.1250); timing numbers NOT in excerpts (unverified). Dolence et al. (2009) cited for direct integration being simpler/quicker for radiative transfer workloads.
- KerrGeoPy paper (arXiv:2406.01413, snippet): claims analytic solutions are more numerically stable than numerical integration over long spans and can be evaluated at any chosen time (no numbers).
- One bound-orbit study (snippet): closed-form fundamental frequencies match direct numerical integration to < 1e-15 relative error in double precision (unnamed).
- Summary: NO modern matched-accuracy wall-clock benchmark found in these results. Most comparative speed claims 2009-2013.

### [20] Plunges / analytic generic plunging (Dyson & van de Meent)
- arXiv:2302.03704 "Kerr-fully Diving into the Abyss: Analytic Solutions to Plunging Geodesics in Kerr" -- Conor Dyson, Maarten van de Meent. v3 30 Oct 2024. CQG DOI 10.1088/1361-6382/acf552 (ADS: Oct 2023; ResearchGate: 8 Sep 2023).
- Closed-form generic plunging geodesics, in Mino time, with elementary + Jacobi elliptic functions. Plunges from innermost precessing stable circular orbit (ISSO): almost entirely elementary, inputs = spin + orbit radius. New inflow equation from ISSO to horizon accounting for inclination. Packaged in BHPT KerrGeodesics (per snippet) and used by KerrGeoPy.
- => "van de Meent 2020" in user list may be a mis-cite; closest van de Meent work found: Dyson & van de Meent 2023 (plunges); van de Meent & Warburton 2018 (isofrequency pairing); van de Meent & Shah 2015 (metric perturbations, PRD). UNVERIFIED.
- Other candidates for limits: arXiv:2309.13832 "Inspiral and Plunging Orbits in Kerr-Newman Spacetimes" (snippet only). arXiv:2606.03558 "A Plunge into the Chasm: Surviving Tidal Effects in Kerr Spacetime" (2026; snippet only, not read).

### [21] FEW generic-Kerr status (snippet-level; arXiv blocked)
- Snippet claim (UNVERIFIED): FEW 2021 paper's relativistic models handled only eccentric Schwarzschild; generic Kerr used an upgraded AAK approximation. => FEW generic Kerr geodesic: NOT exact in 2021 paper.
- Recent: arXiv:2609.37039 "Beyond circular and equatorial orbits: Bayesian constraint on vector charge with the adiabatic inspirals of generic Kerr geodesics" (Sept 2026 per id; snippet) -- orbit-averaged EM Carter-constant flux for eccentric inclined nonresonant geodesics fed into FEW.
- arXiv:2102.02713 "Adiabatic waveforms for EMRIs via multivoice decomposition" (2021; inclined eccentric Kerr inspiral, snippet).
- PRL 128, 231101 (2022) "Adiabatic Waveforms from EMRIs: An Analytical Approach" (generic-orbit model via analytical self-force; snippet).
- arXiv:2411.09700 (Nov 2024): PN fluxes for inclined spherical Kerr orbits (snippet).
- 2026 EPJC "A deep learning framework for amplitude generation of generic EMRIs" (snippet): Teukolsky amplitudes ~1e5 modes over 4 params = bottleneck.

### [5b] GYOTO README (github.com/gyoto/Gyoto) -- ANSWERS
- Lang: C++ core; Yorick is the plug-in/scripting target named in README; Python bindings NOT stated (python/numpy.i present in repo). pip: NOT stated.
- Purpose: "framework for computing orbits and ray-traced images in General relativity". Massive orbits / Kerr timelike: NOT stated in README (ASCL snippet claims Kerr null+timelike; unverified). Generic inclined: NOT stated.
- Observables: ray-traced images (FITS output). Spectra/LTT/redshift/astrometry: NOT stated. Speed: NONE. License GPL-3.0. Citation: Vincent et al. 2011, CQG 28 225011, arXiv:1109.4769; GYOTO 2.0 paper (2024, submitted); Zenodo 2024.
- => Not a pip, not a fast-MCMC tool; general-purpose numerical framework (C++/Yorick).

### [7b] van de Meent 2020 -- CANDIDATE IDENTIFIED (snippet-level)
- M. van de Meent, "Analytic solutions for parallel transport along generic bound geodesics in Kerr spacetime", Class. Quantum Grav. 37, 145007 (2020), DOI 10.1088/1361-6382/ab79d5, arXiv:1906.05090 [gr-qc] (posted 2019). Closed-form (Mino-time) solutions for parallel transport along generic bound Kerr geodesics; spin-spin resonances in quasi-circular EMRIs. Carter/frequency results NOT confirmed in snippet.
- Note: search "Listings ... Classical and Quantum Gravity, not PRD" -> user's "van de Meent 2020" likely = this CQG paper (UNVERIFIED match). Not a timelike-observable code.
- Also: van de Meent earlier self-force generic bound Kerr paper (journal truncated "Phys.") -- unverified.

### [10b] kgeo (achael) -- candidate repo found
- github.com/achael/kgeo ; Andrew Chael. Small Python package for light paths (NULL) in Kerr; follows Gralla & Lupsasca 2019 framework; outputs polarized images + photon-path animations. Deps: SciPy (3rd-kind elliptic integrals), GSL fallback. Not on PyPI per snippet (UNVERIFIED); edit source-file settings to choose backend. Timelike: NO.

### [22] GRay2 (candidate competitor; NOT in user list)
- arXiv:1706.07062 "GRay2: A General Purpose Geodesic Integrator for Kerr Spacetimes" (ApJ 2019, DOI 10.3847/1538-4357/aadfe5 per iopscience link). Repo: github.com/luxsrc/gray (snippet). Integrates particle AND photon geodesics in Kerr. Numerical general-purpose integrator; motivated partly by fast fitting of Galactic Center star orbits (snippet). Lang/pip/speed: TO CHECK (fetch README).

### [23] Gyoto-S2 fitting and S2 orbit-fit landscape (snippet)
- No single public code combines Kerr geodesic + S2 astrometry + Romer delay (search). GC-OrbitFit = public S2 fitter but 1PN (Schwarzschild-type). Elliot Telford S2 orbit simulator (Python, GitHub; Sgr A* spin simulation mention; no Kerr fit confirmed). GRAVITY 2018 (arXiv:1810.12641) and graviton-mass S2 papers (arXiv:2404.09795) = science papers, not code.

### [22b] GRay2 README (github.com/luxsrc/gray) -- ANSWERS
- Lang: NOT stated (built via `lux` framework: `make`, `lux gray`; luxsrc/lux). pip/conda: NOT stated (tarball/GitHub). GPL-3.0 (sidebar).
- "hardware-accelerated geodesic integrator" (GPUs, Xeon Phi). Integrator scheme / tolerances: NOT stated. Speed numbers: NONE.
- Input: reads spacetime + fluid data from HDF5 (g_*, Gamma_* arrays) => NUMERICAL metric input (GRMHD ray-tracing oriented). Kerr / timelike / inclined: NOT stated in README (paper title arXiv:1706.07062 says "for Kerr spacetimes" per search).
- Observables: NOT stated. Citation: not in README.

### [10c] kgeo README (github.com/achael/kgeo) -- ANSWERS
- Lang: Python ("simple python analytic kerr raytracing"); deps SciPy (elliptic integrals of 3rd kind; GSL slower alternative via editing globals). pip: NOT stated (pyproject.toml present; no install command in README). Python 3.8/3.11/3.12.
- NULL ONLY (README: "raytracing null geodesics in the Kerr metric"). Formalism Gralla & Lupsasca 2019 (arXiv:1910.12881). Speed: only "GSL much slower". Outputs: images/polarization/animations NOT explicit; no redshift/time-delay/astrometry stated.
- License GPLv3; Zenodo DOI 10.5281/zenodo.8092968 (badge). CITATION.cff present.

### [24] Gyoto-based / full-GR Kerr S2 modelling papers (snippet-level; codes not confirmed public)
- Grould, Vincent, Paumard & Perrin 2017, A&A 608, A60, "General relativistic effects on the orbit of the S2 star with GRAVITY", arXiv:1709.04492, DOI 10.1051/0004-6361/201731148. Uses GYOTO to produce synthetic S2 astrometry + RVs in Kerr with Lense-Thirring, Shapiro delay, grav. redshift. Simplified lensing (main image only). => Observables + Kerr in one pipeline BUT GYOTO-based, not pip, code not stated public.
- Zhang et al. 2015 (arXiv:1508.06293, "On Testing the Kerr Metric ... via Stellar Orbital Motion: Full General Relativistic Treatment"): full-GR fits of apparent positions + redshifts for mock S0-2/S2; params mass, distance, spin, inclination, eps. Code availability: UNVERIFIED.
- GRAVITY Collab. 2020, "Detection of the Schwarzschild precession in the orbit of the star S2" (arXiv:2004.07187; A&A, aa37813-20): PN model benchmarked vs full-GR GYOTO orbits; precession factor ~1.10 +/- 0.19 (snippet); ~5-6 sigma (snippet). Volume/page not verified.
- => (a) answer so far: no pip-installable package found that does exact-Kerr TIMELIKE + astrometry/RV/redshift + fast. Closest pieces: KerrGeoPy (pip, exact Kerr, no observables), GC-OrbitFit (pip, observables, PN not Kerr), Gyoto/GRay2 (numerical Kerr, not pip, speed unknown), Zhang+2015 & Grould+2017 (papers, Kerr observables, code unknown).

### [25] Speed evidence (snippets; GW/flux models are NOT geodesic evaluations -> not comparable to astrometric forward model)
- FEW 2021 (arXiv:2104.04582, PRD 104 064047): benchmark 2-yr waveform at (p,e)=(10,0.7) on NVIDIA V100; angular-harmonic evaluation < 1 ms; total ms in Fig. 8 NOT in snippet. Mode summation dominant cost.
- FEW Kerr eccentric equatorial, arXiv:2506.09470 (2025, "The Fast and the Frame-Dragging"): ~0.1 s per waveform with GPU (snippet; waveform duration/GPU not given).
- arXiv:2410.17310 (2024): circular equatorial Kerr, 4-yr inspiral: ~30-40 ms/mode on A100 vs ~1.8 s/CPU core; 63 modes ~60 ms GPU vs ~65 s CPU (snippet).
- KerrGeoPy: no timing (README/JOSS/paper-snippet). Paper snippet says: SciPy Legendre elliptic integrals; trajectory() returns 4 callables of Mino time; claims long-span robustness + evaluation at any instant. NO numbers.
- Dexter & Agol 2009 (ApJ 696 1616): title per search result = "A Fast New Public Code for Computing Photon Orbits in a Kerr Spacetime" (iopscience title in search result; page not fetchable). ~100x vs direct integration (Agol 1997 analytic code, photon paths, snippet).
- GRay2 paper (arXiv:1706.07062): benchmarks null-geodesic integration on several CPUs/GPUs (snippet; numbers not retrieved).
- Stein & Warburton (arXiv:1912.07609): new separatrix algorithm ~45x faster (snippet).
- Cieslik (arXiv:2305.07771) Weierstrass closed form: speed NOT retrieved.
- ynogkm (arXiv:1311.4436): timelike Kerr-Newman closed-form code (A&A 2014 per search); timing NOT retrieved.

### [26] Blocked/unreachable hosts this session (record only; no circumvention attempted)
- Proxy 403 policy: arxiv.org, export.arxiv.org, www.aanda.org, iopscience.iop.org, b2find.eudat.eu.
- DNS fail (WebFetch): joss.theoj.org, yachim.cz, kerrgeopy.readthedocs.io, bhptoolkit.org, inspirehep.net.
- Reachable: github.com (READMEs), pypi (not fetched), search snippets.

### [27] DRAFT SYNTHESIS (for final report)
(a) No pip-installable package found combining exact-Kerr TIMELIKE geodesics + astrometry/RV/redshift + published fast speed. pip exact-Kerr codes (KerrGeoPy, EinsteinPy, PyGRO) output geodesic quantities only. Observables codes are PN (GC-OrbitFit) or Keplerian (orbitize!) or non-pip GYOTO/GRay2 (speed unknown). Full-GR Kerr S2 modelling exists in papers (Grould+2017 GYOTO; Zhang+2015) with code not found public.
(b) No retrievable published wall-clock for closed-form Kerr timelike evals. Hints: ~45x (separatrix), ~100x (photon, 1997/2009 analytic vs direct). Limits: KerrGeoPy excludes polar, marginal, a=1; plunges = Dyson&vdM 2023; Fujita&Hikida bound-only; t(lambda) gives coordinate time vs Mino time; LTT inversion = root finding, not in any surveyed tool (own assessment); non-geodesic forces break closed form (EOB-deformed analytic exists only for specific deformation: arXiv:2511.06684 snippet).
(c) Numerics win for non-Kerr/numerical metrics (GRay2 HDF5; Gyoto numerical metrics; PyGRO general metric), added non-geodesic/time-dependent forces, and near degenerate closed-form cases. Closed form wins for many evaluations on fixed (E,L,Q) orbit (KerrGeoPy claim, unverified numerically). Missing datum: same-machine benchmark KerrGeoPy vs RelatiPy at 200 S2 epochs (NOT run).
- Brief corrections: Katz 5th author = Scott A. Hughes; Gralla&Lupsasca = PRD 101 044032 (not 044031); PyGRO = Riccardo Della Monica (single author per README BibTeX); "van de Meent 2020" likely = CQG 37 145007 / arXiv:1906.05090 (parallel transport on generic bound Kerr geodesics; unverified match); KerrGeoPy README: "python implementation of the KerrGeodesics Mathematica library".
- Word budget for final report: <= 700 words.

### [28] STATUS: final report being handed back (<=700 words). Corrections folded in: Agol-1997 100x is the 1997 analytic photon code (not geokerr itself); Zhang+2015 and Grould+2017 are the Kerr S2 modelling papers (only Grould is GYOTO-based per snippet); PyGRO sole author Della Monica; Gralla-Lupsasca PRD 101 044032.
