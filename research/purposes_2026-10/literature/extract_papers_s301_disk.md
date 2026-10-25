# Paper extraction (alphaXiv)

Rules: numbers and formulas copied verbatim with page numbers; "not found" where the paper does not state the item.

Page convention: "p.N" = PDF page index from the tool's <page num="N"> tags (matches the printed page numbers).
Formula convention: copied from the PDF text layer; where the layer scrambles fraction/bracket layout, the reconstructed form is marked "(layout reconstructed)" and needs a check against the PDF.

==============================================================================
## PAPER 1 - arXiv 2607.12664, "Discovery of a star sensitive to the spin of Sgr A*" (GRAVITY Collaboration et al. 2026, Nature)
Source: answer_pdf_queries, first pass. Pages returned: 1, 3-9, 16-18, 20-23, 26.

### Q1. Orbital elements of S301 (Extended Data Table 2, p.20)
The table has two value/error pairs with no A/B labels. Mapping used: A = first pair, B = second pair. Paper 4 (Table 3) adopts the first pair as its truth, which supports this mapping.

| Parameter | A value | A formal error | B value | B formal error |
|---|---|---|---|---|
| a [mas] | 83.0 | 0.7 | 83.0 | 0.7 |
| e | 0.9832 | 0.0010 | 0.9824 | 0.0011 |
| i [deg] | 124.09 | 1.10 | 122.84 | 1.12 |
| Omega [deg] | 73.8 | 3.5 | 256.9 | 3.3 |
| omega [deg] | 293.4 | 2.2 | 115.1 | 2.0 |
| t_P [yr] (time of periastron) | 2023.126 | 0.010 | 2023.125 | 0.011 |
| P [yr] | 8.68 | 0.11 | 8.68 | 0.11 |
| chi^2 | 25.86 | - | 25.53 | - |

- Internal discrepancy: p.6 text says "e = 0.9832/0.9821 (for the two possible orientations)". Table p.20 gives 0.9824 +/- 0.0011 for B. Both reported verbatim.
- p.6 text: "a semi-major axis of a = 83 mas (33% smaller than S2's) and an orbital period of 8.7 years"; "pericenter distance r_p of only 136/142 R_S".
- Mass and distance: p.6 "The potential is defined by the central mass fixed with M_MBH = 4.297 x 10^6 M_sun at a distance of R_0 = 8277 pc." p.23 "a GC distance of R_0 = 8.3 kpc [24]". p.21 "with mass and distance of Sgr A* from [5]".
- Reference-frame conventions: p.20 '(angle conventions follow [4])' = Gillessen et al. 2009. p.6: no coordinate-system offsets fitted (data referenced directly to the NIR counterpart of Sgr A*). p.6, p.20: the two orientations are equally valid; the Romer delay could break the sign ambiguity but does not yet.
- Explicit definitions of i, Omega, omega: none beyond 'angle conventions follow [4]' (p.20). Ref [4] = Gillessen, S. et al. 2009, ApJ 692, 1075 (full-text reference list). The ascending node is not defined in the full text.

### Q2. Data used (p.4-5, p.16, p.20)
- Astrometric epochs: 19. p.5 "Overall, we have 19 astrometric positions of S301 that outline an ellipse on the sky"; p.16 "comprise 19 datasets spanning more than eight years". Dedicated pointings: "eight and five additional astrometric measurements" in 2024 and 2025 (p.4). 2017 and 2021 positions are back-predicted detections, "less precise, as only few exposures ... were taken" (p.5).
- Instrument: GRAVITY only (VLTI). p.4 "we have regularly observed the central arcsecond around Sgr A* with GRAVITY"; p.16 Extended Data Table 1 "Summary of GRAVITY data that yielded an inference of S301". ERIS IFU data searched, S301 not detected (p.5).
- Astrometric precision per epoch: 207 uas per coordinate (RA and Dec), the RMS pixel-grid discretisation error for 0.8 mas pixels; the paper says it 'acts effectively as a noise floor' (Methods, p.19; full text). Forecast for completed GRAVITY+: 100 uas (p.22). p.4 gives ~30 uas as GRAVITY's general precision, not S301-specific.
- Radial velocities: NONE. p.5 "We thus do not yet have any radial velocity information."
- Time span: epoch 2017.2389 (2017-03-26..30) to 2025.5978 (2025-08-06) (p.16). Methods p.16: "more than eight years".
- Extended Data Table 1 (p.16), 19 GRAVITY datasets (date | epoch | N_exp | t_int [h] | p [mas] | program):
  2017-03-26..30 | 2017.2389 | 12 (one pol.) | 1.0 | (-57.974, 35.939) | 60.A-9102(A)
  2021-07-29 | 2021.5759 | 8 | 0.71 | (-45., 45.) | 105.20B2.004
  2023-04-05 | 2023.2604 | 12 | 1.2 | (0.0, 0.0) | 111.24H1.001
  2023-06-01 | 2023.4162 | 13 | 1.3 | (0.0, 0.0) | 111.24H1.002
  2023-06-06 | 2023.4298 | 9 | 0.9 | (0.0, 0.0) | 111.24H1.002
  2023-07-04 | 2023.5063 | 16 | 1.6 | (0.0, 0.0) | 111.24H1.003
  2024-03-25 | 2024.2325 | 3 | 0.3 | (-51.13, 54.59) | 112.25CV.001
  2024-03-27 | 2024.2380 | 2 | 0.2 | (-51.13, 54.59) | 112.25CV.001
  2024-03-28 | 2024.2405 | 6 | 0.6 | (-51.13, 54.59) | 112.25CV.001
  2024-04-29 | 2024.3281 | 4 | 0.4 | (-56.13, 59.59) | 113.268P.001
  2024-05-28 | 2024.4073 | 12 | 1.2 | (-63.08, 63.25) | 113.268P.002
  2024-06-22 | 2024.4756 | 10 | 1.0 | (-63.08, 63.25) | 113.268P.003
  2024-06-25 | 2024.4837 | 9 | 0.9 | (-63.08, 63.25) | 113.268P.003
  2024-07-17 | 2024.5438 | 9 | 0.9 | (-63.08, 63.25) | 113.268P.004
  2025-03-16 | 2025.2070 | 10 | 1.0 | (-61., 51.2) | 114.270V.001
  2025-05-11 | 2025.3602 | 10 | 1.0 | (-95.26, 51.2) | 115.27WT.001
  2025-06-09 | 2025.4393 | 10 | 1.0 | (-95.26, 51.2) | 115.27WT.002
  2025-07-07 | 2025.5159 | 9 | 0.9 | (-69.959, 56.581) | 115.27WT.003
  2025-08-06 | 2025.5978 | 7 | 0.7 | (-72.005, 57.737) | 115.27WT.004

### Q3. Spin forecast (p.8, p.9, p.21-23, p.26)
- Dynamical model: spin enters through Lense-Thirring (LT) precession from the Kerr metric. p.3 "The Kerr metric for a rotating black hole leads to corrections of order beta^3." p.4 "At the beta^3 order, the leading-order term from the Kerr metric contributes". Orbit-averaged LT (Eq. 1, p.21, layout reconstructed): Delta_varpi_LT = -8 pi chi cos(xi) [G M_MBH / (a(1-e^2) c^2)]^(3/2); Delta_Theta_LT = -4 pi chi sin(xi) sin(lambda) [G M_MBH / (a(1-e^2) c^2)]^(3/2). xi = angle between spin axis and orbital angular momentum; lambda = position angle of the spin-axis projection onto the orbital plane (p.22).
- Base model: 1PN Schwarzschild with Romer delay (p.6, p.21). Relativistic Doppler and gravitational redshift "do not matter as they only affect radial velocity measurements" (p.21).
- PN order caveat: p.23 "future robust spin constraints will require modeling up to 2PN order to avoid systematic biases at low spins."
- Assumed future precisions (p.22): astrometric "100 uas (as expected for the final performance of GRAVITY+)"; radial velocity "1 km/s accuracy on the radial velocity as reachable with future ELT/MICADO observations".
- Cadence (p.22-23): "A sampling of 10 data points per year until 2035, with 10 additional data points around pericenter". Simulated observations "between 2026 and 2035" (p.22).
- Timeline: p.8 "to a degree accessible within around a decade, using existing and/or planned instrumentation"; "there is a fair chance to directly measure Sgr A*'s spin from the S301 orbit (Fig. 3 right, Methods section, Extended Data Fig. 9)" within a decade.
- Forecast spin precision: p.8 Fig. 3 caption "These mock data constrain chi with an uncertainty < 0.2, and also the orientation to around +/-30 deg." p.23 "a spin constraint ... chi of < 0.2 (Fig. 3, right), equivalent to a > 5 sigma detection of chi = 1 compared to the non-spinning (chi = 0) case." Fiducial: chi = 1, orientation approx. aligned with stellar orbital angular momentum (p.22).
- Other spin numbers: LT in-plane 0.11 deg chi cos(xi) per revolution (p.8); spin-precession time scale ~58000 yr (p.8).
- Extended mass / disk / perturbers:
  - p.9 "S301's Schwarzschild precession offers the opportunity to tighten constraints on possible deviations from general relativity ... and on a putative extended mass distribution."
  - p.26 "Recent constraints from [30] limit the extended mass within the central 10 mpc to <~ 1200 M_sun."
  - p.26 N-body test: 60 stellar-mass BHs of 20 M_sun within 10 mpc, rho ~ r^-2: "average orbital-plane precession per orbital period of 0.65 (+0.36/-0.56) arcminutes. This contribution is typically subdominant compared to the Lense-Thirring precession for moderate-to-high dimensionless spins and favorable orientations".
  - p.26 disk case: "In an extreme case in which all the mass (allowed by S2 observations) is concentrated within the orbit of S301 in a disk and in which the orientation of the orbit is a few degrees from the disk, the nodal precession due to the disk might be larger, by up to an order of magnitude, than the Lense-Thirring precession. However, for most orientations and for less extreme mass distributions, the Newtonian effect of this mass on the nodal precession is expected to be one or two orders of magnitude smaller than the Lense-Thirring effect."
  - p.26 phase separation: "the Lense-Thirring effect is concentrated near pericenter, whereas perturbations from a granular stellar background tend to produce their largest observable deviations near apocenter."

### Q4. Magnitude, spectral type, pericentre velocity
- K magnitude: m_K = 19.3 +/- 0.3 (p.6; also p.23).
- Spectral type: p.6 "the magnitude matches a main-sequence star of spectral type F1.5, corresponding to a mass <~ 1.5 M_sun". p.23 "a spectral type of F1.5, corresponding to a mass of just below 1.5 M_sun"; also p.23 "Its spectral type then is a late A-type or an early F-type star". Radius ~1.4-1.6 R_sun (p.6-7). L = 5.5 L_sun, absolute V = 2.88 (p.23).
- Pericentre velocity: p.6 "At pericenter, S301 moves with ~ 25 600 / 25 000 km/s or 8.5% / 8.3% of c." Abstract (p.3): "the star's peak velocity reaches 25000 km/s".

### Extras noticed (not asked)
- p.9 "consistent with general relativity at the 1PN level ... resulting PPN bounds are weaker than those derived from S2 [8, 30] with f_SP = 0.94 +/- 0.88."
- p.6 relativistic pericentre advance per orbit 2.0 / 1.9 deg (two solutions); p.8 Schwarzschild advance 1.9 deg per revolution.
- p.6 Hills-mechanism origin, e ~ 0.98 (not needed here).

==============================================================================
## PAPER 2 - arXiv 2607.26713, "Effects of a disk structure on stellar motion at the Galactic Center" (Foschi, Vincent, Paumard, Perrin; A&A)
Source: answer_pdf_queries, first pass. Pages returned: 1, 2, 4, 6, 7, 8, 9, 10. NOT returned in the first pass: p.3 and p.5. Both are now covered by the full-text check (papers 1-3 block).

### Q1. Disk acceleration
- Frames (p.2): BH frame {x_BH, y_BH, z_BH}, harmonic coordinates (t, r, theta, phi), disk in the BH equatorial plane z_BH = 0. Total acceleration (p.2, Eq. 2): "a = -(G M_bh / r^3) r + a^disk + a^1PN".
- Surface density (p.4, Eq. 13): "Sigma(r) = B / r^gamma, where B is a constant." Exponent gamma ~ 2 for the CW stellar disk (p.4, citing Paumard et al. 2006; Bartko et al. 2009); gamma between 1 and 2 considered for generic cases (p.4).
- Mass element (p.4): "dM = 2 pi Sigma(r') r' dr' = 2 pi B r'^(1-gamma) dr'".
- Normalisation to disk mass (p.4, Eq. 16): "M_tot = 2 pi Integral_{r_min}^{r_max} B r'^(1-gamma) dr' = 2 pi { B (r_min^(2-gamma) - r_max^(2-gamma)) / (gamma - 2) for gamma != 2 ; B log(r_max / r_min) for gamma = 2 }".
- Ring acceleration, in-plane (p.4, Eq. 14): same form as Eq. (9) on p.3 with M_ring -> dM and R -> r'. Layer text and closed form: see FULL-TEXT CHECK, papers 1-3 block.
- Ring acceleration, vertical (p.4, Eq. 15): "da^disk_z = - 2 G dM / (pi p q^2) E(k) z".
- Total disk acceleration in the BH frame (p.4, Eq. 17): "a^disk = Integral_{r_min}^{r_max} ( da^disk_d u_d + da^disk_z u_z )".
- Special functions: complete elliptic integrals K(k), E(k) in the PARAMETER convention k = 4Rd/p^2 (Appendix B, p.11, Eqs. B.1-B.2). p, q: Eq. (6), p.3. d = sqrt(x_BH^2 + y_BH^2) (cylindrical radius), z = z_BH. No softening.
- Asymptotics (p.4, Eq. 18): "a_d ~ a_z ~ G M_tot / r^2, if r >> r_max" (point-mass limit). Eq. 19 (r << r_min) is garbled in the layer; p.4 text says a_d, a_z ~ r there.
- Gaussian-frame components (p.4-5): "a^disk = R_disk n + S_disk lambda + W_disk z" (Eq. 20). With psi = omega_obs + f - beta, f = true anomaly, r = p_star / (1 + e cos f) (p.4), layout-reconstructed:
  - Eq. 21: S_disk = [ sin(psi) cos(psi) sin^2(i_BH) / ( cos^2(i_BH) cos^2(psi) + sin^2(psi) )^(1/2) ] * Integral da_d - sin(psi) sin(i_BH) * Integral da_z
  - Eq. 22: W_disk = - [ cos(psi) sin(i_BH) cos(i_BH) / ( cos^2(i_BH) cos^2(psi) + sin^2(psi) )^(1/2) ] * Integral da_d + cos(i_BH) * Integral da_z
  - Eq. 23: R_disk = ( cos^2(i_BH) cos^2(psi) + sin^2(psi) )^(1/2) * Integral da_d + cos(psi) sin(i_BH) * Integral da_z
  (Integrals run from r_min to r_max.)
- Softening: NOT FOUND. The disk is a zero-thickness plane. p.4 "Due to the zero-thickness of the disk, this change of sign ... produces a real discontinuity in S_disk. In a realistic scenario, where the disk has some finite thickness, this change of sign would generally be smoother."

### Q2. Disk parameters and precession per orbit
- Disk mass: M_tot = 10^3 M_sun for the generic-disk figures (p.4, Figs. 3-4, 5, 6, 8, 9, 10). M_tot = 10^4 M_sun for one orientation case (p.7, Fig. 6 caption). CW-disk case: M_extra = 10^3 M_sun (extra, dark) plus luminous M_L = 10^4 M_sun (p.9, Eq. 52 and Table 2 caption).
- Eq. 52 (p.9, layout reconstructed): Sigma = M_extra / (2 pi log(r^L_min / r_min) r^2) + M_L / (2 pi log(r^L_max / r^L_min) r^2). "an extra mass M_extra distributed between r_min and the minimum luminous radius, and a luminous mass M_L = 10^4 M_sun distributed between r^L_min and r^L_max ... with the same power law slope (Paumard et al. 2006)."
- Power-law index: gamma = 1 and gamma = 2 (blue and orange bands in Figs. 5-10). CW disk: gamma = 2 (p.10 "its surface density ... follows a power law with gamma = 2").
- Radial extent: r_max = 0.5 pc for both S2 and S301 (Table 2 caption, p.10). Figs. 3-4 use r_min = 1e-6 pc, r_max = 1 pc (captions, p.4). Table 2 r_min options: 1e-6 pc; r_peri (per star); r^L_min (M_extra = 0). Luminous CW disk (Eq. 51, p.8): r^L_min ~ 0.04 pc, r^L_max ~ 0.2-0.5 pc.
- Plummer comparison (p.4): r0 = 0.012 pc, rho0 = 1.69 x 10^-10 kg/m^3 (values used in Heisel et al. 2022).
- Observer orbit of S2 (p.6, Eqs. 39-41): e = 0.8844, p_star = 0.0272 as; M_bh = 4.29 x 10^6 M_sun, R0 = 8275.9 pc; i_obs = 134.67 deg, omega_obs = 246.28 deg, Omega_obs = 48.21 deg. Eq. 42 (p.6): omega_obs = omega_GRAV + pi, Omega_obs = Omega_GRAV + pi.
- Relative CW-disk orientation: S2: i_BH^CW ~ 86 deg, beta_CW ~ 232 deg (Eq. 53, p.9). S301: i_BH^CW,S301 ~ 21 deg (107 deg), beta_CW,S301 ~ 342 deg (252 deg) (Eq. 54, p.9-10). For S301 the paper takes i_BH = 21 deg and beta = 342 deg (p.10) because they maximise the disk effect.
- Disk orientation (Eq. 50, p.8): i_disk = 127 deg +/- 2 deg, Omega_disk = 279 deg +/- 2 deg. Omega_disk differs by pi from the observed Omega of Paumard et al. (2006), as in Eq. (42).
- Table 2 (p.10), "Maxima shifts per orbit in the observables Delta_varpi and Delta_Theta, using the CW disk configuration and different inner radii, while r_max = 0.5 pc for both S2 and S301. The extra disk mass is M_extra = 10^3 M_sun. When r_min = r^L_min, M_extra = 0". Columns: |Delta X| S2 [deg] | |Delta X| S301 [deg] | r_min
  - Delta_varpi_disk | 0.007 | 0.004 | 10^-6 pc
  - Delta_Theta_disk | 5 x 10^-5 | 0.003 | 10^-6 pc
  - Delta_varpi_disk | 0.024 | 0.012 | r_peri^S2 (r_peri^S301)
  - Delta_Theta_disk | 0.0001 | 0.005 | r_peri^S2 (r_peri^S301)
  - Delta_varpi_disk | 0.0002 | 8 x 10^-6 | r^L_min
  - Delta_Theta_disk | 1 x 10^-6 | 2 x 10^-6 | r^L_min
  - Delta_varpi_LT (maximum LT, from paper 1) | 0.004 | 0.089 | -
  - Delta_Theta_LT (maximum LT) | 0.002 | 0.044 | -
  - Delta_varpi_Sch (in-plane Schwarzschild) | 0.200 | 1.684 | -
  - Delta_varpi_Plummer (S2 only) | 0.044 | - | -
- Text, p.10: "In this configuration, Delta_Theta_disk is roughly nine times smaller than the LT precession for S301 ... we are considering the maximum spin configuration (chi = 1, theta = pi/2), which means that it will likely be smaller". Conclusions p.10: "the LT out-of-plane precession in S301 is roughly nine times larger than that induced by the CW disk."
- Abstract (p.1): S2 "Delta p_star ~ 10 uas"; typical in-plane "Delta omega_obs ~ 0.2 deg"; typical out-of-plane "Delta Theta ~ 0.03 deg"; "For S301, the out-of-plane disk effect is a factor of ten smaller than the spin effect, assuming the most favourable spin configuration."
- Secular semi-latus-rectum change (p.5-6): spherical distributions and PN terms up to order 3 give zero secular change of p_star; the disk does give one (Eq. 48 structure).
- LT for S2: "For S2 the effect of the spin on the orbit is negligible" (p.9, Fig. 8 caption).

### Q3. Method and limitations
- Method: osculating (perturbed two-body) formalism, Gaussian frame {n, lambda, z} (p.4, Poisson & Will 2012). Secular change over one orbit: Delta mu_i = mu_i(f) - mu_i^0 with mu_i(f) = Integral dmu_i/dt dt (Eq. 33, p.6). Observer-independent angles (Eqs. 34-36, p.6): d varpi/dt = d omega_obs/dt + cos(i_obs) d Omega_obs/dt; d Theta/dt = sin(omega_obs) di_obs/dt - sin(i_obs) cos(omega_obs) d Omega_obs/dt; d Xi/dt = -cos(omega_obs) di_obs/dt - sin(i_obs) sin(omega_obs) d Omega_obs/dt.
- Newtonian vs PN: disk treated in a "pure Newtonian framework" (p.2, justified by S-stars being far from Sgr A*). Schwarzschild: first-order PN in harmonic coordinates, via Poisson & Will 2012 and Will 2014 (p.2). Lense-Thirring: 1.5PN equations of motion from Abd El Dayem et al. 2026 (p.2; p.6 Eqs. 37-38 analytic LT variations, layer garbled).
- Limitations and future work (quotes):
  - p.4, p.7: zero-thickness disk; a finite thickness would smooth the sign change of S_disk and the peaks in Delta p_star.
  - p.7: "This means that flattened and spherical distributions can behave differently on the orbital motion, and it is thus necessary to constrain the disk model independently, rather than relying on the current limits obtained for spherically symmetric models. Future work will investigate this point in more depth."
  - p.9: "A proper fit of the non luminous mass contained in the CW disk, or other disk-like structures at the GC, is left for future works."
  - p.10: "the value of M_extra has not been determined yet, and could potentially be larger"; disk mass "can potentially be larger than 10^3 M_sun".
  - p.10: CND (single ring at r_min, M = 10^4 M_sun) and Sgr B2 (point mass 10^7 M_sun at 100 pc): "a 0-th order analysis"; both about two orders of magnitude below the CW disk, "thus are both negligible".
  - Reproducibility: PERSEO Python package (Foschi 2026), GitHub (p.10).

==============================================================================
## PAPER 3 - arXiv 2409.12261, "Improving constraints on the extended mass distribution in the Galactic Center with stellar orbits" (GRAVITY Collaboration 2024, A&A 692 A242)
Source: answer_pdf_queries, first pass. Pages returned: 1, 2, 3, 5, 8, 10. NOT returned in the first pass: p.4, 6, 7, 9 (all now checked in the full text).

### Q1. Orbital elements (Table B.1, p.10; four-star fit that determines f_SP)
Columns in table: a [as], e, i [deg], Omega [deg], omega [deg], t_peri [yr], each with formal error. P (period) is NOT tabulated.

| Star | a [as] | e | i [deg] | Omega [deg] | omega [deg] | t_peri [yr] |
|---|---|---|---|---|---|---|
| S2 | 0.12502 +/- 3e-5 | 0.88444 +/- 6e-5 | 134.67 +/- 0.02 | 228.21 +/- 0.03 | 66.279 +/- 0.029 | 2018.3789 +/- 1e-4 |
| S29 | 0.39025 +/- 9.4e-4 | 0.96880 +/- 9e-5 | 144.24 +/- 0.09 | 4.9259 +/- 0.1590 | 203.68 +/- 0.17 | 2021.4102 +/- [prefactor lost in extraction] e-4 |
| S38 | 0.14249 +/- 4e-5 | 0.81807 +/- 2.2e-4 | 168.69 +/- 0.19 | 122.43 +/- 1.12 | 40.065 +/- 1.118 | 2022.6843 +/- 8e-4 |
| S55 | 0.10424 +/- 5e-5 | 0.72980 +/- 1.8e-4 | 159.59 +/- 0.17 | 319.43 +/- 0.97 | 327.77 +/- 0.93 | 2009.4738 +/- 7.6e-3 |

- Caption (p.10): "The orbital elements are meant in the sense of osculating orbit parameters, using a conversion time close to the respective apocenter times, i.e., 2010.35 for S2, 1977 for S29, 2000 for S38, and 2012 for S55. The reference frame parameters x0 and y0 refer to the epoch 2000."
- S29 t_peri error: the layer reads "2021.4102 x10^-4" with the leading digits missing. Check against the PDF.
- Same table, global parameters (p.10): M_bh = 4.2996 +/- 0.0118 (x10^6 M_sun); R0 = 8275.9 +/- 8.6 pc; x0 = -0.72370 +/- 0.08707 mas (prior -7.79e-4 +/- 7.38e-5); y0 = -0.11921 +/- 0.07391 mas; vx0 = 0.078241 +/- 0.005728 mas/yr; vy0 = 0.040542 +/- 0.005620 mas/yr; vz0 = -1.9983 +/- 1.3325 km/s; f_SP = 1.1350 +/- 0.1104.
- Cross-check (my arithmetic, not in paper): paper 2 Eq. 42 maps S2 to paper 2 Eq. 41: 228.21 + 180 = 408.21 = 48.21 mod 360 (Omega_obs = 48.21) and 66.279 + 180 = 246.28 (omega_obs = 246.28). Consistent.

### Q2. Equations of motion, extended mass, integrator, sampler, frame parameters
- Equations of motion (p.3): "we model their acceleration as given by a first-order PN approximation of GR for a massless test particle (Will 1993) and multiply the 1PN terms by a factor f_SP, such that f_SP = 0 corresponds to a Keplerian closed (non-precessing) orbit and f_SP = 1 to a GR Schwarzschild orbit".
- Spin (Lense-Thirring), not in the fit: p.2 "The Lense-Thirring effect, caused by the spin of the central SMBH, appears at 1.5PN order and gives both an additional contribution to the in-plane precession and a precession of the orbital plane." A_LT = 4 pi chi (R_S / (2 a (1-e^2)))^(3/2), 0.11 arcmin for S2 (p.2). delta_phi_Kerr = delta_phi_Schw - 2 A_LT cos(i); delta_Phi_Kerr = A_LT (p.2).
- Schwarzschild per orbit (p.1): delta_phi_Schw = 3 pi R_S / (a (1-e^2)); 12.1 arcmin per orbit for S2.
- Extended mass (p.5, Eqs. 1-4):
  - Power law: rho(r) = rho0 (r / r0)^s (Eq. 1), r0 = 40 mpc fixed; enclosed mass m(r) = 4 pi rho0 / (3 + s) * r^(3+s) / r0^s (Eq. 2, layout reconstructed). Slopes explored: -3 < s <= 0.
  - Plummer: rho(r) = rho0 [1 + (r/a)^2]^(-5/2) (Eq. 3); enclosed mass m(r) = (4/3) pi rho0 a^3 (r/a)^3 / [1 + (r/a)^2]^(3/2) (Eq. 4, layout reconstructed); total m_tot = (4/3) pi rho0 a^3. Scale radius a in (0, 40] mpc.
  - Parametrisation (p.5): power law rho0 = f_pow M_bh; Plummer rho0 = f_pl M_bh / ((4/3) pi a^3). Priors f_pow >= 0, f_pl >= 0. s and a are fixed, not fitted, in the main runs.
  - Enclosed mass constrained: M_encl,S2 = m(r_peri,S2 < r < r_apo,S2) (p.5), r_peri,S2 ~ 0.6 mpc, r_apo,S2 ~ 9.4 mpc.
- Integrator: NOT STATED (full text checked). p.3 defers to GRAVITY Collaboration 2018a, 2019, 2020.
- Sampler / fitter (p.3, p.5): Levenberg-Marquardt chi^2 minimisation or Metropolis-Hastings MCMC, "using 200 000 realizations".
- Fitted parameters (p.3): f_SP, M_bh, R0, five reference-frame coordinates (x0, y0, vx0, vy0, vz0), six orbital parameters per star (a, e, i, omega, Omega, t0). NACO zero points x0, y0, vx0, vy0 are partially degenerate with f_SP (p.3); seven additional S-stars used to get a prior on them (p.3, Appendix B).

### Q3. Astrometric and RV precision per instrument; epochs per star
- Per-instrument fit precision (NACO, GRAVITY, SINFONI, ERIS): NOT STATED (full text checked; no per-instrument uncertainty appears). Only: p.1 GRAVITY best-case 1 sigma = 30 uas. Forecast assumptions (not fit precisions), p.10: GRAVITY 50 uas, ERIS RV 10 km/s. SHARP astrometry also used for S2 (Fig. 2 caption, p.4).
- Spectral resolution (p.2-3, not precision): SINFONI R ~ 4000; ERIS (SPIFFIER) R ~ 10000.
- Forecast assumptions, NOT fit precisions (p.10, Appendix A): GRAVITY "10 ... data points per year with 50 uas astrometric accuracy"; ERIS "3 radial velocity data points per year with ERIS with 10 km/s accuracy".
- Data per star (p.2):
  - S2: 128 NACO + 82 GRAVITY astrometry; RV: 92 SINFONI, 3 Keck, 2 NACO, 4 GNIRS/GEMINI, 3 ERIS. Span 1992.2-2022.7.
  - S29: 66 NACO (2002.3-2019.7) + 29 GRAVITY (2019.6-2022.7). "We dropped a significant fraction of the available NACO data, which have been affected by confusion events." RV: 17 SINFONI, 2 GNIRS (2007.6-2021.4).
  - S38: 110 NACO (2004.2-2018.6) + 23 GRAVITY (2021.2-2022.7). RV: 8 SINFONI, 2 Keck, 1 ERIS (2008.3-2022.4).
  - S55: 42 NACO (2004.5-2013.6) + 27 GRAVITY (2021.2-2022.6). RV: 2 SINFONI (2014).
  - S24, S42: NACO + GRAVITY astrometry, SINFONI RV. S1, S4, S9, S21, S13: NACO + SINFONI only (p.2).
- Fig. 1 (p.3): 11 S-stars used; S2, S29, S38, S55 are the four most relevant.

### Q4. Final constraints
- Enclosed-mass limit: "upper limit of approximately 1200 M_sun with a 1 sigma confidence level" on mass enclosed within S2's orbit, "consistent with zero" (p.1 abstract; p.8 conclusions). Stellar-cusp prediction "approximately 1210 M_sun" (p.8). 3-sigma limit: stated as 99.7% (p.10 Appendix B), value not given in returned pages.
- f_SP: four-star fit (S2+S29+S38+S55) f_SP = 1.135 +/- 0.110 (MCMC 1.133 +/- 0.113) (p.3); Table B.1 1.1350 +/- 0.1104 (p.10); significance "from ~7 sigma ... to ~10 sigma" (p.5).
  - S2 with NACO prior: f_SP = 0.918 +/- 0.128 (MCMC 0.911 +/- 0.131) (p.3).
  - S2 with GRAVITY-only astrometry + RV: f_SP = 1.016 +/- 0.226 (MCMC 1.024 +/- 0.260), ~4 sigma (p.5).
- Mass and distance of Sgr A*: M_bh = 4.2996 +/- 0.0118 x 10^6 M_sun; R0 = 8275.9 +/- 8.6 pc (Table B.1, p.10).
- Dark-matter spike: "we find no evidence for a significant dark matter spike in the Galactic Center" (p.8).
==============================================================================
## FULL-TEXT CHECK, papers 1-3 (get_paper_content, fullText=true; these lines supersede the first-pass lines they name)

### Paper 1 (2607.12664)
- Astrometric errors, Methods p.19 (verbatim): "The uncertainties of the newly inferred sources are derived from the scatter over the different runs, taking also into account the finite resolution of the pixel grid." and "We derive a discretization error by calculating the root mean square error for a single pixel within the image for both directions, RA. and Dec., independently. Considering a pixel size of 0.8 mas per pixel, this evaluates to ~207 uas and represents the statistical uncertainty in astrometry originating from the pixel grid in GR images. This value acts effectively as a noise floor."
- Fitting-code cross-check, p.19: the two independent fitting codes put S301 "at the same position as with the imaging code, within statistical errors".
- Angle conventions, p.20: "angle conventions follow [4]". Ref [4] = Gillessen, S., et al. 2009, ApJ 692(2), 1075-1109. No further definition of i, Omega, omega in the full text; the ascending node is not defined.
- Extended Data Fig. 5 caption, p.21: posterior of "the six parameter orbit fit in Extended Data Table 2, left column". Left column = set A (first pair), the set paper 4 adopts.
- Mass and distance: p.6 M_MBH = 4.297 x 10^6 M_sun, R_0 = 8277 pc; p.21 "from [5]". Paper 4 (p.2) attributes these values to GRAVITY 2022.

### Paper 2 (2607.26713)
- Eq. (3), p.3 (layer text, PN term in harmonic coordinates): a^1PN = GM/(c^2 r^2) [ (4GM/r - v^2) r/r + 4 rdot v ].
- Eq. (6), p.3 (layer; the radical sign appears as the glyph q): p = sqrt((d+R)^2 + z_BH^2), q = sqrt((d-R)^2 + z_BH^2). Note p^2 - q^2 = 4Rd.
- Eq. (7), p.3 (layer): U_ring = 2GM_ring K(k) / (pi p), k = 4Rd/p^2 (no radical in layer; Appendix B confirms k = 4Rd/p^2).
- Eq. (8), p.3: a_ring = a_d u_d + a_z u_z; u_d = cos(phi_D) u_x + sin(phi_D) u_y (in-plane unit vector along d); u_z = unit vector along z_BH. Ring linear density lambda = M_ring/(2 pi R) (p.3).
- Eq. (9), p.3, LAYER TEXT (verbatim sequence; the layer does not resolve the bracket pairing and the trailing d, which is d(phi) with phi lost): a_d = -GM_ring/(pi p) " K(k) d^2 + 1 - (z^2+R^2)/d^2 ! E(k) q^2 # d,   (9)
- Eq. (10), p.3, LAYER TEXT (clean): a_z = -2GM_ring / (pi p q^2) E(k) z,   (10)
- Eq. (9) closed form, MY DERIVATION from Eq. (7) (a = +grad U, the convention that reproduces Eq. (10); K, E as functions of k = 4Rd/p^2; not a transcription, verify against the PDF):
  a_d = -(GM_ring/(pi p d)) [ K(k) + E(k) (d^2 - R^2 - z^2) / q^2 ].
  Checks: R -> 0 gives -GM_ring d/(p r^2), r^2 = d^2 + z^2, the point-mass limit; the z-component from the same U reproduces Eq. (10).
- Eq. (14), p.4: same structure as Eq. (9) with M_ring -> dM and R -> r'. Eq. (15), p.4 (layer, clean): da_z = -2G dM / (pi p q^2) E(k) z.
- Eq. (16), p.4 (layer): M_tot = 2 pi Integral_{r_min}^{r_max} B r'^(1-gamma) dr' = 2 pi { B (r_min^(2-gamma) - r_max^(2-gamma)) / (gamma - 2) for gamma != 2 ; B log(r_max/r_min) for gamma = 2 }.
- Eq. (17), p.4: a^disk = Integral_{r_min}^{r_max} ( da_d u_d + da_z u_z ). Eq. (18), p.4: a_d ~ a_z ~ GM_tot / r^2 for r >> r_max.
- Singularity: p.3, the acceleration is singular at z = 0 and d = R (K diverges as k -> 1, E finite). Fig. 3 caption, p.4: "the singularity at d = r_max". No softening anywhere in the text.
- Appendix B, p.11: K(k) = Integral_0^{pi/2} dtheta / sqrt(1 - k sin^2 theta) (B.1); E(k) = Integral_0^{pi/2} sqrt(1 - k sin^2 theta) dtheta (B.2); "where k = 4Rd/p^2".
- Appendix C, p.11: Gaussian-frame transformation (C.1)-(C.3), psi_BH = omega_BH + f.
- Eqs. (24)-(26), p.5: R = R_disk + R_1PN; S = S_disk + S_1PN; W = W_disk (the disk gives no W_1PN). R_1PN, S_1PN from Poisson & Will (2012), reproduced in Appendix D.
- Eqs. (27)-(32), p.5: osculating (Gauss) equations after Poisson & Will (2012); di/dt (30) and dOmega/dt (31) are proportional to W; df/dt (32) includes the precession terms. Eq. (27) layer shows "(1 - e cos f)^3", which looks like a layer error; verify.
- Eq. (50), p.8: CW disk "i_disk = 127 deg +/- 2 deg, Omega_disk = 279 deg +/- 2 deg", with Omega_disk differing by pi from the observed Omega of Paumard et al. (2006), as in Eq. (42).
- Eq. (51), p.8: "r^L_min ~ 0.04 pc, r^L_max ~ 0.2 - 0.5 pc" (Paumard 2006; Lu 2009; von Fellenberg 2022; Goldwurm 2026).
- Eq. (52), p.9 (layer): Sigma = M_extra / (2 pi log(r^L_min/r_min) r^2) + M_L / (2 pi log(r^L_max/r^L_min) r^2), M_L = 10^4 M_sun.
- Eqs. (53)-(54), p.9-10 (paper's values): S2 i_BH ~ 86 deg, beta ~ 232 deg; S301 i_BH ~ 21 deg (107 deg), beta ~ 342 deg (252 deg). Used for Table 2: i_BH = 21 deg, beta = 342 deg (p.10).
- MY ARITHMETIC CHECK (not in the paper; uses Eqs. 43-44, Eq. 50 and the paper's own angles):
  - S2: cos i_BH = sin(i_disk) sin(i_obs) cos(dOmega) + cos(i_disk) cos(i_obs) with dOmega = 48.21 - 279 = -230.79 deg gives i_BH = 86.3 deg; cos(beta) = -sin(i_disk) sin(dOmega) gives beta = 231.8 deg. Both match the paper.
  - S301 set B (Omega = 256.9 deg, mapped by +180 deg to 76.9 deg): i = 107.1 deg, beta = 252.4 deg. Both match the paper's parenthetical values.
  - S301 set A (Omega = 73.8 deg mapped to 253.8 deg): i = 20.7 deg matches, but beta = 70.1 deg or 289.9 deg, NOT 342 deg. Unresolved. Check beta = 342 deg before relying on Table 2's S301 disk values.
- CROSS-PAPER CHECKS (my arithmetic on paper 1 numbers; not reconciled in the text):
  - Table 2 in-plane Schwarzschild for S301 = 1.684 deg. Using 3 pi R_S / (a(1-e^2)) with R_S = 10 uas (paper 1 p.4) and set A (a(1-e^2) = 2.77 mas) gives ~1.95 deg, consistent with paper 1's 2.0 deg (p.6).
  - Table 2 'maximum LT' for S301 = 0.089 deg. Paper 1 (p.8) gives 0.11 deg chi cos(xi) per revolution; my evaluation of paper 1's Eq. (1) at chi = 1, aligned gives ~0.11 deg.

### Paper 3 (2409.12261)
- Per-instrument astrometric/RV precision and the integrator are not stated anywhere in the full text (checked by search for uncertainty, accuracy, integrat, solver, numerical).
- S29 t_peri error (Table B.1, p.10): the text layer reads "2021.4102 x10^-4" with the leading digits missing. Recover from the PDF.
- Sec. 4.2, p.6 (verbatim): the 1 sigma upper limit on the enclosed mass within S2's orbit "is lower than ~1200 M_sun for a power-law with slope s <~ -1.2 and for a Plummer profile with scale radius a <~ 8 mpc, which corresponds approximately to the apocenter distance of S2"; retrograde precession "<~ 1 arcmin per orbit", against "~12 arcmin per orbit" Schwarzschild.
- Sec. 4.3, p.6: simulated stellar-cusp enclosed mass within S2's orbit 1210 M_sun (two-component model: 1300 M_sun); "approximately 80 stellar mass black holes within S2's orbit"; "no substantial evidence for a significant enhancement of dark matter density near Sagittarius A*".
- Sec. 4.2, p.6: 1 sigma = 68.3% c.i., 3 sigma = 99.7% c.i. (Appendix B, p.10).
- Fig. 2 caption, p.4: S2 astrometry from SHARP and NACO (black) and GRAVITY (red), 1992 to 2022.

==============================================================================
## PAPER 4 - arXiv 2609.06936, "Catching Up with the Fastest Star in the Galaxy: A Two-Passage GRAVITY+ Campaign to Detect S301's Schwarzschild Precession" (Catanzarite, MNRAS 2026)
Source: answer_pdf_queries (pages 1, 3, 4, 5, 8, 9, 10) and the full text (all 10 pages, returned inline). Page numbers are PDF pages.

### Q1. Planned campaign
- Passages: "the 2031.8 and 2040.5 periapsis passages" (p.1). First periapsis 2031-10-22 (p.3, Sec. 3.4); second at solar elongation 170 deg, "comfortably in season" (p.3); "passage 2 (2040.5) is comfortable" (p.9).
- Epochs: 131 (p.1; Table 2 caption, p.4). Span "between 2028.50 and 2041.69" (p.5, Sec. 4.1).
- Cadence: Design D, constrained Fisher-optimal, "after reserving 0.40 of the budget for the +/-20 d windows around both periapses" (p.3). Season ~March-September from Paranal (p.3). 50 epochs (38.2%) within +/-20 d of a periapsis; 81 in clusters "in 2033 and 2037-2038" (p.5).
- Blackout: 153 d, unobservable 2031-11-06 to 2032-02-01 (p.3).
- Astrometric precision: 207 uas per epoch, "the value achieved on S301 itself in the discovery data" (p.3). Forecast 100 uas for GRAVITY+ (p.3; p.5).
- Radial velocity: none. "Radial velocity was rejected before any data were simulated" (p.3). ERIS-scaled benchmark 1621 km/s per epoch (p.3). MICADO "1 km/s" figure is the discovery paper's (p.9). HARMONI estimate sigma_RV ~ 20.6 km/s, S/N ~741, "a factor ~20 short of the MICADO figure" (p.9). ERIS-precision channel changes omega-dot precision by 4.5% (p.9).
- Results: injected omega-dot 0.2307 deg/yr; recovered 0.2324 +/- 0.0023 deg/yr at 207 uas, "96 sigma" (p.5, Table 3). At 100 uas, 0.2315 +/- 0.0012 deg/yr, "198 sigma" (p.5). Photometry check (not fitted): chi^2 = 137.5 for 131 epochs, reduced 1.05 (p.5).

### Q2. Dynamical model
- Base: Schwarzschild, first post-Newtonian monopole. Eq. (1), p.5 (layer, reconstructed): omega-dot = 6 pi G M_bh / (c^2 a (1-e^2) P), "the Mercury-perihelion formula applied to S301". Used as a check: evaluated on each bootstrap resample's own P and e gives 0.2318 +/- 0.0011 deg/yr, difference 0.16 sigma from the fitted value (p.5).
- Fit (p.5, Sec. 3.5): seven parameters (six Keplerian + omega-dot) applied as omega_eff(t) = omega + omega-dot (t - t_ref); Romer delay in model; emission time by Newton iteration; Levenberg-Marquardt; 1000 bootstrap resamples; M_bh and R0 fixed; "We inject no Lense-Thirring term, so i and Omega stay fixed" (p.5).
- Kerr / spin: "spin-dependent (Kerr) terms, which would require a value for the unmeasured spin" are omitted from the light-curve model (p.2-3, Sec. 3.1). Frame dragging for an aligned spin adds a purely apsidal term, "orbit.lense thirring apsidal rate", validated against the discovery paper's 0.114 deg per orbit at chi = 1 (p.8).
- PN order beyond 1PN: not stated in the text.
- Freeing M_bh and R0 under priors: omega-dot 0.2320 +/- 0.0024 deg/yr; e, i, omega widen by 1.36, 1.56, 1.51x (p.6, Sec. 5.3). Precision cost 1.04x (p.8, Sec. 5.9; p.10).
- Extended-mass check (p.7, Sec. 5.4): post-Newtonian point mass + Plummer (alpha = 5), integrated with the Rubilar & Eckart (2001) Eq. A.7 equation of motion. Bias -0.00013 deg/yr at the 1 sigma bound (<~3000 M_sun; Plummer scale length 0.3 arcsec = 12.04 mpc) and -0.00032 deg/yr at 3 sigma (<~7500 M_sun).

### Q3. Spin contamination
- p.8 (Sec. 5.8): "the possible contamination is 0.01192-0.01325 deg yr^-1: 5.17-5.74 per cent of the Schwarzschild signal, or 5.2-5.7 sigma of the random precision (5.0-5.5 sigma of the reference-marginalized precision)". Inputs: chi = 1 aligned (discovery fiducial) and chi = 0.90 +/- 0.06 (Daly et al. 2024).
- p.9, Table 4: "Sgr A* spin, chi = 0.9-1 | Degeneracy | +0.01192 to +0.01325 (5.0-5.5 sigma) | Unbounded; out of scope."
- Headline (p.8 Sec. 5.9; p.9 Table 4): +/-0.0074 deg/yr, 3.20% of the signal, "a 31 sigma detection". Before the empirical inflation: +/-0.0050 deg/yr (2.17%) with correlated systematics, and +/-0.0024 deg/yr white-noise bottom line.

### Q4. Multi-star fitting or disk degeneracy
- Multi-star: p.1 "a companion paper (Piran et al. 2026) develops the multi-star method that measurement requires". p.3 (Sec. 3.3): Piran's design "synchronizing observations across several stars to subtract Newtonian confusion; with one star and one observable, periapsis coverage is what matters". p.8 (Sec. 5.8): "A measured omega-dot from S301 alone cannot separate pure Schwarzschild precession from Schwarzschild plus aligned frame dragging; Piran et al. (2026)'s multi-star method exists to break that degeneracy."
- Disk: not modelled in this paper. p.7 (Sec. 5.4) is the only disk statement: Piran et al. "model a different quantity, the nodal precession from a flattened disc competing with the spin signal, but likewise find S301 about ten times better protected than S2."
- Extended mass (spherical), p.9 Table 4: -0.00013 to -0.00032 deg/yr (0.06-0.14 sigma); p.10 (v): "extended mass biases omega-dot by at most 0.14 sigma". Not modelled: heavy-tailed encounters "would be dominated by large, non-local jumps" (GRAVITY 2026 Methods, quoted p.7).
- Confusion (p.7, Sec. 5.5): 90th-percentile floor 37-92 uas per epoch; precision x1.00 (p.9 Table 4: "Subdominant").

### Extras
- p.1: "The spin of Sgr A* could shift the rate by 5.17-5.74 per cent and cannot be bounded, so the campaign is a spin-agnostic test of general relativity."
- p.2: S301 at periapsis 136 R_S, 25 600 km/s (fastest known star). Adopted orbit (Table 1, p.2) is "Solution 1" of GRAVITY 2026 Ext. Data Table 2; "a second, nearly degenerate solution is published but not used". This is set A in the first-pass table.
- p.2: effective temperature 7300 K "from standard F1.5V tables".
- p.2: M_sun = (4.297 +/- 0.012) x 10^6, R0 = 8277 +/- 9 pc (GRAVITY 2022), systematics ~40,000 M_sun and ~30 pc.
- p.6 (Sec. 5.1): Romer correction -0.06 to 5.60 d (RMS 3.73 d). Ignoring it biases omega-dot by -0.0023 deg/yr (1.0 sigma) and i, Omega, omega by 0.53, -1.40, -0.75 deg (7-8 sigma).
- p.7 (Sec. 5.2): reference offset 50 uas prior (100 uas conservative); fitted jointly, precision x1.04.
- p.8 (Sec. 5.7): correlated systematics (60 uas per-run calibration; 50 uas OU red noise, 0.5 yr): scatter 2.09x the white-noise value.
- p.9 (Sec. 6.2): HARMONI route changes omega-dot precision by -4.1%. ELT first light 2029; HARMONI science late 2030; passage 1 (Oct 2031) "a stretch goal".
- p.9-10 (Sec. 6.4): a Kepler solver "diverged silently" at e = 0.98; fixed by a bisection-safeguarded Newton iteration, verified on 1.4 x 10^6 random cases to e = 0.9999.
