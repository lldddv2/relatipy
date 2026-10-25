"""Truth elements (GRAVITY S301 discovery paper ED Table 2; GRAVITY 2024 Table B.1) and epochs."""
import numpy as np
import fc

D = fc.DEG
EL = {  # a [mas], e, i, Omega, omega [rad], tp [yr]   (Gillessen+09 / GRAVITY convention)
    'S301A': [83.0, 0.9832, 124.09 * D, 73.8 * D, 293.4 * D, 2023.126],
    'S301B': [83.0, 0.9824, 122.84 * D, 256.9 * D, 115.1 * D, 2023.125],
    'S2':  [125.02, 0.88444, 134.67 * D, 228.21 * D, 66.279 * D, 2018.3789],
    'S29': [390.25, 0.96880, 144.24 * D, 4.9259 * D, 203.68 * D, 2021.4102],
    'S38': [142.49, 0.81807, 168.69 * D, 122.43 * D, 40.065 * D, 2022.6843],
    'S55': [104.24, 0.72980, 159.59 * D, 319.43 * D, 327.77 * D, 2009.4738],
}
# published 1-sigma (same units) -- used only for the realism check of the mock past data
EL_ERR = {
    'S301A': [0.7, 0.0010, 1.10 * D, 3.5 * D, 2.2 * D, 0.010],
    'S2':  [0.03, 6e-5, 0.02 * D, 0.03 * D, 0.029 * D, 1e-4],
    'S29': [0.94, 9e-5, 0.09 * D, 0.159 * D, 0.17 * D, np.nan],
    'S38': [0.04, 2.2e-4, 0.19 * D, 1.12 * D, 1.118 * D, 8e-4],
    'S55': [0.05, 1.8e-4, 0.17 * D, 0.97 * D, 0.93 * D, 7.6e-3],
}
G0 = np.array([fc.M6_FID, fc.R0_FID, 0, 0, 0, 0, 127 * D, 99 * D])   # CW disk (Paumard+06 / Foschi+ Eq. 50)
# 19 GRAVITY epochs of S301 (ED Table 1), 207 uas per coordinate, no RV
S301_PAST = np.array([2017.2389, 2021.5759, 2023.2604, 2023.4162, 2023.4298, 2023.5063, 2024.2325,
                      2024.2380, 2024.2405, 2024.3281, 2024.4073, 2024.4756, 2024.4837, 2024.5438,
                      2025.2070, 2025.3602, 2025.4393, 2025.5159, 2025.5978])


def period_yr(el, M6=fc.M6_FID, R0=fc.R0_FID):
    a_rg = el[0] * 1e-3 * R0 * 1e3 / fc.rg_au(M6)
    return 2 * np.pi * a_rg ** 1.5 * fc.tg_yr(M6)


def peris(el, t0, t1):
    P = period_yr(el); k0 = int(np.floor((t0 - el[5]) / P)) - 1
    tp = el[5] + P * np.arange(k0, k0 + int((t1 - t0) / P) + 4)
    return tp[(tp > t0 - 1) & (tp < t1 + 1)]


def uniform(t0, t1, per_yr):
    n = int(round((t1 - t0) * per_yr))
    return t0 + (np.arange(n) + 0.5) / per_yr


def campaign(el, t0, t1, n, half):
    out = []
    for tp in peris(el, t0, t1):
        tt = np.linspace(tp - half, tp + half, n)
        out += list(tt[(tt >= t0) & (tt <= t1)])
    return np.array(out)


def s301(sol='A', t_end=2041.0, s_ast=100.0, s_rv=1.0, rv_start=2030.0, past=True,
         ast_per_yr=6, n_peri_ast=20, rv_per_yr=4, n_peri_rv=15, half=0.3, nstep=3000):
    el = EL['S301' + sol]
    ta = np.concatenate([uniform(2026.0, t_end, ast_per_yr), campaign(el, 2026.0, t_end, n_peri_ast, half)])
    sa = np.full(ta.shape, s_ast)
    if past:
        ta = np.concatenate([S301_PAST, ta]); sa = np.concatenate([np.full(19, 207.0), sa])
    if s_rv is None or s_rv <= 0:
        tr = None; sr = None
    else:
        tr = np.concatenate([uniform(rv_start, t_end, rv_per_yr), campaign(el, rv_start, t_end, n_peri_rv, half)])
        sr = np.full(tr.shape, s_rv)
    return fc.Star('S301' + sol, el, ta, sa, tr, sr, nstep=nstep)


def friend(name, t_end=2041.0, s_past=50.0, s_fut=30.0, s_rv=10.0, per_yr=4, nstep=1500):
    el = EL[name]
    tp = uniform(2017.0, 2026.0, per_yr); tf = uniform(2026.0, t_end, per_yr)
    if name == 'S2':
        tp = np.concatenate([tp, np.linspace(2018.38 - 0.3, 2018.38 + 0.3, 15)])
    ta = np.concatenate([tp, tf]); sa = np.concatenate([np.full(tp.shape, s_past), np.full(tf.shape, s_fut)])
    tr = uniform(2017.0, t_end, per_yr)
    if name == 'S2':
        tr = np.concatenate([tr, np.linspace(2018.38 - 0.3, 2018.38 + 0.3, 15)])
    return fc.Star(name, el, ta, sa, tr, np.full(tr.shape, s_rv), nstep=nstep)


def spin_dirs(el):
    """6 spin orientations (beta, lambda) w.r.t. the S301 orbital frame (node N, L x N, L)."""
    L = fc.pole(el[2], el[3]); N = np.array([np.cos(el[3]), np.sin(el[3]), 0.0]); M = np.cross(L, N)
    out = {}
    for lab, b, l in [('par', 0, 0), ('perpN', 90, 0), ('perpM', 90, 90), ('b45l45', 45, 45),
                      ('b120l200', 120, 200), ('b70l310', 70, 310)]:
        b *= D; l *= D
        out[lab] = np.sin(b) * np.cos(l) * N + np.sin(b) * np.sin(l) * M + np.cos(b) * L
    return out, (N, M, L)
