"""Shared helpers: initial phase and realistic epochs with every pericentre covered."""
import numpy as np

def f_before_peri(e, P, dt):
    """True anomaly of a Kepler orbit a time dt BEFORE pericentre."""
    Mm = (-2*np.pi*dt/P) % (2*np.pi)
    E = Mm if e < 0.8 else np.pi
    for _ in range(100):
        E -= (E - e*np.sin(E) - Mm)/(1 - e*np.cos(E))
    return 2*np.arctan2(np.sqrt(1+e)*np.sin(E/2), np.sqrt(1-e)*np.cos(E/2))

def epochs(P, baseline, t_first_peri, per_yr=4.0, n_peri=30, half_width_frac=0.04):
    """Uniform monitoring (per_yr per year) plus a dense campaign of n_peri points
    within +-half_width_frac*P of every pericentre t_first_peri + k P."""
    t = list(np.linspace(0.0, baseline, max(int(per_yr*baseline), 20)))
    tp = t_first_peri
    while tp - half_width_frac*P <= baseline:
        lo, hi = max(tp-half_width_frac*P, 0.0), min(tp+half_width_frac*P, baseline)
        if hi > lo: t += list(np.linspace(lo, hi, n_peri))
        tp += P
    t = np.unique(np.round(np.array(t), 12))
    return t[(t > 1e-6) & (t <= baseline)]
