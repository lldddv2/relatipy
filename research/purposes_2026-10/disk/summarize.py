"""Compact Q1 tables from results_*.json (median over the 6 spin orientations)."""
import json, sys, numpy as np
fn = sys.argv[1] if len(sys.argv) > 1 else 'results_A.json'
R = json.load(open(fn))
ast = float(sys.argv[2]) if len(sys.argv) > 2 else 100.0
rows = [('a', 0.0, '-'), ('a_noprior', 0.0, '-'), ('b', 1.0, 'CWD'), ('b', 0.3, 'CWD'), ('b', 1.0, 'P30'), ('b', 0.3, 'P30'),
        ('c0', 0.0, '-'), ('c', 1.0, 'CWD'), ('c', 0.3, 'CWD'), ('c', 1.0, 'P30'), ('c', 0.3, 'P30')]
cols = [(t, rv) for t in (2035.0, 2041.0) for rv in (0, 20.0, 5.0, 1.0)]
print(f"{fn}: S301 future astrometry {ast:.0f} uas; entries sigma_perp(rms of 2 comps)/sigma_par [chi units], median over 6 spins")
print(f"{'case':<16}" + "".join(f"{('%d' % t)+'/RV'+('%g' % rv if rv else '-'):>13}" for t, rv in cols))
for case, md, orient in rows:
    line = f"{case + ' ' + (orient + ' %.0f' % (md * 1e3) if md else ''):<16}"
    for t, rv in cols:
        sel = [r for r in R if r['case'] == case and r['md'] == md and r['orient'] == orient and r['t_end'] == t and r['s_rv'] == rv and r['s_ast'] == ast]
        if not sel: line += f"{'':>13}"; continue
        sp = np.median([np.sqrt((r['sN'] ** 2 + r['sM'] ** 2) / 2) for r in sel]); sl = np.median([r['sL'] for r in sel])
        line += f"{sp:>7.2f}/{sl:<5.2f}"
    print(line)
# spin-orientation dependence and other quantities, fiducial column
print("\nfiducial column 2041/RV1: per-case spread over the 6 spins (max/min of sigma_N), sigma(chi) range, ang err range [deg], sigma(Md) [Msun], sigma(i_d,O_d) [deg]")
for case, md, orient in rows:
    sel = [r for r in R if r['case'] == case and r['md'] == md and r['orient'] == orient and r['t_end'] == 2041.0 and r['s_rv'] == 1.0 and r['s_ast'] == ast]
    if not sel: continue
    sN = [r['sN'] for r in sel]; sc = [r['schi'] for r in sel]; sa = [r['sang_deg'] for r in sel if r['spin'] != 'par'] + [r['sang_deg'] for r in sel if r['spin'] == 'par']
    extra = f" sMd={np.median([r['sMd'] for r in sel]):.0f} sid={np.median([r['sid'] for r in sel]):.1f} sOd={np.median([r['sOd'] for r in sel]):.1f}" if 'sMd' in sel[0] else ''
    print(f"{case + ' ' + orient + ' %.0f' % (md*1e3):<16} sN={np.median(sN):.3f} (max/min {max(sN)/min(sN):.3f}) sM={np.median([r['sM'] for r in sel]):.3f} sL={np.median([r['sL'] for r in sel]):.3f} schi {min(sc):.2f}-{max(sc):.2f} ang {min(sa):.0f}-{max(sa):.0f}deg sM/M={np.median([r['sM6'] for r in sel])/4.297:.1e}{extra}")
