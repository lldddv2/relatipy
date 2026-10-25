"""Aggregate Task A json and Task B per-orientation json into compact tables."""
import json, glob, numpy as np
D = "/home/user/relatipy/research/purposes_2026-10/round2/"

def pc(x):
    x = np.asarray([v for v in x if v is not None and np.isfinite(v)])
    if len(x) == 0: return "   n/a"
    a, b, c = np.percentile(x, [16, 50, 84]); return f"{b:+6.2f} [{a:+6.2f},{c:+6.2f}]"

def pcu(x, fmt="{:.3g}"):
    x = np.asarray([v for v in x if v is not None and np.isfinite(v)])
    a, b, c = np.percentile(x, [16, 50, 84]); return (fmt + " [" + fmt + "," + fmt + "]").format(b, a, c)

print("=== TASK A (S2): bias in units of the fit model's own Fisher sigma; [published-sigma units] ===")
for f in sorted(glob.glob(D + "taskA_spin*.json")):
    for r in json.load(open(f)):
        if r.get('truth', 'full truth') != 'full truth' and 'attrib' not in f: continue
        s = f"{f.split('/')[-1][6:-5]:<16}{r['truth'][:30]:<31}{r['model']:<34} chi2 {r['chi2']:6.3f} rmsGRAV {r['rms_grav_uas']:5.2f} RV {r['rms_rv']:5.2f} |"
        for n in ('M', 'R0', 'fsp', 'f'):
            if 'bias_sig_' + n in r: s += f" {n} {r['bias_sig_' + n]:+.2f}"
        print(s)

fs = sorted(glob.glob(D + "taskB/k*.json"))
if fs:
    R = [json.load(open(f)) for f in fs]
    print(f"\n=== TASK B (S301, a*=0.7, N={len(R)} orientations): median [16,84] ===")
    mods = []
    for r in R:
        for f in r['fits']:
            key = (f['model'], f['sig_ast'])
            if key not in mods: mods.append(key)
    print(f"{'model':<26}{'ast':>4} {'bias a* [sig]':>22} {'bias M [sig]':>22} {'bias R0 [sig]':>22} {'RV rms km/s':>22} {'RV max |res| 2031peri':>26} {'chi2':>20}")
    for m, sa in mods:
        F = [f for r in R for f in r['fits'] if f['model'] == m and f['sig_ast'] == sa]
        print(f"{m:<26}{sa:>4} {pc([f.get('bias_sig_ast') for f in F]):>22} {pc([f.get('bias_sig_M') for f in F]):>22} "
              f"{pc([f.get('bias_sig_R0') for f in F]):>22} {pcu([f['rms_rv'] for f in F]):>22} {pcu([f['rv_peri_max'] for f in F]):>26} {pcu([f['chi2'] for f in F], '{:.3g}'):>20}")
        print(f"{'':<30} abs bias a*: {pcu([f.get('bias_ast') for f in F], '{:+.3f}')}  sig(a*): {pcu([f.get('sig_ast') for f in F])}"
              + (f"  fsp-1: {pcu([f.get('bias_fsp') for f in F], '{:+.4f}')}" if 'fsp' in F[0]['sig'] else ""))
    print("\nFisher (iii) at truth: median [16,84] of sigma(a*), sigma(M)/M, sigma(R0)/R0, spin-axis 1-sigma cone [deg]")
    for key in ("rv_100", "ast_100", "both_100", "ast_30", "both_30"):
        print(f"  {key:<9} s(a*) {pcu([r['fisher_iii'][key]['ast'] for r in R])}  s(M)/M {pcu([r['fisher_iii'][key]['M']/4.297 for r in R])}"
              + (f"  s(R0)/R0 {pcu([r['fisher_iii'][key]['R0']/8.277 for r in R])}" if 'R0' in R[0]['fisher_iii'][key] else "")
              + f"  axis {pcu([r['fisher_iii'][key]['axis_deg'] for r in R])}")
    print("\nfirst-order-T truth vs exact Schwarzschild rays at the 2031 campaign (no refit):")
    psi = [r['peri'][0][2] for r in R]
    print(f"  psi_peri: {pcu(psi, '{:.0f}')} deg ; RV max(+-1 d) {pcu([r['floor_iii']['rv_max_1d'] for r in R])} km/s ;"
          f" RV max(+-15 d) {pcu([r['floor_iii']['rv_max_15d'] for r in R])} ; img max {pcu([r['floor_iii']['ast_max'] for r in R])} uas")
    for r in sorted(R, key=lambda r: r['peri'][0][2]):
        print(f"   k={r['k']:2d} psi_peri={r['peri'][0][2]:5.1f}  1st-order RV err max {r['floor_iii']['rv_max_1d']:8.2f} km/s  img {r['floor_iii']['ast_max']:7.2f} uas")
