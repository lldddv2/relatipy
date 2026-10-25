"""Fisher for S301 (P=8.7 yr, e=0.9825, rp~283 rg), campaigns from 2026, first pericentre 2031.
Orbit orientation unknown here -> isotropic spin axes (60) with S2's orbital angles as placeholder."""
import numpy as np, warnings, sys
warnings.filterwarnings("ignore")
from fisher import jac_star, M0, R00, RAD
from common import epochs, f_before_peri
P=8.7; a=float((P**2*M0*1e6)**(1/3)); e=0.9825
rng=np.random.default_rng(5)
print(f"S301: a={a:.1f} AU e={e} rp={a*(1-e)/0.042414:.0f} rg; a(1-e^2)^.75={a*(1-e*e)**.75/0.084828:.0f} R_S")
for base,lab in [(9.0,"2026-2035 (1 peri)"),(14.0,"2026-2040 (2 peri)")]:
    t1=5.0; t=epochs(P,base,t1); f0=f_before_peri(e,P,t1)
    s=[a,e,134.567*RAD,228.171*RAD,66.263*RAD,f0]
    J=[]
    for k in range(60):
        th=np.arccos(rng.uniform(-1,1)); ph=rng.uniform(0,2*np.pi)
        try: J.append((jac_star(([M0,R00,0.7,th,ph],s,t)),th))
        except Exception: pass
    print(f"\n{lab}: n_ep={len(t)}, a*=0.7, {len(J)} axes")
    print(f"  {'sig_ast':>7} {'sig_rv':>6} | s(a*) p16   p50   p84 | axis p50[deg] | frac<0.1")
    for sa,sv in [(100.,1.),(30.,10.),(30.,1.),(10.,5.),(10.,1.)]:
        sd_a=[]; sd_ax=[]
        for Jk,th in J:
            n=Jk.shape[0]//3; W=np.concatenate([np.full(2*n,1/sa**2),np.full(n,1/sv**2)])
            C=np.linalg.inv(Jk.T@(Jk*W[:,None])); sd=np.sqrt(np.abs(np.diag(C)))
            sd_a.append(sd[2]); sd_ax.append(np.degrees(np.hypot(sd[3],np.sin(th)*sd[4])))
        sd_a=np.array(sd_a)
        print(f"  {sa:7.0f} {sv:6.0f} | {np.percentile(sd_a,16):5.3f} {np.median(sd_a):5.3f} {np.percentile(sd_a,84):5.3f} | {np.median(sd_ax):9.1f} | {np.mean(sd_a<0.1):.2f}",flush=True)
