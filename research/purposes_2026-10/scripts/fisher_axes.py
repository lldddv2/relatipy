import numpy as np, warnings, time
warnings.filterwarnings("ignore")
from fisher import jac_star, setup, M0, R00, RAD
rng=np.random.default_rng(3)
s,t=setup("S62-like")
print("S62-like alone, 20-yr campaign. Distribution over 60 isotropic spin axes.")
print(f"{'a*':>4} {'sig_ast':>7} | s(a*) p16  p50  p84 | axis err[deg] p50 p84 | frac s(a*)<0.1")
for ast in [0.3,0.6,0.9]:
    J=[]; t0=time.perf_counter()
    for k in range(60):
        cth=rng.uniform(-1,1); th=np.arccos(cth); ph=rng.uniform(0,2*np.pi)
        g=[M0,R00,ast,th,ph]
        try: J.append((jac_star((g,s,t)),th))
        except Exception: pass
    for sa,sv in [(30.,10.),(10.,5.)]:
        sd_a=[]; sd_ax=[]
        for Jk,th in J:
            n=Jk.shape[0]//3; W=np.concatenate([np.full(2*n,1/sa**2),np.full(n,1/sv**2)])
            F=Jk.T@(Jk*W[:,None]); C=np.linalg.inv(F); sd=np.sqrt(np.abs(np.diag(C)))
            sd_a.append(sd[2]); sd_ax.append(np.degrees(np.sqrt(sd[3]**2+(np.sin(th)*sd[4])**2)))
        sd_a=np.array(sd_a); sd_ax=np.array(sd_ax)
        print(f"{ast:4.1f} {sa:7.0f} | {np.percentile(sd_a,16):5.3f} {np.median(sd_a):5.3f} {np.percentile(sd_a,84):5.3f} |"
              f" {np.median(sd_ax):8.1f} {np.percentile(sd_ax,84):6.1f} | {np.mean(sd_a<0.1):.2f}   (n={len(J)}, {time.perf_counter()-t0:.0f}s)",flush=True)
