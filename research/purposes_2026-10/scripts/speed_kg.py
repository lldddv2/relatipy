import numpy as np, time, warnings; warnings.filterwarnings("ignore")
import kerrgeopy as kg
# S2-like in r_g: a=24413, e=0.8843 -> p=a(1-e^2); x=cos(60deg)
a=24413.; e=0.8843; p=a*(1-e*e); x=0.5
for spin in [0.9]:
    t0=time.perf_counter(); N=20
    for _ in range(N):
        o=kg.StableOrbit(spin,p,e,x); tf,rf,thf,phf=o.trajectory()
        lam=np.linspace(0, 2*np.pi/o.upsilon_r*2.2, 200)    # ~2 radial periods in Mino time
        r=rf(lam); th=thf(lam); ph=phf(lam); t=tf(lam)
    dt=(time.perf_counter()-t0)/N
    print(f"KerrGeoPy StableOrbit+trajectory+200 Mino-time samples (S2-like p={p:.0f}): {dt*1e3:.1f} ms/eval  -> {1/dt:.0f}/s")
    print("  (no coordinate-time inversion, no Romer inversion, no sky projection, no RV: those add ~10 Newton passes)")
    # with t-inversion at 200 target coordinate times (needed for observations)
    t0=time.perf_counter()
    for _ in range(N):
        o=kg.StableOrbit(spin,p,e,x); tf,rf,thf,phf=o.trajectory()
        T=np.linspace(0,2.2*2*np.pi/kg.fundamental_frequencies(spin,p,e,x)[0],200)
        lam=T/o.gamma
        for _ in range(8):
            h=1e-3*max(1,lam.max()*1e-9); d=(tf(lam+h)-tf(lam-h))/(2*h); lam=lam-(tf(lam)-T)/d
        r=rf(lam); th=thf(lam); ph=phf(lam)
    dt2=(time.perf_counter()-t0)/N
    print(f"  + coordinate-time inversion (8 Newton passes): {dt2*1e3:.1f} ms/eval -> {1/dt2:.0f}/s")
