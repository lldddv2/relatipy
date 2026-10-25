import numpy as np, time
import pn, pn_fast
from common import f_before_peri
R=np.pi/180; a_rg=1035.45/0.042414091583065465; e=0.8843
P=2*np.pi*a_rg**1.5
x0,v0=pn.kepler_to_cart(a_rg,e,134.567*R,228.171*R,66.263*R,f_before_peri(e,1.0,0.125))
t_end=1.25*P; S=np.array([0.,0.,0.9])
ts,Xs,Vs=pn.integrate(2,S,x0,v0,t_end,n_out=4001,rtol=1e-13,atol=1e-15)
for spo in [500,1000,2000,4000,8000]:
    ds=2*np.pi/spo
    t0=time.perf_counter(); T,X,V=pn_fast.integrate(x0,v0,t_end,2,1.0,S,ds); dt=time.perf_counter()-t0
    Xi=np.array([pn_fast._hermite(T,X,V,t)[0] for t in ts[1:-1]])
    err=np.abs(Xi-Xs[1:-1]).max()
    print(f"steps/orbit~{spo:5d}  nsteps={len(T):7d}  {dt*1e3:8.2f} ms  max|dx| vs scipy(1e-13) = {err:.3e} rg = {err*5.1243:.2e} uas")
