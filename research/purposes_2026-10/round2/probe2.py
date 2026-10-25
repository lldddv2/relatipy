import numpy as np, time, warnings
from astropy import units as u
from astropy.constants import c
from relatipy import Kerr
warnings.filterwarnings("ignore")
M=4.297e6*u.M_sun
for spin in (0.0,0.7):
    bh=Kerr(mass=M,spin=spin); rg=bh.r_g; tg=(rg/c).to(u.s)
    a=24414.; e=0.8843; rp=a*(1-e); vp=np.sqrt((1+e)/rp)
    for N in (4000,16000):
        t0=time.perf_counter()
        o=bh.orbit(R=(rp+1)*rg,Theta=1.0*u.rad,Phi=0.3*u.rad,vR=0*c,vTheta=0*u.rad/tg,vPhi=(vp/rp/np.sin(1.0))*u.rad/tg)
        tau=np.linspace(0,1.2e7*(1+np.arange(1)[0]),N)*tg
        s=o.solve(tau_eval=tau,method="dop853",rtol=1e-12)
        t1=time.perf_counter()
        R=s.R.to_value(rg); Th=s.Theta.to_value(u.rad); Ph=s.Phi.to_value(u.rad)
        vR=(s.vR/c).decompose().value; vT=(s.vTheta*tg).to_value(u.rad); vP=(s.vPhi*tg).to_value(u.rad); ut=np.asarray(s.ut); T=(s.t/tg).decompose().value
        t2=time.perf_counter()
        print(f"spin {spin} N {N}: solve {1e3*(t1-t0):.1f} ms, access {1e3*(t2-t1):.1f} ms, t_end/yr {T[-1]*21.1648/3.15576e7:.2f}, ok {s.success}")
