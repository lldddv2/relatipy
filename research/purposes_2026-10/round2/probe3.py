import numpy as np, time, warnings
from astropy import units as u
from astropy.constants import c
from relatipy import Kerr
warnings.filterwarnings("ignore")
for lab,a,e in (("S2",24414.,0.8843),("S301",16210.,0.98273)):
  for rtol in (1e-12,1e-13):
    bh=Kerr(mass=4.297e6*u.M_sun,spin=0.7); rg=bh.r_g; tg=rg/c
    ra=a*(1+e); va=np.sqrt((1-e)/(a*(1+e)))
    o=bh.orbit(R=(ra+1)*rg,Theta=1.0*u.rad,Phi=0.3*u.rad,vR=0*c,vTheta=0*u.rad/tg,vPhi=(va/ra/np.sin(1.0))*u.rad/tg)
    P=2*np.pi*a**1.5
    t0=time.perf_counter(); s=o.solve(tau_span=(0*u.s,1.0*P*tg),method="dop853",rtol=rtol); t1=time.perf_counter()
    tau=(s.tau/tg).decompose().value; R=s.R.to_value(rg); vR=(s.vR/c).decompose().value
    h=np.diff(tau); i=np.argmin(R)
    w=np.sqrt(2/R.min())/R.min()
    print(f"{lab} rtol {rtol}: nsteps {len(tau)} in {1e3*(t1-t0):.0f} ms; h_peri {h[i-1]:.1f} (h*w {h[i-1]*w:.3f}); h_apo {h[0]:.0f}; rmin {R.min():.1f}")
