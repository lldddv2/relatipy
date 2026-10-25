"""Size of the Doppler metric factor omitted by KerrMcmcModel's redshift
1+z = u^t + u_away (flat-space photon momentum).  Static-frame derivation without
light bending gives 1+z = u^t + ghat.u_ON/sqrt(1-2/r), with u_ON = (u^r/sqrt(1-2/r), r u^th, r sin u^ph).
Difference ~ (2/r) u_rad_los + (1/r) u_tan_los.  Evaluated on exact geodesics."""
import numpy as np
from astropy import units as u
from astropy.constants import c
from relatipy import Kerr
import pn
from common import f_before_peri
CK=c.to_value(u.km/u.s); R=np.pi/180; RG_AU=0.042414091583065465
M=4.297e6*u.M_sun; bh=Kerr(mass=M,spin=0.0)
RG_S=(bh.r_g/c).to_value(u.s); RG_PER_YR=(1*u.yr).to_value(u.s)/RG_S
ghat=np.array([0.,np.sin(50*R),np.cos(50*R)])
print(f"{'star':12s}{'rp[rg]':>8} | max|v_los|  max grav+TD  max omitted(lib-vs-static)  [km/s]")
for nm,(a_au,e) in {"S2":(1035.45,0.8843),"S4716-like":(411.,0.756),"S62-like":(749.,0.976),"rp300":(400.,0.968),"rp100":(250.,0.983)}.items():
    a_rg=a_au/RG_AU; P_yr=2*np.pi*a_rg**1.5/RG_PER_YR
    x0,v0=pn.kepler_to_cart(a_rg,e,134.567*R,228.171*R,66.263*R,f_before_peri(e,P_yr,0.1*P_yr))
    orb=bh.orbit(x=x0[0]*bh.r_g,y=x0[1]*bh.r_g,z=x0[2]*bh.r_g,vx=v0[0]*c,vy=v0[1]*c,vz=v0[2]*c)
    tau=np.linspace(0,0.25*P_yr*RG_PER_YR*RG_S,20001)*u.s
    s=orb.solve(tau_span=(tau[0],tau[-1]),tau_eval=tau,method="dop853",rtol=1e-12)
    r=s.r.to_value(u.dimensionless_unscaled) if s.r.unit==u.dimensionless_unscaled else (s.r/bh.r_g).to_value(u.dimensionless_unscaled)
    X=(s.xyz/bh.r_g).to_value(u.dimensionless_unscaled); U=(s.uxyz/c).to_value(u.dimensionless_unscaled)
    ut=s.ut.to_value(u.dimensionless_unscaled)
    rhat=X/np.linalg.norm(X,axis=1)[:,None]
    Urad=(U*rhat).sum(1)[:,None]*rhat; Utan=U-Urad
    f=1/np.sqrt(1-2/r)
    lib = ut + U@ghat                                   # library form (sign aside)
    stat= ut + (Urad@ghat)*f**2 + (Utan@ghat)*f         # static-frame, no bending
    vlos=(U@ghat)/ut
    vv=((U/ut[:,None])**2).sum(1)
    print(f"{nm:12s}{r.min():8.0f} | {np.abs(vlos).max()*CK:10.0f}  {np.max(vv/2+1/r)*CK:11.1f}  {np.abs(stat-lib).max()*CK:14.2f}")
