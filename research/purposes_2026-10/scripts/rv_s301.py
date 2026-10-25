"""S301 pericentre: how wrong are approximate redshift models, in km/s?
Reference: static-frame exact Doppler (no light bending), Schwarzschild-exact for a*=0."""
import numpy as np, warnings; warnings.filterwarnings("ignore")
from astropy import units as u
from astropy.constants import c
from relatipy import Kerr
import pn
from common import f_before_peri
CK=c.to_value(u.km/u.s); R=np.pi/180; RG_AU=0.042414091583065465
M=4.297e6*u.M_sun
for astar in [0.0,0.7]:
    bh=Kerr(mass=M,spin=astar); RG_S=(bh.r_g/c).to_value(u.s); YR=(1*u.yr).to_value(u.s)/RG_S
    P=8.7; a_rg=(P**2*4.297e6)**(1/3)/RG_AU; e=0.9825
    print(f"\na*={astar}  S301-like a={a_rg:.0f} rg rp={a_rg*(1-e):.0f} rg; orbit angles = S2 placeholder; LOS 50 deg from spin axis")
    x0,v0=pn.kepler_to_cart(a_rg,e,134.567*R,228.171*R,66.263*R,f_before_peri(e,P,0.25))
    orb=bh.orbit(x=x0[0]*bh.r_g,y=x0[1]*bh.r_g,z=x0[2]*bh.r_g,vx=v0[0]*c,vy=v0[1]*c,vz=v0[2]*c)
    tau=np.linspace(0,0.5*YR*RG_S,40001)*u.s
    s=orb.solve(tau_span=(tau[0],tau[-1]),tau_eval=tau,method="dop853",rtol=1e-13)
    r=(s.R/bh.r_g).to_value(u.dimensionless_unscaled) if s.R.unit!=u.dimensionless_unscaled else s.R.value
    X=(s.xyz/bh.r_g).to_value(u.dimensionless_unscaled); U=(s.uxyz/c).to_value(u.dimensionless_unscaled)
    ut=s.ut.to_value(u.dimensionless_unscaled)
    for los_name,ghat in [("50deg",np.array([0.,np.sin(50*R),np.cos(50*R)])),("edge",np.array([0.,1.,0.]))]:
        rhat=X/np.linalg.norm(X,axis=1)[:,None]; Urad=(U*rhat).sum(1)[:,None]*rhat; Utan=U-Urad
        f=1/np.sqrt(1-2/r)
        z_ref = ut + (Urad@ghat)*f**2 + (Utan@ghat)*f - 1          # static frame, no bending (sign: +ghat = away)
        z_lib = ut + U@ghat - 1                                      # relatipy KerrMcmcModel
        v=U/ut[:,None]; vlos=v@ghat; vv=(v**2).sum(1)
        z_mult = (1+vlos)/np.sqrt(1-vv)/np.sqrt(1-2/r) - 1          # SR Doppler x GR redshift, flat photon
        z_1pn  = vlos + vv/2 + 1/r                                   # additive 1PN
        z_1pn5 = vlos + vv/2 + 1/r + vlos*(vv/2+1/r)                 # + leading cross term
        d=lambda z: np.abs(z-z_ref).max()*CK
        print(f"  LOS {los_name:5s}: max|v_los|={np.abs(vlos).max()*CK:7.0f}  max(z_ref)*c={z_ref.max()*CK:7.0f} km/s | "
              f"err lib={d(z_lib):6.1f}  mult={d(z_mult):6.1f}  1PN={d(z_1pn):6.1f}  1PN+cross={d(z_1pn5):6.1f} km/s")
