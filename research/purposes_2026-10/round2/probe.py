import numpy as np, time, warnings
from astropy import units as u
from astropy.constants import c, G, M_sun
from relatipy import Kerr, KerrMcmcModel
warnings.filterwarnings("ignore")
M=4.297e6*u.M_sun; bh=Kerr(mass=M,spin=0.0); rg=bh.r_g; tg=(rg/c).to(u.s)
print("rg AU",rg.to(u.au), "tg s",tg)
m=KerrMcmcModel()
# Cartesian obs-frame state: x=Dec, y=RA, z=away
x0=np.array([1000.,2000.,3000.]); v0=np.array([0.01,-0.02,0.03])
orb=dict(x=x0[0]*rg,y=x0[1]*rg,z=x0[2]*rg,vx=v0[0]*c,vy=v0[1]*c,vz=v0[2]*c,t=0*u.s)
D=8.277*u.kpc; Dg=(D/rg).decompose().value
al,de,vr=m.get_ra_dec_vr(kerr=dict(mass=M,spin=0.0,vec=(0,0,1)),orbit=orb,sol={"t_eval":np.array([0.0,1.0])*tg,"method":"dop853","rtol":1e-12},distance=D)
print("alpha*D/rg",al[0]/206265*Dg, " delta*D/rg",de[0]/206265*Dg)
ut=1/np.sqrt(1-2/np.linalg.norm(x0)-v0@v0)  # rough (Schw. not exact)
print("v_los km/s",vr[0]," (u^t(1+vz)-1)c approx",(ut*(1+v0[2])-1)*c.to_value(u.km/u.s))
# Romer convention: t_obs at t_e=0 should be Z=3000 if t_obs=t_e+Z
al2,de2,vr2=m.get_ra_dec_vr(kerr=dict(mass=M,spin=0.0,vec=(0,0,1)),orbit=orb,sol={"t_obs":np.array([3000.0,-3000.0])*tg,"method":"dop853","rtol":1e-12},distance=D)
print("t_obs=+Z -> alpha diff",(al2[0]-al[0])/206265*Dg," t_obs=-Z -> ",(al2[1]-al[0])/206265*Dg)
# backward integration test with Kerr.orbit
o=bh.orbit(x=x0[0]*rg,y=x0[1]*rg,z=x0[2]*rg,vx=v0[0]*c,vy=v0[1]*c,vz=v0[2]*c)
try:
    s=o.solve(tau_span=(0*u.s,-1e5*tg),method="dop853",rtol=1e-12); print("backward ok",s.success,len(s.t), s.t[-1]/tg)
except Exception as e: print("backward fail:",repr(e)[:200])
# BL input
ob=bh.orbit(R=3000*rg,Theta=1.0*u.rad,Phi=0.5*u.rad,vR=0.001*c,vTheta=1e-6*u.rad/tg,vPhi=2e-6*u.rad/tg)
s=ob.solve(tau_eval=np.linspace(0,1e3,3)*tg,method="dop853",rtol=1e-12)
print("BL ok R",(s.R/rg).decompose()[:2], "vR",(s.vR/c).decompose()[:2],"vTh",(s.vTheta*tg).to(u.rad)[:2],"ut",s.ut[:2])
