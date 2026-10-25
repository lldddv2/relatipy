import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
M=1.0
def u_at(b,psi):
    f=lambda ph,y:[y[1],-y[0]+3*M*y[0]**2]
    s=solve_ivp(f,[0,psi],[0,1/b],rtol=1e-12,atol=1e-15)
    return s.y[0,-1],s.y[1,-1]
print("astrometry: exact b vs BL-proj, harmonic-proj, weak-field formula (units M)")
for re in [280.,2800.]:
  for psi in np.radians([30,60,90,120,150]):
    g=lambda b:u_at(b,psi)[0]-1/re
    b=brentq(g,0.3*re*np.sin(psi),3*re*np.sin(psi)+50,xtol=1e-12)
    bl=re*np.sin(psi); hm=(re-M)*np.sin(psi); wf=hm+2*M*np.tan(psi/2)
    print(f"r={re:6.0f} psi={np.degrees(psi):4.0f}  b-BL={b-bl:7.3f}  b-harm={b-hm:7.3f}  b-weakfield={b-wf:8.4f}")
# redshift at pericentre, S301-like, velocity along t-hat (in photon plane)
print("RV at pericentre r=280, e=0.982 (km/s), velocity in photon plane")
c=299792.458; rp=280.; e=0.982; p=rp*(1+e)
L=np.sqrt(p**2/(p-3-e**2)); E=np.sqrt((p-2-2*e)*(p-2+2*e)/(p*(p-3-e**2)))
ut=E/(1-2/rp); uph=L/rp**2; lap=np.sqrt(1-2/rp)
vloc=rp*uph/(ut*lap); gam=1/np.sqrt(1-vloc**2)
print("v_loc=",vloc*c)
for psi in np.radians([30,60,90,120,150]):
    g=lambda b:u_at(b,psi)[0]-1/rp
    b=brentq(g,0.3*rp*np.sin(psi),3*rp*np.sin(psi)+50,xtol=1e-12)
    up=u_at(b,psi)[1]
    sa=b*lap/rp; ca=np.sign(up)*np.sqrt(1-sa**2)
    # star velocity along t-hat (toward observer side) magnitude vloc; n_gamma=(ca,sa)
    exact=gam*(1-vloc*sa)/lap
    stat=gam*(1-vloc*np.sin(psi))/lap
    rela=ut - rp*uph*np.sin(psi)
    pred_rel=-(2/rp)*gam*vloc*np.tan(psi/2)*np.sin(psi/2)**2
    print(f"psi={np.degrees(psi):4.0f} static-exact={c*(stat-exact):7.1f} relatipy-exact={c*(rela-exact):7.1f} relatipy-static={c*(rela-stat):7.1f}  pred(exact-relatipy)={c*pred_rel:7.1f}")
# 2PN precession check
for (p,e,ra_uas,lab) in [(579.,0.982,1363*120.8,'S301'),(5330.,0.8847,1958*120.8,'S2'),(288.,0.2,None,'QPE')]:
    ex=quad(lambda x:np.sqrt(p/(p-6-2*e*np.cos(x))),0,2*np.pi,epsabs=1e-14,epsrel=1e-14)[0]-2*np.pi
    pn1=6*np.pi/p; pn2=1.5*np.pi*(18+e**2)/p**2
    print(lab,"dphi exact",ex," 1PN",pn1," 2PN",pn2," resid",ex-pn1-pn2, " 2PN uas@apo", None if ra_uas is None else pn2*ra_uas)
