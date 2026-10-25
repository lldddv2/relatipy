"""Pick the correct 1.5PN spin-orbit sign/ordering by comparing to the EXACT Kerr
geodesic from relatipy in the weak field, and to analytic Lense-Thirring rates."""
import numpy as np
from astropy import units as u
from astropy.constants import G, c
from relatipy import Kerr
import pn

M = 4.297e6*u.M_sun
bh0 = Kerr(mass=M, spin=0.0)
RG = bh0.r_g.to(u.m).value
T0 = (bh0.r_g/c).to(u.s).value

a_rg, e, inc, Om, om = 3000.0, 0.3, 1.0, 0.7, 0.4
astar = 0.9
S = astar*np.array([0.0,0.0,1.0])       # spin along +z, J = a* in G=M=c=1

# ---- analytic Lense-Thirring per-orbit rates (Will 2008 / Merritt+2010 form)
P_rg = 2*np.pi*a_rg**1.5
dOm_per_orbit = 4*np.pi*astar/(a_rg**1.5*(1-e**2)**1.5)
dom_per_orbit = -12*np.pi*astar*np.cos(inc)/(a_rg**1.5*(1-e**2)**1.5)
dom_schw = 6*np.pi/(a_rg*(1-e**2))
print(f"a={a_rg} rg  e={e}  inc={inc:.3f}  a*={astar}  P={P_rg:.4g} rg/c")
print(f"analytic per orbit: dOmega_LT={dOm_per_orbit:+.6e} rad   "
      f"domega_LT={dom_per_orbit:+.6e} rad   domega_Schw={dom_schw:+.6e} rad")

x0, v0 = pn.kepler_to_cart(a_rg, e, inc, Om, om, 0.0)
N_ORB = 20
t_end = N_ORB*P_rg

# ---- EXACT Kerr via relatipy
bh = Kerr(mass=M, spin=astar)
orb = bh.orbit(x=x0[0]*bh.r_g, y=x0[1]*bh.r_g, z=x0[2]*bh.r_g,
               vx=v0[0]*c, vy=v0[1]*c, vz=v0[2]*c)
tau = np.linspace(0, t_end*T0, 6000)*u.s
sol = orb.solve(tau_span=(tau[0],tau[-1]), tau_eval=tau, method="dop853", rtol=1e-13)
print("relatipy status", sol.status, "success", sol.success)
Xk = (sol.xyz/bh.r_g).to_value(u.dimensionless_unscaled)
Vk = (sol.vxyz/c).to_value(u.dimensionless_unscaled)
ak,ek,ick,Omk,omk = pn.elements_from_state(Xk,Vk)
def secular(ang, t):
    aa = np.unwrap(ang)
    A = np.polyfit(t, aa, 1)
    return A[0]        # rad per unit t
tk = sol.tau.to_value(u.s)/T0
rate_Om_k = secular(Omk, tk); rate_om_k = secular(omk, tk)
print(f"EXACT KERR  : dOmega/dt={rate_Om_k:+.6e}  -> per orbit {rate_Om_k*P_rg:+.6e}")
print(f"              domega/dt={rate_om_k:+.6e}  -> per orbit {rate_om_k*P_rg:+.6e}"
      f"   (expect {dom_schw+dom_per_orbit:+.6e})")

for var in [1,2,3,4,5]:
    t,X,V = pn.integrate(2, S, x0, v0, t_end, n_out=6000, variant=var)
    a_,e_,i_,O_,w_ = pn.elements_from_state(X,V)
    rO = secular(O_,t)*P_rg; rw = secular(w_,t)*P_rg
    print(f"PN variant {var}: dOmega/orbit={rO:+.6e} (err {rO-dOm_per_orbit:+.2e})"
          f"  domega/orbit={rw:+.6e} (err {rw-(dom_schw+dom_per_orbit):+.2e})")

t,X,V = pn.integrate(1, [0,0,0], x0, v0, t_end, n_out=6000)
a_,e_,i_,O_,w_ = pn.elements_from_state(X,V)
print(f"PN 1PN only : dOmega/orbit={secular(O_,t)*P_rg:+.6e}  domega/orbit={secular(w_,t)*P_rg:+.6e}"
      f"  (expect 0 and {dom_schw:+.6e})")
