"""Failure-rate census of KerrMcmcModel under small random parameter perturbations."""
import numpy as np, warnings, sys
warnings.filterwarnings("ignore")
from astropy import units as u
from relatipy import KerrMcmcModel
from fisher import setup, STARS, M0, R00, RAD
m=KerrMcmcModel(); rng=np.random.default_rng(1)
def ev(g,s,t,tk,rt,meth="dop853"):
    Mm,R0,ast,th,ph=g; aa,e,inc,Om,om,f0=s
    vec=(np.sin(th)*np.cos(ph), np.sin(th)*np.sin(ph), np.cos(th))
    try:
        m.get_ra_dec_vr(kerr=dict(mass=Mm*1e6*u.M_sun,spin=ast,vec=vec),
          orbit=dict(a=aa*u.au,e=e,inc=inc*u.rad,Omega=Om*u.rad,omega=om*u.rad,f=f0*u.rad),
          sol={tk:t*u.yr,"method":meth,"rtol":rt}, distance=R0*u.kpc); return 0
    except Exception: return 1
N=int(sys.argv[1]) if len(sys.argv)>1 else 150
print(f"{'star':12s}{'mode':8s}{'meth':8s}{'rtol':>7}  fail/N   (random spin axis, a* U[0,1], params x(1+1e-3 N(0,1)))")
for n in STARS:
    s0,t=setup(n)
    for tk in ["t_obs","t_eval"]:
        for meth,rt in [("dop853",1e-12),("dop853",1e-10),("radau",1e-10)]:
            fails=0
            rng=np.random.default_rng(7)
            for k in range(N):
                g=[M0*(1+1e-3*rng.normal()),R00*(1+1e-3*rng.normal()),rng.uniform(0,1),
                   np.arccos(rng.uniform(-1,1)),rng.uniform(0,2*np.pi)]
                s=[s0[0]*(1+1e-3*rng.normal()),min(s0[1]+1e-3*rng.normal(),0.999)]+[x+1e-3*rng.normal() for x in s0[2:]]
                fails+=ev(g,s,t,tk,rt,meth)
            print(f"{n:12s}{tk:8s}{meth:8s}{rt:7.0e}  {fails:3d}/{N}",flush=True)
