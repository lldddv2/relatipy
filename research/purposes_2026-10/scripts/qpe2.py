"""QPE-like timing: orbiter around a spinning SMBH crosses the (spin-aligned) disk twice
per orbit.  Truth: exact Kerr (relatipy).  Model: 1PN+1.5PN (numba), a* free + all
orbital params + mass scale; BL gauge applied (r->r+1).  Observable: crossing times
with Romer delay.  Output: irreducible timing residual and a* bias vs number of orbits."""
import numpy as np, sys, time, warnings
warnings.filterwarnings("ignore")
from astropy import units as u
from astropy.constants import c
from relatipy import Kerr
from scipy.optimize import least_squares
from numba import njit
import pn, pn_fast
R=np.pi/180
A_RG, E, INC, ASTAR = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
NORB=int(sys.argv[5])
MSUN=4e5; RG_S=MSUN*4.925490947e-6   # r_g/c in seconds for 4e5 Msun
ghat=np.array([0., np.sin(60*R), np.cos(60*R)])

@njit(cache=True)
def crossings(T,X,V,ghat):
    out=[]
    for i in range(len(T)-1):
        z0=X[i,2]; z1=X[i+1,2]
        if z0==0.0 or z0*z1<0:
            h=T[i+1]-T[i]; s=0.5
            for it in range(50):   # Newton on cubic Hermite z(s)
                h00=2*s**3-3*s**2+1; h10=s**3-2*s**2+s; h01=-2*s**3+3*s**2; h11=s**3-s**2
                z=h00*z0+h10*h*V[i,2]+h01*z1+h11*h*V[i+1,2]
                d=((6*s**2-6*s)*z0+(3*s**2-4*s+1)*h*V[i,2]+(-6*s**2+6*s)*z1+(3*s**2-2*s)*h*V[i+1,2])
                ds=z/d; s-=ds
                if abs(ds)<1e-13: break
            h00=2*s**3-3*s**2+1; h10=s**3-2*s**2+s; h01=-2*s**3+3*s**2; h11=s**3-s**2
            los=0.0
            for k in range(3):
                xk=h00*X[i,k]+h10*h*V[i,k]+h01*X[i+1,k]+h11*h*V[i+1,k]; los+=xk*ghat[k]
            out.append(T[i]+s*h+los)
    return np.array(out)

def exact_times(x0,v0,t_end,rtol=1e-12):
    bh=Kerr(mass=MSUN*u.M_sun,spin=ASTAR); T0=(bh.r_g/c).to(u.s)
    orb=bh.orbit(x=x0[0]*bh.r_g,y=x0[1]*bh.r_g,z=x0[2]*bh.r_g,vx=v0[0]*c,vy=v0[1]*c,vz=v0[2]*c)
    t0=time.perf_counter()
    s=orb.solve(tau_span=(0*u.s,t_end*T0),method="dop853",rtol=rtol)
    w=time.perf_counter()-t0
    assert s.success, s.message
    T=(s.t/T0).to_value(u.dimensionless_unscaled); X=np.ascontiguousarray((s.xyz/bh.r_g).to_value(u.dimensionless_unscaled))
    V=np.ascontiguousarray((s.vxyz/c).to_value(u.dimensionless_unscaled))
    return crossings(T,X,V,ghat), w, s.integration.n_steps

x0,v0=pn.kepler_to_cart(A_RG,E,INC*R,0.3,0.5,0.0)
P=2*np.pi*A_RG**1.5; t_end=NORB*P*1.02
tc_true,w,nst=exact_times(x0,v0,t_end)
tc_true=tc_true[:2*NORB]
print(f"a={A_RG} rg e={E} inc={INC}deg a*={ASTAR} N_orb={NORB}  P={P*RG_S/3600:.2f} h (M=4e5 Msun)  "
      f"exact Kerr cost {w:.2f}s, {nst} steps, {len(tc_true)} crossings",flush=True)
spo=600; ds=2*np.pi/spo
P_ex=float(np.mean(np.diff(tc_true[::2]))); A0=(P_ex/(2*np.pi))**(2/3)
print(f"  measured crossing period {P_ex:.1f} rg/c vs Kepler {P:.1f} -> a0={A0:.3f} (input {A_RG})",flush=True)
def model(th, n):
    ar,ee,i_,O_,w_,f_,ast,mu,dt=th
    xx,vv=pn.kepler_to_cart(ar,ee,i_,O_,w_,f_)
    tn=(n/2+3)*P_ex*1.3/mu
    Tm,Xm,Vm=pn_fast.integrate_pre(xx,vv,tn,2,1.0,np.array([0.,0.,ast]),ds,int(spo*(n/2+5)*1.8))
    rr=np.sqrt((Xm**2).sum(1)); fac=((rr+1.0)/rr)[:,None]
    tc=crossings(Tm*mu,np.ascontiguousarray(Xm*fac*mu),np.ascontiguousarray(Vm*fac),ghat)+dt
    if len(tc)<n: return np.full(n,1e9)
    return tc[:n]
for nfit in [int(x) for x in sys.argv[6].split(',')]:
    n=2*nfit; y=tc_true[:n]
    p0=np.array([A0,E,INC*R,0.3,0.5,0.0,ASTAR,1.0,0.0])
    lo=[A0*.85,max(E-.05,0),(INC-3)*R,0.3-.1,0.5-.3,-.3,-2.,0.85,-2000.]
    hi=[A0*1.15,E+.05,(INC+3)*R,0.3+.1,0.5+.3,.3,2.,1.15,2000.]
    t0=time.perf_counter()
    r=least_squares(lambda th: model(th,n)-y, p0, bounds=(lo,hi), x_scale=[1e-3,1e-6,1e-6,1e-6,1e-6,1e-6,1e-3,1e-7,1e-2],
                    xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=400)
    res=r.fun
    print(f"  fit first {nfit:5d} orbits: rms={np.sqrt((res**2).mean()):.3e} rg/c = {np.sqrt((res**2).mean())*RG_S:9.3f} s,"
          f" max={np.abs(res).max()*RG_S:9.3f} s | a*_fit={r.x[6]:.5f} (true {ASTAR}) | nfev={r.nfev} {time.perf_counter()-t0:.0f}s",flush=True)
# unfitted model drift for scale: 1PN without spin vs truth after best fit of spin-free model on first 10%
