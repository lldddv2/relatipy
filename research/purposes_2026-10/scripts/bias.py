"""Truth = exact Kerr (relatipy).  Model = PN (numba).  Same projection + Romer.
Fit ALL nuisance params (6 elements, mu, kappa, dt) and optionally f_SP or a*.
Reports residual floor and BIAS on the physics parameter.  Astrometry only."""
import numpy as np, sys, time, json
from astropy import units as u
from astropy.constants import c
from relatipy import Kerr
from scipy.optimize import least_squares
import pn, pn_fast
from common import epochs, f_before_peri
UAS=5.124331470709855; RG_AU=0.042414091583065465
M=4.297e6*u.M_sun
RG_S=(Kerr(mass=M,spin=0.0).r_g/c).to_value(u.s); RG_PER_YR=(1*u.yr).to_value(u.s)/RG_S
R=np.pi/180

def exact(astar, x0, v0, t_end, n=40000):
    bh=Kerr(mass=M,spin=astar)
    orb=bh.orbit(x=x0[0]*bh.r_g,y=x0[1]*bh.r_g,z=x0[2]*bh.r_g,vx=v0[0]*c,vy=v0[1]*c,vz=v0[2]*c)
    tau=np.linspace(0,t_end*1.05*RG_S,n)*u.s
    s=orb.solve(tau_span=(tau[0],tau[-1]),tau_eval=tau,method="dop853",rtol=1e-13)
    assert s.success
    return (s.t.to_value(u.s)/RG_S, np.ascontiguousarray((s.xyz/bh.r_g).to_value(u.dimensionless_unscaled)),
            np.ascontiguousarray((s.vxyz/c).to_value(u.dimensionless_unscaled)))

def basis(ghat):
    ref=np.array([1.,0.,0.]) if abs(ghat[0])<0.9 else np.array([0.,1.,0.])
    e1=ref-(ref@ghat)*ghat; e1/=np.linalg.norm(e1); return e1, np.cross(ghat,e1)

def run(name, a_au, e, astar, fits, inc=134.567, Om=228.171, om=66.263, baseline=20.0, los=50.0, spo=1500):
    a_rg=a_au/RG_AU; P_yr=2*np.pi*a_rg**1.5/RG_PER_YR; t1=min(2.0,0.5*P_yr)
    f0=f_before_peri(e,P_yr,t1)
    x0,v0=pn.kepler_to_cart(a_rg,e,inc*R,Om*R,om*R,f0)
    t_end=baseline*RG_PER_YR; tobs=epochs(P_yr,baseline*0.97,t1)*RG_PER_YR
    ghat=np.array([0.,np.sin(los*R),np.cos(los*R)]); e1,e2=basis(ghat)
    T,X,V=exact(astar,x0,v0,t_end)
    ytrue=pn_fast.sky(T,X,V,ghat,e1,e2,tobs)[0].ravel()*UAS
    ds=2*np.pi/spo
    print(f"\n==== {name}: a={a_au} AU ({a_rg:.0f} rg) e={e} rp={a_rg*(1-e):.0f} rg P={P_yr:.2f} yr  truth a*={astar}  n_ep={len(tobs)}",flush=True)
    print(f"  {'model':<26}{'rms[uas]':>9}{'max[uas]':>9}{'phys-param':>22}{'nfev':>6}",flush=True)
    out={}
    for lab,order,free in fits:
        base=[a_rg,e,inc*R,Om*R,om*R,f0,1.0,1.0,0.0]
        lo=[a_rg*.9,max(e-.05,.01),(inc-5)*R,(Om-5)*R,(om-5)*R,f0-.3,.8,.8,-.02*t_end]
        hi=[a_rg*1.1,min(e+.05,.9995),(inc+5)*R,(Om+5)*R,(om+5)*R,f0+.3,1.2,1.2,.02*t_end]
        xs=[a_rg*1e-5,1e-7,1e-7,1e-7,1e-7,1e-7,1e-6,1e-7,t_end*1e-7]
        if free=="fsp": base+= [1.0]; lo+=[-2.0]; hi+=[4.0]; xs+=[1e-3]
        if free=="spin": base+=[astar if astar>0 else 0.5]; lo+=[-3.0]; hi+=[3.0]; xs+=[1e-2]
        Sfix=np.array([0.,0.,astar if order>=2 else 0.0])
        def resid(th):
            ar,ee,i_,O_,w_,ff,mu,kap,dt=th[:9]
            fsp=th[9] if free=="fsp" else 1.0
            S=np.array([0.,0.,th[9]]) if free=="spin" else Sfix
            xx,vv=pn.kepler_to_cart(ar,ee,i_,O_,w_,ff)
            Tm,Xm,Vm=pn_fast.integrate(xx,vv,t_end*1.03/mu,order,fsp,S,ds)
            y=pn_fast.sky(Tm*mu,Xm*mu,Vm,ghat,e1,e2,tobs+dt)[0].ravel()*UAS*kap
            return y-ytrue
        t0=time.perf_counter()
        r=least_squares(resid,np.array(base),bounds=(lo,hi),xtol=1e-15,ftol=1e-15,gtol=1e-15,x_scale=xs,max_nfev=3000)
        rms=float(np.sqrt((r.fun**2).mean())); mx=float(np.abs(r.fun).max())
        pp=""
        if free=="fsp": pp=f"f_SP={r.x[9]:.5f}"
        if free=="spin": pp=f"a*={r.x[9]:.4f} (true {astar})"
        out[lab]=dict(rms=rms,mx=mx,x=r.x.tolist())
        print(f"  {lab:<26}{rms:9.3f}{mx:9.3f}{pp:>22}{r.nfev:6d}  [{time.perf_counter()-t0:.0f}s]",flush=True)
    return out

if __name__=="__main__":
    STARS={"S2":(1035.45,0.8843),"S62-like":(749.0,0.976),"S4716-like":(411.0,0.756),"rp300":(400.0,0.968)}
    allr={}
    for nm in sys.argv[1:]:
        a,e=STARS[nm]
        allr[(nm,0.0)]=run(nm,a,e,0.0,[("Newtonian",0,None),("1PN",1,None),("1PN x f_SP (free)",1,"fsp")])
        allr[(nm,0.9)]=run(nm,a,e,0.9,[("1PN (no spin)",1,None),("1PN+1.5PN a*=0.9 fixed",2,None),("1PN+1.5PN a* free",2,"spin")])
        json.dump({f"{k[0]}|{k[1]}":v for k,v in allr.items()},open(f"bias_{'_'.join(sys.argv[1:])}.json","w"),indent=1)
