"""Fisher forecast for the Sgr A* spin VECTOR (a*, theta_s, phi_s) from exact Kerr
geodesics, star by star and jointly.  Globals: M, R0, a*, theta_s, phi_s.
Per star: a, e, inc, Omega, omega, f0.  Central differences on the exact model."""
import numpy as np, json, sys, warnings, time
from multiprocessing import Pool
from astropy import units as u
from relatipy import KerrMcmcModel
from common import epochs, f_before_peri
warnings.filterwarnings("ignore")
RAD=np.pi/180; M0,R00=4.297,8.277
_m=None
def mdl():
    global _m
    if _m is None: _m=KerrMcmcModel()
    return _m

def obs(g, s, t):
    Mm,R0,ast,th,ph = g
    aa,e,inc,Om,om,f0 = s
    vec=(np.sin(th)*np.cos(ph), np.sin(th)*np.sin(ph), np.cos(th))
    al,de,vr = mdl().get_ra_dec_vr(kerr=dict(mass=Mm*1e6*u.M_sun, spin=ast, vec=vec),
        orbit=dict(a=aa*u.au,e=e,inc=inc*u.rad,Omega=Om*u.rad,omega=om*u.rad,f=f0*u.rad),
        sol={"t_obs":t*u.yr,"method":"dop853","rtol":RTOL}, distance=R0*u.kpc)
    return np.concatenate([al*1e6, de*1e6, vr])

RTOL=1e-12
HG=[1e-5,1e-5,2e-4,1e-4,1e-4]; HS=[1e-6,1e-6,1e-6,1e-6,1e-6,1e-6]
def dcol(fp, fm, hh):
    for mult in (1.0, 2.7, 1/2.7, 7.3):
        try: return (fp(hh*mult)-fm(hh*mult))/(2*hh*mult)
        except Exception: pass
    raise RuntimeError("derivative failed at all steps")

def jac_star(args):
    g,s,t = args
    cols=[]
    for i,h in enumerate(HG):
        gp=list(g); gm=list(g)
        hh = h*(abs(g[i]) if i<2 else 1.0)
        def fp(d,i=i): q=list(g); q[i]+=d; return obs(q,s,t)
        def fm(d,i=i): q=list(g); q[i]-=d; return obs(q,s,t)
        cols.append(dcol(fp,fm,hh))
    for i,h in enumerate(HS):
        sp=list(s); sm=list(s)
        hh = h*(abs(s[i]) if i==0 else 1.0)
        def fp(d,i=i): q=list(s); q[i]+=d; return obs(g,q,t)
        def fm(d,i=i): q=list(s); q[i]-=d; return obs(g,q,t)
        cols.append(dcol(fp,fm,hh))
    return np.array(cols).T     # (3n, 11)

STARS={  # a[AU], e, inc, Omega, omega (deg)
 "S2":        (1035.45,0.8843,134.567,228.171,66.263),
 "S62-like":  (749.0,  0.976, 72.8,  122.6,  42.6),
 "S4716-like":(411.0,  0.756, 161.2,  10.0,   0.0),
 "S29":       (3500.0, 0.97,  144.0,   7.0,  205.0),
}
def setup(name, baseline=20.0):
    a,e,i,O,w = STARS[name]
    P=float(np.sqrt(a**3/(M0*1e6))); t1=min(2.0,0.5*P)
    t=epochs(P,baseline,t1); f0=f_before_peri(e,P,t1)
    return [a,e,i*RAD,O*RAD,w*RAD,f0], t

if __name__=="__main__":
    ast_true=float(sys.argv[1]) if len(sys.argv)>1 else 0.9
    th_s, ph_s = 60*RAD, 30*RAD
    g=[M0,R00,ast_true,th_s,ph_s]
    names=list(STARS)
    t0=time.perf_counter()
    with Pool(2) as pool:
        setups=[setup(n) for n in names]
        J=pool.map(jac_star,[(g,s,t) for s,t in setups])
    print(f"[jacobians {time.perf_counter()-t0:.0f}s]  a*_true={ast_true}  axis th=60 ph=30 deg",flush=True)
    res={}
    for sa,sv in [(50.,10.),(30.,10.),(10.,5.),(3.,2.)]:
        print(f"\n### sigma_ast={sa} uas  sigma_rv={sv} km/s  (20-yr campaign)")
        print(f"{'stars':<28}{'s(a*)':>9}{'s(th)deg':>10}{'s(ph)deg':>10}{'s(M)/M':>10}{'s(R0)/R0':>10}")
        combos=[[n] for n in names]+[["S2","S62-like"],["S2","S4716-like","S29"],names]
        for combo in combos:
            ng=5; ns=6*len(combo); F=np.zeros((ng+ns,ng+ns))
            for k,nm in enumerate(combo):
                Jk=J[names.index(nm)]; n=Jk.shape[0]//3
                W=np.concatenate([np.full(2*n,1/sa**2),np.full(n,1/sv**2)])
                Fk=Jk.T@(Jk*W[:,None])
                idx=list(range(5))+list(range(ng+6*k,ng+6*k+6))
                F[np.ix_(idx,idx)]+=Fk
            try:
                C=np.linalg.inv(F); sd=np.sqrt(np.clip(np.diag(C),0,None))
            except np.linalg.LinAlgError:
                sd=np.full(ng+ns,np.inf)
            lab="+".join(combo)
            res[(sa,lab)]=sd[:5].tolist()
            print(f"{lab:<28}{sd[2]:9.3f}{sd[3]/RAD:10.1f}{sd[4]/RAD:10.1f}{sd[0]/M0:10.2e}{sd[1]/R00:10.2e}")
    json.dump({f"{k[0]}|{k[1]}":v for k,v in res.items()},open(f"fisher_{ast_true}.json","w"),indent=1)
