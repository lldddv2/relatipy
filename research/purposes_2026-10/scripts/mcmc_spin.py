"""Exact-Kerr MCMC forecast: S2 + S62-like star, sigma_ast=10 uas, sigma_rv=5 km/s,
20-yr campaign.  17 params: M, R0, a*, cos(th_s), ph_s + 6 elements per star.
Truth a*=0.9, spin axis th=60 ph=30 deg (same as Fisher run)."""
import numpy as np, emcee, time, sys, json, warnings
from multiprocessing import Pool
from astropy import units as u
from relatipy import KerrMcmcModel
from fisher import setup, M0, R00, RAD
warnings.filterwarnings("ignore")
SA, SV = 10.0, 5.0
STAR_NAMES = sys.argv[1].split(",") if len(sys.argv)>1 else ["S2","S62-like"]
NSTEP = int(sys.argv[2]) if len(sys.argv)>2 else 6000
TAG = sys.argv[3] if len(sys.argv)>3 else "S2_S62"
_m=None
def obs(g,s,t):
    global _m
    if _m is None: _m=KerrMcmcModel()
    Mm,R0,ast,cth,ph=g; aa,e,inc,Om,om,f0=s
    sth=np.sqrt(max(0.0,1-cth*cth)); vec=(sth*np.cos(ph),sth*np.sin(ph),cth)
    al,de,vr=_m.get_ra_dec_vr(kerr=dict(mass=Mm*1e6*u.M_sun,spin=ast,vec=vec),
        orbit=dict(a=aa*u.au,e=e,inc=inc*u.rad,Omega=Om*u.rad,omega=om*u.rad,f=f0*u.rad),
        sol={"t_obs":t*u.yr,"method":"dop853","rtol":1e-11},distance=R0*u.kpc)
    return np.concatenate([al*1e6,de*1e6,vr])

rng=np.random.default_rng(42)
g_true=[M0,R00,0.9,np.cos(60*RAD),30*RAD]
SETUPS=[setup(n) for n in STAR_NAMES]
DATA=[]
for s,t in SETUPS:
    y=obs(g_true,s,t); n=len(t)
    sig=np.concatenate([np.full(2*n,SA),np.full(n,SV)])
    DATA.append((y+rng.normal(0,1,y.size)*sig, sig, t))
NS=len(SETUPS)
def logp(th):
    Mm,R0,ast,cth,ph=th[:5]
    if not (0<=ast<=1 and -1<=cth<=1 and 2<Mm<7 and 7<R0<9.5): return -np.inf
    ph=ph%(2*np.pi); lp=0.0
    for k,(yd,sig,t) in enumerate(DATA):
        s=th[5+6*k:11+6*k]
        if not (0<s[1]<1 and s[0]>0): return -np.inf
        try: y=obs([Mm,R0,ast,cth,ph],s,t)
        except Exception: return -np.inf
        lp+=-0.5*np.sum(((y-yd)/sig)**2)
    return lp

if __name__=="__main__":
    p_true=np.array(g_true+sum([list(s) for s,_ in SETUPS],[]))
    ndim=len(p_true); nw=max(2*ndim+4,40)
    scale=np.array([1e-4,1e-4,1e-3,1e-3,1e-3]+[1e-5]*6*NS)*np.maximum(np.abs(p_true),1.0)
    p0=p_true+scale*rng.normal(size=(nw,ndim)); p0[:,2]=np.clip(p0[:,2],0,0.999)
    t0=time.perf_counter(); print(f"ndim={ndim} walkers={nw} steps={NSTEP} logp(true)={logp(p_true):.1f}",flush=True)
    backend=emcee.backends.HDFBackend(f"chain_{TAG}.h5") if False else None
    with Pool(3) as pool:
        sam=emcee.EnsembleSampler(nw,ndim,logp,pool=pool,moves=[(emcee.moves.DEMove(),0.8),(emcee.moves.DESnookerMove(),0.2)])
        for i,_ in enumerate(sam.sample(p0,iterations=NSTEP)):
            if (i+1)%500==0:
                ch=sam.get_chain()
                try: tau=emcee.autocorr.integrated_time(ch[ch.shape[0]//2:],quiet=True)[:5]
                except Exception: tau=[np.nan]*5
                c=ch[ch.shape[0]//2:,:,:].reshape(-1,ndim)
                print(f"step {i+1}  acc={sam.acceptance_fraction.mean():.2f}  {time.perf_counter()-t0:.0f}s  "
                      f"a*={np.median(c[:,2]):.3f}+-{np.std(c[:,2]):.3f}  cos_th={np.median(c[:,3]):.3f}+-{np.std(c[:,3]):.3f}  "
                      f"ph={np.median(c[:,4])/RAD:.1f}+-{np.std(c[:,4])/RAD:.1f}deg  tau(glob)={np.round(tau,0)}",flush=True)
                np.save(f"chain_{TAG}.npy",ch); np.save(f"lnp_{TAG}.npy",sam.get_log_prob())
