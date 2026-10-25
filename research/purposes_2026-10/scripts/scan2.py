"""Irreducible (marginalised) frame-dragging signal over (a, e).  Inject exact Kerr
a*=0.9; refit exact Kerr a*=0 with all 8 params free.  20-yr campaign, first
pericentre at t=2 yr, 4 epochs/yr + 30-epoch dense campaign at every pericentre."""
import numpy as np, time, json, warnings, sys
from astropy import units as u
from relatipy import KerrMcmcModel
from scipy.optimize import least_squares
from common import epochs, f_before_peri
warnings.filterwarnings("ignore")
M0, R00 = 4.297, 8.277
RG_AU = 0.042414091583065465
m = KerrMcmcModel(); RAD = np.pi/180

def model(th, t, spin, vec, rtol=1e-12):
    Mm,R0,aa,e,inc,Om,om,f0 = th
    try:
        al,de,vr = m.get_ra_dec_vr(kerr=dict(mass=Mm*1e6*u.M_sun, spin=spin, vec=tuple(vec)),
            orbit=dict(a=aa*u.au, e=e, inc=inc*u.rad, Omega=Om*u.rad, omega=om*u.rad, f=f0*u.rad),
            sol={"t_obs": t*u.yr, "method":"dop853", "rtol":rtol}, distance=R0*u.kpc)
    except Exception: return None
    return np.concatenate([al*1e6, de*1e6, vr])

def irreducible(a_au, e, spin=0.9, vec=(0.,0.707,0.707), baseline=20.0, fit_spin=0.0,
                inc=134.567, Om=228.171, om=66.263):
    P = float(np.sqrt(a_au**3/(M0*1e6)))
    t1 = min(2.0, 0.5*P)
    t = epochs(P, baseline, t1)
    f0 = f_before_peri(e, P, t1)
    p_true=[M0,R00,a_au,e,inc*RAD,Om*RAD,om*RAD,f0]
    truth = model(p_true,t,spin,vec)
    if truth is None: return None
    n=len(t); w=np.concatenate([np.full(2*n,1/30.),np.full(n,1/10.)])
    def resid(th):
        y=model(th,t,fit_spin,(0.,0.,1.))
        return np.full(3*n,1e7) if y is None else (y-truth)*w
    lo=[2.,7.,a_au*.5,max(e-.3,.01),-np.pi,-2*np.pi,-2*np.pi,-2*np.pi]
    hi=[7.,9.5,a_au*1.5,min(e+.3,.995),np.pi,2*np.pi,2*np.pi,2*np.pi]
    r=least_squares(resid,np.array(p_true),bounds=(lo,hi),xtol=1e-15,ftol=1e-15,gtol=1e-15,
        x_scale=[1e-3,1e-3,a_au*1e-5,1e-6,1e-6,1e-6,1e-6,1e-6])
    ast=r.fun[:2*n]*30.; rv=r.fun[2*n:]*10.
    return dict(a_au=a_au,e=e,P=P,rperi_rg=a_au*(1-e)/RG_AU,n=n,
        rms=float(np.sqrt((ast**2).mean())),mx=float(np.abs(ast).max()),
        rv_rms=float(np.sqrt((rv**2).mean())),rv_mx=float(np.abs(rv).max()),
        dM=float(r.x[0]-M0), dR0=float(r.x[1]-R00), nfev=int(r.nfev))

if __name__=="__main__":
    print(f"{'a[AU]':>7}{'e':>6}{'P[yr]':>8}{'rp[rg]':>8}{'n':>5} |{'ast_rms':>9}{'ast_max':>9}{'rv_rms':>8}{'rv_max':>8}{'dM/1e6':>10}{'nfev':>5}",flush=True)
    rows=[]; t0=time.perf_counter()
    for a_au in [1035.45, 600., 400., 250., 150., 80., 40.]:
        for e in [0.3,0.6,0.8,0.9,0.95,0.976]:
            if a_au*(1-e)/RG_AU < 30: continue
            d=irreducible(a_au,e)
            if d is None: print(f"{a_au:7.1f}{e:6.2f} FAILED",flush=True); continue
            rows.append(d)
            print(f"{d['a_au']:7.1f}{d['e']:6.3f}{d['P']:8.2f}{d['rperi_rg']:8.0f}{d['n']:5d} |"
                  f"{d['rms']:9.3f}{d['mx']:9.3f}{d['rv_rms']:8.3f}{d['rv_mx']:8.3f}{d['dM']:+10.2e}{d['nfev']:5d}",flush=True)
            json.dump(rows,open("scan2.json","w"),indent=1)
    print(f"[{time.perf_counter()-t0:.0f}s]")
