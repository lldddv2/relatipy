"""Same benchmark on ACCEPTED steps only (no tau_eval interpolation), plus failure messages."""
import numpy as np, time, sys, json, warnings
warnings.filterwarnings("ignore")
from astropy import units as u
from relatipy import Kerr
from relatipy.coordinates import KerrOrbitalElements
import kerrgeopy as kg
A,P,E,X = float(sys.argv[1]),float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]); NORB=int(sys.argv[5])
meths = sys.argv[6].split(",") if len(sys.argv)>6 else ["dop853","radau","projection_radau","dp45"]
rtols = [float(x) for x in sys.argv[7].split(",")] if len(sys.argv)>7 else [1e-10,1e-12]
bh=Kerr(mass=1*u.Msun, spin=A); rg=bh.r_g; T0=(bh.r_g/(299792458*u.m/u.s)).to(u.s)
el=KerrOrbitalElements(p=P*rg,e=E,x=X,q_r0=0*u.rad,q_theta0=0*u.rad,q_phi0=0*u.rad)
ko=kg.StableOrbit(A,P,E,X); tf,rf,thf,phf=ko.trajectory()
Tr=2*np.pi/kg.fundamental_frequencies(A,P,E,X)[0]
s0=bh.orbit(elements=el).solve(tau_span=(0*u.s,50*Tr*T0),method="dop853",rtol=1e-12)
ratio=(s0.t[-1]/s0.tau[-1]).to_value(u.dimensionless_unscaled); tau_end=NORB*Tr/ratio
def kg_at_t(tt):
    lam=tt/ko.gamma
    for _ in range(14):
        h=1e-3; d=(tf(lam+h)-tf(lam-h))/(2*h); lam=lam-(tf(lam)-tt)/d
    return rf(lam),thf(lam),phf(lam)
print(f"a*={A} p={P} e={E} x={X} N_orb={NORB} (accepted steps only)")
print(f"{'method':<17}{'rtol':>7}{'wall[s]':>8}{'steps':>9}{'dE/E':>9}{'dLz/Lz':>9}{'dQ/Q':>9}{'|dr|max':>9}{'|dphi|max':>10}{'|dth|max':>9}")
for meth in meths:
    for rt in rtols:
        t0=time.perf_counter()
        s=bh.orbit(elements=el).solve(tau_span=(0*u.s,tau_end*T0),method=meth,rtol=rt)
        w=time.perf_counter()-t0
        if not s.success:
            print(f"{meth:<17}{rt:7.0e}{w:8.2f}  status={s.status} msg='{s.message}' reached tau={float((s.tau[-1]/T0).decompose()):.3e}/{tau_end:.3e}",flush=True); continue
        rel=lambda q: float(np.max(np.abs(np.asarray(q)-np.asarray(q)[0]))/abs(np.asarray(q)[0]))
        dE,dL,dQ=rel(s.get_E()),rel(s.get_Lz()),rel(s.get_Q())
        idx=np.linspace(0,len(s.tau)-1,4000).astype(int)
        tt=(s.t[idx]/T0).to_value(u.dimensionless_unscaled)
        Rr=(s.R[idx]/rg).to_value(u.dimensionless_unscaled) if s.R.unit!=u.dimensionless_unscaled else s.R.value[idx]
        rk,thk,phk=kg_at_t(tt)
        dr=np.max(np.abs(Rr-rk)); dth=np.max(np.abs(s.Theta.to_value(u.rad)[idx]-thk)); dph=np.max(np.abs(s.Phi.to_value(u.rad)[idx]-phk))
        print(f"{meth:<17}{rt:7.0e}{w:8.2f}{s.integration.n_steps:9d}{dE:9.1e}{dL:9.1e}{dQ:9.1e}{dr:9.1e}{dph:10.1e}{dth:9.1e}",flush=True)
