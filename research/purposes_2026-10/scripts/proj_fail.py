import numpy as np, warnings; warnings.filterwarnings("ignore")
from astropy import units as u
from relatipy import Kerr
from relatipy.coordinates import KerrOrbitalElements
import kerrgeopy as kg
A,P,E,X=0.9,10,0.5,0.5
bh=Kerr(mass=1*u.Msun,spin=A); rg=bh.r_g; T0=(bh.r_g/(299792458*u.m/u.s)).to(u.s)
el=KerrOrbitalElements(p=P*rg,e=E,x=X,q_r0=0*u.rad,q_theta0=0*u.rad,q_phi0=0*u.rad)
rp,ra=P/(1+E),P/(1-E); thmin=np.arccos(np.sqrt(1-X**2)); thmax=np.pi-thmin
print(f"allowed r in [{rp:.6f},{ra:.6f}]  theta in [{thmin:.6f},{thmax:.6f}]")
for rt in [1e-10,1e-9,1e-8,1e-11]:
    s=bh.orbit(elements=el).solve(tau_span=(0*u.s,3.48e5*T0),method="projection_radau",rtol=rt)
    R=(s.R/rg).to_value(u.dimensionless_unscaled) if s.R.unit!=u.dimensionless_unscaled else s.R.value
    th=s.Theta.to_value(u.rad)
    print(f"rtol={rt:.0e} status={s.status} last r={R[-1]:.6f} (r-rp={R[-1]-rp:+.2e}, ra-r={ra-R[-1]:+.2e})  "
          f"th={th[-1]:.6f} (th-thmin={th[-1]-thmin:+.2e}, thmax-th={thmax-th[-1]:+.2e})  n={len(R)}  tau_end={float((s.tau[-1]/T0).decompose()):.3e}")
