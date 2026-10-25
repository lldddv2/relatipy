import numpy as np, time
from astropy import units as u
from astropy.constants import G, c
from relatipy import KerrMcmcModel

M = 4.297e6*u.M_sun; R0 = 8277*u.pc
a_S2 = 1035.45*u.au
m = KerrMcmcModel()
base = dict(
    kerr=dict(mass=M, spin=0.5, vec=(0.,0.,1.)),
    orbit=dict(a=a_S2, e=0.8843, inc=134.567*u.deg, Omega=228.171*u.deg, omega=66.263*u.deg, f=0*u.deg),
    distance=R0,
)
def bench(n_ep, rtol, method, span=16.0, reps=5, tkind="t_eval"):
    t = np.linspace(0.0, span, n_ep)*u.yr
    sol = {tkind: t, "method": method, "rtol": rtol}
    try:
        m.get_ra_dec_vr(**base, sol=sol)
    except Exception as e:
        return f"FAIL {type(e).__name__}: {str(e)[:70]}"
    t0=time.perf_counter()
    for _ in range(reps): m.get_ra_dec_vr(**base, sol=sol)
    dt=(time.perf_counter()-t0)/reps
    return f"{dt*1e3:8.2f} ms   {1/dt:8.1f} eval/s"

print("=== n_epochs scan (rtol=1e-10, dop853, 1 S2 orbit=16yr) ===")
for n in [10,50,100,200,500,1000]:
    print(f"n={n:5d}  {bench(n,1e-10,'dop853')}")
print("\n=== rtol scan (n=200) ===")
for rt in [1e-13,1e-11,1e-10,1e-8,1e-6,1e-4]:
    print(f"rtol={rt:8.0e} dop853 {bench(200,rt,'dop853')}")
    print(f"rtol={rt:8.0e} radau  {bench(200,rt,'radau')}")
print("\n=== span scan (n=200, rtol=1e-10) — S2 has ~2 orbits of data (1992-2026) ===")
for sp in [1.0,4.0,16.0,34.0,100.0]:
    print(f"span={sp:6.1f} yr  {bench(200,1e-10,'dop853',span=sp)}")
print("\n=== t_obs (Romer inversion) vs t_eval, n=200 ===")
print("t_eval", bench(200,1e-10,'dop853',tkind='t_eval'))
print("t_obs ", bench(200,1e-10,'dop853',tkind='t_obs'))
