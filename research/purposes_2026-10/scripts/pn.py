"""Post-Newtonian reference integrators (geometric units G=c=M=1; lengths in r_g).
Validated against relatipy's exact Kerr geodesic in the weak field."""
import numpy as np
from scipy.integrate import solve_ivp

def kepler_to_cart(a, e, inc, Om, om, f, mu=1.0):
    """Classical elements -> Cartesian position/velocity (same frame convention
    as relatipy: right-handed, z along spin axis)."""
    p = a*(1-e**2)
    r = p/(1+e*np.cos(f))
    # perifocal
    xp, yp = r*np.cos(f), r*np.sin(f)
    h = np.sqrt(mu*p)
    vxp, vyp = -mu/h*np.sin(f), mu/h*(e+np.cos(f))
    cO,sO = np.cos(Om),np.sin(Om); ci,si = np.cos(inc),np.sin(inc)
    cw,sw = np.cos(om),np.sin(om)
    R = np.array([
      [cO*cw-sO*sw*ci, -cO*sw-sO*cw*ci,  sO*si],
      [sO*cw+cO*sw*ci, -sO*sw+cO*cw*ci, -cO*si],
      [sw*si,            cw*si,           ci   ]])
    return R@np.array([xp,yp,0.0]), R@np.array([vxp,vyp,0.0])

def rhs(order, Svec):
    S = np.asarray(Svec, float)
    def f(t, y):
        x, v = y[:3], y[3:]
        r = np.linalg.norm(x); n = x/r
        a = -n/r**2
        if order >= 1:
            vv = v@v; rdot = v@n
            a = a + (1.0/r**2)*(n*(4.0/r - vv) + 4.0*rdot*v)
        if order >= 2 and np.any(S):
            a = a + (2.0/r**3)*(3.0*(n@S)*np.cross(n, v) + np.cross(v, S))
        return np.concatenate([v, a])
    return f

def rhs_alt(order, Svec, variant):
    S = np.asarray(Svec, float)
    def f(t, y):
        x, v = y[:3], y[3:]
        r = np.linalg.norm(x); n = x/r
        a = -n/r**2
        if order >= 1:
            vv = v@v; rdot = v@n
            a = a + (1.0/r**2)*(n*(4.0/r - vv) + 4.0*rdot*v)
        if np.any(S):
            nxv = np.cross(n, v); vxS = np.cross(v, S)
            if   variant==1: dv = 3.0*(n@S)*nxv + vxS
            elif variant==2: dv = 3.0*(n@S)*nxv - vxS
            elif variant==3: dv = 3.0*n*(nxv@S) - vxS
            elif variant==4: dv = 3.0*n*(nxv@S) + vxS
            elif variant==5: dv = 3.0*n*(nxv@S) + 3.0*(n@S)*nxv - vxS
            a = a + (2.0/r**3)*dv
        return np.concatenate([v, a])
    return f

def integrate(order, Svec, x0, v0, t_end, n_out=4000, rtol=1e-12, atol=1e-14, variant=None):
    f = rhs(order,Svec) if variant is None else rhs_alt(order,Svec,variant)
    t = np.linspace(0, t_end, n_out)
    s = solve_ivp(f, (0,t_end), np.concatenate([x0,v0]), t_eval=t,
                  method="DOP853", rtol=rtol, atol=atol)
    return s.t, s.y[:3].T, s.y[3:].T

def elements_from_state(x, v, mu=1.0):
    """osculating a,e,inc,Omega,omega from Cartesian (Newtonian)."""
    r = np.linalg.norm(x,axis=-1); vv=np.einsum('...i,...i->...',v,v)
    h = np.cross(x,v); hn = np.linalg.norm(h,axis=-1)
    E = vv/2 - mu/r; a = -mu/(2*E)
    ev = np.cross(v,h)/mu - x/r[...,None]; e = np.linalg.norm(ev,axis=-1)
    inc = np.arccos(h[...,2]/hn)
    nvec = np.stack([-h[...,1], h[...,0], np.zeros_like(h[...,0])],axis=-1)
    nn = np.linalg.norm(nvec,axis=-1)
    Om = np.arctan2(nvec[...,1], nvec[...,0])
    cosw = np.einsum('...i,...i->...',nvec,ev)/(nn*e)
    om = np.arccos(np.clip(cosw,-1,1))
    om = np.where(ev[...,2] < 0, 2*np.pi-om, om)
    return a,e,inc,np.mod(Om,2*np.pi),np.mod(om,2*np.pi)
