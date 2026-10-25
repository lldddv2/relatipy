import numpy as np, warnings; warnings.filterwarnings("ignore")
from fisher import jac_star, setup, M0, R00, RAD
ch=np.load("chain_S2_S62.npy"); c=ch[ch.shape[0]//2:].reshape(-1,ch.shape[2]); med=np.median(c,axis=0)
print("posterior median globals:", np.round(med[:5],4), " th=%.1f deg"%np.degrees(np.arccos(med[3])))
def sig(g, names=("S2","S62-like"), stars=None, sa=10., sv=5.):
    F=np.zeros((5+6*len(names),)*2)
    for k,nm in enumerate(names):
        s,t=setup(nm)
        if stars is not None: s=list(stars[k])
        J=jac_star((g,s,t)); n=J.shape[0]//3
        W=np.concatenate([np.full(2*n,1/sa**2),np.full(n,1/sv**2)])
        idx=list(range(5))+list(range(5+6*k,11+6*k)); F[np.ix_(idx,idx)]+=J.T@(J*W[:,None])
    return np.sqrt(np.diag(np.linalg.inv(F)))
g_true=[M0,R00,0.9,60*RAD,30*RAD]
sd=sig(g_true); print(f"Fisher @truth : s(a*)={sd[2]:.4f} s(th)={np.degrees(sd[3]):.2f} s(ph)={np.degrees(sd[4]):.2f} deg")
g_med=[med[0],med[1],med[2],np.arccos(med[3]),med[4]]
stars=[med[5:11],med[11:17]]
sd=sig(g_med,stars=stars); print(f"Fisher @median: s(a*)={sd[2]:.4f} s(th)={np.degrees(sd[3]):.2f} s(ph)={np.degrees(sd[4]):.2f} deg")
print(f"MCMC          : s(a*)={np.std(c[:,2]):.4f} s(th)={np.degrees(np.std(np.arccos(c[:,3]))):.2f} s(ph)={np.degrees(np.std(c[:,4])):.2f} deg")
