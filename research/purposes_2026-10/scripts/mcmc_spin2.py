"""Same posterior as mcmc_spin.py, but 96 walkers initialised from the Fisher covariance."""
import numpy as np, emcee, time, sys, warnings
from multiprocessing import Pool
warnings.filterwarnings("ignore")
NSTEP=int(sys.argv[1]) if len(sys.argv)>1 else 4000
sys.argv=[sys.argv[0],"S2,S62-like","0","S2_S62"]
import mcmc_spin as M
from fisher import jac_star, RAD
if __name__=="__main__":
    g=[M.g_true[0],M.g_true[1],M.g_true[2],np.arccos(M.g_true[3]),M.g_true[4]]
    ns=len(M.SETUPS); F=np.zeros((5+6*ns,)*2)
    for k,(s,t) in enumerate(M.SETUPS):
        J=jac_star((g,s,t)); n=J.shape[0]//3
        J=J.copy(); J[:,3]=J[:,3]/(-np.sin(g[3]))        # d/dth -> d/dcos(th)
        W=np.concatenate([np.full(2*n,1/M.SA**2),np.full(n,1/M.SV**2)])
        idx=list(range(5))+list(range(5+6*k,11+6*k)); F[np.ix_(idx,idx)]+=J.T@(J*W[:,None])
    C=np.linalg.inv(F); C=0.5*(C+C.T)
    p_true=np.array(M.g_true+sum([list(s) for s,_ in M.SETUPS],[]))
    rng=np.random.default_rng(11); nw=96; ndim=len(p_true)
    p0=rng.multivariate_normal(p_true,C,size=nw); p0[:,2]=np.clip(p0[:,2],0.01,0.999); p0[:,3]=np.clip(p0[:,3],-0.999,0.999)
    print(f"Fisher sigmas (a*, cos_th, ph[deg]) = {np.sqrt(C[2,2]):.4f} {np.sqrt(C[3,3]):.4f} {np.degrees(np.sqrt(C[4,4])):.2f}; nw={nw}",flush=True)
    t0=time.perf_counter()
    with Pool(4) as pool:
        sam=emcee.EnsembleSampler(nw,ndim,M.logp,pool=pool,moves=[(emcee.moves.StretchMove(),0.5),(emcee.moves.DEMove(),0.5)])
        for i,_ in enumerate(sam.sample(p0,iterations=NSTEP)):
            if (i+1)%250==0:
                ch=sam.get_chain(); c=ch[-250:].reshape(-1,ndim)
                print(f"step {i+1} acc={sam.acceptance_fraction.mean():.2f} {time.perf_counter()-t0:.0f}s  last250: "
                      f"a*={np.median(c[:,2]):.3f}+-{c[:,2].std():.4f} cos_th={np.median(c[:,3]):.3f}+-{c[:,3].std():.4f} "
                      f"ph={np.degrees(np.median(c[:,4])):.1f}+-{np.degrees(c[:,4].std()):.2f}",flush=True)
                np.save("chain2_S2_S62.npy",ch)
