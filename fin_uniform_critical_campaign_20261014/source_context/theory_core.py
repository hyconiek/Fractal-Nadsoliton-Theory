from __future__ import annotations
import math, json
import numpy as np
from scipy.optimize import root

Q=12
OMEGA=0.18575
PHI=0.1625
ETA=1.8
GF_CENTER=3.5156447168395917
GF_INTERVAL=(3.5156447068395917,3.5156447268395917)
A_INTERVAL=(-0.2125716400840138741,-0.2125716260204590083)
B_INTERVAL=(0.1189825274504987752,0.1189826671032236644)


def build_rank7():
    j=np.arange(Q)
    W=np.zeros((Q,Q),float)
    for i in range(Q):
        for k in range(Q):
            if i==k: continue
            d=min(abs(i-k),Q-abs(i-k))
            W[i,k]=math.cos(OMEGA*d+PHI)/(1+d**ETA)
    L=np.diag(W.sum(1))-W
    lap=np.fft.fft(L[0]).real[:7]
    cols=[]
    for k in (3,4,5):
        cols += [np.sqrt(lap[k]/6)*np.cos(2*np.pi*k*j/Q),
                 np.sqrt(lap[k]/6)*np.sin(2*np.pi*k*j/Q)]
    cols += [np.sqrt(lap[6]/12)*(-1.)**j]
    X=np.column_stack(cols)
    A=X@X.T
    return W,L,lap,X,A

W,LAPL,lap,X,A=build_rank7()
EVEN_IDX=np.array([0,2,4,6])
Xe=X[:,EVEN_IDX]


def softmax_vec(z):
    z=np.asarray(z,float)
    zz=z-z.max()
    e=np.exp(zz)
    return e/e.sum()


def phi(theta,g):
    theta=np.asarray(theta,float)
    z=X@theta
    m=z.max()
    return float(theta@theta/(2*g) - (m+math.log(np.exp(z-m).mean())))


def grad_phi(theta,g):
    theta=np.asarray(theta,float)
    return theta/g - X.T@softmax_vec(X@theta)


def hess_phi(theta,g):
    theta=np.asarray(theta,float)
    p=softmax_vec(X@theta)
    S=np.diag(p)-np.outer(p,p)
    return np.eye(7)/g-X.T@S@X


def even_grad(theta4,g):
    p=softmax_vec(Xe@theta4)
    return theta4/g-Xe.T@p


def even_hess(theta4,g):
    p=softmax_vec(Xe@theta4)
    S=np.diag(p)-np.outer(p,p)
    return np.eye(4)/g-Xe.T@S@Xe


def solve_fold():
    y0=np.array([1.3643114281782407,1.4331027014285425,1.4080457012550374,1.0056820773201627,GF_CENTER])
    def eq(y):
        return np.r_[even_grad(y[:4],y[4]),np.linalg.det(even_hess(y[:4],y[4]))]
    sol=root(eq,y0,tol=1e-12)
    if not sol.success and np.linalg.norm(eq(sol.x))>1e-10:
        raise RuntimeError(sol.message)
    th4=sol.x[:4]; g=float(sol.x[4])
    H=even_hess(th4,g)
    evals,V=np.linalg.eigh(H)
    v=V[:,0]
    if v[0]<0: v=-v
    p=softmax_vec(Xe@th4)
    a=float(-v@th4/g**2)
    zv=Xe@v; centered=zv-p@zv
    b=float(-np.sum(p*centered**3))
    th7=np.zeros(7); th7[EVEN_IDX]=th4
    v7=np.zeros(7); v7[EVEN_IDX]=v
    return dict(theta4=th4,g=g,evals4=evals,v4=v,p=p,a=a,b=b,theta7=th7,v7=v7)


def finite_loo_rate(p,i,j,g):
    p=np.asarray(p,float)
    Ninv=None
    raise RuntimeError('Use finite_loo_rate_from_counts for a concrete N')


def finite_loo_rate_from_counts(n,i,j,g):
    n=np.asarray(n,int); N=int(n.sum())
    m=n.copy(); m[i]-=1
    s=softmax_vec((g/N)*(A@m))
    return float(n[i]*s[j])


def limiting_edge_rate(p,i,j,g):
    p=np.asarray(p,float)
    s=softmax_vec(g*(A@p))
    return float(p[i]*s[j])


def hamiltonian(p,xi,g):
    p=np.asarray(p,float); xi=np.asarray(xi,float)
    s=softmax_vec(g*(A@p))
    return float(np.sum(p[:,None]*s[None,:]*(np.exp(xi[None,:]-xi[:,None])-1)))


def free_energy_grad(p,g):
    p=np.asarray(p,float)
    return np.log(12*p)+1-g*(A@p)


def drift(p,g):
    return softmax_vec(g*(A@p))-np.asarray(p,float)


def diffusion_D(p,g):
    p=np.asarray(p,float); s=softmax_vec(g*(A@p))
    return np.diag(p)+np.diag(s)-np.outer(p,s)-np.outer(s,p)


def balanced_counts(p,N):
    x=np.asarray(p,float)*N
    n=np.floor(x).astype(int)
    rem=N-int(n.sum())
    if rem:
        order=np.argsort(-(x-n)); n[order[:rem]]+=1
    return n


def fold_generator_metrics(N,lamb):
    f=solve_fold(); gf=f['g']; pf=f['p']; thf=f['theta7']; v=f['v7']; a=f['a']; b=f['b']
    gn=gf+lamb*N**(-2/3)
    n=balanced_counts(pf,N); pp=n/N
    theta=gn*(X.T@pp)
    Z=N**(1/3)*float(v@(theta-thf))
    drift_z=var_z=third_abs=0.0; total=0.0
    for i in range(Q):
        if n[i]==0: continue
        m=n.copy(); m[i]-=1
        probs=softmax_vec((gn/N)*(A@m))
        for jj in range(Q):
            if jj==i: continue
            rate=float(n[i]*probs[jj])
            dz=N**(1/3)*(gn/N)*float(v@(X[jj]-X[i]))
            total+=rate; drift_z+=rate*dz; var_z+=rate*dz*dz; third_abs+=rate*abs(dz)**3
    fac=N**(1/3)
    pred=-gf*(a*lamb+(b/2)*Z*Z)
    return dict(N=N,lambda_=lamb,Z_actual=Z,scaled_drift=fac*drift_z,predicted_drift=pred,
                drift_error=fac*drift_z-pred,scaled_variance=fac*var_z,target_variance=2*gf,
                variance_error=fac*var_z-2*gf,scaled_abs_third_jump_rate=fac*third_abs,total_rate=total)


def moment_lift_matrix():
    j=np.arange(Q)
    c={k:np.cos(2*np.pi*k*j/Q) for k in range(1,7)}
    s={k:np.sin(2*np.pi*k*j/Q) for k in range(1,6)}
    # normalization + seven retained feature means + four raw products
    M=np.vstack([np.ones(Q), *[X[:,r] for r in range(7)], c[3]*c[4], s[3]*c[4], c[3]*c[5], s[3]*c[5]])
    return M


def hidden_witness_counts(sign=1,N=120):
    if N!=120: raise ValueError('witness frozen at N=120')
    j=np.arange(Q)
    h=np.cos(2*np.pi*2*j/Q) # k=2 hidden cosine; 2h is integer-valued
    n=np.rint(10+sign*2*h).astype(int)
    assert n.sum()==120 and np.all(n>=0)
    return n


def generator_of_mu5_squared_g0(n):
    n=np.asarray(n,int); N=int(n.sum())
    y=X[:,4] # c5 feature column
    mu=float((n/N)@y)
    val=0.0
    for i in range(Q):
        if n[i]==0: continue
        for j in range(Q):
            if j==i: continue
            rate=n[i]/Q
            mun=mu+(y[j]-y[i])/N
            val += rate*(mun*mun-mu*mu)
    return val


def soft_jump_cumulants_at_fold():
    f=solve_fold(); g=f['g']; p=f['p']; v=f['v4']; z=Xe@v
    vals=[]
    for i in range(Q):
        for j in range(Q):
            vals.append((p[i]*p[j],g*(z[j]-z[i])))
    ms={k:sum(w*d**k for w,d in vals) for k in range(1,5)}
    return ms


def normalized_sde_constants():
    f=solve_fold(); a=f['a']; b=f['b']; g=f['g']
    c=(4/b)**(1/3)
    d=c*c/(2*g)
    kappa=-a*c/2
    return dict(Z_scale=c,s_time_per_tau=d,lambda_to_canonical=kappa,
                canonical='dY=(Lambda-Y^2)d_tau+dW, Lambda=kappa*lambda')

if __name__=='__main__':
    f=solve_fold()
    out={
        'laplacian_eigenvalues_0_6':lap.tolist(),
        'A_diag':float(A[0,0]),'A_rank':int(np.linalg.matrix_rank(A,1e-10)),
        'fold':{k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in f.items() if k not in ('p','theta7','v7')},
        'fold_p':f['p'].tolist(),'fold_theta7':f['theta7'].tolist(),'fold_v7':f['v7'].tolist(),
        'moment_lift_rank':int(np.linalg.matrix_rank(moment_lift_matrix(),1e-11)),
        'hidden_witness_difference':generator_of_mu5_squared_g0(hidden_witness_counts(+1))-generator_of_mu5_squared_g0(hidden_witness_counts(-1)),
        'hidden_witness_formula':float(lap[5]/7200),
        'soft_cumulants':soft_jump_cumulants_at_fold(),
        'normalized_sde':normalized_sde_constants(),
    }
    print(json.dumps(out,indent=2))
