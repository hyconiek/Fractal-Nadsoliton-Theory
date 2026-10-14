"""Reproducible finite q12 FIN parent with explicit reversible resolved/hidden projections."""
import json, hashlib, math
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix, diags
from scipy.special import logsumexp, softmax
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply
Q=12; G=5.145228719489142; N=1728; Y=144
HERE=Path(__file__).resolve().parent

def write(rel,obj):
    p=HERE/rel;p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False,default=lambda o: o.item() if isinstance(o,np.generic) else str(o))+'\n')
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def A7():
    d=np.arange(Q)
    W=np.zeros(Q)
    for k in range(1,Q):
        r=min(k,Q-k)
        W[k]=math.cos(.18575*r+.1625)/(1+r**1.8)
    eig=np.array([sum(W[r]*(1-math.cos(2*math.pi*k*r/Q)) for r in range(1,Q)) for k in range(Q)])
    # formula used by versioned safe-seed: Fourier real rank-seven projection, modes 3,4,5,6
    a=np.array([sum(eig[k]/6*math.cos(2*math.pi*k*r/Q) for k in (3,4,5))+(eig[6]/12)*(-1)**r for r in range(Q)])
    return np.array([np.roll(a,i) for i in range(Q)])
def models():
    A=A7();trace=np.trace(A);P=np.eye(Q)-np.ones((Q,Q))/Q
    eig=np.zeros(Q);eig[[3,4,5,6,7,8,9]]=trace/7
    a=np.fft.ifft(eig).real
    return {'FIN_A7':A,'FULL_POTTS_TRACE':P*trace/11,'FLAT_P7_TRACE':np.array([np.roll(a,i) for i in range(Q)])}

def build(A,g=G):
    states=np.indices((Q,Q,Q)).reshape(3,-1).T
    n=np.zeros((N,Q),int)
    for k in range(3): n[np.arange(N),states[:,k]]+=1
    E=(g/6)*np.einsum('ni,ij,nj->n',n,A,n)
    pi=np.exp(E-logsumexp(E));yi=states[:,0]*Q+states[:,1]
    piY=np.bincount(yi,weights=pi,minlength=Y)
    # R[other 1, other 2, target label] is full categorical refresh including self-events
    R=np.empty((Q,Q,Q))
    for i in range(Q):
      for j in range(Q): R[i,j,:]=softmax((g/3)*(A[:,i]+A[:,j]))
    rows=[];cols=[];dat=[]
    for ix,(a,b,c) in enumerate(states):
      for j in range(Q):
        rows.extend((ix,ix,ix))
        cols.extend((((j*Q+b)*Q+c),((a*Q+j)*Q+c),((a*Q+b)*Q+j)))
        dat.extend((R[b,c,j]/3,R[a,c,j]/3,R[a,b,j]/3))
    P=csr_matrix((dat,(rows,cols)),shape=(N,N));L=3*(P-csr_matrix(np.eye(N)))
    spi=np.sqrt(pi);U=csr_matrix((np.sqrt(pi/piY[yi]),(np.arange(N),yi)),shape=(N,Y))
    S=diags(spi)@L@diags(1/spi)
    SU=S@U;AA=(U.T@SU).toarray(); SQ=SU-U@csr_matrix(AA)
    # parent reversible transform and exact 2-copy equilibrium marginal
    return {'A7':A,'g':g,'states':states,'pi':pi,'piY':piY,'P':P,'L':L,'S':S.tocsr(),'U':U,'A':AA,'SU':SU,'C':SQ,'R':R,'yi':yi}
def make_prep(m,kind):
    p=np.zeros(N)
    if kind=='conditional_00': p[:12]=m['pi'][:12]/m['piY'][0]
    elif kind.startswith('delta_'):
        k=int(kind.split('_')[-1]);p[k]=1.
    else: raise ValueError(kind)
    return p

def project_hidden(m,v): return v-m['U']@(m['U'].T@v)
def z_from_p(m,p): return p/np.sqrt(m['pi'])
def observe_z(m,z): return np.sqrt(m['piY'])*(m['U'].T@z)
def markov(m,t):
    init=np.zeros(Y);init[0]=1
    return np.sqrt(m['piY'])*(expm(m['A']*t)@(init/np.sqrt(m['piY'])))
def poisson_tail(z,m):
    from scipy.special import gammainc
    # Poisson Pr[X>m] = regularized lower incomplete Gamma(m+1,z)
    return float(gammainc(m+1,z))
def uniform_p(m,p0,t,cutoff=65):
    # P refresh includes 'self' outcomes (not a CTMC jump); e^(3t(P-I)) uniformization
    z=3*t
    p=p0.copy();w=math.exp(-z);ans=w*p.copy()
    for k in range(1,cutoff+1):
        p=m['P'].T@p
        w*=z/k
        ans+=w*p
    return ans

def uniform_tensor(m,p0,t,cutoff=65):
    # INDEPENDENT tensor-contraction method for reference rather than csr matrix product.
    z=3*t;R=m['R'];a,b,c=np.indices((Q,Q,Q));
    r1=R[b,c,a];r2=R[a,c,b];r3=R[a,b,c]
    p=p0.reshape(Q,Q,Q).copy();w=math.exp(-z);out=w*p.copy()
    for k in range(1,cutoff+1):
        p=(r1*p.sum(axis=0)[None,:,:]+r2*p.sum(axis=1)[:,None,:]+r3*p.sum(axis=2)[:,:,None])/3
        w*=z/k;out+=w*p
    return out.ravel()
def prep_names(): return ['conditional_00','delta_0','delta_1']
def observed(m,p): return p.reshape(Q,Q,Q).sum(axis=2).ravel()
def tv(a,b): return float(.5*np.abs(a-b).sum())

def model_diagnostics(m):
    L=m['L'];P=m['P'];pi=m['pi'];S=m['S'];U=m['U'];A=m['A'];C=m['C']
    d=(S-S.T).tocoo();sym=np.max(np.abs(d.data)) if len(d.data) else 0.
    ident=U.T@(S@(S@U));semigroup=np.max(np.abs(ident.toarray()-A@A-(C.T@C).toarray()))
    return {'P_rowsum_inf':float(np.max(np.abs(np.asarray(P.sum(axis=1)).ravel()-1))), 'L_rowsum_inf':float(np.max(np.abs(np.asarray(L.sum(axis=1)).ravel()))),'stationarity_inf':float(np.max(np.abs(L.T@pi))), 'S_symmetry_inf':float(sym),'U_orthogonality_inf':float(np.max(np.abs((U.T@U).toarray()-np.eye(Y)))),'semigroup_defect_identity_inf':float(semigroup),'min_P':float(min(P.data)), 'piY_mass':float(m['piY'].sum()),'A_top_eigenvalue':float(np.linalg.eigvalsh(A)[-1]), 'S_norm_upper_exact':3.}

def orthogonal_basis(m,order, all_hidden=False):
    """Hidden Krylov vectors Q S^k z_p, training preparations fixed; optionally
    all 11 independent hidden-coordinate directions for Y=00, irrespective of
    later holdout hidden label. Cap 32 total directions. Modified Gram-Schmidt x2.
    """
    initial=[z_from_p(m,make_prep(m,n)) for n in prep_names()]
    vectors=[];labels=[]
    if all_hidden:
        for h in range(Q):
            v=z_from_p(m,make_prep(m,'delta_'+str(h)))
            vectors.append(project_hidden(m,v));labels.append('preparation_complete_hidden_'+str(h))
    for j in range(order+1):
        for i,v in enumerate(initial):
            vectors.append(project_hidden(m,v));labels.append(f'Q_S^{j}_'+prep_names()[i])
        initial=[m['S']@v for v in initial]
    basis=[];accepted=[]
    for v,name in zip(vectors,labels):
        v=np.asarray(v).ravel().copy()
        for _ in range(2):
            for b in basis: v-=b*np.dot(b,v)
        nrm=np.linalg.norm(v)
        if nrm>1e-9:
            basis.append(v/nrm);accepted.append(name)
        if len(basis)==32: break
    B=np.column_stack(basis) if basis else np.empty((N,0))
    return B,accepted

def reduced(m,B):
    U=m['U'];S=m['S'];A=m['A'];r=B.shape[1]
    if r:
        Z=np.asarray(U.T@(S@B));H=B.T@(S@B)
        T=np.block([[A,Z],[Z.T,H]])
    else:T=A.copy()
    def w_init(kind):
        z=z_from_p(m,make_prep(m,kind));return np.r_[np.asarray(U.T@z),B.T@z]
    def w_to_y(w):return np.sqrt(m['piY'])*w[:Y]
    def value(kind,t):return w_to_y(expm(T*t)@w_init(kind))
    return T,w_init,w_to_y,value

def memory_operator_stats(m,B):
    U=m['U']; S=m['S'];C=m['C'];A=m['A'];r=B.shape[1]
    if r==0:return {}
    Cp=np.asarray(B.T@C)
    H=-B.T@(S@B)
    eig=np.linalg.eigvalsh(H)
    v= C[:,0].toarray().ravel() if hasattr(C[:,0],'toarray') else np.asarray(C[:,0]).ravel()
    K0=float(np.dot(v,v));Kprime=float(v@(S@project_hidden(m,v)))
    H_e=np.linalg.eigvalsh((H+H.T)/2)
    # Khat memory is a partial Galerkin approximation; true H has zero eigenvalue? coupled relevant range only
    k0red=float((Cp[:,0]@Cp[:,0]));
    return {'k0_00':K0,'kprime_00':Kprime,'reduced_k0_00':k0red,'reduced_hidden_H_eig_min':float(eig.min()), 'reduced_hidden_H_eig_max':float(eig.max()), 'full_S_2norm_upper_analytic':3.,'Khat0_psd':bool(k0red>=-1e-10)}

def universal_hidden_basis(m,cap=32,degree=4):
    """Preparation-complete block-Krylov response compression without selecting a hidden label.
    First keep every independent hidden initial vector at Y=00. For each next Krylov
    degree, keep full residual span if fits, else choose leading SVD response directions.
    SVD uses source-operator vectors only, never reference trajectory outputs.
    """
    originals=[z_from_p(m,make_prep(m,'delta_'+str(k))) for k in range(Q)]
    Z=np.column_stack(originals)
    basis=[];labels=[]
    for j in range(degree+1):
        Zp=Z-np.asarray(m['U']@(m['U'].T@Z))
        if basis:
            B=np.column_stack(basis)
            Zp=Zp-B@(B.T@Zp)
            Zp=Zp-B@(B.T@Zp)
        W,s,vt=np.linalg.svd(Zp,full_matrices=False)
        new_count=min(int(sum(s>1e-8)),cap-len(basis))
        for k in range(new_count):
            vv=np.asarray(W[:,k]).ravel()
            # double-pass to preserve orthonormality under repeated transforms
            for _ in range(2):
              for b in basis:vv-=b*np.dot(b,vv)
            vv=vv/np.linalg.norm(vv)
            basis.append(vv);labels.append(f'ALL_HIDDEN_krylov_degree_{j}_SVD_rank_{k+1}')
        if len(basis)>=cap:break
        Z=m['S']@Z
    return np.column_stack(basis),labels

def preparation_conditioned_basis(m,prep,order):
    """Frozen rule: derive each prediction's response space from supplied initial
    distribution, without fitting any trajectory values, plus all resolved U.
    Classifies as preparation-conditioned ROM, not global reduced Markov process.
    """
    z=z_from_p(m,make_prep(m,prep))
    b=[];names=[]
    for k in range(order+1):
        v=project_hidden(m,z).copy()
        for _ in range(2):
            for b0 in b:v-=np.dot(b0,v)*b0
        n=np.linalg.norm(v)
        if n>1e-9: b.append(v/n);names.append('Q_S^'+str(k)+'_initial_density')
        z=m['S']@z
        if len(b)==32:break
    return np.column_stack(b) if b else np.zeros((N,0)),names

def channel_conditioned_basis(m,prep,depth=5,cap=32):
    """Frozen channel-aware Krylov/Galerkin construction. Source and initial
    preparation determine every direction. Readout fixed at first-stage
    event Y=(0,0) versus complement; no fitted rates or trajectories.
    Sorted total polynomial degree for first/second propagation and filtering.
    """
    z=z_from_p(m,make_prep(m,prep));B=[];labels=[]
    Dz=[];v=z.copy()
    for j in range(depth+1):
        Dz.append(v.copy());v=m['S']@v
    candidates=[]
    for k in range(depth+1):
        candidates.append((k,0,k,False,Dz[k]))
    for k in range(depth+1):
        v=np.zeros(N);v[:12]=Dz[k][:12]
        for ell in range(depth+1-k):
            candidates.append((k+ell,ell,k,True,v.copy()));v=m['S']@v
    candidates.sort(key=lambda row:(row[0],row[3],row[1],row[2]))
    for degree,ell,k,filtered,v in candidates:
        v=project_hidden(m,v).copy()
        for _ in range(2):
            for b in B:v-=np.dot(v,b)*b
        nrm=np.linalg.norm(v)
        if nrm>1e-8:
            B.append(v/nrm);labels.append(f'Q_S^{ell}_'+('D_E_' if filtered else '')+f'S^{k}_z0')
        if len(B)>=cap:break
    return np.column_stack(B) if B else np.empty((N,0)),labels
