import numpy as np, math, json
import scipy.sparse as sp, scipy.sparse.linalg as spla
from scipy.special import softmax
G=5.145228719489142; N=9; Q=12; theta=2.0; kap=18.0
j=np.arange(Q)
W=np.array([[0 if i==k else math.cos(.18575*min(abs(i-k),Q-abs(i-k))+.1625)/(1+min(abs(i-k),Q-abs(i-k))**1.8) for k in range(Q)] for i in range(Q)])
Lap=np.diag(W.sum(1))-W; lam=np.fft.fft(Lap[0]).real[:7]; cols=[]
for k in (3,4,5): cols += [np.sqrt(lam[k]/6)*np.cos(2*np.pi*k*j/Q),np.sqrt(lam[k]/6)*np.sin(2*np.pi*k*j/Q)]
cols += [np.sqrt(lam[6]/12)*(-1.)**j]; X=np.column_stack(cols); A=X@X.T
z=np.load(f'/mnt/data/defect_research_inputs/base_N{N}.npz',allow_pickle=True);st=z['states'];lab=z['labels'];pi=z['pi']; idx={tuple(map(int,s)):i for i,s in enumerate(st)}
j0=np.where(lab==0)[0]; loc={int(g):a for a,g in enumerate(j0)}; parent=np.arange(len(j0)); sz=np.ones(len(j0),int)
def find(a):
 while parent[a]!=a: parent[a]=parent[parent[a]]; a=parent[a]
 return a
def union(a,b):
 a=find(a);b=find(b)
 if a==b:return
 if sz[a]<sz[b]:a,b=b,a
 parent[b]=a;sz[a]+=sz[b]
for aa,gidx in enumerate(j0):
 n=st[gidx]
 for i in np.where(n>0)[0]:
  for jj in range(Q):
   if jj==i:continue
   n2=n.copy(); n2[i]-=1;n2[jj]+=1;h=idx[tuple(map(int,n2))]
   if lab[h]==0:union(aa,loc[int(h)])
seedg=np.where((st[:,0]==N)&(st[:,1:].sum(1)==0))[0][0]; sr=find(loc[int(seedg)])
ids=j0[np.array([find(i)==sr for i in range(len(j0))])]; pos={int(g):a for a,g in enumerate(ids)}
lw=np.log(pi[ids])+kap*st[ids,0]/N;lw-=lw.max();mu=np.exp(lw);mu/=mu.sum()
rows=[];cc=[];dat=[];diag=np.zeros(len(ids))
for aa,gidx in enumerate(ids):
 n=st[gidx]
 for i in np.where(n>0)[0]:
  ni=int(n[i]); m=n.astype(float).copy();m[i]-=1; field=(G/N)*(A@m);field[0]+=theta;q=softmax(field)
  for jj in range(Q):
   if jj==i:continue
   n2=n.copy();n2[i]-=1;n2[jj]+=1;h=idx[tuple(map(int,n2))];bb=pos.get(int(h))
   if bb is None:continue
   rate=ni*q[jj]; rows.append(aa);cc.append(bb);dat.append(rate);diag[aa]-=rate
rows.extend(range(len(ids)));cc.extend(range(len(ids)));dat.extend(diag.tolist())
Qm=sp.csr_matrix((dat,(rows,cc)),shape=(len(ids),len(ids)))
sq=np.sqrt(mu);S=sp.diags(sq)@Qm@sp.diags(1/sq);S=(S+S.T)*.5
vals,vecs=spla.eigsh(S,k=5,which='LA',tol=1e-11,maxiter=200000);ord=np.argsort(vals)[::-1];vals=vals[ord];vecs=vecs[:,ord]
# Eigenfunctions in L2(mu): f = v/sqrt(mu), normalize centered automatically for nonzero modes.
f=vecs[:,1:3]/sq[:,None]
pop=st[ids]/N
# natural Fourier count observables (unscaled), real and imaginary k=1..6 on label cycle
obs={}
for k in range(1,7):
 c=np.cos(2*np.pi*k*j/Q); s=np.sin(2*np.pi*k*j/Q)
 obs[f'cos{k}']=pop@c; obs[f'sin{k}']=pop@s
obs['D']=(N-st[ids,0]).astype(float); obs['n0']=st[ids,0].astype(float)
def wcorr(x,y):
 xm=np.sum(mu*x); ym=np.sum(mu*y); dx=x-xm;dy=y-ym
 den=np.sqrt(np.sum(mu*dx*dx)*np.sum(mu*dy*dy));return float(np.sum(mu*dx*dy)/den) if den else 0.
corr={name:[wcorr(x,f[:,a]) for a in range(2)] for name,x in obs.items()}
# rotation by 3 labels permutation if stays within component; transform eigenfunctions and estimate 2x2 action
idpos={tuple(map(int,st[g])):a for a,g in enumerate(ids)}
perm=[]; valid=True
for g in ids:
 n=st[g]; nr=np.roll(n,3); b=idpos.get(tuple(map(int,nr)))
 if b is None:valid=False;break
 perm.append(b)
action=None
if valid:
 perm=np.array(perm); # f_rot(x)=f(R x) or perm
 # L2 projection matrix <f_a, f_b after rotation>
 action=np.array([[np.sum(mu*f[:,a]*f[perm,b]) for b in range(2)] for a in range(2)])
# sign cut conductance for first slow mode
cuts=[]
for a in range(2):
 mask=f[:,a]>=0; mass=float(mu[mask].sum());
 # flow from mask to complement: sum mu_i q_ij
 coo=Qm.tocoo();sel=(mask[coo.row]) & (~mask[coo.col]) & (coo.row!=coo.col);flow=float(np.sum(mu[coo.row[sel]]*coo.data[sel]));phi=flow/min(mass,1-mass)
 cuts.append({'mode':a+1,'positive_mass':mass,'flow':flow,'conductance_sign_cut':phi})
seedpos=pos[seedg]
seed_slow_eigenfunctions=f[seedpos,:].tolist()
out={'eigenvalues':vals.tolist(),'seed_slow_eigenfunctions':seed_slow_eigenfunctions,'correlations':corr,'rotation_by_3_valid':valid,'slow_subspace_rotation_action':action.tolist() if action is not None else None,'sign_cuts':cuts}
open('/mnt/data/defect_research_inputs/PREP_N9_SLOW_MODE_330.json','w').write(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
