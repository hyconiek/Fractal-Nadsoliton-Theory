import numpy as np, math, json, time
import scipy.sparse as sp, scipy.sparse.linalg as spla
from scipy.special import softmax
G=5.145228719489142;Q=12;theta=2.0
j=np.arange(Q);W=np.array([[0 if i==k else math.cos(.18575*min(abs(i-k),Q-abs(i-k))+.1625)/(1+min(abs(i-k),Q-abs(i-k))**1.8) for k in range(Q)] for i in range(Q)]);Lap=np.diag(W.sum(1))-W;lam=np.fft.fft(Lap[0]).real[:7];cols=[]
for k in (3,4,5):cols += [np.sqrt(lam[k]/6)*np.cos(2*np.pi*k*j/Q),np.sqrt(lam[k]/6)*np.sin(2*np.pi*k*j/Q)]
cols += [np.sqrt(lam[6]/12)*(-1.)**j];X=np.column_stack(cols);A=X@X.T

def maincomp(st,lab,N):
 idx={tuple(map(int,s)):i for i,s in enumerate(st)};j0=np.where(lab==0)[0];loc={int(g):a for a,g in enumerate(j0)};parent=np.arange(len(j0));sz=np.ones(len(j0),int)
 def find(a):
  while parent[a]!=a:parent[a]=parent[parent[a]];a=parent[a]
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
    n2=n.copy();n2[i]-=1;n2[jj]+=1;h=idx[tuple(map(int,n2))]
    if lab[h]==0:union(aa,loc[int(h)])
 seedg=np.where((st[:,0]==N)&(st[:,1:].sum(1)==0))[0][0];sr=find(loc[int(seedg)]);ids=j0[np.array([find(i)==sr for i in range(len(j0))])];return ids,idx

def analyze(N):
 z=np.load(f'/mnt/data/defect_research_inputs/base_N{N}.npz',allow_pickle=True);st=z['states'];lab=z['labels'];pi=z['pi']; ids,idx=maincomp(st,lab,N);pos={int(g):a for a,g in enumerate(ids)};kap=theta*N
 # stationary weights on component
 lw=np.log(pi[ids])+kap*st[ids,0]/N;lw-=lw.max();mu=np.exp(lw);mu/=mu.sum()
 rows=[];cols=[];dat=[];diag=np.zeros(len(ids));maxdb=0.
 for aa,gidx in enumerate(ids):
  n=st[gidx]
  for i in np.where(n>0)[0]:
   ni=int(n[i]);m=n.astype(float).copy();m[i]-=1;field=(G/N)*(A@m);field[0]+=theta;q=softmax(field)
   for jj in range(Q):
    if jj==i:continue
    n2=n.copy();n2[i]-=1;n2[jj]+=1;h=idx[tuple(map(int,n2))]
    bb=pos.get(int(h))
    if bb is None:continue
    rate=ni*q[jj];rows.append(aa);cols.append(bb);dat.append(rate);diag[aa]-=rate
 rows.extend(range(len(ids)));cols.extend(range(len(ids)));dat.extend(diag.tolist());Qm=sp.csr_matrix((dat,(rows,cols)),shape=(len(ids),len(ids)))
 # detailed-balance residual on sparse entries via symmetric transform
 sq=np.sqrt(mu); inv=1/sq; S=sp.diags(sq)@Qm@sp.diags(inv);asym=S-S.T;asymres=float(np.max(np.abs(asym.data))) if asym.nnz else 0.
 vals=spla.eigsh((S+S.T)*.5,k=min(4,len(ids)-1),which='LA',return_eigenvectors=False,tol=1e-10,maxiter=100000);vals=np.sort(vals)[::-1];gap=float(-vals[1]);
 pmin=float(mu.min()); # crude TV mixing time upper from L2: tv<=.5/sqrt(pmin)*exp(-gap t)
 t1=float(math.log(.5/math.sqrt(pmin)/.01)/gap) if gap>0 else math.inf
 return {'N':N,'kappa':kap,'theta':theta,'component_states':len(ids),'nnz':int(Qm.nnz),'detailed_balance_sym_residual':asymres,'eigs_LA':vals.tolist(),'spectral_gap':gap,'relaxation_time_1_over_gap':1/gap,'crude_time_to_TV_1pct_from_worst_state':t1,'pi_min_component':pmin}
out=[]
for N in range(3,11):
 t=time.time();r=analyze(N);r['seconds']=time.time()-t;out.append(r);print(json.dumps(r),flush=True)
open('/mnt/data/defect_research_inputs/PREP_MIXING_330.json','w').write(json.dumps({'rows':out},indent=2))
