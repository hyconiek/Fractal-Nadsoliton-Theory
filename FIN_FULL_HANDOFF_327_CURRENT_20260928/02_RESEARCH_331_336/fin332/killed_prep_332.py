import numpy as np, math, json, scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply
from scipy.special import softmax
G=5.145228719489142;Q=12;theta=2.;times=[1,2,4,8]
j=np.arange(Q);W=np.array([[0 if i==k else math.cos(.18575*min(abs(i-k),Q-abs(i-k))+.1625)/(1+min(abs(i-k),Q-abs(i-k))**1.8) for k in range(Q)] for i in range(Q)]);Lap=np.diag(W.sum(1))-W;lam=np.fft.fft(Lap[0]).real[:7];cols=[]
for k in (3,4,5):cols += [np.sqrt(lam[k]/6)*np.cos(2*np.pi*k*j/Q),np.sqrt(lam[k]/6)*np.sin(2*np.pi*k*j/Q)]
cols += [np.sqrt(lam[6]/12)*(-1.)**j];X=np.column_stack(cols);A=X@X.T
base='/mnt/data/work327330/FIN_RESEARCH_CONTINUATION_327_330_20260928/fin330'
out=[]
for N in range(3,11):
 z=np.load(f'{base}/base_N{N}.npz');st=z['states'];lab=z['labels'];pi=z['pi'];idx={tuple(map(int,s)):ii for ii,s in enumerate(st)};j0=np.where(lab==0)[0];pos={int(g):a for a,g in enumerate(j0)}
 lw=np.log(pi[j0])+theta*st[j0,0];lw-=lw.max();target=np.exp(lw);target/=target.sum()
 rows=[];cols=[];dat=[];diag=np.zeros(len(j0));out_rate=np.zeros(len(j0))
 for aa,gidx in enumerate(j0):
  n=st[gidx];tot=0.;outg=0.
  for i in np.where(n>0)[0]:
   ni=int(n[i]);m=n.astype(float).copy();m[i]-=1;field=(G/N)*(A@m);field[0]+=theta;q=softmax(field)
   for jj in range(Q):
    if jj==i:continue
    rate=ni*q[jj];tot+=rate;n2=n.copy();n2[i]-=1;n2[jj]+=1;h=idx[tuple(map(int,n2))];bb=pos.get(int(h))
    if bb is not None: rows.append(aa);cols.append(bb);dat.append(rate)
    else: outg+=rate
  diag[aa]=-tot;out_rate[aa]=outg
 rows.extend(range(len(j0)));cols.extend(range(len(j0)));dat.extend(diag.tolist());K=sp.csr_matrix((dat,(rows,cols)),shape=(len(j0),len(j0)))
 seedg=np.where((st[:,0]==N)&(st[:,1:].sum(1)==0))[0][0];v=np.zeros(len(j0));v[pos[int(seedg)]]=1
 rec={'N':N,'theta':theta,'j0_states':len(j0),'by_t':{}}
 # compute individually to avoid interpolation issues
 for t in times:
  p=expm_multiply(K.T*t,v,traceA=float(K.diagonal().sum()));surv=float(p.sum());cond=p/surv;tv=.5*np.abs(cond-target).sum()
  # derivative of survival = - p @ out_rate
  hazard=float((p@out_rate)/surv)
  rec['by_t'][str(t)]={'survival_no_exit':surv,'exit_by_t':1-surv,'TV_survivor_to_hardwall_target':float(tv),'conditional_hazard_at_t':hazard}
 out.append(rec);print(json.dumps(rec),flush=True)
json.dump({'rows':out},open('/mnt/data/work331/KILLED_PREP_332.json','w'),indent=2)
