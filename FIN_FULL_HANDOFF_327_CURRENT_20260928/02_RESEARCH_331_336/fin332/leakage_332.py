import numpy as np, math, json
from scipy.special import softmax
G=5.145228719489142;Q=12;theta=2.0
j=np.arange(Q);W=np.array([[0 if i==k else math.cos(.18575*min(abs(i-k),Q-abs(i-k))+.1625)/(1+min(abs(i-k),Q-abs(i-k))**1.8) for k in range(Q)] for i in range(Q)]);Lap=np.diag(W.sum(1))-W;lam=np.fft.fft(Lap[0]).real[:7];cols=[]
for k in (3,4,5):cols += [np.sqrt(lam[k]/6)*np.cos(2*np.pi*k*j/Q),np.sqrt(lam[k]/6)*np.sin(2*np.pi*k*j/Q)]
cols += [np.sqrt(lam[6]/12)*(-1.)**j];X=np.column_stack(cols);A=X@X.T
base='/mnt/data/work327330/FIN_RESEARCH_CONTINUATION_327_330_20260928/fin330'
out=[]
for N in range(3,11):
 z=np.load(f'{base}/base_N{N}.npz');st=z['states'];lab=z['labels'];pi=z['pi'];idx={tuple(map(int,s)):ii for ii,s in enumerate(st)};j0=np.where(lab==0)[0]
 lw=np.log(pi[j0])+theta*st[j0,0];lw-=lw.max();mu=np.exp(lw);mu/=mu.sum()
 rates=[]
 for gidx in j0:
  n=st[gidx]; rr=0.
  for i in np.where(n>0)[0]:
   ni=int(n[i]);m=n.astype(float).copy();m[i]-=1;field=(G/N)*(A@m);field[0]+=theta;q=softmax(field)
   for jj in range(Q):
    if jj==i:continue
    n2=n.copy();n2[i]-=1;n2[jj]+=1;h=idx[tuple(map(int,n2))]
    if lab[h]!=0: rr+=ni*q[jj]
  rates.append(rr)
 rates=np.array(rates); seed=np.where((st[:,0]==N)&(st[:,1:].sum(1)==0))[0][0];seedpos=np.where(j0==seed)[0][0]
 out.append({'N':N,'theta':theta,'seed_exit_rate':float(rates[seedpos]),'mean_exit_flux_under_conditional_target':float(mu@rates),'max_exit_rate_in_J0':float(rates.max()),'target_boundary_mass_rate_positive':float(mu[rates>0].sum()),'naive_mean_escape_time_from_stationary_flux':float(1/(mu@rates)) if mu@rates>0 else None})
 print(out[-1],flush=True)
json.dump({'rows':out},open('/mnt/data/work331/UNRESTRICTED_LEAKAGE_332.json','w'),indent=2)
