import numpy as np, math, json
from pathlib import Path
G=5.145228719489142
Q=12
# exact repository finite-rank construction
j=np.arange(Q)
W=np.array([[0.0 if i==k else math.cos(0.18575*min(abs(i-k),Q-abs(i-k))+0.1625)/(1+min(abs(i-k),Q-abs(i-k))**1.8) for k in range(Q)] for i in range(Q)])
A=np.diag(W.sum(axis=1))-W
L=np.fft.fft(A[0]).real[:7]
cols=[]
for k in (3,4,5):
    cols += [np.sqrt(L[k]/6)*np.cos(2*np.pi*k*j/Q), np.sqrt(L[k]/6)*np.sin(2*np.pi*k*j/Q)]
cols += [np.sqrt(L[6]/12)*(-1.0)**j]
X=np.column_stack(cols); A7=X@X.T
# defect coefficients
idx=np.arange(1,12)
delta=A7[0,0]-A7[0,idx]
B=A7[np.ix_(idx,idx)]-A7[idx,0,None]-A7[0,idx][None,:]+A7[0,0]

def log_multinom_ratio(N,m):
    D=int(m.sum())
    # log[N!/(N-D)! / prod m!]
    return math.lgamma(N+1)-math.lgamma(N-D+1)-sum(math.lgamma(int(x)+1) for x in m)

def log_def_ratio(N,m,kappa=0.0):
    D=int(m.sum())
    return log_multinom_ratio(N,m)-G*float(delta@m)+(G/(2*N))*float(m@B@m)-kappa*D/N

rows=[]; validation=[]
for N in range(3,11):
    z=np.load(f'/mnt/data/defect_research_inputs/base_N{N}.npz',allow_pickle=True)
    states=z['states']; labels=z['labels']; pi=z['pi']
    j0=np.where(labels==0)[0]
    D=N-states[j0,0]
    base=pi[j0]
    # Seed state index [N,0,...]
    seed=np.where((states[:,0]==N)&(states[:,1:].sum(axis=1)==0))[0][0]
    lp0=math.log(pi[seed])
    # validate defect formula on every J0 state with nonzero pi
    errs=[]
    for ii in j0:
        m=states[ii,1:].astype(int)
        obs=math.log(pi[ii])-lp0
        pred=log_def_ratio(N,m,0)
        errs.append(abs(obs-pred))
    validation.append({'N':N,'max_abs_log_ratio_error':max(errs),'n_j0':len(j0)})
    for kap in [0,3,6,9,12]:
        lw=np.log(base)+kap*states[j0,0]/N
        lw-=lw.max(); ww=np.exp(lw); ww/=ww.sum()
        meanD=float(ww@D); varD=float(ww@((D-meanD)**2))
        dist={int(d):float(ww[D==d].sum()) for d in sorted(set(D.tolist()))}
        tails={str(K):float(ww[D>K].sum()) for K in range(0,min(10,N)+1)}
        rows.append({'N':N,'kappa':kap,'mean_D':meanD,'var_D':varD,'dist_D':dist,'tail':tails})
out={'g':G,'laplacian_eigenvalues':L.tolist(),'delta':delta.tolist(),'B':B.tolist(),'validation':validation,'rows':rows,
     'notes':['Tail P(D>K) is monotone nonincreasing in kappa because d/dk E[1_{D>K}]=-(1/N)Cov(1_{D>K},D)<=0. Therefore kappa=0 is the exact worst case over kappa>=0.']}
Path('/mnt/data/defect_research_inputs/DEFECT_PREPARATION_327.json').write_text(json.dumps(out,indent=2))
print('validation')
for v in validation: print(v)
print('\nWorst-case kappa=0 tail masses:')
for N in range(3,11):
 r=next(r for r in rows if r['N']==N and r['kappa']==0)
 print(N,'meanD',r['mean_D'],'K2',r['tail'].get('2'),'K3',r['tail'].get('3'),'K4',r['tail'].get('4'),'K6',r['tail'].get('6'))
