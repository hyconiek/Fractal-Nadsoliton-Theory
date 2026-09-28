import json, math, numpy as np
from pathlib import Path
from scipy.linalg import expm
ROOT=Path('/mnt/data/work331');BURN=24.;ETA=.05;KGRID=np.arange(13.)
vals=np.round(np.cos(2*np.pi*np.arange(12)/12),12);yv=np.array(sorted(set(vals)));l2b=np.array([np.where(yv==vals[j])[0][0] for j in range(12)],int);sizes=np.bincount(l2b,minlength=7)
B=np.zeros((12,7));B[np.arange(12),l2b]=1;C=np.zeros((7,7))
for a in range(7):
 for b in range(7): C[a,b]=(1-ETA)+(ETA*(sizes[b]-1)/11) if a==b else ETA*sizes[b]/11
E=B@C

def tv(a,b):return float(.5*np.abs(np.asarray(a)-np.asarray(b)).sum())
def Q12(R,rho):
 MM=np.zeros((6,6))
 for k in range(1,7):
  for d in range(1,6): MM[k-1,d-1]=2*(math.cos(2*math.pi*k*d/12)-1)
  MM[k-1,5]=math.cos(math.pi*k)-1
 q=np.linalg.solve(MM,-rho*np.asarray(R));Q=np.zeros((12,12))
 for a in range(12):
  for d in range(1,6): Q[a,(a+d)%12]+=q[d-1];Q[a,(a-d)%12]+=q[d-1]
  Q[a,(a+6)%12]+=q[5]
 np.fill_diagonal(Q,-Q.sum(1));return Q

def params(N):
 if N==7:
  A=json.load(open('/mnt/data/ex306/fin306/slip_training_rows_306.json'))['allrows'];rhoD={3:.131439786195646,4:.071692277481914,5:.0401867559160965,6:.0226112751871144};tr=[3,4,5,6]
  rho=float(np.exp(np.polyval(np.polyfit(tr,np.log([rhoD[n] for n in tr]),1),7)));R=[]
  for ki in range(6):R.append(1. if ki==3 else float(np.polyval(np.polyfit(tr,[A[str(n)]['R_eff'][ki] for n in tr],1),7)))
  return rho,.5/rho,Q12(R,rho)
 if N==8:
  d=json.load(open('/mnt/data/ex312/fin312/PHASE_A_FROZEN_BEFORE_N8_312.json'))['N8_prediction'];rho=d['rho_pred'];return rho,.5/rho,Q12(d['R_pred'],rho)
 if N==9:
  d=json.load(open('/mnt/data/ex315/fin315/PHASE_A_FROZEN_BEFORE_N9_315.json'));return d['rho_model']['rho9_pred'],d['transition']['dt_pred'],np.array(d['transition']['Qpred'])
 if N==10:
  d=json.load(open('/mnt/data/ex317/fin317/PHASE_A_FROZEN_BEFORE_N10_317.json'));return d['rho_model']['rho10_pred'],d['transition']['dt_pred'],np.array(d['transition']['Qpred'])

def get_spec(N,rank):
 if N==7:p=ROOT/f'N7_SPEC_R{rank}.npz'
 elif N==8:p=ROOT/'fin312/N8_SPECTRAL_CACHE_312.npz'
 elif N==9:p=ROOT/'fin315/N9_SPECTRAL_CORE_315.npz'
 else:p=ROOT/'N10_SPEC_R16.npz'
 z=np.load(p);return z['eigvals'][:rank],z['eigvecs'][:,:rank]

def run_rank(N,rank):
 z=np.load(ROOT/f'base_N{N}.npz');states=z['states'];labels=z['labels'];pi=z['pi'];D=N-states[:,0];j0=np.where(labels==0)[0] if 'j0' not in z.files else z['j0'];n=len(states)
 lam,U=get_spec(N,rank);sp=np.sqrt(pi);rho,dt,Qp=params(N)
 M=[];Mt=[];tails=[]
 for k in KGRID:
  x=np.log(pi[j0])+k*states[j0,0]/N;x-=x.max();w=np.exp(x);w/=w.sum();mu=np.zeros(n);mu[j0]=w;M.append(mu);tails.append(float(w[D[j0]>6].sum()))
  keep=j0[D[j0]<=6];x=np.log(pi[keep])+k*states[keep,0]/N;x-=x.max();w=np.exp(x);w/=w.sum();mu=np.zeros(n);mu[keep]=w;Mt.append(mu)
 M=np.array(M);Mt=np.array(Mt)
 # O: 8 observed outputs, U explicit
 O=np.zeros((n,8));good=labels>=0;O[good,:7]=E[labels[good]];O[~good,7]=1
 # F latent 12 + U t1
 F=np.zeros((n,13));F[good,labels[good]]=1;F[~good,12]=1
 # spectral coefficients and propagations
 def ptime(M0,t): return (((M0/sp[None,:])@U*np.exp(lam*t)[None,:])@U.T)*sp[None,:]
 P1=ptime(M,BURN);P1t=ptime(Mt,BURN)
 Bc=U.T@(sp[:,None]*O);G=((U*np.exp(lam*dt)[None,:])@Bc)/sp[:,None]
 def joint(P):
  J=np.zeros((13,8,8))
  for a in range(8):J[:,a,:]=(P*O[:,a][None,:])@G
  # rank trunc may lose min mass, renorm each joint to probability
  J/=J.sum((1,2))[:,None,None];return J
 J=joint(P1);Jt=joint(P1t)
 prior=P1@F;priort=P1t@F;prior/=prior.sum(1)[:,None];priort/=priort.sum(1)[:,None]
 # normalized localized 12 priors and same effective transition
 pf=prior[:,:12]/prior[:,:12].sum(1)[:,None];pt=priort[:,:12]/priort[:,:12].sum(1)[:,None]
 Pe=expm(Qp*dt)
 def eff(p):
  J12=np.einsum('ni,ij->nij',p,Pe);return np.einsum('ia,nij,jb->nab',E,J12,E,optimize=True)
 Efull=eff(pf);Etr=eff(pt)
 # explicit U hybrid: t1 U absorbing, localized prior unnormalized from trunc
 Jhy=np.zeros((13,8,8))
 for c in range(13):
  p12=priort[c,:12];J12=np.einsum('i,ij->ij',p12,Pe);Jhy[c,:7,:7]=np.einsum('ia,ij,jb->ab',E,J12,E,optimize=True);Jhy[c,7,7]=priort[c,12]
  Jhy[c]/=Jhy[c].sum()
 return {'rank':rank,'last_lambda':float(lam[-1]),'tail_max':max(tails),'tail_by_kappa':tails,
  'defect_micro_joint_max':max(tv(J[i],Jt[i]) for i in range(13)),'defect_micro_joint':[tv(J[i],Jt[i]) for i in range(13)],
  'prior13_max':max(tv(prior[i],priort[i]) for i in range(13)),'prior_loc12_max':max(tv(pf[i],pt[i]) for i in range(13)),
  'effective_prep_component_max':max(tv(Efull[i],Etr[i]) for i in range(13)),
  'U_t1_full_max':float(prior[:,12].max()),'U_t1_trunc_max':float(priort[:,12].max()),
  'hybrid_explicitU_vs_micro_max':max(tv(J[i],Jhy[i]) for i in range(13)), 'hybrid_by_kappa':[tv(J[i],Jhy[i]) for i in range(13)]}

def main():
 out={'task':331,'rows':[]}
 ranks={7:[32,48],8:[24,32],9:[16,24]}
 for N in [7,8,9]:
  print('N',N,flush=True); rr=[]
  for r in ranks[N]:
   x=run_rank(N,r);rr.append(x);print(r,json.dumps(x,indent=2),flush=True)
  out['rows'].append({'N':N,'low':rr[0],'high':rr[1],
    'rank_convergence':{k:abs(rr[1][k]-rr[0][k]) for k in ['defect_micro_joint_max','prior13_max','prior_loc12_max','effective_prep_component_max','hybrid_explicitU_vs_micro_max']}})
 json.dump(out,open(ROOT/'DEFECT_OBSERVED_331_N7_N9.json','w'),indent=2)
main()
