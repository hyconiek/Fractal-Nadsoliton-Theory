import numpy as np, json, math
from pathlib import Path
ROOT=Path('/mnt/data/work331');BURN=24.;KGRID=np.arange(13.)
# use highest available ranks
specs={7:(ROOT/'N7_SPEC_R48.npz',48),8:(ROOT/'fin312/N8_SPECTRAL_CACHE_312.npz',32),9:(ROOT/'fin315/N9_SPECTRAL_CORE_315.npz',24),10:(ROOT/'N10_SPEC_R16.npz',16)}
out={'task':336,'burn':BURN,'rows':[]}
for N,(spfile,rank) in specs.items():
 z=np.load(ROOT/f'base_N{N}.npz');states=z['states'];labels=z['labels'];pi=z['pi'];sqrtpi=np.sqrt(pi);D=N-states[:,0]
 s=np.load(spfile);lam=s['eigvals'][:rank];U=s['eigvecs'][:,:rank]
 F=np.zeros((len(labels),13));good=labels>=0;F[good,labels[good]]=1;F[~good,12]=1
 # backward response to outcomes at burn: G=e^{Q t}F
 Bc=U.T@(sqrtpi[:,None]*F);G=((U*np.exp(lam*BURN)[None,:])@Bc)/sqrtpi[:,None]
 ids=np.where((labels==0)&(D<=6))[0];resp=G[ids]; # may tiny neg due spectral trunc; retain raw and normalized clipped diagnostic
 # derive predicted priors for all kappa by exact preparation weights on ids
 pri=[]
 for k in KGRID:
  lw=np.log(pi[ids])+k*states[ids,0]/N;lw-=lw.max();w=np.exp(lw);w/=w.sum();p=w@resp;p=np.maximum(p,0);p/=p.sum();pri.append(p)
 pri=np.array(pri)
 np.savez_compressed(ROOT/f'DEFECT_RESPONSE_N{N}_336.npz',state_ids=ids,states=states[ids],responses=resp,pi=pi[ids])
 out['rows'].append({'N':N,'rank':rank,'retained_states':len(ids),'response_min':float(resp.min()),'response_rowsum_maxdev':float(np.max(np.abs(resp.sum(1)-1))), 'last_lambda':float(lam[-1]), 'priors':pri.tolist()})
 print(N,len(ids),resp.min(),np.max(np.abs(resp.sum(1)-1)),flush=True)
json.dump(out,open(ROOT/'DEFECT_RESPONSE_336.json','w'),indent=2)
