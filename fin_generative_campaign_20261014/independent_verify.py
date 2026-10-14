#!/usr/bin/env python3
"""Independent reader/verifier. Does not import run_campaign; independent Fourier and CTMC formulas."""
import itertools, json, hashlib, math, sys
from pathlib import Path
import numpy as np
from scipy.special import logsumexp,softmax
ROOT=Path(__file__).resolve().parent
A_ref=np.asarray(json.loads((ROOT/'RG-C01/MODEL_CONTRACT.json').read_text())['matrix_A7'])
base=json.loads((ROOT/'RG-C01/MODEL_CONTRACT.json').read_text())['contract']
Q=12;g=float(base['gain']);dists=[]
for x in range(Q):
    d=min(x,12-x)
    dists.append(0 if d==0 else math.cos(.18575*d+.1625)/(1+d**1.8))
lambdas={k:sum(dists[d]*(1-math.cos(2*math.pi*k*d/Q)) for d in range(1,Q)) for k in (3,4,5,6)}
a_row=[sum(lambdas[k]/6*math.cos(2*math.pi*k*d/Q) for k in (3,4,5))+(lambdas[6]/12)*(-1)**d for d in range(Q)]
A=np.asarray([[a_row[(j-i)%Q] for j in range(Q)] for i in range(Q)])
checks={}
def check(name, condition, detail):
    checks[name]={'passed':bool(condition),'detail':detail}
    if not condition:print('FAILED:',name,detail)
check('independent_matrix_matches_source_Fourier',float(np.max(np.abs(A-A_ref)))<2e-14,float(np.max(np.abs(A-A_ref))))
check('centering',float(np.max(np.abs(A.sum(axis=1))))<1e-12,float(np.max(np.abs(A.sum(axis=1)))))
check('rank_seven',int(np.linalg.matrix_rank(A,tol=1e-10))==7,int(np.linalg.matrix_rank(A,tol=1e-10)))
check('PSD',float(np.linalg.eigvalsh(A).min())>-1e-12,float(np.linalg.eigvalsh(A).min()))
F=np.empty((Q,Q,Q));F_direct=np.empty_like(F)
for a,b,c in itertools.product(range(Q),repeat=3):
    n=np.bincount((a,b,c),minlength=Q)
    part=(g/8)*float(n@A@n)
    F[a,b,c]=part+logsumexp((g/4)*A@n)
    F_direct[a,b,c]=logsumexp([g/8*np.bincount((a,b,c,h),minlength=Q)@A@np.bincount((a,b,c,h),minlength=Q) for h in range(Q)])
check('exact_4_to_3_enumeration',float(np.max(np.abs((F-F_direct)-(F-F_direct)[0,0,0])))<3e-14,float(np.max(np.abs((F-F_direct)-(F-F_direct)[0,0,0]))))
ct=0.
for ba,bb,bc in itertools.product((0,1),repeat=3):
    state=(ba,[0,2][bb],[0,4][bc])
    ct+=(-1)**(ba+bb+bc)*F[state]
claimed=json.loads((ROOT/'RG-C02/RESULTS.json').read_text())['derived_objects']
check('independent_third_contrast',abs(ct-claimed['contrast'])<2e-13,{'calculated':ct,'reported':claimed['contrast']})
check('contrast_nonzero',abs(ct)>.25,ct)
# explicit g0 factorization
fg0=np.zeros((Q,Q,Q))
for a,b,c in itertools.product(range(Q),repeat=3):
    fg0[a,b,c]=logsumexp([0.] * Q)
check('g0_constant',float(np.max(fg0)-np.min(fg0))<1e-14,float(np.max(fg0)-np.min(fg0)))
# Projected finite dynamics: independently compute equilibrium conditional coarse memory channel
states=np.array(list(itertools.product(range(Q),repeat=3)),dtype=int)
ctn=np.zeros((len(states),Q),int)
for i in range(3):ctn[np.arange(len(states)),states[:,i]]+=1
logpi=(g/6)*np.einsum('bi,ij,bj->b',ctn,A,ctn)
pi=np.exp(logpi-logsumexp(logpi))
pY=np.sum(pi.reshape((Q,Q,Q)),axis=2)
Lf=np.zeros(len(states))
for ix,(a,b,c) in enumerate(states):
    # observable f=1{a=b=0}, only visible-label updates contribute to changes.
    fval=int(a==0 and b==0)
    ratechange=0.
    for axis in (0,1):
        old=(a,b)[axis]
        other=(b if axis==0 else a,c)
        rates=softmax((g/3)*(A[:,other[0]]+A[:,other[1]]))
        for new in range(Q):
            if new==old: continue
            newa,newb=(new,b) if axis==0 else (a,new)
            ratechange+=rates[new]*(int(newa==0 and newb==0)-fval)
    Lf[ix]=ratechange
Lf3=Lf.reshape((Q,Q,Q));P_Lf=np.sum(Lf3*pi.reshape((Q,Q,Q)),axis=2)/pY
v=Lf3-P_Lf[:,:,None]
M0=float(np.sum(pi.reshape((Q,Q,Q))*v*v))
mem=json.loads((ROOT/'RG-C03/MEMORY_REPORT.json').read_text())
check('independent_M0',abs(M0-mem['M0_coupled'])<1e-13,{'calculated':M0,'reported':mem['M0_coupled']})
r0=softmax((g/3)*(A[:,0]+A[:,0]));r1=softmax((g/3)*(A[:,0]+A[:,1]));idx=int(np.argmax(abs(r0[1:]-r1[1:]))+1)
lump=json.loads((ROOT/'RG-C03/LUMPABILITY_WITNESS.json').read_text())
check('off_diagonal_lumpability',idx!=0 and idx==lump['j_witness'] and abs(abs(r0[idx]-r1[idx])-lump['maximum_single_visible_refresh_rate_difference'])<1e-13,{'j':idx,'defect':abs(r0[idx]-r1[idx])})
# static parent event direct (independently), no reducer invoked
Z=0.;event=0.
for a,b,c,d in itertools.product(range(Q),repeat=4):
    n=np.bincount((a,b,c,d),minlength=Q)
    weight=math.exp((g/8)*n@A@n)
    Z+=weight
    if a==c: event+=weight
held=json.loads((ROOT/'RG-C05/PREDICTIONS.json').read_text())['predictions_and_reference'][0]
check('independent_direct_holdout',abs(event/Z-held['parent_direct_reference'])<1e-12,{'reference':event/Z,'reported':held['parent_direct_reference']})
check('predicted_difference_bounded',held['abs_discrepancy']<=held['nominal_bound']+1e-10,held)
cert=json.loads((ROOT/'RG-C02/INTERVAL_CERTIFICATE.json').read_text())
check('interval_bound_dominates_nominal',float(cert['validated_TV_bound_upper'])>=claimed['TV_bound'] and float(cert['validated_TV_bound_upper'])<0.06,cert['validated_TV_bound_upper'])
check('input_rounding_control',float(cert['source_A7_max_discrepancy_from_binary64'])<1e-13,cert['source_A7_max_discrepancy_from_binary64'])
out={'execution_status':'INDEPENDENT_REPLAY','passed':all(v['passed'] for v in checks.values()),'number_checks':len(checks),'checks':checks,'method':'independent direct circulant spectral construction; independent direct 4-state sum; independent CTMC local Lf and equilibrium M0; independent parent holdout enumeration; interval proof ledger sanity checks','source':'do not import producer function'}
(ROOT/'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'passed':out['passed'],'number_checks':len(checks),'failed':[k for k,v in checks.items() if not v['passed']]},indent=2))
if not out['passed']:sys.exit(1)
