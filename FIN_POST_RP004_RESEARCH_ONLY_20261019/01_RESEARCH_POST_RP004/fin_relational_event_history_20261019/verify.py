#!/usr/bin/env python3
"""Exact algebra and bounded finite-matrix controls; NOT a formal proof assistant."""
from pathlib import Path
import itertools, json, hashlib
import sympy as S
import numpy as np
from scipy.linalg import expm
R=S.sqrt(3); z=S.Rational
B=S.Matrix([
 [2-R,1,2,3,2+R,2], [1,3,4,3,1,0], [2,4,2,0,2,2],
 [3,3,0,3,3,0], [2+R,1,2,3,2-R,2], [4,0,4,0,4,0]])
Bin=S.Matrix([
 [-R/12,-z(1,12),0,z(1,12),R/12,z(1,12)],
 [-z(1,12),z(1,12),z(1,6),z(1,12),-z(1,12),-z(1,12)],
 [0,z(1,6),0,-z(1,6),0,z(1,12)],
 [z(1,12),z(1,12),-z(1,6),z(1,12),z(1,12),-z(1,12)],
 [R/12,-z(1,12),0,z(1,12),-R/12,z(1,12)],
 [z(1,6),-z(1,6),z(1,6),-z(1,6),z(1,6),-z(1,12)]])
assert all(S.expand(x)==0 for x in (B*Bin-S.eye(6)))
assert all(S.expand(x)==0 for x in (Bin*B-S.eye(6)))
assert S.simplify(B.det()+3456*R)==0
l=S.Matrix([12+z(k,10) for k in range(1,7)])
p=S.Matrix([l[1],l[0],l[2],l[3],l[4],l[5]])
w=Bin*l; v=Bin*p
assert all(S.simplify(x-y)==0 for x,y in zip(w,[z(16,15)+R/30,1,z(61,60),1,z(16,15)-R/30,1]))
assert all(S.simplify(x-y)==0 for x,y in zip(v,[z(43,40)+R/40,z(59,60),1,1,z(43,40)-R/40,z(31,30)]))
assert min(float(x) for x in list(w)+list(v))>0
s=2*sum(w[:5])+w[5]; s2=2*sum(v[:5])+v[5]
assert s==z(113,10) and s2==s
assert all(S.expand(x)==0 for x in B*w-l)
assert all(S.expand(x)==0 for x in B*v-p)

def cubic(x): return 2*sum(q**3 for q in x[:5])+x[5]**3
cubic_difference=S.expand(cubic(w)-cubic(v))
assert cubic_difference==z(11,24000), cubic_difference
assert 12*cubic_difference==z(11,2000)

def makeQ(weights):
  weights=np.array([float(x) for x in weights]); q=np.zeros((12,12))
  for i in range(12):
    for j in range(12):
      if i!=j:
        q[i,j]=weights[min((i-j)%12,(j-i)%12)-1]
    q[i,i]=-q[i,:].sum()
  return q
qa=makeQ(w); qb=makeQ(v)
sa=float(s)
Ja=np.eye(12)+qa/sa; Jb=np.eye(12)+qb/sa
assert np.max(np.abs(np.diag(Ja)))<1e-12
assert np.max(np.abs(Ja.sum(1)-1))<1e-13
assert np.min(Ja)>-1e-13
assert np.max(np.abs(Jb.sum(1)-1))<1e-13
max_cont=0;max_event=0;max_poly=0;max_modal=0
for t in [0,.01,.05,.1,.3,1,2,5]:
    pa=expm(t*qa);pb=expm(t*qb)
    ea=np.trace(pa)/12;eb=np.trace(pb)/12
    exact=(1+2*sum(np.exp(-float(x)*t) for x in l[:5])+np.exp(-float(l[5])*t))/12
    max_cont=max(max_cont,abs(ea-eb),abs(ea-exact),abs(eb-exact))
for n in list(range(0,9)) +[14,22]:
    ma=np.linalg.matrix_power(Ja,n);mb=np.linalg.matrix_power(Jb,n)
    ea=np.trace(ma)/12;eb=np.trace(mb)/12
    exact=(1+2*sum((1-float(x)/sa)**n for x in l[:5])+(1-float(l[5])/sa)**n)/12
    max_event=max(max_event,abs(ea-eb),abs(ea-exact),abs(eb-exact))
for t in [.07,.9]:
    # bounded Poisson-series reconstruction
    import math
    poiss=sum(np.exp(-sa*t)*(sa*t)**n/math.factorial(n)*np.linalg.matrix_power(Ja,n) for n in range(70))
    max_poly=max(max_poly,np.max(abs(poiss-expm(qa*t))))
for i in range(1,7):
    j=np.arange(12)
    f=np.sqrt(2)*np.cos(2*np.pi*i*j/12) if i<6 else (-1.)**j
    eig=1-float(l[i-1])/sa
    max_modal=max(max_modal,float(np.linalg.norm(Ja@f-eig*f,np.inf)))
assert max(max_cont,max_event,max_poly,max_modal)<3e-11
# Diagnostic only: permutations in decimal-frozen strict FIN radial tuple
strict=np.array([.46998567264502017,.19204355169010282,.09142861427792497,.0470291687456504,.02413122336363006,.011070817321442113])
Bn=np.array(B.evalf(),float); ln=Bn@strict
pos=[]
for perm in itertools.permutations(range(5)):
  ll=np.array([ln[i] for i in perm]+[ln[5]])
  ww=np.linalg.solve(Bn,ll)
  if min(ww)>1e-11: pos.append(perm)
assert pos==[(0,1,2,3,4),(4,1,2,3,0)],pos
results={
 'scope':'12-state single-label radial symmetric Markov comparison family; not FIN native HB',
 'proof_checks':{'B_inverse_exact':True,'det_B':'-3456*sqrt(3)','counterpair_positive':True,
 'isospectral_weight_permutation':'swap_lambda_1_lambda_2','counterpair_nonisomorphic_cubic_difference':'11/24000',
 'signed_cubic_full_directed_edge_difference':'11/2000','common_exit':'113/10'},
 'counterpair_w':[str(x) for x in w], 'counterpair_v':[str(x) for x in v],
 'sampled_numeric_control':{'max_heat_trace_error':max_cont,'max_jump_return_error':max_event,
 'max_poissonization_error':max_poly,'max_modal_eigen_error':max_modal},
 'decimal_strict_diagnostic_positive_permutations':[list(x) for x in pos],
 'limitations':['No inference on fundamental state/time generator','Unlabelled two-time return statistics only, not complete unlabelled multi-time law','Strict-decimal permutation enumeration is not a certificate for true transcendental source weights']
}
Path(__file__).with_name('RESULTS.json').write_text(json.dumps(results,indent=2,ensure_ascii=False)+'\n')
print('PASS exact symbolic + bounded finite matrix tests; counterpair cubic edge invariant =',cubic_difference)
print(json.dumps(results['sampled_numeric_control'],indent=2))
