#!/usr/bin/env python3
"""Independent exact finite-state checks for FIN Z12 anonymous event histories.
Mathematical comparison model only; not a validation of physical FIN.
"""
import json, math
from pathlib import Path
import sympy as sp
import numpy as np

D=Path(__file__).resolve().parent
z=sp.Rational
r=sp.sqrt(3)
w=[z(16,15)+r/30,z(1),z(61,60),z(1),z(16,15)-r/30,z(1)]
v=[z(43,40)+r/40,z(59,60),z(1),z(1),z(43,40)-r/40,z(31,30)]
s=z(113,10)

# The generator comparison class is inherited; reconstruct it independently.
def jump(q):
    den=2*sum(q[:5])+q[5]
    assert den>0
    return sp.Matrix(12,12, lambda i,j:
                     z(0) if i==j else q[min((i-j)%12,(j-i)%12)-1]/den)

J=jump(w); K=jump(v)
assert all(sp.simplify(x)>0 for x in w+v)
assert all(sum(J[i,j] for j in range(12))==1 for i in range(12))
assert all(sum(K[i,j] for j in range(12))==1 for i in range(12))
assert J==J.T and K==K.T

# Return data for the previously frozen exact witness pair.
r2a=sp.trace(J**2)/12; r2b=sp.trace(K**2)/12
r3a=sp.trace(J**3)/12; r3b=sp.trace(K**3)/12
assert sp.simplify(r2a-r2b)==0
assert sp.simplify(r3a-r3b)==0

def alternating_p(q,n):
    den=2*sum(q[:5])+q[5]
    return sp.expand((2*sum(x**n for x in q[:5])+q[5]**n)/den**n)

A=alternating_p(w,3);B=alternating_p(v,3)
diff=sp.factor(A-B)
assert diff==sp.Rational(11,24*113**3)
assert sp.simplify(A-sum(J[0,j]**3 for j in range(12)))==0
assert sp.simplify(B-sum(K[0,j]**3 for j in range(12)))==0

# With mandatory unequal adjacent event states, the four-index equality
# partition has exactly five possibilities.
def partition_probs(m3,r2,r3):
    return {'ABAB':m3, 'ABAC':r2-m3, 'ABCB':r2-m3,
            'ABCA':r3, 'ABCD':1-2*r2-r3+m3}
pa=partition_probs(A,r2a,r3a); pb=partition_probs(B,r2b,r3b)
assert sum(pa.values())==sum(pb.values())==1
assert all(float(x)>0 for x in list(pa.values())+list(pb.values()))
tv=sp.simplify(sum(abs(pa[x]-pb[x]) for x in pa)/2)
assert tv==2*diff and tv==sp.Rational(11,12*113**3)

# Independent trajectory enumeration using a completely different route.
def brute_numeric(M):
    m=np.array(M.tolist(),dtype=float)
    pp={k:[] for k in pa}
    for i in range(12):
        for j in range(12):
            if i==j:continue
            for k in range(12):
                if k==j:continue
                for l in range(12):
                    if l==k:continue
                    val=m[i,j]*m[j,k]*m[k,l]/12
                    if i==k and j==l: nm='ABAB'
                    elif i==k: nm='ABAC'
                    elif j==l: nm='ABCB'
                    elif i==l: nm='ABCA'
                    else: nm='ABCD'
                    pp[nm].append(val)
    return {k:math.fsum(v) for k,v in pp.items()}
na=brute_numeric(J); nb=brute_numeric(K)
assert max(abs(na[k]-float(pa[k])) for k in pa)<4e-15
assert max(abs(nb[k]-float(pb[k])) for k in pb)<4e-15

# Newton's identities recover the full unordered 11-entry row of transition
# probabilities from consecutive alternating equality words through length 12.
def rowpoly_from_moments(q):
    p=[0]+[alternating_p(q,j) for j in range(1,12)]
    e=[sp.Integer(1)]
    for k in range(1,12):
        e.append(sp.expand(sum((-1)**(j-1)*e[k-j]*p[j] for j in range(1,k+1))/k))
    x=sp.Symbol('x')
    recovered=sp.Poly(sum((-1)**k*e[k]*x**(11-k) for k in range(12)),x)
    den=2*sum(q[:5])+q[5]
    actual=sp.Poly((x-q[5]/den)*sp.prod((x-q[i]/den)**2 for i in range(5)),x)
    assert sp.expand(recovered.as_expr()-actual.as_expr())==0
    return recovered

# Check with fully rational exact weights including a duplicate, avoiding
# inference based on generic distinctness.
q=[z(1),z(2),z(3),z(4),z(5),z(6)]
poly=rowpoly_from_moments(q)
poly_w=rowpoly_from_moments(w)

# All alternating moments can coincide but circular placement differs.
q2=[z(2),z(1),z(3),z(4),z(5),z(6)]
Jq=jump(q); Jq2=jump(q2)
assert all(alternating_p(q,n)==alternating_p(q2,n) for n in range(1,12))
nonalt=sp.factor(sp.trace(Jq**3-Jq2**3)/12)
assert nonalt==z(-11,7776)

# Fixed Fourier-position eigenvalues and exact Isospectrality of original
# witness are independently checked through trace powers n=1..12.
for n in range(1,13):
    assert sp.simplify(sp.trace(J**n)-sp.trace(K**n))==0

# Continuous-time alternating same-two-label history:
# probability T(t,s,u) = sum_d P_t(d)P_s(d)P_u(d), with cubic Taylor leading
# term m3*(s_rate)^3 t*s*u. The return-only d=0 term cancels exactly for
# equal heat traces; strict coefficient difference is 11/24000.
ct_coefficient=sp.simplify((s**3)*diff)
assert ct_coefficient==z(11,24000)

result={
 'scope':'pre-supplied single-label D12-symmetric positive jump Markov comparison model, NOT native interacting FIN or sourced time',
 'theorem_E3':{
   'three_index_equalities_indistinguishable':True,
   'four_index_ABAB_probability_A':str(A),
   'four_index_ABAB_probability_B':str(B),
   'difference_exact':str(diff),
   'difference_float':float(diff),
   'full_equality_record_TV_exact':str(tv),
   'full_equality_record_TV_float':float(tv),
   'partitions_A':{k:str(x) for k,x in pa.items()},
   'partitions_B':{k:str(x) for k,x in pb.items()},
   'enumerated_max_abs_error_A':max(abs(na[k]-float(pa[k])) for k in pa),
   'enumerated_max_abs_error_B':max(abs(nb[k]-float(pb[k])) for k in pb),
 },
 'theorem_E4':{
   'max_alternating_length_for_all_row_moments':12,
   'row_polynomial_rational_control':str(poly.as_expr()),
   'strict_ring_shell_identification':'normalised w6 as unique odd-multiplicity value; multiset of normalised w1..w5; no cyclic assignment',
   'isoperiodic_alternating_control_rational':True,
   'nonalternating_three_step_return_difference':str(nonalt)
 },
 'continuous_time_corollary':{
   'leading_t_s_u_coefficient_difference':str(ct_coefficient),
   'requires_given_common_model_time':True
 },
 'not_claimed':['intrinsic source of generator','gauge-free circle geometry','physical clock','native 144-state N2 HB','experimental validation','ToE']
}
(D/'RESULTS.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
print('PASS: four-event equality TV =',tv,'; alternating word moment polynomial and counterexample.')
print('ΔP(ABAB)=',diff,' ~',float(diff))
print('Brute enumeration max error',result['theorem_E3']['enumerated_max_abs_error_A'],result['theorem_E3']['enumerated_max_abs_error_B'])
