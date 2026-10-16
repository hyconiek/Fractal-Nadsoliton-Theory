#!/usr/bin/env python3
"""Standalone exact rational bounded controls for KR, not a stochastic-limit proof."""
import itertools, json, hashlib
from fractions import Fraction as F
from pathlib import Path

q,N=4,2
sg=(1,1,-1,-1)
states=[s for s in itertools.product(range(N+1),repeat=q) if sum(s)==N]
ids={s:i for i,s in enumerate(states)}
m=lambda s: sum(a*b for a,b in zip(s,sg))*F(1,N)
weight=lambda s: F(2 if max(s)==1 else 1,16) # N=2, iid uniform multinomial
assert sum(map(weight,states))==1
edges={}
for i,s in enumerate(states):
 for a in range(q):
  if not s[a]:continue
  for b in range(q):
   if a==b:continue
   ss=list(s);ss[a]-=1;ss[b]+=1;j=ids[tuple(ss)]
   edges[tuple(sorted((i,j)))] = True
# base iid-uniform independent-copy refresh rate n_a/q
base={}
for (i,j) in edges:
 s,t=states[i],states[j]
 aa=next(a for a in range(q) if s[a]>t[a]);bb=next(b for b in range(q) if s[b]<t[b]);
 base[(i,j)]=weight(s)*F(s[aa],q)
 assert base[(i,j)] == weight(t)*F(t[bb],q)
# A rational, divergence-free parity drift wheel on reversible CONDUCTANCES.
wheel={(0,6):F(1),(0,3):F(-1),(1,6):F(-1),(1,4):F(1),(3,8):F(-1),(4,8):F(1)}
assert set(wheel)<=set(base)
epsilon=min(base[e] for e in wheel)/4
pert={e:base[e]+epsilon*wheel.get(e,0) for e in base}
assert all(x>0 for x in pert.values())

def measures(conductances):
 drift=[F(0) for _ in states];gamma=[F(0) for _ in states];capacity=[F(0) for _ in states]
 for (i,j),c in conductances.items():
  delta=m(states[j])-m(states[i]);
  drift[i]+=c/weight(states[i])*delta;drift[j]-=c/weight(states[j])*delta
  gamma[i]+=c/weight(states[i])*delta**2;gamma[j]+=c/weight(states[j])*delta**2
  capacity[i]+=c/weight(states[i]);capacity[j]+=c/weight(states[j])
 return drift,gamma,capacity

b0,G0,cap0=measures(base)
b1,G1,cap1=measures(pert)
assert b0==b1
assert G0!=G1
assert all(x<=N for x in cap1) # independently controlled clock capacity
assert sum(weight(states[i])*G1[i] for i in range(len(states))) == sum(weight(states[i])*G0[i] for i in range(len(states)))
assert sum(weight(states[i])*G1[i] for i in range(len(states))) == -2*sum(weight(states[i])*m(states[i])*b1[i] for i in range(len(states)))
# Same-parity thinning (non-wheel) retains parity Gamma and drift exactly.
eta=F(1,3)
thin={e:(c*eta if m(states[e[0]])==m(states[e[1]]) else c) for e,c in base.items()}
b2,G2,cap2=measures(thin)
assert b0==b2 and G0==G2 and base!=thin
assert all(x<=N for x in cap2)
# Generator symmetry check every oriented edge for wheel + thinning.
for conductances in [base,pert,thin]:
 for (i,j),c in conductances.items():
  qij=c/weight(states[i]);qji=c/weight(states[j]);assert weight(states[i])*qij == weight(states[j])*qji
# Jump and bracket scaling: delta Z=2 N^(-3/4), integrated accelerated qv= N Gamma(m).
statuses={
 'state_count':len(states),'edge_count':len(edges),'wheel_support':len(wheel),
 'wheel_epsilon':str(epsilon),'stationary_gamma':str(sum(weight(states[i])*G1[i] for i in range(len(states)))),
 'max_absolute_gamma_gap':str(max(abs(a-b) for a,b in zip(G0,G1))),
 'wheel_gamma_nonzero_states':sum(a!=b for a,b in zip(G0,G1)),
 'drift_exact_equality':True,'averaged_gamma_exact_equality':True,
 'same_parity_thinning_drift_and_gamma_equality':True,'rate_cap_per_state':str(N),
 'source_scope':'q4/N2 rational diagnostic only; NOT q12/A7 critical countermodel',
 'conclusions':'Statewise bracket not determined by stationary Gibbs+exact parity drift. Different critical limit not established.'
}
root=Path(__file__).resolve().parent
(root/'FINITE_EXACT_CHECKS.json').write_text(json.dumps(statuses,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(statuses,indent=2))
