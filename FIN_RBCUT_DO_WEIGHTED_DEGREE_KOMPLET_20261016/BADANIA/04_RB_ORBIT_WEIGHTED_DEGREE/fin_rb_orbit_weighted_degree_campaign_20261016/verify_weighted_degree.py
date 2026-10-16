#!/usr/bin/env python3
"""Exact finite rational checks of the abstract weighted-cut reduction.
Artificial symmetric graph only; DOES NOT instantiate true MT q12.
"""
from fractions import Fraction as F
from collections import defaultdict
from itertools import product
import json
from pathlib import Path
root=Path(__file__).resolve().parent
O=[(level,phase) for level in range(3) for phase in 'AB']
V=[(level,phase,j) for level,phase in O for j in range(6)]
H={};U={}
def add(x,y,w,allowed):
    assert x[0]+1==y[0]
    assert (x,y) not in U
    U[x,y]=F(w)
    if allowed:H[x,y]=F(w)
for p in 'AB':
    other='B' if p=='A' else 'A'
    for j in range(6):
        for off,w in [(0,F(1)),(1,F(1,4)),(-1,F(1,4))]:
            add((0,p,j),(1,p,(j+off)%6),w,True)
        for off,w in [(0,F(1,3)),(1,F(1,6)),(-1,F(1,6))]:
            add((0,p,j),(1,other,(j+off)%6),w,False)
        for target in 'AB':add((1,p,j),(2,target,j),F(1),True)
W={x:F(1) for x in V}
s=defaultdict(F)
for (x,y),w in U.items():s[x]+=w;s[y]-=w
N=12
p_min=(F(1)-F(1,100))/12-F(1,10000)
p_max=(F(1)+F(1,100))/12+F(1,10000)
assert p_min==F(103,1250) and p_max==F(158,1875)
# The real FIN reference has the independently derived conservative bound
# a_0=(103/15000)*exp(-24); not an exact rational quantity.
a=min(U.values())/N
# All native q^0>=a N, native total rate <= N.
assert all(w>=a*N*W[x] for (x,y),w in U.items())
assert all(sum((w for (u,v),w in U.items() if x in (u,v)),F())<=N*W[x] for x in V)
assert all(s[(0,p,j)]==F(13,6) for p in 'AB' for j in range(6))
assert all(s[(1,p,j)]==F(-1,6) for p in 'AB' for j in range(6))
assert all(s[(2,p,j)]==F(-2) for p in 'AB' for j in range(6))

def orbit(x):return x[:2]
So={o:sum((s[x] for x in V if orbit(x)==o),F()) for o in O}
Uh=defaultdict(F);edgecount=defaultdict(int)
for (x,y),w in H.items():
    Uh[orbit(x),orbit(y)]+=w
    edgecount[orbit(x),orbit(y)]+=1

def cuts():
    data=[]
    for mask in range(1<<len(O)):
        J={o for i,o in enumerate(O) if mask & (1<<i)}
        supply=sum((So[o] for o in J),F())
        native_cap=sum((w for (a,b),w in Uh.items() if a in J and b not in J),F())
        weighted_out=sum((F(cnt) for (a,b),cnt in edgecount.items() if a in J and b not in J),F())
        data.append((J,supply,native_cap,weighted_out))
    return data

def rat(num,den):
    if num<=0:return F(0)
    if den==0:return None
    return num/den
ratiosK=[rat(s,uh) for J,s,uh,b in cuts()]
ratiosQ=[rat(s,N*b) for J,s,uh,b in cuts()]
assert all(v is not None for v in ratiosK+ratiosQ)
K=max(ratiosK);Q=max(ratiosQ)
assert K==F(13,9)
# Exact primal flow from prior independent symmetric toy.
f={}
for p in 'AB':
    f[(0,p),(1,p)]=F(13)
    for tgt in 'AB':f[(1,p),(2,tgt)]=F(6)
assert all(f[key]<=K*Uh[key] for key in f)
ce={(x,y):f[orbit(x),orbit(y)]*w/Uh[orbit(x),orbit(y)] for (x,y),w in H.items()}
ss=defaultdict(F)
for (x,y),c in ce.items():ss[x]+=c;ss[y]-=c
assert all(ss[x]==s[x] for x in V)
C=max(sum((val for e,val in ce.items() if x in e),F())/(N*W[x]) for x in V)
assert C==F(25,72)
# Dual phi=1-level and beta=1 at middle gives C*>=25/72.
phi={x:F(1-x[0]) for x in V};beta={x:F(x[0]==1) for x in V}
assert all(phi[x]-phi[y]<=beta[x]+beta[y] for x,y in H)
Cdual=sum(phi[x]*s[x] for x in V)/(N*sum(beta[x]*W[x] for x in V))
assert Cdual==C
assert Q<=C<=K<=Q/a, (Q,C,K,a,Q/a)
assert a*K<=C
# The fundamental pointwise inequality c_e<=C N W_x holds for this flow.
assert all(c<=C*N*W[x] and c<=C*N*W[y] for (x,y),c in ce.items())
# Opposite example: no uniform ellipticity, unbounded K but uniformly capped rates.
# Four vertices 0->2->3 native references t=1/N, deleted 0->1->3 refs 1.
nonuniform=[]
for n in (2,12,120,1200):
    t=F(1,n);native=[(0,2,t),(2,3,t),(0,1,F(1)),(1,3,F(1))]
    source=[F(0)]*4
    for x,y,c in native:source[x]+=c;source[y]-=c
    assert source==[1+t,F(0),F(0),-1-t]
    Kt=(1+t)/t
    Ct=(1+t)/n
    nonuniform.append({'N':n,'K_star':str(Kt),'C_feasible':str(Ct),'a_local':str(t/n)})
    assert Kt==n+1 and Ct<=1
results={'scope':'ARTIFICIAL_D6_36_VERTEX_TOY_AND_NONUNIFORM_4_VERTEX_NEGATIVE_CONTROL_NOT_REAL_Q12_MT',
 'checks':{
    'native_ellipticity':True,'exact_MT_occupation_floor':True,'native_clock_bound':True,'full_orbit_cuts_64':True,
    'sharp_relative_cut_K':True,'weighted_incidence_cut_Q':True,
    'Q_le_C_le_K_le_Q_over_a':True,'primal_exact_vertexwise':True,
    'C_dual_optimal':True,'nonuniform_lower_bound_needed':True},
 'exact':{'p_min_true_MT':str(p_min),'p_max_true_MT':str(p_max),'a':str(a),'Q':str(Q),'Cstar':str(C),'Kstar':str(K),'Q_over_a':str(Q/a),'alpha_K':str(a*K),'number_of_orbit_cuts':2**len(O)},
 'negative_control':nonuniform}
assert all(results['checks'].values())
(root/'FINITE_CHECKS.json').write_text(json.dumps(results,indent=2,ensure_ascii=False)+'\n')
print('PASS',len(results['checks']),'test families; a',a,'Q',Q,'Cstar',C,'Kstar',K)
