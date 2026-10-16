#!/usr/bin/env python3
"""Exact rational orbit flow / dual fixtures (NOT the actual q12 MT domain)."""
from fractions import Fraction as F
from collections import defaultdict
from itertools import product
import json
from pathlib import Path

P=Path(__file__).resolve().parent
phases=('A','B')
V=[(l,p,j) for l in range(3) for p in phases for j in range(6)]

def edge_key(x,y):
    assert x[0]+1==y[0]
    return (x,y)

C0={}
H={}

def put(x,y,w,allowed):
    k=edge_key(x,y); assert k not in C0
    C0[k]=F(w)
    if allowed: H[k]=F(w)

for p in phases:
    other='B' if p=='A' else 'A'
    for j in range(6):
        for step,w in [(0,F(1)),(1,F(1,4)),(-1,F(1,4))]:
            put((0,p,j),(1,p,(j+step)%6),w,True)
        for step,w in [(0,F(1,3)),(1,F(1,6)),(-1,F(1,6))]:
            put((0,p,j),(1,other,(j+step)%6),w,False)
        for target in phases:
            put((1,p,j),(2,target,j),F(1),True)

def div(C):
    d=defaultdict(F)
    for (x,y),v in C.items():d[x]+=v; d[y]-=v
    return {x:d[x] for x in V}

s=div(C0)
assert all(s[(0,p,j)]==F(13,6) for p in phases for j in range(6))
assert all(s[(1,p,j)]==-F(1,6) for p in phases for j in range(6))
assert all(s[(2,p,j)]==-F(2) for p in phases for j in range(6))
assert sum(s.values())==0

# D6 action on six copy indices is the analogue of the parity-preserving subgroup.
for sign,shift in product((1,-1),range(6)):
    g=lambda x:(x[0],x[1],(sign*x[2]+shift)%6)
    assert all(C0.get((g(x),g(y)))==c for (x,y),c in C0.items())
    assert all(H.get((g(x),g(y)))==c for (x,y),c in H.items())
    assert all(s[g(x)]==s[x] for x in V)

orbs=[(l,p) for l in range(3) for p in phases]
def orb(x):return (x[0],x[1])
S={o:sum((s[x] for x in V if orb(x)==o),F(0)) for o in orbs}
Hbar=defaultdict(F)
for (x,y),v in H.items():Hbar[(orb(x),orb(y))]+=v

def cut(Sset):
    demand=sum((S[x] for x in Sset),F(0))
    capacity=sum((v for (a,b),v in Hbar.items() if a in Sset and b not in Sset),F(0))
    return demand, capacity
max_ratio=F(0)
max_sets=[]
for mask in range(1<<len(orbs)):
    pick={o for i,o in enumerate(orbs) if mask>>i &1}
    demand,capacity=cut(pick)
    assert not (demand>0 and capacity==0),('invalid ideal',pick)
    if demand>0:
        r=demand/capacity
        if r>max_ratio:max_ratio=r;max_sets=[pick]
        elif r==max_ratio:max_sets.append(pick)
assert max_ratio==F(13,9),max_ratio

# A quotient flow saturating the worst lower-layer cut.
fbar={}
for p in phases:
    fbar[((0,p),(1,p))]=F(13) # 6 copies of 13/6
    for target in phases:fbar[((1,p),(2,target))]=F(6)
assert all(0<=v<=max_ratio*Hbar[e] for e,v in fbar.items())
# The exact lifting formula c_e=f_OO' (c0_e/C_OO').
c={}
for e,cap in H.items():
    x,y=e
    c[e]=fbar[(orb(x),orb(y))]*cap/Hbar[(orb(x),orb(y))]
assert div(c)==s
assert all(0<=c[e]<=max_ratio*H[e] for e in H)
assert all(sum((v for e,v in c.items() if x in e),F(0))<=F(25,6) for x in V)
assert any(sum((v for e,v in c.items() if x in e),F(0))==F(25,6) for x in V)

# Vertex rate cap dual certificate phi(level)=(1,0,-1), beta(mid)=1;
# same inequality on each allowed directed edge.
phi={x:F(1-x[0]) for x in V}
beta={x:F(int(x[0]==1)) for x in V}
assert all(phi[x]-phi[y]<=beta[x]+beta[y] for x,y in H)
numer=sum((phi[x]*s[x] for x in V),F(0))
weighted=sum(beta.values(),F(0))
N=12
assert numer==F(50) and weighted==F(12)
Cstar=numer/(N*weighted)
assert Cstar==F(25,72)
assert max(sum((v for e,v in c.items() if x in e),F(0)) for x in V)/N==Cstar

# Even a source induced by a reversible native flow can be made impossible
# after support restriction: fully block one source orbit.
Hbad={e:v for e,v in H.items() if orb(e[0])!=(0,'A')}
cut_source=sum((s[x] for x in V if orb(x)==(0,'A')),F(0))
cut_capacity=sum((v for (x,y),v in Hbad.items() if orb(x)==(0,'A')),F(0))
assert cut_source==F(13) and cut_capacity==0

flags={
 'native_D6_invariance_12_maps':True,
 'native_divergence_from_full_graph':div(C0)==s,
 'exact_64_orbit_cuts':True,
 'sharp_relative_K':max_ratio==F(13,9),
 'exact_quotient_primal':True,
 'exact_lift_all_36_vertices':div(c)==s,
 'native_relative_cap':True,
 'exact_vertex_C_dual':Cstar==F(25,72),
 'exact_vertex_C_primal':True,
 'zero_capacity_negative_control':cut_capacity==0 and cut_source>0,
}
def q(x):return str(x)
result={
 'scope':'TOY_D6_SYMMETRIC_36_VERTEX_RATIONAL_FIXTURE_NOT_REAL_MT_D_N',
 'n_vertices':len(V),'n_native_directed_edges':len(C0),
 'n_allowed_directed_edges':len(H),'n_orbits':len(orbs),
 'native_K_star':q(max_ratio),'vertex_C_star_N12':q(Cstar),
 'cut_witness':{'orbit':'0A','demand':q(F(13)),'out_capacity':q(F(9))},
 'vertex_dual':{'phi':'1-level','beta':'indicator(level==1)','numerator':q(numer),'sum_beta_W':q(weighted)},
 'checks':flags,'all_pass':all(flags.values())}
(P/'EXACT_CHECKS.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
print(json.dumps(result,indent=2))

# Independently enumerate EVERY full cut in a smaller Z2-symmetric graph.
smallV=[(lev,j) for lev in range(3) for j in range(2)]
smallC0={}
smallH={}
for j in range(2):
    for k in range(2):
        e=((0,j),(1,k));smallC0[e]=F(1)
        if j==k:smallH[e]=F(1)
        e2=((1,j),(2,k));smallC0[e2]=F(1);smallH[e2]=F(1)
smallSource=defaultdict(F)
for (x,y),v in smallC0.items():smallSource[x]+=v;smallSource[y]-=v
assert all(smallSource[(l,j)]==[F(2),F(0),-F(2)][l] for l,j in smallV)

def small_ratio(Sset):
    d=sum((smallSource[x] for x in Sset),F(0))
    cc=sum((v for (x,y),v in smallH.items() if x in Sset and y not in Sset),F(0))
    if d>0 and cc==0:return None
    return d/cc if d>0 else F(0)
fullMax=F(0)
for mask in range(1<<len(smallV)):
    sel={smallV[i] for i in range(len(smallV)) if mask>>i&1}
    ratio=small_ratio(sel)
    assert ratio is not None
    fullMax=max(fullMax,ratio)
orbMax=F(0)
for mask in range(1<<3):
    sel={(l,j) for l in range(3) for j in range(2) if mask>>l &1}
    orbMax=max(orbMax,small_ratio(sel))
assert fullMax==orbMax==F(2)
result['extra_full_cut_exhaustive_control']={'full_cuts':64,'orbit_cuts':8,'K_star':'2','pass':True}
(P/'EXACT_CHECKS.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
print('PASS independent six-vertex exhaustive full-cut/orbit-cut comparison')
