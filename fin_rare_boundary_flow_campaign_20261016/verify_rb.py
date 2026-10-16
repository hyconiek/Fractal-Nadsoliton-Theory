#!/usr/bin/env python3
"""Independent, bounded exact checks for the NEW RB handoff. NOT an all-N q12 proof."""
from fractions import Fraction as F
from itertools import product
import json, pathlib
P=pathlib.Path(__file__).resolve().parent
checks={}
def check(label, predicate):
    if not predicate: raise AssertionError(label)
    checks[label]=True

# Fourier normalization and D12 action are exact integer-valued for k=3.
c=[1,0,-1,0]*3
s=[0,1,0,-1]*3
sigma=[(-1)**j for j in range(12)]
for a,b in product(range(12), repeat=2):
    if a==b: continue
    dc,ds=c[b]-c[a],s[b]-s[a]
    check(f"fourier_jump_{a}_{b}",dc in (-2,-1,0,1,2) and ds in (-2,-1,0,1,2))
    check(f"parity_jump_{a}_{b}",sigma[b]-sigma[a] in (-2,0,2))
for a in range(12):
    check(f"rotation_{a}",c[(a+1)%12]==-s[a] and s[(a+1)%12]==c[a] and sigma[(a+1)%12]==-sigma[a])

# Exact rational toy: 12 vertices arranged in A: a/c branches and B: b branch.
# Departure layer E includes a3,b3. Only cross-phase edge not touching E is b0->a1,
# which is illegal in H, but present in a reversible, positive native reference.
eps=F(1,10)
phase={f'{q}{i}':('B' if q=='b' else 'A') for q in 'abc' for i in range(4)}
level={f'{q}{i}':i for q in 'abc' for i in range(4)}
E={'a3','b3'}
edge={}
def put(x,y,v):
    assert level[y]==level[x]+1
    assert (x,y) not in edge
    edge[x,y]=v
# Intra-phase chains all have strictly positive rates; c branch provides A parity-cut traffic.
for q in 'abc':
    for i in range(3):
        value=(F(1)+eps if q=='a' and i in (0,2) else F(2)+eps if q=='c' else F(1)+eps if q=='b' and i in (0,1) else eps)
        put(f'{q}{i}',f'{q}{i+1}',value)
# One added intra A link and two allowed E incident cross-phase links.
put('c2','a3',eps)
put('a2','b3',F(1)+eps)
put('b2','a3',eps)
# Cross-phase interior forbidden by H.
put('b0','a1',F(1))
H={(x,y):v for (x,y),v in edge.items() if phase[x]==phase[y] or x in E or y in E}
W={x:F(1) for x in phase}
source={x:F(0) for x in phase}
for (x,y),value in edge.items():
    source[x]+=value;source[y]-=value
L={'a0','a1'}
def downset(L):return all(x in L for x,y in H if y in L)
check('toy_total_zero',sum(source.values())==0)
check('toy_all_native_positive',all(v>0 for v in edge.values()))
check('toy_cross_edge_rejected',('b0','a1') not in H)
check('toy_boundary_cross_edges_allowed',('a2','b3') in H and ('b2','a3') in H)
check('toy_lower_ideal',downset(L))
check('toy_negative_cut',sum(source[x] for x in L)==-F(9,10))
for q in ('A','B'):
    check(f'toy_phase_balance_{q}',sum(source[x] for x in phase if phase[x]==q)==0)
    for i in range(3):
        prefix=sum(source[x] for x in phase if phase[x]==q and level[x]<=i)
        check(f'toy_phase_level_cut_{q}_{i}',prefix>=0)
for i in range(3):
    check(f'toy_global_level_cut_{i}',sum(source[x] for x in phase if level[x]<=i)>=0)
for x,b in source.items():
    if b>0:check(f'toy_positive_drift_move_{x}',any(u==x for u,v in H))
    if b<0:check(f'toy_negative_drift_move_{x}',any(v==x for u,v in H))
# Exact native cut defect formula for a H-lower ideal; native forbidden edge carries negative cut.
up=sum(v for (x,y),v in edge.items() if x in L and y not in L)
down=sum(v for (x,y),v in edge.items() if x not in L and y in L)
allowed_up=sum(v for (x,y),v in H.items() if x in L and y not in L)
forbidden_up=sum(v for (x,y),v in edge.items() if (x,y) not in H and x in L and y not in L)
forbidden_down=sum(v for (x,y),v in edge.items() if (x,y) not in H and x not in L and y in L)
check('toy_defect_ledger',sum(source[x] for x in L)==up-down==allowed_up+forbidden_up-forbidden_down)
check('toy_forbidden_down',forbidden_down==F(1) and forbidden_up==0 and allowed_up==eps)
# Exact dual phi=1_L (where incidence is + lower and - upper), beta=0.
check('toy_dual_edge',all((int(x in L)-int(y in L))>=0 for x,y in H))
check('toy_dual_negative',sum(F(int(x in L))*v for x,v in source.items())<0)
# Coarse phase/level cuts cannot see this H-ideal violation.
summary={
 'scope':'exact synthetic 12-vertex directed-level network; NOT native q12/A7 nor MT D_N',
 'vertex_count':len(phase),'legal_reference_edges':len(edge),'allowed_edges':len(H),
 'forbidden_edges':[list(t) for t in edge if t not in H],
 'departure_layer':sorted(E), 'negative_lower_ideal':sorted(L),
 'negative_cut':str(sum(source[x] for x in L)),
 'allowed_exit_flow':str(allowed_up),'forbidden_inward_flow':str(forbidden_down),
 'per_phase_cuts':{q:[str(sum(source[x] for x in phase if phase[x]==q and level[x]<=i)) for i in range(3)] for q in ('A','B')},
 'checks_passed':len(checks),'failed':0
}
(P/'EXACT_DIAGNOSTICS.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
print('PASS',len(checks),'exact finite checks;',len(phase),'toy vertices;',len(edge),'native edges;',len(H),'allowed')
