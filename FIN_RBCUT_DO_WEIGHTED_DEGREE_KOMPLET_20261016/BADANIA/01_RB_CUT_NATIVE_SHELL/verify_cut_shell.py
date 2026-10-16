#!/usr/bin/env python3
"""Independent bounded rational q12 Fourier and cut checks; NOT an all-N MT ideal sign proof."""
from fractions import Fraction as F
from random import Random
from itertools import product
import json
from pathlib import Path

out = Path(__file__).resolve().parent
checks = {}
def ck(key, yes):
    if not yes:
        raise AssertionError(key)
    checks[key] = True

C = [1,0,-1,0]*3
S = [0,1,0,-1]*3
P = [(-1)**j for j in range(12)]
# Identity for every possible (ordered) opposite-parity label replacement.
for a,b in product(range(12),repeat=2):
    if P[a]==P[b]:
        continue
    dc = C[b]-C[a]
    ds = S[b]-S[a]
    ck(f'opposite_pair_{a}_{b}', abs(dc)==abs(ds)==1 and dc*dc+ds*ds==2 and P[b]-P[a] in (-2,2))
    for cc,ss in [(0,0),(11,-7),(-31,17)]:
        # From R=2(C^2+S^2)/N, exact with N symbolic:
        exact_n_delta_R = 2*((cc+dc)**2+(ss+ds)**2-cc*cc-ss*ss)
        formula = 4*(cc*dc+ss*ds)+4
        ck(f'exact_delta_{a}_{b}_{cc}_{ss}', exact_n_delta_R==formula)

rng=Random(20261016)
for N in [12,24,120,1200]:
    counts = [N//12]*12
    for i in range(96):
        a=rng.choice([j for j in range(12) if counts[j]>0]);b=rng.randrange(12)
        if a==b: continue
        old_C=sum(x*y for x,y in zip(C,counts));old_S=sum(x*y for x,y in zip(S,counts))
        old_m=sum(x*y for x,y in zip(P,counts))
        old_R=F(2*(old_C**2+old_S**2),N)
        counts[a]-=1;counts[b]+=1
        new_C=sum(x*y for x,y in zip(C,counts));new_S=sum(x*y for x,y in zip(S,counts))
        new_m=sum(x*y for x,y in zip(P,counts))
        new_R=F(2*(new_C**2+new_S**2),N)
        ck(f'counts_valid_N{N}_{i}', sum(counts)==N and all(v>=0 for v in counts))
        ck(f'count_step_delta_N{N}_{i}',N*(new_R-old_R)==4*(old_C*(C[b]-C[a])+old_S*(S[b]-S[a]))+2*((C[b]-C[a])**2+(S[b]-S[a])**2))
        ck(f'parity_step_N{N}_{i}', new_m-old_m==P[b]-P[a])

# Small exact transshipment toy with a forbidden cross-phase incoming edge.
# It illustrates why a shell-supported *small* negative cut is still negative.
e = {('a0','a1'):F(1,10),('a1','a2'):F(1,10),('b0','a1'):F(1)}
phase = {'a0':'A','a1':'A','a2':'A','b0':'B'}
H = {xy:w for xy,w in e.items() if phase[xy[0]]==phase[xy[1]]}
L={'a0','a1'}
src={x:F(0) for x in phase}
for (a,b),w in e.items():src[a]+=w;src[b]-=w
ck('toy_source_total_zero',sum(src.values())==0)
ck('toy_legal_lower_ideal',all(a in L for (a,b) in H if b in L))
ck('toy_bad_ideal',sum(src[x] for x in L)==-F(9,10))
ck('toy_forbidden_down_equals_one',sum(w for (a,b),w in e.items() if (a,b) not in H and b in L and a not in L)==1)
ck('toy_allowed_up_one_tenth',sum(w for (a,b),w in H.items() if a in L and b not in L)==F(1,10))
ck('toy_cut_identity',sum(src[x] for x in L)==sum(w for (a,b),w in e.items() if a in L and b not in L)-sum(w for (a,b),w in e.items() if b in L and a not in L))
for i,epsilon in enumerate([F(1),F(1,10),F(1,1000),F(1,1000000)]):
    ck(f'toy_rescale_stays_negative_{i}',epsilon*sum(src[x] for x in L)<0)

result = {'scientific_scope':'exact q12 Fourier one-copy arithmetic and abstract finite-flow negative controls; no all-N MT graph solved',
          'opposite_parity_ordered_pairs':72,'diagnostic_count_N':[12,24,120,1200],
          'toy_negative_cut':'-9/10','checks':len(checks),'pass':all(checks.values()),'failures':0}
(out/'verification.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(f"PASS {len(checks)} bounded exact integer/rational assertions. NO native all-N ideal certification.")
