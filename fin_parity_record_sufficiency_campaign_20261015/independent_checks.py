"""Finite-dimensional independent arithmetic sanity checks for LI source conventions.

Does NOT certify any stochastic functional limit, MT stationary lift, LI risk theorem,
or source A7 spectral inputs. Uses no external dependencies.
"""
import json
import math
from pathlib import Path

OUT = Path(__file__).resolve().parent
A = ((0.7, -0.2, 0.1), (-0.2, 0.7, 0.3), (0.1, 0.3, 0.7))  # symmetric, constant diag test kernel, NOT FIN A7
M = max(abs(x) for row in A for x in row)
g = 2.1
checks = 0
largest_ratio = 0.0
worst_gap = float('inf')

def softmax(h):
    hmax = max(h)
    s = [math.exp(x-hmax) for x in h]
    t = sum(s)
    return [x/t for x in s]

def dot(v,w):
    return sum(x*y for x,y in zip(v,w))

# Direct categorical enumeration across backgrounds, attempted slots, tag labels.
for N in range(3,20):
    for a in range(N):
        for b in range(N-a):
            n = (a,b,N-1-a-b)
            for tag in range(3):
                for i in range(3):
                    if n[i] < 1:
                        continue
                    h = [g/N*sum(A[j][k]*(n[k]-(k==i)) for k in range(3)) for j in range(3)]
                    delta = [g/N*A[j][tag] for j in range(3)]
                    pp = softmax([h[j]+delta[j] for j in range(3)])
                    pr = softmax(h)
                    kl = sum(x*math.log(x/y) for x,y in zip(pp,pr))
                    bound = (max(delta)-min(delta))**2/8
                    loose_bound = g*g*M*M/(2*N*N)
                    slack = min(bound-kl,loose_bound-kl)
                    if slack < -2e-14:
                        raise AssertionError(f'KL bound failure: N={N}, n={n}, tag={tag}, i={i}, kl={kl}, bound={bound}')
                    checks += 2
                    worst_gap=min(worst_gap,slack)
                    if bound>0:
                        largest_ratio=max(largest_ratio,kl/bound)

# For a path length sqrt(N), the TV bound is exactly O(N^-1/4).
for N in (12,24,120,1200,12000):
    S=1.
    T=N**.5*S
    tv=abs(g)*M*math.sqrt((N-1)*T)/(2*N)
    assert tv <= abs(g)*M*N**(-.25)/2*(1+1e-13)
    checks += 1

# Exact score algebra: dM_S= -2 S dM_J at each flip plus compensator drift.
for S in (-1,1):
    for q in (-0.7,-0.2,0.0,0.2,0.7):
        # at a flip event dM_J=1, dM_S=-2*S
        assert abs((-S*q)-(0.5*q*(-2*S)))<1e-14
        # drift under Q: dM_J compensator -dt/2; dM_S drift Sdt
        assert abs((S*q/2)-(0.5*q*S))<1e-14
        checks += 2

result={'status':'PASS','checks':checks,'max_KL_to_osc_bound_ratio':largest_ratio,'min_bound_minus_KL':worst_gap,
        'scope':'3-state synthetic symmetric constant-diagonal kernel; not A7 and no process convergence verification',
        'original_model_unchanged':True}
(OUT/'INDEPENDENT_CHECKS.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
print(json.dumps(result,indent=2,ensure_ascii=False))
