#!/usr/bin/env python3
"""Small, reproducible diagnostics for AH-001/002; not a proof or empirical measurement."""
from fractions import Fraction
from math import exp, log, ceil
from random import Random
import json
from pathlib import Path

N=12
weights=[1,2,3,4,5,6]
s=2*sum(weights[:5])+weights[5]
J=[[Fraction(0) for _ in range(N)] for __ in range(N)]
for i in range(N):
    for j in range(N):
        if i!=j:
            d=min((i-j)%N,(j-i)%N)
            J[i][j]=Fraction(weights[d-1],s)
assert s==36
assert all(sum(J[i])==1 for i in range(N))
a=min(J[i][j] for i in range(N) for j in range(N) if i!=j)
b=10*a*a
J2=[[sum(J[i][k]*J[k][j] for k in range(N)) for j in range(N)] for i in range(N)]
assert min(x for row in J2 for x in row)>=b
assert min(J2[i][j] for i in range(N) for j in range(N) if i!=j)>b # Lower bound need not be attained

# Enumeration without 12! permutations: one nontrivial permutation of names.
sigma=[5,1,11,6,0,7,3,10,4,8,2,9]
K=[[J[sigma[i]][sigma[j]] for j in range(N)] for i in range(N)]
assert all(sum(K[i])==1 for i in range(N))

# For a path, standardize to first-appearance names and collect rows; reconstruct
# the original matrix up to the random first-visit permutation.
rng=Random(20261020)
steps=3200
current=rng.randrange(N)
path=[current]
for _ in range(steps):
    u=rng.random()
    q=0.
    next_state=None
    for j in range(N):
        q+=float(J[current][j])
        if u<q:
            next_state=j
            break
    if next_state is None: next_state=N-1
    assert next_state!=current
    path.append(next_state)
    current=next_state
name={}
code=[]
for i in path:
    if i not in name: name[i]=len(name)
    code.append(name[i])
assert len(name)==N
assert set(name.values())==set(range(N))
rows=[0]*N
edges=[[0]*N for _ in range(N)]
for x,y in zip(code[:-1],code[1:]):
    rows[x]+=1
    edges[x][y]+=1
inv={v:k for k,v in name.items()}
max_error=max(abs(edges[i][j]/rows[i]-float(J[inv[i]][inv[j]]))
              for i in range(N) for j in range(N))
assert max_error < 0.07, (max_error, rows, name, edges[0], [float(J[inv[0]][inv[j]]) for j in range(N)])

# Explicit analytic event moment controls: ABAB and cyclic three-step loop.
alt=sum(2*Fraction(weights[d-1],s)**3 for d in range(1,6))+Fraction(weights[5],s)**3
J2diag=sum(J[0][k]*J[k][0] for k in range(N))
assert alt>0 and J2diag>0

# Quantitative AH bound for demonstrative admissible parameters. It is a 
# deliberately large sufficient number, not a practical sample target.
a_example=Fraction(1,20)
b_example=10*a_example*a_example
delta=Fraction(1,10)
eta=Fraction(1,20)
m_needed=max(ceil(2/float(b_example)),ceil(8/float(b_example*b_example)*log(24/float(eta))),
             ceil(384/(float(delta*delta*b_example*b_example))*log(528/float(eta))))
L_needed=2*m_needed+2
assert 12*exp(-float(b_example*b_example)*m_needed/8)+264*exp(-float(delta*delta*b_example*b_example)*m_needed/384)<=float(eta)

out={"method":"exact rational source checks and short seeded diagnostic only", "states":N,
     "example_weights":weights,"example_exit_rate":s,"example_min_offdiag_a":str(a),
     "two_step_offdiag_min":str(min(J2[i][j] for i in range(N) for j in range(N) if i!=j)),
     "two_step_minorization_b":str(b),"sample_seed":20261020,"sample_transition_count":steps,
     "first_visit_code_covered":len(name),"sample_max_entry_reconstruction_error":max_error,
     "example_ABAB_probability":str(alt),"example_two_step_return_probability":str(J2diag),
     "sample_size_formula_check":{"a":"1/20","delta":"1/10","eta":"1/20","m":m_needed,"L":L_needed},
     "notes":["Not an independent proof of the infinite path theorem or Hoeffding-Azuma.",
       "No physical source, clock, FIN HB/count-process or experimental inference."]}
Path(__file__).with_name("RESULTS.json").write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n")
print(json.dumps(out,indent=2,ensure_ascii=False))
