from pathlib import Path
import json
import numpy as np
import sympy as s
H=Path(__file__).resolve().parent
x=s.symbols('x',real=True);f=1/s.cosh(x)**2
assert s.simplify(s.expand_trig(-s.diff(f,x,2)+(4-6/s.cosh(x)**2)*f))==0
# Factorization B*B, B=d/dx+2 tanh(x), proves nonnegativity on H1.
a=2*s.tanh(x)
assert s.simplify(a*a-s.diff(a,x)-(4-6/s.cosh(x)**2))==0
k=np.array([np.cos(.18575*d+.1625)/(1+d**1.8) for d in range(1,12)])
lam=np.array([sum(k*np.cos(2*np.pi*j*np.arange(1,12)/12)) for j in range(12)])
K=np.array([[0 if i==j else k[(j-i)%12-1] for j in range(12)] for i in range(12)])
assert np.max(abs(np.sort(np.linalg.eigvalsh((K+K.T)/2))-np.sort(lam)))<1e-12
threshold=lam[0]-lam
rows=[{'mode':j,'multiplicity':1 if j==6 else 2,'threshold':float(threshold[j]),'lowest_eigenvalue_at_w_1_25':float(1.25-threshold[j]),'lowest_eigenvalue_at_w_1_75':float(1.75-threshold[j])} for j in range(1,7)]
assert sum(r['multiplicity'] for r in rows)==11
assert min(r['lowest_eigenvalue_at_w_1_25'] for r in rows)>.098
# Conservative analytic bound: |lambda0-lambdaj| <= 2 sum |k_d|.
# Actual classification uses the explicitly stated floating kernel evaluation.
assert abs(threshold[6]-max(threshold))<1e-12
out={'status':'CONDITIONAL_LINEAR_TRANSVERSE_STABILITY','rows':rows,'minimum_gap':float(1.25-max(threshold)),'critical_w':float(max(threshold)),'numeric_precision_note':'Double precision evaluation, not an outward-rounded interval certificate. Exact criterion and factorization are separate.'}
(H/'results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
