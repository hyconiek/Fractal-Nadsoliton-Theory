from __future__ import annotations
import math,sys,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
H=ROOT/'inputs/FR223_20260916'
sys.path.insert(0,str(H/'src'))
import off_face

def weights_formula(r,s,t,y, mutate=None):
    a=math.sqrt(3.0)
    sr=math.sqrt(r)
    if mutate=='missing_sqrt_r': sr=1.0
    ap=a if mutate!='sqrt3_sign' else -a
    w=np.array([
      1.0,
      sr*s*t**(2+ap)*y,
      r*s*t,
      sr*t**2*y,
      s*t**3,
      sr*s*t**(2-ap)*y,
      r*t**4,
      sr*s*t**(2-ap)*y,
      s*t**3,
      sr*t**2*y,
      r*s*t,
      sr*s*t**(2+ap)*y,
    ],float)
    if mutate=='y_even': w[2]*=y
    return w

def p_formula(r,s,t,y,mutate=None):
    w=weights_formula(r,s,t,y,mutate); return w/w.sum()

def fields(r,s,t,y):
    return [-.5*math.log(r),-(2/3)*math.log(s),-2*math.log(t),-.5*math.log(y)]

def run(seed=17009,n=128):
    rng=np.random.default_rng(seed)
    maxerr=0.; worst=None
    mutation_min={'sqrt3_sign':float('inf'),'missing_sqrt_r':float('inf'),'y_even':float('inf')}
    for i in range(n):
        # broad interior points including near residual-hull edges
        lr=rng.uniform(math.log(1/900),0); ls=rng.uniform(math.log(1/128),0)
        lt=rng.uniform(math.log(1/9),0); ly=rng.uniform(math.log(1e-6),0)
        r,s,t,y=map(math.exp,(lr,ls,lt,ly))
        p=p_formula(r,s,t,y); q=off_face.p_from_fields(fields(r,s,t,y))
        e=float(np.max(np.abs(p-q)))
        if e>maxerr:maxerr=e;worst=[r,s,t,y]
        for m in mutation_min:
            mutation_min[m]=min(mutation_min[m],float(np.max(np.abs(p_formula(r,s,t,y,m)-q))))
    # limiting cases and anchor positivity
    limits={}
    for name,z in {'uniform':(1,1,1,1),'deep_all':(1/900,1/128,1/9,1e-6),'odd_suppressed':(.2,.3,.4,1e-12)}.items():
        p=p_formula(*z); limits[name]={'sum':float(p.sum()),'min':float(p.min()),'p0':float(p[0])}
    out={'task':'R7N-009','status':'NUMERICAL_DIAGNOSTIC_PLUS_DIRECT_FORMULA_IMPLEMENTATION','seed':seed,'samples':n,
         'max_abs_probability_discrepancy_vs_original_field_definition':maxerr,'worst_compact_point':worst,
         'mutation_min_detected_discrepancy':mutation_min,'limits':limits,
         'formula_source':'FIN_Rank7_Next_Campaign_Plan_EN_20260919.md Section 4.1',
         'note':'Agreement is a diagnostic of the independently transcribed compact formula; it is not by itself an exact symbolic proof.'}
    return out
if __name__=='__main__':
    out=run(); print(json.dumps(out,indent=2));
    (ROOT/'results').mkdir(exist_ok=True);json.dump(out,open(ROOT/'results/R7N-009_compact_formula_diagnostics.json','w'),indent=2)
