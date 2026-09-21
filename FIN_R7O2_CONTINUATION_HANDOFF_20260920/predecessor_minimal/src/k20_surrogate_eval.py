from __future__ import annotations
from pathlib import Path
from fractions import Fraction as F
import sys,json,math
ROOT=Path(__file__).resolve().parents[1]; IR=ROOT/'inputs/intake_review_20260919';sys.path.insert(0,str(IR));import scientific_rechecks as sr
DATA=json.load(open(ROOT/'results/R7N-043_K20_surrogate.json'))
EPS=F(DATA['uniform_full_gradient_error_bound']); TWOPI=2*sr.iv.pi
TERMS=[]
for r in DATA['terms']:
    lo,hi=map(F,r['coefficient_interval']); TERMS.append((sr.iv.mpf([sr.I(lo),sr.I(hi)]),tuple(r['frequency'])))

def grad_box(box):
    z=[sr.iv.mpf([sr.I(F(a)),sr.I(F(b))]) for a,b in box]
    g=[sr.I(0),sr.I(0),sr.I(0)]
    for coef,e in TERMS:
        ang=TWOPI*sum((sr.I(e[k])*z[k] for k in range(3)),sr.I(0)); sn=sr.iv.sin(ang)
        for k in range(3):
            if e[k]: g[k]+=-sr.I(2*e[k])*coef*sn
    return g

def classify(box):
    g=grad_box(box)
    for k,v in enumerate(g):
        lo,hi=sr.bounds(v)
        if lo>EPS or hi<-EPS: return True,k,[str(lo),str(hi)]
    return False,None,None

def numeric(ph):
    out=[0.,0.,0.]
    for coef,e in TERMS:
        lo,hi=sr.bounds(coef); c=float((lo+hi)/2); ang=sum(e[k]*ph[k] for k in range(3)); sn=math.sin(ang)
        for k in range(3): out[k]+=-2*c*e[k]*sn
    return out
