from __future__ import annotations
from pathlib import Path
from fractions import Fraction as F
import sys,json,math,time
ROOT=Path('/mnt/data/fin_rank7_next_campaign');IR=ROOT/'inputs/intake_review_20260919';sys.path.insert(0,str(IR))
import scientific_rechecks as sr
AMP=list(map(sr.I,['0.1131879146','0.1698528641','0.2269339093'])); Z6=sr.I('-0.3380663037');RT3=sr.iv.sqrt(3);RT12=sr.iv.sqrt(12);TWOPI=2*sr.iv.pi

def grad_box(box):
    z=[sr.iv.mpf([sr.I(F(a)),sr.I(F(b))]) for a,b in box]
    phi=[TWOPI*x for x in z]
    hs=[]; ds=[]
    for j in range(12):
        h=Z6*((-1)**j)/RT12; row=[]
        for a,k,p in zip(AMP,(3,4,5),phi):
            ang=TWOPI*sr.I(F(k*j,12))+p; aa=a/RT3
            h+=aa*sr.iv.cos(ang); row.append(-aa*sr.iv.sin(ang))
        hs.append(h); ds.append(row)
    # fixed shift chosen from interval midpoint upper proxy; any scalar is exact algebraically
    mids=[sr.mid(h) for h in hs]; m=sr.I(repr(max(mids)))
    w=[sr.iv.exp(h-m) for h in hs]; s=sum(w,sr.I(0)); p=[x/s for x in w]
    return [sum((p[j]*ds[j][a] for j in range(12)),sr.I(0)) for a in range(3)]

def classify(box):
    g=grad_box(box)
    for k,v in enumerate(g):
        lo,hi=sr.bounds(v)
        if lo>0 or hi<0:return True,k,[str(lo),str(hi)]
    return False,None,None

if __name__=='__main__':
    import random
    d=json.load(open(ROOT/'checkpoints/R7N-037_quartic_cover_normalized.json'))
    cells=d['safe_leaves'][:100]+d['root_leaves'][:100]
    st=time.time(); n=0
    for rec in cells:
        ok,*_=classify(rec['box']); n+=ok
    print(json.dumps({'cells':len(cells),'safe':n,'elapsed':time.time()-st},indent=2))
