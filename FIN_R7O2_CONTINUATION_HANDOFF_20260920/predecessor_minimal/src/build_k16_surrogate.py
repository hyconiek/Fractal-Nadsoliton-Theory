from __future__ import annotations
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
from math import comb, factorial
import sys,json,time,random,math
ROOT=Path('/mnt/data/fin_rank7_next_campaign'); IR=ROOT/'inputs/intake_review_20260919';sys.path.insert(0,str(IR))
import scientific_rechecks as sr
sr.iv.dps=55

def add_iv(a,b): return a+b

def bounds(x): return sr.bounds(x)
def maxabs(x):
    lo,hi=bounds(x); return max(abs(lo),abs(hi))

def serialize_exactish(x): return sr.serialize(x)

r3,r4,r5=map(sr.I,['0.1131879146','0.1698528641','0.2269339093']); z6=sr.I('-0.3380663037')
a=r3/(2*sr.iv.sqrt(3)); b=r4/(2*sr.iv.sqrt(3)); c=r5/(2*sr.iv.sqrt(3)); d=z6/sr.iv.sqrt(12)
BASE=[(3,(1,0,0),a),(9,(-1,0,0),a),(4,(0,1,0),b),(8,(0,-1,0),b),(5,(0,0,1),c),(7,(0,0,-1),c),(6,(0,0,0),d)]
ZERO=sr.I(0)

def poly_mul(A,B):
    out={}
    for e,x in A.items():
        for f,y in B.items():
            g=(e[0]+f[0],e[1]+f[1],e[2]+f[2])
            out[g]=out.get(g,ZERO)+x*y
    return out

def build(N=16):
    st=time.time(); state={(0,(0,0,0)):sr.I(1)}; moments={0:{(0,0,0):sr.I(1)}}
    for n in range(1,N+1):
        new={}
        for (lab,e),v in state.items():
            for l,f,coef in BASE:
                g=((lab+l)%12,(e[0]+f[0],e[1]+f[1],e[2]+f[2]))
                new[g]=new.get(g,ZERO)+v*coef
        state=new; moments[n]={e:v for (lab,e),v in state.items() if lab==0}
        print('moment',n,len(moments[n]),flush=True)
    kappas={1:moments[1]}
    for n in range(2,N+1):
        out=dict(moments[n])
        for j in range(1,n):
            prod=poly_mul(kappas[j],moments[n-j]); scale=sr.I(comb(n-1,j-1))
            for e,v in prod.items(): out[e]=out.get(e,ZERO)-scale*v
        kappas[n]=out
        print('kappa',n,len(out),flush=True)
    agg={}
    for n in range(2,N+1):
        scale=sr.I(F(1,factorial(n)))
        for e,v in kappas[n].items(): agg[e]=agg.get(e,ZERO)+scale*v
    # Pair +/- frequencies. Exact real symmetry means coefficient belongs to both intervals; intersect.
    terms=[]; seen=set()
    for e,v in agg.items():
        if e==(0,0,0) or e in seen: continue
        ne=(-e[0],-e[1],-e[2]); w=agg.get(ne)
        if w is None: raise RuntimeError(('missing conjugate',e))
        vlo,vhi=bounds(v); wlo,whi=bounds(w); lo=max(vlo,wlo); hi=min(vhi,whi)
        if lo>hi: raise RuntimeError(('symmetry intervals disjoint',e,vlo,vhi,wlo,whi))
        ci=sr.iv.mpf([sr.I(lo),sr.I(hi)])
        # canonical sign for stable output
        ce=e if next(x for x in e if x!=0)>0 else ne
        if ce!=e: ci=ci # same real coefficient
        contrib=[2*maxabs(ci)*abs(ce[k]) for k in range(3)]
        terms.append((ce,ci,contrib));seen.add(e);seen.add(ne)
    # Keep terms until rigorous omitted componentwise L1 <= 3.5e-9.
    terms.sort(key=lambda x:max(x[2]),reverse=True)
    keep=[]; drop=[F(0),F(0),F(0)]
    # Start by keeping terms with max contribution >=1e-9, then compute exact drop.
    for e,ci,cc in terms:
        if max(cc)>=F(1,10**9): keep.append((e,ci,cc))
        else:
            for k in range(3): drop[k]+=cc[k]
    # if drop too big, promote largest omitted until target
    omitted=[x for x in terms if x not in keep]
    while max(drop)>F(35,10**10):
        e,ci,cc=omitted.pop(0);keep.append((e,ci,cc))
        for k in range(3): drop[k]-=cc[k]
    keep.sort(key=lambda x:x[0])
    # Cauchy tail at R=3.7
    R=sr.I(F(37,10));H=(r3+r4+r5)/sr.iv.sqrt(3)+abs(z6)/sr.iv.sqrt(12);Amax=r5/sr.iv.sqrt(3);RH=R*H
    assert bounds(RH-sr.iv.pi/2)[1]<0
    M=R*Amax*sr.iv.exp(2*RH)/sr.iv.cos(RH)
    tail=M*(R**(-(N+1)))/(1-1/R); tail_hi=bounds(tail)[1]
    total_eps=max(drop)+tail_hi
    rows=[]
    for e,ci,cc in keep:
        rows.append({'frequency':list(e),'coefficient_interval':serialize_exactish(ci),'gradient_component_abs_bounds':[str(x) for x in cc]})
    out={'task':'R7N-042/043-surrogate','N':N,'fixture':{'r3':'0.1131879146','r4':'0.1698528641','r5':'0.2269339093','z6':'-0.3380663037'},
         'derivation':'K_N=sum_{n=1}^N kappa_n/n! from exact Fourier moment convolution on Z12; retained real cosine resonances are paired +/- phase frequencies.',
         'retained_terms':len(rows),'all_nonzero_frequency_pairs':len(terms),'terms':rows,
         'dropped_gradient_component_bounds':[str(x) for x in drop],
         'complex_analytic_tail':{'R':'37/10','H_interval':sr.serialize(H),'RH_interval':sr.serialize(RH),'M_interval':sr.serialize(M),'tail_interval':sr.serialize(tail),
          'proof':'For |alpha|<=R and RH<pi/2, Re <exp(alpha h)> >= exp(-RH) cos(RH)>0. Thus G_a=partial_phi_a log<exp(alpha h)> is analytic and |G_a|<=R A_a exp(2RH)/cos(RH); Cauchy bounds the coefficient tail after N.'},
         'uniform_full_gradient_error_bound':str(total_eps),'elapsed_seconds':time.time()-st}
    (ROOT/'results/R7N-043_K16_surrogate.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='terms'},indent=2))
    return out
if __name__=='__main__':build()
