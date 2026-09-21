from __future__ import annotations
from pathlib import Path
from fractions import Fraction as F
from math import comb, factorial
import sys,json,time
ROOT=Path('/mnt/data/fin_rank7_next_campaign'); IR=ROOT/'inputs/intake_review_20260919';sys.path.insert(0,str(IR))
import scientific_rechecks as sr
sr.iv.dps=60

N=20
DROP_TARGET=F(1,10**10)
R=F(19,5)  # 3.8; near-optimal simple rational for N=20
r3,r4,r5=map(sr.I,['0.1131879146','0.1698528641','0.2269339093']); z6=sr.I('-0.3380663037')
a=r3/(2*sr.iv.sqrt(3)); b=r4/(2*sr.iv.sqrt(3)); c=r5/(2*sr.iv.sqrt(3)); d=z6/sr.iv.sqrt(12)
BASE=[(3,(1,0,0),a),(9,(-1,0,0),a),(4,(0,1,0),b),(8,(0,-1,0),b),(5,(0,0,1),c),(7,(0,0,-1),c),(6,(0,0,0),d)]
ZERO=sr.I(0)

def bounds(x): return sr.bounds(x)
def maxabs(x):
    lo,hi=bounds(x); return max(abs(lo),abs(hi))

def poly_mul(A,B):
    out={}
    for e,x in A.items():
        for f,y in B.items():
            g=(e[0]+f[0],e[1]+f[1],e[2]+f[2])
            out[g]=out.get(g,ZERO)+x*y
    return out

def build():
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
    terms=[]; seen=set()
    for e,v in agg.items():
        if e==(0,0,0) or e in seen: continue
        ne=(-e[0],-e[1],-e[2]); w=agg.get(ne)
        if w is None: raise RuntimeError(('missing conjugate',e))
        vlo,vhi=bounds(v); wlo,whi=bounds(w); lo=max(vlo,wlo); hi=min(vhi,whi)
        if lo>hi: raise RuntimeError(('symmetry intervals disjoint',e))
        ci=sr.iv.mpf([sr.I(lo),sr.I(hi)])
        ce=e if next(x for x in e if x!=0)>0 else ne
        cc=[2*maxabs(ci)*abs(ce[k]) for k in range(3)]
        terms.append((ce,ci,cc)); seen.add(e); seen.add(ne)
    terms.sort(key=lambda x:max(x[2]),reverse=True)
    # Greedy retain: start empty, drop all, promote largest until each componentwise omitted L1 <= target.
    drop=[sum((x[2][k] for x in terms),F(0)) for k in range(3)]
    keep=[]; idx=0
    while max(drop)>DROP_TARGET:
        e,ci,cc=terms[idx]; idx+=1; keep.append((e,ci,cc))
        for k in range(3): drop[k]-=cc[k]
    keep.sort(key=lambda x:x[0])
    RI=sr.I(R); H=(r3+r4+r5)/sr.iv.sqrt(3)+abs(z6)/sr.iv.sqrt(12); Amax=r5/sr.iv.sqrt(3); RH=RI*H
    assert bounds(RH-sr.iv.pi/2)[1]<0
    M=RI*Amax*sr.iv.exp(2*RH)/sr.iv.cos(RH)
    tail=M*(RI**(-(N+1)))/(1-1/RI); tail_hi=bounds(tail)[1]
    total_eps=max(drop)+tail_hi
    rows=[]
    for e,ci,cc in keep:
        rows.append({'frequency':list(e),'coefficient_interval':sr.serialize(ci),'gradient_component_abs_bounds':[str(x) for x in cc]})
    out={
      'task':'R7N-043-K20-surrogate','N':N,
      'fixture':{'r3':'0.1131879146','r4':'0.1698528641','r5':'0.2269339093','z6':'-0.3380663037'},
      'derivation':'K_N=sum_{n=1}^N kappa_n/n! from exact Fourier moment convolution on Z12; retained real cosine resonances paired +/- phase frequencies.',
      'retained_terms':len(rows),'all_nonzero_frequency_pairs':len(terms),'terms':rows,
      'drop_target':str(DROP_TARGET),'dropped_gradient_component_bounds':[str(x) for x in drop],
      'complex_analytic_tail':{'R':str(R),'H_interval':sr.serialize(H),'RH_interval':sr.serialize(RH),'M_interval':sr.serialize(M),'tail_interval':sr.serialize(tail),
        'proof':'For |alpha|<=R and RH<pi/2, Re <exp(alpha h)> >= exp(-RH) cos(RH)>0. Thus each phase-gradient component of log<exp(alpha h)> is analytic; Cauchy bounds the coefficient tail after N.'},
      'uniform_full_gradient_error_bound':str(total_eps),'elapsed_seconds':time.time()-st
    }
    p=ROOT/'results/R7N-043_K20_surrogate.json'; p.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='terms'},indent=2))
if __name__=='__main__': build()
