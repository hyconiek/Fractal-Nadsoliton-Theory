#!/usr/bin/env python3
from __future__ import annotations
import os
import json, math, sys
from pathlib import Path
from fractions import Fraction
import numpy as np
import mpmath as mp

WORK=Path(os.environ.get('MP7_WORK_ROOT','/mnt/data/fin_rank7_mathphysics_next'))
R7P=Path(os.environ.get('R7P_ROOT','/mnt/data/r7p_source/extracted/fin_rank7_followup'))
sys.path.insert(0,str(R7P))
from src import coexistence_certificate as cc

mp.mp.dps=80; mp.iv.dps=60; iv=mp.iv
coex=json.loads((R7P/'certificates/R7P-026_equal_energy_event.json').read_text())
sadd_src=json.loads((R7P/'certificates/R7P-029_barrier_saddle.json').read_text())
trans=json.loads((R7P/'certificates/R7P-028_crossing_transversality.json').read_text())
raw={int(k):v for k,v in coex['spectral_intervals'].items()}
GLO=float(coex['root_box'][4][0]); GHI=float(coex['root_box'][4][1])
GMID=(mp.mpf(coex['root_box'][4][0])+mp.mpf(coex['root_box'][4][1]))/2

def rat(s): return float(Fraction(s))
Llo=np.array([rat(raw[k][0]) for k in (3,4,5,6)])
Lhi=np.array([rat(raw[k][1]) for k in (3,4,5,6)])
dlo=np.nextafter(np.array([Llo[0]/6,Llo[1]/6,Llo[2]/6,Llo[3]/12]),-np.inf)
dhi=np.nextafter(np.array([Lhi[0]/6,Lhi[1]/6,Lhi[2]/6,Lhi[3]/12]), np.inf)
a_lo=np.nextafter(np.sqrt(dlo),-np.inf); a_hi=np.nextafter(np.sqrt(dhi),np.inf)
B=np.nextafter(GHI*dhi,np.inf)

# exact unscaled character table
rt3=iv.sqrt(iv.mpf(3)); half=iv.mpf('0.5')
c3=[1,0,-1,0]*3; c4=[iv.mpf(1),-half,-half]*4
c5=[iv.mpf(1),-rt3/2,half,iv.mpf(0),-half,rt3/2,-iv.mpf(1),rt3/2,-half,iv.mpf(0),half,-rt3/2]
alt=[1 if j%2==0 else -1 for j in range(12)]
F=[[iv.mpf(c3[j]),c4[j],c5[j],iv.mpf(alt[j])] for j in range(12)]

def ivpoint(x): return iv.mpf([repr(float(np.nextafter(x,-np.inf))),repr(float(np.nextafter(x,np.inf)))])
def lo_float(x): return float(np.nextafter(float(x.a),-np.inf))
def hi_float(x): return float(np.nextafter(float(x.b), np.inf))

def mean_at_point(J):
    X=[ivpoint(x) for x in J]; h=[sum(F[j][i]*X[i] for i in range(4)) for j in range(12)]
    e=[iv.exp(x) for x in h]; Z=sum(e,iv.mpf(0)); p=[x/Z for x in e]
    return [sum(p[j]*F[j][i] for j in range(12)) for i in range(4)]

def contract(lo,hi,maxiter=30):
    lo=lo.copy();hi=hi.copy()
    for _ in range(maxiter):
        ml=mean_at_point(lo);mu=mean_at_point(hi);nl=lo.copy();nh=hi.copy()
        for i in range(4):
            mlb=max(0.0,lo_float(ml[i])); mub=max(0.0,hi_float(mu[i]))
            tl=np.nextafter(GLO*dlo[i]*mlb,-np.inf); tu=np.nextafter(GHI*dhi[i]*mub,np.inf)
            nl[i]=max(lo[i],tl);nh[i]=min(hi[i],tu)
        if np.any(nl>nh): return None
        delta=max(float(np.max(np.abs(nl-lo))),float(np.max(np.abs(nh-hi))))
        lo,hi=nl,nh
        if delta<2e-15:break
    return lo,hi

# parametric 4D Krawczyk branch tube over the whole event g-box
_,Cp=cc.point_features();_,_,Cnorm,_=cc.interval_features(); I=cc._I

def stat_point(s,g):
    FF,JJ,_,_=cc.point_eval(list(s)+[g],Cp)
    return FF[:4],mp.matrix([[JJ[i,j] for j in range(4)] for i in range(4)])

def param_branch(seed,name,radius='1e-8'):
    def f(*s): return tuple(stat_point(s,GMID)[0])
    center=list(mp.findroot(f,tuple(mp.mpf(str(x)) for x in seed),tol=mp.mpf('1e-60'),maxsteps=60))
    _,H=stat_point(center,GMID); R=H**-1; rad=mp.mpf(radius)
    X=[I([mp.nstr(x-rad,90),mp.nstr(x+rad,90)]) for x in center]; G=I(coex['root_box'][4])
    x0=[I(mp.nstr(x,90)) for x in center]
    F0,J0,_,_=cc.interval_eval(x0+[G],Cnorm); FX,JX,_,_=cc.interval_eval(X+[G],Cnorm)
    JB=[[JX[i][j] for j in range(4)] for i in range(4)]
    K0=[x0[i]-sum(I(mp.nstr(R[i,j],90))*F0[j] for j in range(4)) for i in range(4)]
    E=[[I(1 if i==j else 0)-sum(I(mp.nstr(R[i,k],90))*JB[k][j] for k in range(4)) for j in range(4)] for i in range(4)]
    dx=[I([mp.nstr(-rad,90),mp.nstr(rad,90)]) for _ in range(4)]
    K=[K0[i]+sum(E[i][j]*dx[j] for j in range(4)) for i in range(4)]
    inc=[float(X[i].a)<float(K[i].a) and float(K[i].b)<float(X[i].b) for i in range(4)]
    assert all(inc)
    return {'name':name,'center':[str(x) for x in center], 'sbox':[[lo_float(x),hi_float(x)] for x in X],
            'krawczyk':[[lo_float(x),hi_float(x)] for x in K], 'strict_inclusion':True,
            'phi_over_tube':[lo_float(FX[4]),hi_float(FX[4])]}

locseed=[float(x) for x in coex['root_center'][:4]]
sadseed=[(float(x[0])+float(x[1]))/2 for x in sadd_src['saddle_root_box']]
loc=param_branch(locseed,'localized'); sad=param_branch(sadseed,'saddle')
branch_boxes={x['name']:(np.array([z[0] for z in x['sbox']]),np.array([z[1] for z in x['sbox']])) for x in (loc,sad)}

def spreimage(lo,hi): return np.nextafter(lo/a_hi,-np.inf),np.nextafter(hi/a_lo,np.inf)
def contained(name,lo,hi):
    sl,sh=spreimage(lo,hi); L,U=branch_boxes[name]; return bool(np.all(sl>=L) and np.all(sh<=U))
def uniform_contained(lo,hi): return bool(np.all(lo>=-1e-18) and np.all(hi<=1e-8))

def uniform_banach():
    J=[iv.mpf([0,'1e-8']) for _ in range(4)]
    h=[sum(F[j][i]*J[i] for i in range(4)) for j in range(12)];e=[iv.exp(x) for x in h];Z=sum(e,iv.mpf(0));p=[x/Z for x in e]
    mu=[sum(p[j]*F[j][i] for j in range(12)) for i in range(4)]
    cov=[[sum(p[j]*F[j][i]*F[j][k] for j in range(12))-mu[i]*mu[k] for k in range(4)] for i in range(4)]
    rows=[]
    for i in range(4):
        rs=0.0
        for k in range(4):
            x=iv.mpf(repr(float(GHI*dhi[i])))*cov[i][k];rs+=max(abs(lo_float(x)),abs(hi_float(x)))
        rows.append(float(np.nextafter(rs,np.inf)))
    return rows,max(rows)

initial=contract(np.zeros(4),B,50); assert initial is not None
stack=[('R',initial[0],initial[1],0)];records=[];terms=[]
while stack:
    path,ilo,ihi,depth=stack.pop(); c=contract(ilo,ihi,30)
    if c is None:
        records.append({'path':path,'depth':depth,'status':'EXCLUDED','input_lo':ilo.tolist(),'input_hi':ihi.tolist()});continue
    lo,hi=c; term=None
    if uniform_contained(lo,hi):term='uniform'
    elif contained('localized',lo,hi):term='localized'
    elif contained('saddle',lo,hi):term='saddle'
    if term:
        rec={'path':path,'depth':depth,'status':'TERMINAL_'+term.upper(),'contracted_lo':lo.tolist(),'contracted_hi':hi.tolist()};records.append(rec);terms.append(rec);continue
    if depth>=100: raise RuntimeError(('depth',depth,lo,hi))
    k=int(np.argmax((hi-lo)/B));mid=float((lo[k]+hi[k])/2)
    records.append({'path':path,'depth':depth,'status':'SPLIT','contracted_lo':lo.tolist(),'contracted_hi':hi.tolist(),'dimension':k,'midpoint':mid})
    llo=lo.copy();lhi=hi.copy();lhi[k]=mid;rlo=lo.copy();rhi=hi.copy();rlo[k]=mid
    stack.append((path+'1',rlo,rhi,depth+1));stack.append((path+'0',llo,lhi,depth+1))

rows,q=uniform_banach(); assert q<1
names=sorted(set(r['status'].replace('TERMINAL_','').lower() for r in terms));assert names==['localized','saddle','uniform']
# Saddle is strictly above zero through the whole event parameter box.
assert sad['phi_over_tube'][0]>0
# Event certificate gives localized Phi=0 simultaneously with stationarity for each supplied tuple.
# Its s-box must lie inside the localized stationary branch tube.
coex_slo=np.array([float(x[0]) for x in coex['root_box'][:4]]);coex_shi=np.array([float(x[1]) for x in coex['root_box'][:4]])
assert np.all(coex_slo>=branch_boxes['localized'][0]) and np.all(coex_shi<=branch_boxes['localized'][1])

counts={}
for r in records:counts[r['status']]=counts.get(r['status'],0)+1
out={
 'task':'MP7-017',
 'scientific_state':'PROVED_INTERVAL_ASSISTED_FIRST_GLOBAL_TRANSITION',
 'event_g_box':coex['root_box'][4],
 'event_stationary_energy_certificate':'R7P-026 unique localized stationary root with Phi=0 in 5D (s,g) box',
 'stationary_branch_tubes':[loc,sad],
 'global_exhaustion_records':len(records),'status_counts':counts,'max_depth':max(r['depth'] for r in records),'terminals':terms,
 'uniform_banach_q_upper':q,'uniform_banach_rows':rows,
 'saddle_energy_lower_over_event_box':sad['phi_over_tube'][0],
 'event_global_minimizers_aligned':['uniform','localized'],
 'localized_full_D12_orbit_size':12,
 'first_transition_argument':'at g_eq all primal energies are >=0; for any g<g_eq, V_g(p)=V_geq(p)+(g_eq-g)||X7^T p||^2/2, so every nonuniform p has positive energy; R7P-028 gives negative localized branch slope through the event',
 'conclusion':'for each supplied strict spectral tuple, its certified local equal-energy event is the first global transition: uniform is unique global below it, and at it the uniform state and the 12 localized D12 images are global minimizers',
 'scope':'finite supplied FIN rank-seven model; no physical gain provenance or selector'
}
(WORK/'results/MP7-017_first_global_transition.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ('terminals',)},indent=2))
