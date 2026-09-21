#!/usr/bin/env python3
from __future__ import annotations
import os
import json, math, sys, hashlib
from pathlib import Path
from fractions import Fraction
import numpy as np
import mpmath as mp

WORK=Path(os.environ.get('MP7_WORK_ROOT','/mnt/data/fin_rank7_mathphysics_next'))
R7P=Path(os.environ.get('R7P_ROOT','/mnt/data/r7p_source/extracted/fin_rank7_followup'))
sys.path.insert(0,str(R7P))
from src import coexistence_certificate as cc
mp.iv.dps=60
iv=mp.iv
G=3.7

# accepted spectral rational intervals from the R7P equal-energy certificate (same strict tuple enclosure)
coex=json.loads((R7P/'certificates/R7P-026_equal_energy_event.json').read_text())
raw={int(k):v for k,v in coex['spectral_intervals'].items()}

def rat(s): return float(Fraction(s))
Llo=np.array([rat(raw[k][0]) for k in (3,4,5,6)])
Lhi=np.array([rat(raw[k][1]) for k in (3,4,5,6)])
dlo=np.array([Llo[0]/6,Llo[1]/6,Llo[2]/6,Llo[3]/12])
dhi=np.array([Lhi[0]/6,Lhi[1]/6,Lhi[2]/6,Lhi[3]/12])
a_lo=np.sqrt(dlo); a_hi=np.sqrt(dhi)
# outward inflate float parameter endpoints
for arr,side in [(dlo,-1),(a_lo,-1)]:
    arr[:] = np.nextafter(arr, -np.inf)
for arr,side in [(dhi,1),(a_hi,1)]:
    arr[:] = np.nextafter(arr, np.inf)
B=np.nextafter(G*dhi,np.inf)

# exact unscaled character table using interval sqrt(3)
rt3=iv.sqrt(iv.mpf(3)); half=iv.mpf('0.5')
c3=[1,0,-1,0]*3
c4=[iv.mpf(1),-half,-half]*4
c5=[iv.mpf(1),-rt3/2,half,iv.mpf(0),-half,rt3/2,-iv.mpf(1),rt3/2,-half,iv.mpf(0),half,-rt3/2]
alt=[1 if j%2==0 else -1 for j in range(12)]
F=[[iv.mpf(c3[j]),c4[j],c5[j],iv.mpf(alt[j])] for j in range(12)]

def ivpoint(x):
    # Float x denotes a covering endpoint. Inflate one ulp both ways before interval evaluation.
    return iv.mpf([repr(float(np.nextafter(x,-np.inf))),repr(float(np.nextafter(x,np.inf)))])

def lo_float(x): return float(np.nextafter(float(x.a),-np.inf))
def hi_float(x): return float(np.nextafter(float(x.b), np.inf))

def mean_at_point(J):
    X=[ivpoint(x) for x in J]
    h=[sum(F[j][i]*X[i] for i in range(4)) for j in range(12)]
    hmax=None  # interval exp is safe without shifting at this bounded domain
    e=[iv.exp(x) for x in h]; Z=sum(e,iv.mpf(0)); p=[x/Z for x in e]
    mu=[sum(p[j]*F[j][i] for j in range(12)) for i in range(4)]
    return mu

def contract(lo,hi,maxiter=30):
    lo=lo.copy(); hi=hi.copy()
    for _ in range(maxiter):
        ml=mean_at_point(lo); mu=mean_at_point(hi)
        nl=lo.copy(); nh=hi.copy()
        for i in range(4):
            # Ginibre/cooperativity theorem gives m_i >=0 exactly.
            mlb=max(0.0,lo_float(ml[i]))
            mub=max(0.0,hi_float(mu[i]))
            tl=np.nextafter(G*dlo[i]*mlb,-np.inf)
            tu=np.nextafter(G*dhi[i]*mub, np.inf)
            nl[i]=max(lo[i],tl); nh[i]=min(hi[i],tu)
        if np.any(nl>nh): return None
        delta=max(float(np.max(np.abs(nl-lo))),float(np.max(np.abs(nh-hi))))
        lo,hi=nl,nh
        if delta<2e-15: break
    return lo,hi

roots=json.loads((WORK/'results/MP7-015_g37_local_interval_roots.json').read_text())['certificates']
root_sboxes={r['name']:(np.array([float(x[0]) for x in r['box']]),np.array([float(x[1]) for x in r['box']])) for r in roots}

def spreimage(lo,hi):
    slo=np.nextafter(lo/a_hi,-np.inf); shi=np.nextafter(hi/a_lo,np.inf)
    return slo,shi

def contained_in_source(name,lo,hi):
    slo,shi=spreimage(lo,hi); L,U=root_sboxes[name]
    return bool(np.all(slo>=L) and np.all(shi<=U))

def uniform_contained(lo,hi):
    # global terminal is in the explicit J cube on which Banach contraction is paid below
    return bool(np.all(lo>=-1e-18) and np.all(hi<=1e-8))

# certify the local uniform Banach cube with interval Jacobian in J variables.
def uniform_banach():
    J=[iv.mpf([0,'1e-8']) for _ in range(4)]
    h=[sum(F[j][i]*J[i] for i in range(4)) for j in range(12)]
    e=[iv.exp(x) for x in h]; Z=sum(e,iv.mpf(0)); p=[x/Z for x in e]
    mu=[sum(p[j]*F[j][i] for j in range(12)) for i in range(4)]
    cov=[[sum(p[j]*F[j][i]*F[j][k] for j in range(12))-mu[i]*mu[k] for k in range(4)] for i in range(4)]
    rows=[]
    for i in range(4):
        rs=0.0
        for k in range(4):
            x=iv.mpf(repr(float(G*dhi[i])))*cov[i][k]
            a=lo_float(x); b=hi_float(x); rs += max(abs(a),abs(b))
        rows.append(float(np.nextafter(rs,np.inf)))
    # T(0)=0 and ||DT||_inf <=q<1 implies cube invariance and unique fixed point 0.
    return rows,max(rows)

# Root-preserving binary tree. A contracted shell is excluded by the isotone fixed-point lemma.
initial=contract(np.zeros(4),B,50)
assert initial is not None
stack=[('R',initial[0],initial[1],0)]
records=[]; terminals=[]
while stack:
    path,ilo,ihi,depth=stack.pop()
    c=contract(ilo,ihi,30)
    if c is None:
        records.append({'path':path,'depth':depth,'input_lo':ilo.tolist(),'input_hi':ihi.tolist(),'status':'EXCLUDED_EMPTY_AFTER_ISOTONE_CONTRACTION'})
        continue
    lo,hi=c
    term=None
    if uniform_contained(lo,hi): term='uniform'
    elif contained_in_source('localized',lo,hi): term='localized'
    elif contained_in_source('saddle',lo,hi): term='saddle'
    if term:
        rec={'path':path,'depth':depth,'input_lo':ilo.tolist(),'input_hi':ihi.tolist(),'contracted_lo':lo.tolist(),'contracted_hi':hi.tolist(),'status':'TERMINAL_'+term.upper()}
        records.append(rec); terminals.append(rec); continue
    if depth>=90:
        raise RuntimeError(f'depth cap with unresolved box {path}: {lo} {hi}')
    widths=(hi-lo)/B
    k=int(np.argmax(widths)); mid=float((lo[k]+hi[k])/2)
    # overlapping split at exact stored float midpoint covers the contracted parent.
    llo=lo.copy(); lhi=hi.copy(); lhi[k]=mid
    rlo=lo.copy(); rhi=hi.copy(); rlo[k]=mid
    records.append({'path':path,'depth':depth,'input_lo':ilo.tolist(),'input_hi':ihi.tolist(),'contracted_lo':lo.tolist(),'contracted_hi':hi.tolist(),'status':'SPLIT','dimension':k,'midpoint':mid})
    stack.append((path+'1',rlo,rhi,depth+1)); stack.append((path+'0',llo,lhi,depth+1))

rows,q=uniform_banach()
counts={}
for r in records: counts[r['status']]=counts.get(r['status'],0)+1
terminal_names=[r['status'].replace('TERMINAL_','').lower() for r in terminals]
assert sorted(set(terminal_names))==['localized','saddle','uniform']
assert q<1

# Rigorous energies on the two nonzero source Krawczyk boxes.
_,_,Cnorm,_=cc.interval_features()
energy_intervals={}
for r in roots:
    XX=[iv.mpf(x) for x in r['box']]+[iv.mpf('3.7')]
    FF,_,_,_=cc.interval_eval(XX,Cnorm)
    energy_intervals[r['name']]=[lo_float(FF[4]),hi_float(FF[4])]
assert energy_intervals['localized'][0]>0 and energy_intervals['saddle'][0]>0

out={
 'task':'MP7-016',
 'scientific_state':'PROVED_INTERVAL_ASSISTED_GLOBAL_STATIONARY_EXHAUSTION_AT_G_37_10',
 'gain':'37/10',
 'coordinate':'unscaled aligned fields J=(J3,J4,J5,J6)',
 'stationary_equation':'J_i = g d_i m_i(J), d=(lambda3/6,lambda4/6,lambda5/6,lambda6/12)',
 'initial_box_lo':[0,0,0,0],
 'initial_box_hi':B.tolist(),
 'spectral_intervals':raw,
 'contractor_basis':'accepted nonnegative covariances imply every m_i is isotone; on [l,u], m_i(l)<=m_i(J)<=m_i(u)',
 'processed_records':len(records),
 'status_counts':counts,
 'max_depth':max(r['depth'] for r in records),
 'terminals':terminals,
 'uniform_banach_cube_J_upper':1e-8,
 'uniform_DTout_inf_row_upper':rows,
 'uniform_contraction_q_upper':q,
 'nonzero_local_uniqueness_source':'MP7-015 parametric Krawczyk s-boxes; every terminal J box maps wholly into its source s-box for all accepted feature scalings',
 'root_count_aligned_nonnegative':3,
 'root_names':['uniform','saddle','localized'],
 'energy_uniform_exact':0.0,
 'energy_nonzero_root_intervals':energy_intervals,
 'global_minimum_conclusion':'uniform theta=0 is the unique global minimizer at g=37/10 after MP7-011 phase alignment; both other stationary roots have rigorously positive energy',
 'tree_records':records,
 'scope':'all aligned nonnegative stationary points at exact g=37/10 over the accepted strict spectral enclosure; full-X7 global minimizers reduce to this chart by MP7-011, boundary nonzero roots excluded by MP7-014'
}
(WORK/'results/MP7-016_g37_global_exhaustion.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ('tree_records',)},indent=2))
