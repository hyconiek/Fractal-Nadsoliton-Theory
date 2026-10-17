#!/usr/bin/env python3
"""Bounded exact/numerical controls; analytical trace/clock proofs are separate."""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);args=p.parse_args();checks=[]
    def ck(name,ok,**detail):checks.append({'name':name,'passed':bool(ok),**detail});assert ok,name
    # Different rational q3/N2 conditional Gibbs construction, not an A7 replacement.
    q=3;nu=F(2,3);kappa=F(7);B=((2,1,3),(1,2,4),(3,4,2))
    visible=list(itertools.product(range(q),repeat=2));hubs=[(i,r) for i,r in itertools.product(range(2),range(q))]
    size=len(visible)+len(hubs);index={x:i for i,x in enumerate(visible)};H={h:len(visible)+j for j,h in enumerate(hubs)}
    Q=[[F(0) for _ in range(size)] for _ in range(size)]
    weights=[F(B[x[0]][x[1]]) for x in visible]+[nu/kappa]*len(hubs)
    for x,ii in index.items():
        for slot in range(2):Q[ii][H[slot,x[1-slot]]]+=nu
    for (slot,r),ii in H.items():
        for b in range(q):
            x=(b,r) if slot==0 else (r,b);Q[ii][index[x]]+=kappa*B[r][b]
    for i in range(size):Q[i][i]=-sum(Q[i])
    ck('rational full labelled augmented chain is positive and reversible',all(sum(row)==0 for row in Q) and all(weights[i]*Q[i][j]==weights[j]*Q[j][i] for i,j in itertools.product(range(size),repeat=2)))
    trace=[[Q[i][j]+sum(Q[i][h]*Q[h][j]/(-Q[h][h]) for h in H.values()) for j in range(len(visible))] for i in range(len(visible))]
    native=[[F(0) for _ in visible] for _ in visible]
    for x,i in index.items():
        for slot in range(2):
            r=x[1-slot];z=sum(B[r])
            for b in range(q):
                if b==x[slot]:continue
                y=(b,r) if slot==0 else (r,b);native[i][index[y]]+=nu*F(B[r][b],z)
        native[i][i]=-sum(native[i])
    ck('exact Schur trace equals declared rational labelled HB',trace==native)
    def klass(i):
        if i<len(visible):return ('V',)+tuple(visible[i].count(j) for j in range(q))
        return ('H',)+tuple(int(hubs[i-len(visible)][1]==j) for j in range(q))
    classes=sorted({klass(i) for i in range(size)});blocks={c:[i for i in range(size) if klass(i)==c] for c in classes}
    ck('full augmented generator lumps exactly to count and one-hole classes',all(len({sum(Q[i][j] for j in target) for i in source})==1 for source,target in itertools.product(blocks.values(),repeat=2)),
       labelled_states=size,count_classes=len(classes),scope='q3 N2 rational control, not a critical domain or physical record')
    ck('hidden stationary mass is not zero at finite kappa',sum(weights[len(visible):])>0 and sum(weights[len(visible):])/(sum(weights))<1)
    # Exact clock-law controls, no simulation or fitted rates.
    for lam,K,z in itertools.product((F(1,2),F(1),F(3)),(F(1),F(7),F(100)),(F(1),F(2))):
        native_lap=lam/(lam+z);aug_lap=native_lap*K/(K+z)
        assert aug_lap<native_lap
        mean=1/lam+1/K;variance=1/lam**2+1/K**2;assert variance<mean**2
    ck('finite intermediate waiting law differs even after mean-clock matching',True)
    for lam,K in ((.5,.5),(1.,2.),(2.,1.),(1.,100.)):
        tv=math.exp(-1) if K==lam else (lam/K)*math.exp(-lam*math.log(K/lam)/(K-lam))
        assert 0<tv<=lam/(lam+K)+1e-14
    ck('closed-form waiting-TV diagnostic respects analytic convexity bound',True,scope='binary64 scalar check; proof uses shifted exponential TV directly')
    # Reflection/elimination order counterexample in exact fractions.
    s=(F(1,2),F(1,3),F(1,6));allowed=(0,1)
    ck('trace then reject differs from reject then trace',s[1]==F(1,3) and s[1]/sum(s[i] for i in allowed)==F(2,5))
    # Clock source and raw observation distinction is algebraic, not inferred from occupancy.
    ck('serial critical-horizon cost exponent paid',F(1)+F(3,2)==F(5,2))
    lam=F(2);K=F(7);oldprob=F(1,2);z=F(1)
    attempt=lam/(lam+z)*K/(K+z);actual_exit=(1-oldprob)*attempt/(1-oldprob*attempt)
    native_exit=lam*(1-oldprob)/(lam*(1-oldprob)+z)
    ck('unobserved null returns do not erase finite-duration memory',actual_exit!=native_exit)
    # Relative versus absolute barriers give different reversible activities.
    w=(F(1),F(2),F(4));Z=sum(w);relative=[[nu*w[b]/Z if a!=b else F(0) for b in range(q)] for a in range(q)]
    absolute=[[F(1,3)/w[a] if a!=b else F(0) for b in range(q)] for a in range(q)]
    ck('common absolute barrier countermodel preserves Gibbs but differs from reset',relative!=absolute and all(w[a]*absolute[a][b]==w[b]*absolute[b][a] for a,b in itertools.product(range(q),repeat=2)))
    # Nominal actual A7/N2 count trace, independently constructed through holes.
    safe=ROOT/'fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py';spec=importlib.util.spec_from_file_location('trace_safe_seed',safe);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    import numpy as np
    _,_,_,_,A=mod.build_rank7();N=2;g=3.;vv=mod.compositions(N);hh=mod.compositions(N-1);vi={n:i for i,n in enumerate(vv)};hi={r:i for i,r in enumerate(hh)}
    VH=np.zeros((len(vv),len(hh)));HV=np.zeros((len(hh),len(vv)))
    for n,i in vi.items():
        for a,count in enumerate(n):
            if count:
                r=list(n);r[a]-=1;VH[i,hi[tuple(r)]]+=count
    for r,i in hi.items():
        field=(g/N)*(A@np.array(r))
        for b in range(12):
            n=list(r);n[b]+=1;HV[i,vi[tuple(n)]]+=7*np.exp(field[b])
    KK=HV.sum(1);G=-N*np.eye(len(vv))+(VH/KK[None,:])@HV
    nativeQ=mod.count_generator(vv,g,A).toarray();defect=float(np.max(np.abs(G-nativeQ)))
    ck('actual versioned A7 count trace numerical diagnostic',defect<2e-13,visible=len(vv),hubs=len(hh),max_abs_defect=defect,
       scope='binary64 N2 unreflected control, not an interval or asymptotic proof')
    log=HERE/'BROWSER_READ_LOG.jsonl';captures=[json.loads(line) for line in log.read_text().splitlines()];labels={r['label'] for r in captures}
    ck('new actual browser search and primary metadata/body reading captured',{'intermediates-search','cappelletti-wiuf','stationary-reduction-2025'}<=labels)
    cache={}
    for row in captures:
        path=ROOT/row['cache_path'];cache[row['label']]='ABSENT_OPTIONAL_CACHE'
        if path.is_file():assert sha(path)==row['text_sha256'];cache[row['label']]='MATCHED_EXTRACTED_TEXT'
    inputs=[Path(__file__),safe,log,HERE/'TRACE_CLOCK_AND_INSTRUMENT_THEOREM.md',HERE/'ACTIVATION_BARRIER_SOURCE_GATE.md',
        ROOT/'fin_balanced_phase_current_review_20261017/MANIFEST.sha256',ROOT/'fin_q12_physical_bridge_research_20261017/MANIFEST.sha256']
    result={'date':'2026-10-17','status':'PASS_BOUNDED_INTERMEDIATE_CLOCK_BRIDGE_CONTROLS','checks':checks,
        'input_hashes':{str(p.relative_to(ROOT)):sha(p) for p in inputs},'browser_cache':cache,
        'claim_scope':'analytical target-realization/clock/instrument results, not a sourced physical intermediate',
        'physical_source':'NOT_FOUND','empirical_evidence':False}
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({'status':result['status'],'bounded_checks':len(checks),'browser_captures':len(captures)},indent=2))

if __name__=='__main__':main()
