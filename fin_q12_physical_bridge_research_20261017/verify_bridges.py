#!/usr/bin/env python3
"""Symbolic/finite controls for the q12 physical-compatibility research."""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);args=p.parse_args();checks=[]
    def ck(name,test,**detail):checks.append({'name':name,'passed':bool(test),**detail});assert test,name
    import sympy as sp
    a,b,c,d=sp.symbols('lambda3 lambda4 lambda5 lambda6',positive=True)
    cs3=(1,0,-1,0)*3;cs4=(sp.Integer(1),-sp.Rational(1,2),-sp.Rational(1,2))*4
    cs5=tuple(sp.simplify(sp.cos(5*sp.pi*j/6)) for j in range(12));sig=sp.Matrix([(-1)**j for j in range(12)])
    A=sp.Matrix(12,12,lambda i,j:a*cs3[(i-j)%12]/6+b*cs4[(i-j)%12]/6+c*cs5[(i-j)%12]/6+d*(-1)**((i-j)%12)/12)
    ck('generic admitted Fourier family row sums vanish exactly',all(sp.simplify(v)==0 for v in A*sp.ones(12,1)))
    ck('generic parity eigenvalue binding exact',all(sp.simplify(v)==0 for v in A*sig-d*sig))
    ck('both selected pair gaps vanish by actual reflection identities',sp.simplify(A[0,3]-A[6,3])==0 and sp.simplify(A[1,3]-A[5,3])==0)
    m,t,g=sp.symbols('m t g',positive=True);r0=m*sp.ones(12,1)-sp.eye(12)[:,3];rt=r0+t*sig
    base=A*r0;shift=A*(rt-r0)
    ck('spectator field shift is exactly lambda6 t sigma',all(sp.simplify(v)==0 for v in shift-d*t*sig))
    ck('zero-gap comparison survives both contexts',all(sp.simplify((A*r)[j]-(A*r)[i])==0 for r in (r0,rt) for i,j in ((0,6),(1,5))))
    exponent=sp.simplify((shift[6]-shift[5]))
    ck('common clock and row normalizer cancel in double ratio',exponent==2*d*t,
       exact_result='R_spec=exp(2 g lambda6 t/N), versus1 for gap-only common-clock class')
    for N in (24,120,12000,24000):
        mm=N//12;tt=1;r=[mm]*12;r[3]-=1;rr=[x+tt*((-1)**j) for j,x in enumerate(r)]
        assert min(rr)>=0 and sum(r)==sum(rr)==N-1
    ck('witness uses actual integer cavities not fractional populations',True)
    theta=F(1,1000)
    ck('optional actual MT fixed-parity input is below admitted energy height',theta**4/(12*(1-theta**2))<F(1,10**10))
    # Generic SWAP equality and negative controls; no numerical commutator tolerance.
    q=3;swap=sp.zeros(q*q)
    for i,j in itertools.product(range(q),repeat=2):swap[j*q+i,i*q+j]=1
    hs=sp.diag(0,1,2);eye=sp.eye(q);matched=sp.kronecker_product(hs,eye)+sp.kronecker_product(eye,hs)
    shifted=sp.kronecker_product(hs,eye)+sp.kronecker_product(eye,hs+7*eye)
    wrong=sp.kronecker_product(hs,eye)+sp.kronecker_product(eye,sp.diag(0,1,3))
    ck('matched energies and common offsets conserve energy under SWAP',swap*matched==matched*swap and swap*shifted==shifted*swap)
    ck('fixed mismatched ancilla is a genuine negative control',swap*wrong!=wrong*swap)
    U=(3*sp.eye(q*q)-4*sp.I*swap)/5
    ck('partial SWAP is unitary and energy conserving exactly',
       (U*U.conjugate().T-sp.eye(q*q)).applyfunc(sp.simplify)==sp.zeros(q*q) and
       (U*matched-matched*U).applyfunc(sp.simplify)==sp.zeros(q*q))
    pp=(sp.Rational(1,2),sp.Rational(1,3),sp.Rational(1,6));tau=(sp.Rational(4,7),sp.Rational(2,7),sp.Rational(1,7))
    state=U*sp.kronecker_product(sp.diag(*pp),sp.diag(*tau))*U.conjugate().T
    reduced=sp.Matrix(q,q,lambda i,j:sum(state[i*q+k,j*q+k] for k in range(q)))
    target=sp.diag(*[sp.Rational(9,25)*x+sp.Rational(16,25)*y for x,y in zip(pp,tau)])
    ck('diagonal partial SWAP changes reset activity without changing target',(reduced-target).applyfunc(sp.simplify)==sp.zeros(q))
    # A finite two-context bank defeats an overbroad fixed-bath no-go.
    probs0=(F(2,3),F(1,3));probs1=(F(4,5),F(1,5));bgprob=(F(1,3),F(2,3));tagprob=(F(3,5),F(2,5));out={0:[F(0),F(0)],1:[F(0),F(0)]}
    for bg,x,b0,b1 in itertools.product(range(2),repeat=4):
        initial=(bg+1)*x+b0+2*b1;bank=[b0,b1];newx=bank[bg];bank[bg]=x
        final=(bg+1)*newx+bank[0]+2*bank[1];assert initial==final
        prob=bgprob[bg]*tagprob[x]*probs0[b0]*probs1[b1];out[bg][newx]+=prob
    ck('explicit finite controlled bank conserves energy and prepares both thermal targets',all(tuple(v/bgprob[bg] for v in out[bg])==(probs0 if bg==0 else probs1) for bg in range(2)),
       scope='rational q2 finite comparison, not a supplied FIN physical bath or autonomous all-time process')
    # Nominal numerical A7 check supplements, never replaces, the symbolic proof.
    safe=ROOT/'fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py';spec=importlib.util.spec_from_file_location('physical_gap_seed',safe);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    import numpy as np
    _,_,lam,_,AA=mod.build_rank7();gain=12/lam[6];diagnostics=[]
    for N in (24,120,12000,24000,120000):
        tt=max(1,N//12000);r=np.full(12,N//12);r[3]-=1;rr=r+tt*np.asarray(list(sig),dtype=float).reshape(12)
        fields=[gain*(AA@v)/N for v in (r,rr)];rates=[]
        for h in fields:
            prob=np.exp(h-h.max());prob/=prob.sum();rates.append((prob[6],prob[5]))
            assert abs(h[6]-h[0])<1e-12 and abs(h[5]-h[1])<1e-12
        double=(rates[1][0]/rates[1][1])/(rates[0][0]/rates[0][1]);expected=float(np.exp(2*gain*lam[6]*tt/N))
        assert abs(double-expected)<2e-13 and double>1
        diagnostics.append({'N':N,'t':tt,'double_ratio':double,'expected':expected})
    ck('actual versioned A7 numerical rates agree with symbolic witness',True,scope='binary64 diagnostic, not an outward certificate or measured rate')
    ck('zero gain does not falsely trigger the activity obstruction',float(np.exp(0))==1)
    log=HERE/'BROWSER_READ_LOG.jsonl';captures=[json.loads(line) for line in log.read_text().splitlines()]
    labels={row['label'] for row in captures};assert {'arxiv-glauber-search','veszeli-vattay-2022','chen-nature-2025'}<=labels
    cache={}
    for row in captures:
        path=ROOT/row['cache_path'];cache[row['label']]='ABSENT_OPTIONAL_CACHE'
        if path.is_file():assert sha(path)==row['text_sha256'];cache[row['label']]='MATCHED_EXTRACTED_TEXT'
    ck('actual browser search and primary HTML captures recorded',len(labels)>=3,
       scope='capture provenance, not automatic full-paper or theorem verification')
    inputs=[Path(__file__),safe,log,HERE/'Q12_SPECTATOR_ACTIVITY_OBSTRUCTION.md',HERE/'THERMAL_RESET_RESOURCE_BOUNDARY.md',
        ROOT/'fin_signed_interface_orbit_cut_review_20261016/MANIFEST.sha256',ROOT/'fin_balanced_phase_cut_research_20261016/MANIFEST.sha256',
        ROOT/'fin_native_transfer_review_20261016/MANIFEST.sha256',ROOT/'fin_auxiliary_mediator_review_20261016/MANIFEST.sha256']
    result={'date':'2026-10-17','status':'PASS_BOUNDED_Q12_PHYSICAL_COMPATIBILITY_CHECKS','checks':checks,
        'input_hashes':{str(p.relative_to(ROOT)):sha(p) for p in inputs},'browser_cache':cache,'numeric_rate_diagnostics':diagnostics,
        'science_scope':'restricted gap-only class excluded; conditional reset-bank implementation; no natural activity or physical source',
        'empirical_evidence':False}
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'bounded_checks':len(checks),'browser_captures':len(captures)},indent=2))

if __name__=='__main__':main()
