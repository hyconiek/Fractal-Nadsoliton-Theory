#!/usr/bin/env python3
"""Exact bounded source algebra/rate checks, not a PDE proof assistant."""
from __future__ import annotations
import argparse
import ast
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
NV=9


def constant(a=0):
    return {(0,)*NV:F(a)} if a else {}


def var(i):
    return {tuple(int(j==i) for j in range(NV)):F(1)}


def add(*ps):
    out={}
    for p in ps:
        for e,a in p.items():
            out[e]=out.get(e,F(0))+a
    return {e:a for e,a in out.items() if a}


def scale(p,a):
    return {e:F(a)*v for e,v in p.items() if a*v}


def mul(p,q):
    out={}
    for e,a in p.items():
        for f,b in q.items():
            key=tuple(x+y for x,y in zip(e,f))
            out[key]=out.get(key,F(0))+a*b
    return {e:a for e,a in out.items() if a}


def power(p,n):
    out=constant(1)
    for _ in range(n):
        out=mul(out,p)
    return out


def diff(p,i):
    out={}
    for e,a in p.items():
        if e[i]:
            k=list(e)
            k[i]-=1
            out[tuple(k)]=a*e[i]
    return out


def bracket(p,q):
    return add(mul(diff(p,0),diff(q,6)),scale(mul(diff(p,6),diff(q,0)),-1),
               mul(diff(p,1),diff(q,7)),scale(mul(diff(p,7),diff(q,1)),-1))


def science():
    checks={}
    def ck(name,condition):
        assert name not in checks
        checks[name]=bool(condition)
        assert condition,name
    # Threshold phase matrix. Variables stand for cos(theta),sin(theta),m.
    co,si,m=[var(i) for i in range(3)]
    mat=((mul(m,co),si),(scale(mul(power(m,2),si),-1),mul(m,co)))
    det=add(mul(mat[0][0],mat[1][1]),scale(mul(mat[0][1],mat[1][0]),-1))
    circle=add(power(co,2),power(si,2))
    ck('threshold_phase_determinant_symbolic',det==mul(power(m,2),circle))
    n11=add(mul(power(m,2),power(mat[0][0],2)),power(mat[1][0],2))
    n22=add(mul(power(m,2),power(mat[0][1],2)),power(mat[1][1],2))
    n12=add(mul(power(m,2),mul(mat[0][0],mat[0][1])),mul(mat[1][0],mat[1][1]))
    ck('threshold_phase_energy_identity',n11==mul(power(m,4),circle)
       and n22==mul(power(m,2),circle) and not n12)
    ck('uniform_high_frequency_cone_gap',F(64)/(64+F(23,4))>F(9,16))
    ck('source_weight_HS_integrable',2*(6-3)>1)
    ck('regularized_symbol_boundary_error_vanishes',8-1-1>0)
    # Exact mixed density (K is factored out); state vars Z,C,U,V,P,R.
    Z,C,U,V,f,g,P,R,h4=[var(i) for i in range(9)]
    eta=add(mul(C,g),V)
    psi=add(mul(Z,f),U)
    e3=add(scale(power(eta,3),2),scale(mul(eta,power(psi,2)),6))
    pure3=add(scale(power(mul(C,g),3),2),scale(mul(mul(C,g),power(mul(Z,f),2)),6))
    mixed3=add(e3,scale(pure3,-1))
    dz=scale(add(mul(mul(mul(C,g),f),U),mul(mul(Z,power(f,2)),V),mul(mul(f,V),U)),12)
    dc=add(scale(mul(mul(C,power(g,2)),V),12),scale(mul(g,power(V,2)),6),
           scale(mul(mul(mul(Z,g),f),U),12),scale(mul(g,power(U,2)),6))
    ck('exact_mixed_cubic_Z_derivative',diff(mixed3,0)==dz)
    ck('exact_mixed_cubic_C_derivative',diff(mixed3,1)==dc)
    ck('exact_mixed_second_Z_derivative',diff(mixed3,0) and
       diff(diff(mixed3,0),0)==scale(mul(power(f,2),V),12))
    ck('exact_mixed_second_C_derivative',
       diff(diff(mixed3,1),1)==scale(mul(power(g,2),V),12))
    ck('exact_mixed_cross_derivative',
       diff(diff(mixed3,0),1)==scale(mul(mul(g,f),U),12))
    e4=add(scale(power(eta,4),F(1,2)),scale(mul(power(eta,2),power(psi,2)),3),
           scale(mul(h4,power(psi,4)),F(1,4)))
    pure4=add(scale(power(mul(C,g),4),F(1,2)),
              scale(mul(power(mul(C,g),2),power(mul(Z,f),2)),3),
              scale(mul(h4,power(mul(Z,f),4)),F(1,4)))
    mixed4=add(e4,scale(pure4,-1))
    def local_gradient_support(poly,degree):
        for j in (0,1,6,7):
            for e in diff(poly,j):
                if not (e[2]+e[3]>=1 and e[4]+e[5]>=1
                        and sum(e[i] for i in (0,1,2,3,6,7))==degree):
                    return False
        return True
    ck('quartic_mixed_gradient_local_support_all_g4',local_gradient_support(mixed4,3))
    # A structural fixture at the allowed w=5/4, d treated as formal factor1.
    G3=add(scale(mul(power(Z,2),R),F(1,2)),scale(mul(mul(Z,P),C),-6),
           scale(mul(power(P,2),R),-2),scale(mul(power(C,2),R),F(-4,3)),
           scale(power(R,3),F(-8,27)))
    mixed_bracket=bracket(mixed3,G3)
    ck('quartic_bracket_gradient_support_fixture',local_gradient_support(mixed_bracket,3))
    # Exact phase-to-position/momentum forcing, not position-velocity guessing.
    qa=lambda a,b:(a[0]+b[0],a[1]+b[1])
    qs=lambda a,c:(a[0]*c,a[1]*c)
    qm=lambda a,b:(a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
    ii=(F(0),F(1))
    phase_cases=0
    for om in (F(2),F(7,3)):
        a,p,r,s,Vv=F(1,5),F(2,7),F(3,11),F(4,13),F(3)
        up,um=(a,p/om),(a,-p/om)
        dp=qa(qa(qs(qm(ii,up),-om),qs(ii,Vv*a/om)),(r,s))
        dm=qa(qa(qs(qm(ii,um),om),qs(ii,-Vv*a/om)),(r,-s))
        adot=qs(qa(dp,dm),F(1,2))
        pdot=qs(qm((F(0),F(-1)),qa(dp,qs(dm,-1))),om/2)
        assert adot==(p+r,0) and pdot==(-(om*om-Vv)*a+om*s,0)
        assert adot!=(p,0)
        phase_cases+=1
    ck('actual_QP_phase_Duhamel_forcing_conversion',phase_cases==2)
    # Exact affine rate endpoints: slopes negative, so these certify open ranges.
    rates=lambda alpha:(1-alpha/2,2-F(3,2)*alpha,3-2*alpha)
    endpoint=rates(F(4,3))
    ck('full_transfer_open_rate_interval',endpoint==(F(1,3),F(0),F(1,3)))
    ck('radiation_open_rate_interval',1-F(2)/2==0 and 2-F(2)==0)
    samples=(F(1,2),F(1),F(5,4),F(13,10))
    ck('full_rate_fixture_positivity',all(min(rates(a))>0 for a in samples))
    ck('endpoint_not_falsely_admitted',min(rates(F(4,3)))==0 and min(rates(F(7,5)))<0)
    return {'status':'PASS_BOUNDED_SOURCE_ALGEBRA_NOT_FULL_PROOF_CERTIFICATE',
            'checks':checks,'check_count':len(checks),'phase_conversion_cases':phase_cases,
            'rate_fixtures':[str(a) for a in samples],
            'arithmetic':'standard-library Fraction sparse polynomials and exact affine rates',
            'limitations':'No machine proof of oscillatory kernel estimates, nonlinear energy bounds or actual trajectory limits.'}


INPUTS=[
    'AGENTS.md','research_tasks/RP004_SHARED_CONTRACT.json',
    'fin_quadratic_localization_research_20261021/PROOF_OR_OPEN_LEMMA.md',
    'fin_quadratic_localization_research_20261021/CHART_AND_REMAINING_GATE.md',
    'fin_quadratic_localization_research_20261021/MANIFEST.sha256',
    'fin_canonical_modal_state_map_20261021/PROOF_OR_OPEN_LEMMA.md',
    'fin_canonical_modal_state_map_20261021/MANIFEST.sha256',
    'fin_spatial_cascade_research_20261021/PROOF_OR_OPEN_LEMMA.md',
    'fin_spatial_cascade_research_20261021/MANIFEST.sha256',
    'fin_even_continuum_threshold_research_20261021/PROOF_OR_OPEN_LEMMA.md',
    'fin_even_continuum_threshold_research_20261021/MANIFEST.sha256',
    'fin_coherent_mode_radiation_review_20261019/PROOF_COMPLETIONS_AND_SCOPE.md',
    'fin_coherent_mode_radiation_review_20261019/MANIFEST.sha256',
]
PROTECTED=[p for p in INPUTS if p.endswith('/MANIFEST.sha256')]


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def check_manifest(path):
    n=0
    for row in path.read_text().splitlines():
        if not row.strip():
            continue
        h,name=row.split(None,1)
        name=name.lstrip('*')
        assert any(p.is_file() and digest(p)==h for p in (path.parent/name,ROOT/name)),(str(path),name)
        n+=1
    return n


def main(record):
    result=science()
    result['protected_manifests']={p:check_manifest(ROOT/p) for p in PROTECTED}
    result['protected_entries']=sum(result['protected_manifests'].values())
    result['input_count']=len(INPUTS)
    for p in HERE.glob('*.py'):
        ast.parse(p.read_text())
    for p in HERE.glob('*.json'):
        json.loads(p.read_text())
    for p in HERE.glob('*.md'):
        txt=p.read_text()
        assert all(line==line.rstrip() for line in txt.splitlines()),str(p)
        prose=re.sub(r'```.*?```','',txt,flags=re.S)
        prose=re.sub(r'`[^`]*`','',prose)
        prose=re.sub(r'(?m)^(?: {4}|\t).*\n?','',prose)
        for link in re.findall(r'\[[^\]]*\]\(([^)]+)\)',prose):
            if not re.match(r'^[a-z]+://',link) and not link.startswith('#'):
                assert (p.parent/link.split('#')[0]).is_file(),(str(p),link)
    inputs=''.join(f'{digest(ROOT/p)}  {p}\n' for p in INPUTS)
    payload=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if record:
        (HERE/'RESULTS.json').write_text(payload)
        (HERE/'INPUTS.sha256').write_text(inputs)
        own=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.sha256')
        (HERE/'MANIFEST.sha256').write_text(''.join(f'{digest(p)}  {p.name}\n' for p in own))
    else:
        assert (HERE/'RESULTS.json').read_text()==payload
        assert (HERE/'INPUTS.sha256').read_text()==inputs
        check_manifest(HERE/'MANIFEST.sha256')
        names={r.split(None,1)[1].lstrip('*') for r in (HERE/'MANIFEST.sha256').read_text().splitlines()}
        assert names=={p.name for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.sha256'}
    print(json.dumps({'status':result['status'],'bounded_checks':result['check_count'],
                      'inputs':len(INPUTS),'protected_entries':result['protected_entries'],
                      'mode':'own_metadata_record' if record else 'read_only'},sort_keys=True))


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    group=ap.add_mutually_exclusive_group(required=True)
    group.add_argument('--record',action='store_true')
    group.add_argument('--check',action='store_true')
    main(ap.parse_args().record)
