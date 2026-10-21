#!/usr/bin/env python3
"""Exact rational identities; not a proof of the unclosed trace budget."""
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


def pt(a):
    a=list(map(F,a))
    while a and not a[-1]:
        a.pop()
    return a


def pa(*ps):
    out=[F(0)]*max((len(p) for p in ps),default=0)
    for p in ps:
        for j,a in enumerate(p):
            out[j]+=a
    return pt(out)


def ps(p,a):
    return pt([v*F(a) for v in p])


def pm(p,q):
    out=[F(0)]*max(0,len(p)+len(q)-1)
    for j,a in enumerate(p):
        for k,b in enumerate(q):
            out[j+k]+=a*b
    return pt(out)


def pp(p,n):
    out=[F(1)]
    for _ in range(n):
        out=pm(out,p)
    return out


def pd(p):
    return pt([j*p[j] for j in range(1,len(p))])


class RF:
    def __init__(self,n,d=(1,)):
        self.n,self.d=pt(n),pt(d)
        assert self.d
    def deriv(self):
        return RF(pa(pm(pd(self.n),self.d),ps(pm(self.n,pd(self.d)),-1)),pm(self.d,self.d))
    def __eq__(self,o):
        return pm(self.n,o.d)==pm(o.n,self.d)


NV=17
def v(i):
    return {tuple(int(j==i) for j in range(NV)):F(1)}


def aa(*ps_):
    out={}
    for p in ps_:
        for e,a in p.items():
            out[e]=out.get(e,F(0))+a
    return {e:a for e,a in out.items() if a}


def ss(p,a):
    return {e:F(a)*z for e,z in p.items() if a*z}


def mm(p,q):
    out={}
    for e,a in p.items():
        for f,b in q.items():
            g=tuple(x+y for x,y in zip(e,f))
            out[g]=out.get(g,F(0))+a*b
    return {e:a for e,a in out.items() if a}


def dx(p):
    # phi0..3, a0..2, w0..2, pi0..1, F10..1, F20, mass2, spare
    mapping={0:1,1:2,2:3,4:5,5:6,7:8,8:9,10:11,12:13}
    out={}
    for e,a in p.items():
        for i,j in mapping.items():
            if e[i]:
                k=list(e); k[i]-=1; k[j]+=1; k=tuple(k)
                out[k]=out.get(k,F(0))+a*e[i]
    return {e:a for e,a in out.items() if a}


def science():
    checks={}
    def ck(name,condition):
        assert name not in checks
        checks[name]=bool(condition)
        assert condition,name
    r=[F(1),F(0),F(1)]
    phi_prime=RF([1],r)
    ck('actual_arctan_third_derivative',phi_prime.deriv().deriv()==RF([-2,0,6],pp(r,3)))
    log_zeta_prime=RF([0,-1],r)
    ck('actual_conjugated_potential',log_zeta_prime.deriv()==RF([-1,0,1],pp(r,2)))
    a=RF([1],pp(r,3))
    ck('actual_auxiliary_weight_second_derivative',a.deriv().deriv()==RF([-6,0,42],pp(r,5)))
    upper=pa(ps(pp(r,2),12),[6,0,-42])
    sos=pa(ps(pp([F(-3,4),0,1],2),12),[F(45,4)])
    ck('auxiliary_weight_upper_bound_SOS',upper==sos)
    lower=pa(ps(pp(r,2),12),[-6,0,42])
    ck('auxiliary_weight_lower_bound_positive_coefficients',lower==[6,0,66,0,12])
    ck('weight_log_derivative_bound_SOS',
       pa(ps(pp(r,2),F(9,4)),[0,0,-9])==ps(pp([-1,0,1],2),F(9,4)))
    # The interval variance/Fubini kernel, exactly for arbitrary s.
    integral_xy=pa(ps(pm([1,0,-1],[1,1]),F(1,2)),ps(pm([1,-1],[-1,0,1]),F(-1,2)))
    ck('interval_mean_Fubini_kernel_exact',integral_xy==[1,0,-1])
    # Exact forced-wave local identities modulo one spatial derivative.
    ph0,ph1,ph2,ph3,a0,a1,a2,w0,w1,w2,p0,p1,f10,f11,f20,m2,_=[v(i) for i in range(NV)]
    position_flux=aa(mm(ph0,w1),ss(mm(ph1,w0),F(1,2)))
    timeI=aa(mm(aa(mm(ph0,aa(p1,f11)),ss(mm(ph1,aa(p0,f10)),F(1,2))),p0),
             mm(position_flux,aa(w2,ss(mm(m2,w0),-1),f20)))
    B=aa(mm(ph1,mm(w1,w1)),ss(mm(ph3,mm(w0,w0)),F(-1,4)))
    EI=aa(mm(aa(mm(ph0,f11),ss(mm(ph1,f10),F(1,2))),p0),mm(position_flux,f20))
    fluxI=aa(ss(mm(ph0,aa(mm(p0,p0),mm(w1,w1),ss(mm(m2,mm(w0,w0)),-1))),F(1,2)),
             ss(mm(ph1,mm(w0,w1)),F(1,2)),ss(mm(ph2,mm(w0,w0)),F(-1,4)))
    ck('forced_virial_local_identity_exact',aa(timeI,B,ss(EI,-1))==dx(fluxI))
    timeJ=mm(a0,aa(mm(aa(p0,f10),p0),mm(w0,aa(w2,ss(mm(m2,w0),-1),f20))))
    mainJ=aa(mm(a0,mm(p0,p0)),ss(mm(a0,mm(w1,w1)),-1),
             ss(mm(a2,mm(w0,w0)),F(1,2)),ss(mm(m2,mm(a0,mm(w0,w0))),-1),
             mm(a0,aa(mm(f10,p0),mm(w0,f20))))
    fluxJ=aa(mm(a0,mm(w0,w1)),ss(mm(a1,mm(w0,w0)),F(-1,2)))
    ck('forced_momentum_auxiliary_identity_exact',aa(timeJ,ss(mainJ,-1))==dx(fluxJ))
    # Darboux threshold normalization, with t=tanh(x) and Dx=(1-t^2)d_t.
    Dt=lambda p:pm([1,0,-1],pd(p))
    r0=pt([F(-1,2),0,F(3,2)])
    Tstar_r0=pa(Dt(Dt(r0)),ps(pm([0,1],Dt(r0)),3),ps(r0,2))
    ck('actual_Darboux_threshold_image',Tstar_r0==[2])
    ck('actual_inverse_threshold_image',ps([-1,0,3],F(1,2))==r0)
    # Rational certificate, with the actual pi upper bound already paid.
    c_grad,c_q=F(122,7),F(33,14)
    cpos=F(1,24)
    grad_margin=F(3,4)-cpos*c_grad
    mean_margin=F(1,4)-cpos*c_q
    ck('position_coercivity_margins',(grad_margin,mean_margin)==(F(1,42),F(17,112)))
    aux=F(23,4)+6+F(9,2)
    ck('uniform_mass_auxiliary_constant',aux==F(65,4)<17)
    delta=cpos/(2*17)
    cphase=min(delta,cpos-17*delta)
    cyoung=cphase/2
    trace=F(1,2)/cyoung
    ck('phase_coercivity_constants',(delta,cphase,cyoung,trace)==(F(1,816),F(1,816),F(1,1632),F(816)))
    return {'status':'PASS_BOUNDED_VIRIAL_ALGEBRA_WITH_UNPAID_TRACE_BUDGET',
            'checks':checks,'check_count':len(checks),
            'certificate':{'position_coercivity':str(cpos),'gradient_margin':str(grad_margin),
                           'mean_margin':str(mean_margin),'delta':str(delta),
                           'phase_coercivity':str(cphase),'after_Young':str(cyoung),
                           'free_frame_trace_coefficient':str(trace)},
            'arithmetic':'Fraction rational functions and exact differential-polynomial flux identities',
            'limitations':'No bound on the actual integrated trace, no full modal FGR coercivity, no actual damping ODE or physical source theorem.'}


INPUTS=[
    'AGENTS.md','research_tasks/RP004_SHARED_CONTRACT.json',
    'fin_quadratic_localization_research_20261021/PROOF_OR_OPEN_LEMMA.md',
    'fin_quadratic_localization_research_20261021/MANIFEST.sha256',
    'fin_subcritical_field_transfer_research_20261021/LINEAR_WEIGHTED_PROPAGATION.md',
    'fin_subcritical_field_transfer_research_20261021/MANIFEST.sha256',
    'fin_two_mode_dissipation_research_20261021/VIRIAL_SCOPE_AND_THRESHOLD_GATE.md',
    'fin_two_mode_dissipation_research_20261021/NEXT_JOINT_VIRIAL_GATE.md',
    'fin_two_mode_dissipation_research_20261021/MANIFEST.sha256',
    'fin_spatial_cascade_research_20261021/PROOF_OR_OPEN_LEMMA.md',
    'fin_spatial_cascade_research_20261021/MANIFEST.sha256',
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
        (HERE/'CONTINUUM_COERCIVITY_CERTIFICATE.json').write_text(json.dumps({
            'scope':'exact continuum coercivity with rank defect; NOT full joint modal FGR coercivity',
            'values':result['certificate']},indent=2,sort_keys=True)+'\n')
        (HERE/'INPUTS.sha256').write_text(inputs)
        own=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.sha256')
        (HERE/'MANIFEST.sha256').write_text(''.join(f'{digest(p)}  {p.name}\n' for p in own))
    else:
        assert (HERE/'RESULTS.json').read_text()==payload
        assert (HERE/'INPUTS.sha256').read_text()==inputs
        assert json.loads((HERE/'CONTINUUM_COERCIVITY_CERTIFICATE.json').read_text())['values']==result['certificate']
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
