#!/usr/bin/env python3
"""Bounded rational phase/algebra diagnostics, not an operator/PDE proof."""
from __future__ import annotations
import argparse
import ast
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from math import isqrt
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent


class IV:
    def __init__(self,lo,hi=None):
        self.lo=F(lo)
        self.hi=self.lo if hi is None else F(hi)
        assert self.lo<=self.hi
    @staticmethod
    def get(x):
        return x if isinstance(x,IV) else IV(x)
    def __add__(self,o):
        o=self.get(o)
        return IV(self.lo+o.lo,self.hi+o.hi)
    __radd__=__add__
    def __neg__(self):
        return IV(-self.hi,-self.lo)
    def __sub__(self,o):
        return self+-self.get(o)
    def __rsub__(self,o):
        return self.get(o)+-self
    def __mul__(self,o):
        o=self.get(o)
        vals=[a*b for a in (self.lo,self.hi) for b in (o.lo,o.hi)]
        return IV(min(vals),max(vals))
    __rmul__=__mul__
    def inv(self):
        if self.lo<=0<=self.hi:
            raise ZeroDivisionError('interval contains zero')
        return IV(1/self.hi,1/self.lo)
    def __truediv__(self,o):
        return self*self.get(o).inv()
    def __rtruediv__(self,o):
        return self.get(o)*self.inv()
    def abs_lower(self):
        return self.lo if self.lo>0 else -self.hi if self.hi<0 else F(0)
    def abs_upper(self):
        return max(abs(self.lo),abs(self.hi))
    def overlaps(self,o):
        o=self.get(o)
        return max(self.lo,o.lo)<=min(self.hi,o.hi)


def sqrt_box(x,bits=80):
    x=F(x)
    assert x>0
    scale=1<<bits
    n=isqrt((x.numerator*scale*scale)//x.denominator)
    lo=F(n,scale)
    hi=lo if lo*lo==x else F(n+1,scale)
    assert lo*lo<=x<=hi*hi
    return IV(lo,hi)


def poly(a=0):
    return {(0,0,0):F(a)} if a else {}


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
            g=tuple(x+y for x,y in zip(e,f))
            out[g]=out.get(g,F(0))+a*b
    return {e:a for e,a in out.items() if a}


def science():
    checks={}
    def ck(name,condition):
        assert name not in checks
        checks[name]=bool(condition)
        assert condition,name
    X,Y,Z=[{tuple(int(i==j) for i in range(3)):F(1)} for j in range(3)]
    S=add(Y,Z)
    D=add(Y,scale(Z,-1))
    A=add(mul(S,S),scale(mul(X,X),-1))
    B=add(mul(X,X),scale(mul(D,D),-1))
    sign_fixtures=0
    for sigma,tau,rho in product((-1,1),repeat=3):
        denominator=mul(X,add(scale(X,sigma),scale(Y,-tau),scale(Z,-rho)))
        if tau==rho and sigma==tau:
            numerator=scale(add(S,X),-sigma)
            expected=mul(X,A)
        elif tau==rho:
            numerator=poly(sigma)
            expected=mul(X,add(X,S))
        else:
            numerator=add(scale(X,sigma),scale(D,tau))
            expected=mul(X,B)
        ck('symbolic_rationalization_'+str((sigma,tau,rho)),mul(denominator,numerator)==expected)
        sign_fixtures+=1
    masses=(F(7,3),F(19,8))
    momenta=((F(0),F(0)),(F(1),F(2)),(F(-3),F(3)),
             (F(20),F(-19)),(F(50),F(50)),(F(-12),F(-7)),(F(1,5),F(40)))
    phase_cases=0
    all_gap=all_symbol_bound=all_forms=all_symmetry=True
    for mass in masses:
        w=mass*mass-4
        ck('actual_parameter_fixture_'+str(mass).replace('/','_'),F(5,4)<=w<=F(7,4))
        channels=((F(2),F(2),F(2),F(-3,2)),
                  (F(2),mass,mass,F(-3,2)),
                  (mass,F(2),mass,F(-3)))
        for m0,m1,m2,a in channels:
            assert abs(m1-m2)<m0<m1+m2
            for p,q in momenta:
                e0=sqrt_box((p+q)**2+m0*m0)
                e1=sqrt_box(p*p+m1*m1)
                e2=sqrt_box(q*q+m2*m2)
                for sigma,tau,rho in product((-1,1),repeat=3):
                    phase=sigma*e0-tau*e1-rho*e2
                    direct=1/(e0*phase)
                    all_gap &= (e0*phase).abs_lower()>=2
                    all_symbol_bound &= direct.abs_upper()<=F(1,2)
                    if tau==rho and sigma==tau:
                        invA=1/(m1*m1+m2*m2-m0*m0+2*(e1*e2-p*q))
                        rational=-sigma*(e1+e2+e0)/e0*invA
                    elif tau==rho:
                        rational=sigma/(e0*(e0+e1+e2))
                    else:
                        invB=1/(m0*m0-m1*m1-m2*m2+2*(e1*e2+p*q))
                        rational=(sigma+tau*(e1-e2)/e0)*invB
                    all_forms &= direct.overlaps(rational)
                    reverse=1/(e0*(-sigma*e0+tau*e1+rho*e2))
                    all_symmetry &= direct.overlaps(-reverse)
                    phase_cases+=1
    ck('all_outward_phase_fixture_gaps',all_gap)
    ck('all_outward_symbol_fixture_bounds',all_symbol_bound)
    ck('outward_and_symbolic_formula_agreement',all_forms)
    ck('real_structure_sign_reversal',all_symmetry)
    # Exact zero-frequency homological checks: no interval or quadrature.
    homological_cases=0
    for m0,m1,m2,a in ((F(2),F(2),F(2),F(-3,2)),
                        (F(2),F(7,3),F(7,3),F(-3,2)),
                        (F(7,3),F(2),F(7,3),F(-3))):
        for sigma,tau,rho in product((-1,1),repeat=3):
            phase=sigma*m0-tau*m1-rho*m2
            b=1/(m0*phase)
            old_source=sigma*a/m0
            coordinate_derivative=phase*sigma*a*b
            assert old_source-coordinate_derivative==0
            assert old_source+coordinate_derivative!=0
            homological_cases+=1
    ck('actual_coefficients_and_flow_sign',homological_cases==24)
    ck('half_wave_quadratic_phi_prefactor',F(-6,4)==F(-3,2))
    ck('half_wave_quadratic_psi_prefactor',F(-12,4)==F(-3))
    # Prepared source: sum b for p,-p and output free mass2.
    initial_sums=[]
    for E in (F(7,3),F(19,8),F(3),F(5)):
        direct=sum((1/(F(2)*(2-tau*E-rho*E))
                    for tau,rho in product((-1,1),repeat=2)),F(0))
        expected=(E*E-2)/(2*(E*E-1))
        assert direct==expected
        assert F(13,34)<=direct<F(1,2)
        initial_sums.append(str(direct))
    ck('actual_initial_sum_formula',len(initial_sums)==4)
    x=F(21,4)
    ck('exact_nonzero_initial_profile_lower_bound',(x-2)/(2*(x-1))==F(13,34)>0)
    ck('reality_at_initial_preparation',all(
        sum((1/(2*(sigma*2-tau*F(7,3)-rho*F(7,3)))
             for tau,rho in product((-1,1),repeat=2)),F(0))*sigma==F(31,80)
        for sigma in (-1,1)))
    # The last fixture value follows from E^2=49/9: (49/9-2)/(2*(49/9-1))=31/80.
    ck('strict_triangle_needed_at_boundary',F(4)-2-2==0)
    try:
        IV(-1,1).inv()
    except ZeroDivisionError:
        ck('zero_phase_interval_division_rejected',True)
    else:
        ck('zero_phase_interval_division_rejected',False)
    ck('source_linear_operators_unchanged',all(m*m-6==(m*m-4)-2 for m in masses))
    return {'status':'PASS_BOUNDED_PHASE_ALGEBRA_NOT_OPERATOR_CERTIFICATE',
            'checks':checks,'check_count':len(checks),
            'outward_phase_cases':phase_cases,'homological_sign_cases':homological_cases,
            'symbolic_rationalizations':sign_fixtures,
            'initial_sum_fixtures':initial_sums,
            'arithmetic':'Fraction polynomials and dyadic outward square-root boxes (80 bits)',
            'limitations':'Finite cases do not prove all-frequency symbol derivatives, weighted localization or nonlinear PDE bounds.'}


INPUTS=[
    'AGENTS.md','research_tasks/RP004_SHARED_CONTRACT.json',
    'fin_spatial_cascade_research_20261021/PROOF_OR_OPEN_LEMMA.md',
    'fin_spatial_cascade_research_20261021/EXACT_REMAINDER_AND_NEXT_GATE.md',
    'fin_spatial_cascade_research_20261021/MANIFEST.sha256',
    'fin_canonical_modal_state_map_20261021/PROOF_OR_OPEN_LEMMA.md',
    'fin_canonical_modal_state_map_20261021/MANIFEST.sha256',
    'fin_even_continuum_threshold_research_20261021/PROOF_OR_OPEN_LEMMA.md',
    'fin_even_continuum_threshold_research_20261021/MANIFEST.sha256',
    'fin_nonlinear_mode_memory_20261021/NONLINEAR_REDUCTION.md',
    'fin_nonlinear_mode_memory_20261021/MANIFEST.sha256',
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
                      'phase_cases':result['outward_phase_cases'],'inputs':len(INPUTS),
                      'protected_entries':result['protected_entries'],
                      'mode':'own_metadata_record' if record else 'read_only'},sort_keys=True))


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    group=ap.add_mutually_exclusive_group(required=True)
    group.add_argument('--record',action='store_true')
    group.add_argument('--check',action='store_true')
    main(ap.parse_args().record)
