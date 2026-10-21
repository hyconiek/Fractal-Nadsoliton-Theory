#!/usr/bin/env python3
"""Exact finite algebra for ET, not a whole-PDE or asymptotic certificate.

pi, sqrt3, and sqrt2 are kept as named formal factors. No floating value
is used to decide a cancellation. Metadata authorship is own-packet only.
"""
from __future__ import annotations
import argparse
import ast
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent


def trim(a):
    a=list(map(F,a))
    while a and not a[-1]:
        a.pop()
    return a


def add(*ps):
    out=[F(0)]*max((len(p) for p in ps),default=0)
    for p in ps:
        for i,a in enumerate(p):
            out[i]+=a
    return trim(out)


def scale(p,a):
    return trim([F(a)*x for x in p])


def mul(p,q):
    out=[F(0)]*max(0,len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q):
            out[i+j]+=a*b
    return trim(out)


def power(p,n):
    out=[F(1)]
    for _ in range(n):
        out=mul(out,p)
    return out


def diff(p):
    return trim([i*p[i] for i in range(1,len(p))])


def integral(p):
    return sum((a*F(2,i+1) for i,a in enumerate(p) if i%2==0),F(0))


def val(p,t):
    return sum((a*F(t)**i for i,a in enumerate(p)),F(0))


def science():
    checks={}
    def ck(name,condition):
        assert name not in checks
        checks[name]=bool(condition)
        assert condition,name
    one=[F(1)]
    x=[F(0),F(1)]
    # In this block x=k^2. Fourier recurrence ratios remain rational.
    for n in (2,3,4,5,6,8):
        ratio=[F(n,n+1),F(1,n*(n+1))]
        even=add([2,F(-1)+F(3,n)],scale(ratio,-3))
        even_expected=scale([1,1],F(2-n,n+1))
        ck('Jost_even_reduction_p'+str(n),even==even_expected)
        odd=add([F(-2,n)+3,F(1,n)],scale(ratio,F(-3*(n+1),n+2)))
        odd_expected=scale([4,1],F(n-1,n*(n+2)))
        ck('Jost_odd_reduction_p'+str(n),odd==odd_expected)
    ck('bound_f_continuum_transform_zero',scale([1,1],F(2-2,3))==[])
    Bnum=mul([1,1],add([F(-1,4)],scale([F(3,4),F(1,12)],F(1,2))))
    ck('actual_B_profile_transform',Bnum==scale(mul([1,1],[3,1]),F(1,24)))
    ck('actual_J_extra_k_squared_factor',
       add(one,scale([F(2,3),F(1,6)],F(-3,2)))==scale(x,F(-1,4)))
    # Now t=tanh(x); ordinary polynomial integration is exact.
    t=[F(0),F(1)]
    s2=[F(1),F(0),F(-1)]
    resonance=[F(-1,2),F(0),F(3,2)]
    Dx=lambda a:mul(s2,diff(a))
    H=lambda a:add(scale(Dx(Dx(a)),-1),mul([-2,0,6],a))
    ck('actual_even_threshold_ODE',H(resonance)==scale(resonance,4))
    ck('threshold_state_nonzero_at_both_ends',val(resonance,-1)==val(resonance,1)==1)
    R=lambda n:integral(mul(resonance,power(s2,(n-2)//2)))
    ck('threshold_orthogonal_to_bound_f',R(2)==0)
    ck('threshold_C_unnormalized',R(6)==F(-32,105))
    ck('threshold_early_cancellation_R6_R8',R(6)==R(8))
    ck('threshold_R10',R(10)==F(-1024,3465))
    # int t^(2j)*sqrt(1-t^2) dt / pi = Catalan(j)/(2*4^j).
    half=lambda j:F(comb(2*j,j), (j+1)*2*4**j)
    B0pi=F(1,2)*(3*half(2)-half(1))
    ck('threshold_B_unnormalized_pi_factor',B0pi==F(1,32))
    ck('normalized_B_diagonal_factor',F(9,8)*B0pi**2==F(9,8192))
    ck('normalized_C_diagonal_factor',F(27,64)*R(6)**2==F(48,1225))
    ck('rank_one_threshold_matrix',F(9,8192)*F(48,1225)==F(54,1120**2))
    # q=pi^2 below. Compute the null-transform numerator at k^2.
    numerator=add(scale([F(5,16),F(-1,24)],F(-64,105)),
                  [F(256,315),F(-8,105)])
    unnormalized=scale(numerator,F(1,2))  # D(0)=2.
    ck('null_profile_constant_cancellation',F(-32,105)+F(1024,105)*B0pi==0)
    ck('null_profile_second_coefficient',unnormalized==[F(98,315),F(-8,315)])
    ck('normalized_null_second_coefficient_without_sqrt3',
       scale(unnormalized,F(3,8))==[F(49,420),F(-1,105)])
    ck('null_second_coefficient_positive_pi_upper_bound',
       F(49)-4*F(22,7)**2==F(465,49)>0)
    # Actual mixed moment, retaining BOTH shape and continuum contributions.
    Fprofile=scale(mul(t,power(s2,2)),F(3,4))
    HF=H(Fprofile)
    ck('actual_HF_profile',HF==add(scale(mul(t,power(s2,2)),-9),
                                  scale(mul(t,power(s2,3)),18)))
    J_over_sqrt3_after_jacobian=scale(mul(resonance,t),F(1,2))
    JF=integral(mul(J_over_sqrt3_after_jacobian,Fprofile))
    JHF=integral(mul(J_over_sqrt3_after_jacobian,HF))
    ck('actual_JF_zero_not_shape_only',JF==0)
    ck('actual_JHF_nonzero',JHF==F(-32,385))
    ck('actual_g4_zero_source_fourth_derivative',-72*JHF==F(2304,385)>0)
    ck('actual_g4_nonzero_source_initial_coefficient',-F(3,8)*R(6)==F(4,35))
    ck('mixed_k4_spectral_density_without_sqrt3',F(-1,16)*F(1,4)==F(-1,64))
    ck('mixed_gaussian_tail_prefactor_without_sqrt3_over_sqrtpi',F(-1,64)*6==F(-3,32))
    # Rational uniform carrier separation; no rounding at band boundaries.
    wl,wh=F(5,4),F(7,4)
    ck('uniform_omega_brackets',F(11,10)**2<wl and wh<F(4,3)**2)
    ck('uniform_shape_brackets',F(17,10)**2<3<F(7,4)**2)
    ck('uniform_continuum_mass_brackets',F(9,4)**2<4+wl and 4+wh<F(12,5)**2)
    gaps=[F(9,4)-F(4,3),3*F(11,10)-F(12,5),
          F(11,10)+F(17,10)-F(12,5),F(9,4)-(F(7,4)-F(11,10))]
    ck('all_declared_carrier_gaps',min(gaps)==F(2,5)>F(1,4))
    return {'status':'PASS_BOUNDED_EXACT_ALGEBRA_NOT_FULL_PDE_PROOF',
            'checks':checks,'check_count':len(checks),
            'arithmetic':'Fraction; pi and radicals treated as formal factors',
            'limitations':'No stationary-phase machine certificate, actual slow envelope, or full growing-window PDE bound'}


INPUTS=[
    'AGENTS.md',
    'research_tasks/RP004_SHARED_CONTRACT.json',
    'fin_canonical_modal_state_map_20261021/PROOF_OR_OPEN_LEMMA.md',
    'fin_canonical_modal_state_map_20261021/STATE_MAP_AND_REMAINING_GATE.md',
    'fin_canonical_modal_state_map_20261021/MANIFEST.sha256',
    'fin_nonlinear_mode_memory_20261021/NONLINEAR_REDUCTION.md',
    'fin_nonlinear_mode_memory_20261021/SPECTRAL_MEMORY_KERNEL.md',
    'fin_nonlinear_mode_memory_20261021/MANIFEST.sha256',
    'fin_slow_memory_modulation_20261021/PROOF_OR_OPEN_LEMMA.md',
    'fin_slow_memory_modulation_20261021/MANIFEST.sha256',
    'fin_coherent_mode_radiation_review_20261019/PROOF_COMPLETIONS_AND_SCOPE.md',
    'fin_coherent_mode_radiation_review_20261019/MANIFEST.sha256',
    'fin_rp004_flux_review_20261019/BOOTSTRAP_COMPLETION.md',
    'fin_rp004_flux_review_20261019/MANIFEST.sha256',
]
PROTECTED=[x for x in INPUTS if x.endswith('/MANIFEST.sha256')]


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def check_manifest(path):
    count=0
    for row in path.read_text().splitlines():
        if not row.strip():
            continue
        h,name=row.split(None,1)
        name=name.lstrip('*')
        assert any(p.is_file() and digest(p)==h for p in (path.parent/name,ROOT/name)),(str(path),name)
        count+=1
    return count


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
        names={row.split(None,1)[1].lstrip('*') for row in (HERE/'MANIFEST.sha256').read_text().splitlines()}
        assert names=={p.name for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.sha256'}
    print(json.dumps({'status':result['status'],'bounded_checks':result['check_count'],
                      'inputs':len(INPUTS),'protected_entries':result['protected_entries'],
                      'mode':'record_own_packet_only' if record else 'read_only'},sort_keys=True))


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    group=ap.add_mutually_exclusive_group(required=True)
    group.add_argument('--record',action='store_true')
    group.add_argument('--check',action='store_true')
    main(ap.parse_args().record)
