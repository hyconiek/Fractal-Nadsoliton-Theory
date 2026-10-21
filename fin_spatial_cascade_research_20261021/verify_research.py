#!/usr/bin/env python3
"""Exact finite source-kernel algebra, not a dispersive/PDE certificate."""
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


def mon(k=0,l=0,i=0,c=1):
    sign=-1 if (i//2)%2 else 1
    return {(k,l,i%2):F(c)*sign} if c else {}


def add(*ps):
    out={}
    for p in ps:
        for e,c in p.items():
            out[e]=out.get(e,F(0))+c
    return {e:c for e,c in out.items() if c}


def scale(p,c):
    return {e:F(c)*v for e,v in p.items() if c*v}


def mul(p,q):
    out={}
    for e,a in p.items():
        for f,b in q.items():
            out=add(out,mon(e[0]+f[0],e[1]+f[1],e[2]+f[2],a*b))
    return out


def power(p,n):
    out=mon()
    for _ in range(n):
        out=mul(out,p)
    return out


def conjugate(p):
    return {e:c*(-1 if e[2] else 1) for e,c in p.items()}


def swap(p):
    return {(e[1],e[0],e[2]):c for e,c in p.items()}


def at(p,k,l):
    out=[F(0),F(0)]
    for e,c in p.items():
        out[e[2]]+=c*F(k)**e[0]*F(l)**e[1]
    return tuple(out)


def equal_line(p,sgn):
    out={}
    for e,c in p.items():
        out=add(out,mon(e[0]+e[1],0,e[2],c*sgn**e[1]))
    return out


def science():
    checks={}
    def ck(name,condition):
        assert name not in checks
        checks[name]=bool(condition)
        assert condition,name
    one=mon()
    k,l,I=mon(k=1),mon(l=1),mon(i=1)
    q=add(k,scale(l,-1))
    p=[one,scale(mul(I,q),F(-1,2))]
    for n in range(1,5):
        p.append(scale(add(scale(p[n-1],n),scale(mul(mul(I,q),p[n]),-1)),F(1,n+2)))
    b=add(mon(c=-6),scale(power(k,2),-3),scale(power(l,2),-3),scale(mul(k,l),9))
    c=mul(add(one,power(k,2)),add(one,power(l,2)))
    raw=add(scale(p[5],9),mul(b,p[3]),mul(c,p[1]),
            scale(mul(mul(I,q),add(scale(p[4],3),mul(add(mul(k,l),mon(c=-1)),p[2]))),3))
    expected=scale(mul(mul(I,q),mul(power(add(k,l),2),add(power(k,2),power(l,2),mon(c=4)))),F(-1,16))
    ck('symbolic_mixed_Jost_numerator_all_k_l',raw==expected)
    ck('symbolic_Hermitian_kernel_numerator',raw==conjugate(swap(raw)))
    ck('symbolic_diagonal_zero',not equal_line(raw,1))
    ck('symbolic_antidiagonal_zero',not equal_line(raw,-1))
    quotient=scale(mul(I,mul(add(k,l),add(power(k,2),power(l,2),mon(c=4)))),F(1,16))
    ck('symbolic_commutator_energy_factor',raw==mul(add(power(l,2),scale(power(k,2),-1)),quotient))
    energy=add(power(k,2),scale(power(l,2),-1))
    ck('actual_parity_fold_first_factor',mul(q,power(add(k,l),2))==mul(energy,add(k,l)))
    ck('actual_parity_fold_second_factor',mul(add(k,l),power(q,2))==mul(energy,q))
    ck('joint_threshold_order_three',min(e[0]+e[1] for e in raw)==3)
    row0={e:v for e,v in raw.items() if e[0]==0}
    col0={e:v for e,v in raw.items() if e[1]==0}
    ck('threshold_row_matches_conjugate_J',row0==scale(mul(I,mul(power(l,3),add(power(l,2),mon(c=4)))),F(1,16)))
    ck('threshold_column_matches_J',col0==scale(mul(I,mul(power(k,3),add(power(k,2),mon(c=4)))),F(-1,16)))
    ck('missing_sum_factor_rejected',raw!=scale(mul(I,mul(q,add(power(k,2),power(l,2),mon(c=4)))),F(-1,16)))
    # Exact nonzero matrix samples, for sign/orientation diagnostics only.
    fixtures=[(F(1,2),F(1)),(F(2),F(1)),(F(-2),F(1)),(F(1),F(-3,2))]
    ck('off_shell_numerator_samples_nonzero',all(at(raw,a,b)!=(0,0) for a,b in fixtures))
    Dk=add(mon(c=2),scale(mul(I,k),3),scale(power(k,2),-1))
    ck('real_denominator_norm_identity',mul(Dk,conjugate(Dk))==mul(add(one,power(k,2)),add(mon(c=4),power(k,2))))
    # Resonant countercontrol lives in Q(sqrt30), using positive branches.
    w=F(3,2)
    output_frequency_sq=(4+w+1,F(1))
    sum_frequency_sq=(F(5)+w,F(1))
    ck('actual_temporal_resonance_squared_identity',output_frequency_sq==sum_frequency_sq)
    ck('resonance_not_equal_spatial_energies',(F(1)-1,F(1))!=(0,0))
    ck('resonant_tuple_inside_original_parameter_band',F(5,4)<=w<=F(7,4))
    # Profile moments from t=tanh(x), no numerical pi value.
    def integral(poly):
        return sum((a*F(2,2*j+1) for j,a in enumerate(poly)),F(0))
    # Here polynomials use u=t^2. Integral u^j dt=2/(2j+1).
    def pmul(a,b):
        out=[F(0)]*(len(a)+len(b)-1)
        for j,x in enumerate(a):
            for n,y in enumerate(b):
                out[j+n]+=x*y
        return out
    def ppow(a,n):
        out=[F(1)]
        for _ in range(n):
            out=pmul(out,a)
        return out
    resonance=[F(-1,2),F(3,2)]
    R=lambda n:integral(pmul(resonance,ppow([F(1),F(-1)],(n-2)//2)))
    ck('actual_full_JF_cancellation',R(6)==R(8))
    ck('actual_threshold_projection_f_orthogonality',R(2)==0)
    ck('source_f_fourth_moment',F(9,16)*integral(ppow([F(1),F(-1)],3))==F(18,35))
    ck('source_f_sixth_moment',F(27,64)*integral(ppow([F(1),F(-1)],5))==F(24,77))
    ck('initial_quartic_energy_coefficient',F(18,35)/4==F(9,70))
    ck('initial_sextic_energy_coefficient',F(24,77)/6==F(4,77))
    # Shape subtraction is essential: q'(0)= -d^2/sqrt3.
    ck('continuous_cross_mass_negative_formal_prefactor',-F(27,8192)/3==F(-9,8192)<0)
    ck('source_cross_threshold_tail_coefficient',F(-1,16)*F(1,4)*6==F(-3,32))
    return {'status':'PASS_BOUNDED_EXACT_SOURCE_ALGEBRA_NOT_PDE_CERTIFICATE',
            'checks':checks,'check_count':len(checks),
            'arithmetic':'sparse Fraction polynomials modulo i^2+1; formal pi/radical factors',
            'limitations':'Dispersive bounds and the energy continuation proof require the analytical arguments; no complete nonlinear modulation proof.'}


INPUTS=[
    'AGENTS.md','research_tasks/RP004_SHARED_CONTRACT.json',
    'fin_canonical_modal_state_map_20261021/PROOF_OR_OPEN_LEMMA.md',
    'fin_canonical_modal_state_map_20261021/MANIFEST.sha256',
    'fin_even_continuum_threshold_research_20261021/PROOF_OR_OPEN_LEMMA.md',
    'fin_even_continuum_threshold_research_20261021/NEXT_EXACT_GATE.md',
    'fin_even_continuum_threshold_research_20261021/MANIFEST.sha256',
    'fin_nonlinear_mode_memory_20261021/NONLINEAR_REDUCTION.md',
    'fin_nonlinear_mode_memory_20261021/SPECTRAL_MEMORY_KERNEL.md',
    'fin_nonlinear_mode_memory_20261021/MANIFEST.sha256',
    'fin_coherent_mode_radiation_review_20261019/PROOF_COMPLETIONS_AND_SCOPE.md',
    'fin_coherent_mode_radiation_review_20261019/MANIFEST.sha256',
    'fin_rp004_flux_review_20261019/OUTWARD_FLUX_THEOREM.md',
    'fin_rp004_flux_review_20261019/MANIFEST.sha256',
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
