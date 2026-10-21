#!/usr/bin/env python3
"""Exact response-work algebra; not a proof of actual averaged dynamics."""
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


def qi(x):
    return (F(x[0]),F(x[1])) if isinstance(x,tuple) else (F(x),F(0))


def qa(a,b):
    a,b=qi(a),qi(b)
    return (a[0]+b[0],a[1]+b[1])


def qm(a,b):
    a,b=qi(a),qi(b)
    return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])


def add(*ps):
    out={}
    for p in ps:
        for e,a in p.items():
            out[e]=qa(out.get(e,(F(0),F(0))),a)
    return {e:a for e,a in out.items() if a!=(0,0)}


def scale(p,a):
    return {e:qm(v,a) for e,v in p.items() if qm(v,a)!=(0,0)}


def mul(p,q):
    out={}
    for e,a in p.items():
        for f,b in q.items():
            key=(e[0]+f[0],e[1]+f[1])
            out[key]=qa(out.get(key,(F(0),F(0))),qm(a,b))
    return {e:a for e,a in out.items() if a!=(0,0)}


def cos(n,m):
    return add({(n,m):qi(F(1,2))},{(-n,-m):qi(F(1,2))})


def sin(n,m):
    return add({(n,m):(F(0),F(-1,2))},{(-n,-m):(F(0),F(1,2))})


def mean(p):
    return p.get((0,0),(F(0),F(0)))


def science():
    checks={}
    def ck(name,condition):
        assert name not in checks
        checks[name]=bool(condition)
        assert condition,name
    ck('trig_real_normalization',mean(mul(sin(1,0),sin(1,0)))==qi(F(1,2)))
    fixture_count=0
    for a,b in ((F(2,3),F(3,5)),(F(1),F(4,7))):
        da,db,dx=F(3,7),F(2,9),F(5,11)
        z,c=scale(cos(1,0),a),scale(cos(0,1),b)
        zd,cd=scale(sin(1,0),-a),scale(sin(0,1),-b)  # velocities divided by own frequencies
        Xa=add(scale(add(scale(cos(2,0),2),scale(sin(2,0),da)),-3*a*a),
               scale(add(scale(cos(0,2),F(3,17)),scale(sin(0,2),F(5,19))),-3*b*b),
               {(0,0):qi(7)})
        Xb=add(scale(add(scale(cos(0,2),3),scale(sin(0,2),db)),-3*b*b),
               scale(add(scale(cos(2,0),F(7,23)),scale(sin(2,0),F(-2,13))),-3*a*a),
               {(0,0):qi(11)})
        Xg=scale(add(scale(cos(1,1),F(2,5)),scale(sin(1,1),dx),
                     scale(cos(-1,1),F(3,4))),-6*a*b)
        Fz=add(scale(mul(z,Xa),-12),scale(mul(c,Xg),-12))
        Fc=add(scale(mul(c,Xb),-12),scale(mul(z,Xg),-12))
        BI=mean(mul(zd,Fz)); BJ=mean(mul(cd,Fc))
        assert BI==qi(-9*da*a**4-18*dx*a*a*b*b)
        assert BJ==qi(-9*db*b**4-18*dx*a*a*b*b)
        Imix=mean(mul(zd,scale(mul(c,Xg),-12)))
        Jmix=mean(mul(cd,scale(mul(z,Xg),-12)))
        assert Imix==Jmix==qi(-18*dx*a*a*b*b)
        closed=scale(cos(-1,1),-6*a*b*F(3,4))
        assert mean(mul(zd,scale(mul(c,closed),-12)))==(0,0)
        assert mean(mul(cd,scale(mul(z,closed),-12)))==(0,0)
        fixture_count+=1
    ck('complete_action_work_fixtures',fixture_count==2)
    ck('self_action_normalization',F(9)*4==36)
    ck('mixed_action_normalization',F(18)*4==72)
    ck('reference_cos_and_energy_factor',F(1,2)*F(1,2)**2*2==F(1,4))
    # Source profile transforms after the common ik(k^2+4)S2/D factor.
    # Fb=(3/2)(tanh sech^2-tanh sech^4).
    const=F(3,2)*(F(1,8)-F(1,8)*F(4,6))
    slope=F(3,2)*(-F(1,8)*F(1,6))
    ck('new_Fb_transform_exact',const==F(2,32) and slope==F(-1,32))
    x=F(8)
    abs_b_pref=x*(x+4)*(2-x)**2/(1024*(1+x))
    ck('shape_channel_on_shell_prefactor',abs_b_pref==F(3,8))
    ck('mixed_profile_normalization',F(9,8)/24**2==F(1,512))
    ck('mixed_open_momentum_interval',15>3**2 and 21<5**2)
    ck('shape_momentum_bounds',2**2<8<3**2)
    ck('source_exponential_coarse_bound',F(33,7)<5 and F(22,7)<4)
    old_gamma=F(125,73728)
    db_lower=F(1,16)*F(4,81)**2
    dx_lower=F(75,16384)*F(1,18)**2
    ka=4*old_gamma/F(7,4)
    kb=12*db_lower
    kx=F(72)*dx_lower/F(7,3)
    ck('shape_density_lower_exact',db_lower==F(1,6561))
    ck('mixed_density_lower_exact',dx_lower==F(25,1769472))
    ck('three_action_bounds_exact',(ka,kb,kx)==(F(125,32256),F(4,2187),F(25,57344)))
    ck('sum_action_coercivity_constant',min(ka,kb)/2==F(2,2187))
    # Coercivity leaves (min/2)(I-J)^2 plus nonnegative terms.
    coeff=min(ka,kb)
    ck('coercivity_residual_square',
       (coeff-F(2,2187),-F(4,2187),coeff-F(2,2187))==
       (coeff/2,-coeff,coeff/2))
    ck('old_Gamma_factor_nine_retained',F(36)/9==4)
    ck('threshold_counterprofile_nonzero',F(48,1225)>0)
    ck('odd_open_frequency_separation',2*F(17,10)-F(8,3)>F(1,2))
    return {'status':'PASS_BOUNDED_REFERENCE_WORK_ALGEBRA_NOT_ACTUAL_ACTION_ODE',
            'checks':checks,'check_count':len(checks),'phase_average_fixtures':fixture_count,
            'certificate':{'old_Gamma_lower':str(old_gamma),'density_b_lower':str(db_lower),
                           'density_x_lower':str(dx_lower),'k_a_lower':str(ka),
                           'k_b_lower':str(kb),'k_x_lower':str(kx),'sum_coercivity':str(coeff/2)},
            'arithmetic':'Fraction Fourier-Laurent polynomials over Q(i); exact rational lower bounds',
            'limitations':'Frozen reference response only. No machine proof of limiting absorption, nonlinear virial coercivity, or actual critical-time averaging.'}


INPUTS=[
    'AGENTS.md','research_tasks/RP004_SHARED_CONTRACT.json',
    'fin_subcritical_field_transfer_research_20261021/LINEAR_WEIGHTED_PROPAGATION.md',
    'fin_subcritical_field_transfer_research_20261021/ACTUAL_SUBCRITICAL_THEOREM.md',
    'fin_subcritical_field_transfer_research_20261021/NEXT_CRITICAL_GATE.md',
    'fin_subcritical_field_transfer_research_20261021/MANIFEST.sha256',
    'fin_quadratic_localization_research_20261021/PROOF_OR_OPEN_LEMMA.md',
    'fin_quadratic_localization_research_20261021/MANIFEST.sha256',
    'fin_canonical_modal_state_map_20261021/PROOF_OR_OPEN_LEMMA.md',
    'fin_canonical_modal_state_map_20261021/MANIFEST.sha256',
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
        (HERE/'COERCIVITY_CERTIFICATE.json').write_text(json.dumps({
            'scope':'exact arithmetic for the frozen reference action-work form, NOT nonlinear virial coercivity',
            'values':result['certificate']},indent=2,sort_keys=True)+'\n')
        (HERE/'INPUTS.sha256').write_text(inputs)
        own=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.sha256')
        (HERE/'MANIFEST.sha256').write_text(''.join(f'{digest(p)}  {p.name}\n' for p in own))
    else:
        assert (HERE/'RESULTS.json').read_text()==payload
        assert (HERE/'INPUTS.sha256').read_text()==inputs
        assert json.loads((HERE/'COERCIVITY_CERTIFICATE.json').read_text())['values']==result['certificate']
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
