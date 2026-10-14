#!/usr/bin/env python3
"""Fail-closed composite FIN/flat decision helper for READY-001.

The helper validates all statistical inputs, binds the frozen model source hash,
and uses mpmath interval arithmetic for grid-point likelihoods/confidence radius.
Grid-to-continuum payment uses the admitted conservative softmax Lipschitz bound.
This is conditional mathematical tooling, not empirical calibration coverage.
"""
from __future__ import annotations
import hashlib, importlib.util, json, math
from decimal import Decimal, InvalidOperation
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
import numpy as np
import mpmath as mp

ROOT=Path(__file__).resolve().parent
CTX=ROOT/'source_context'
MODELS=('FIN_A7','FLAT_P7_TRACE')
K=13
mp.iv.dps=60

class DecisionInputError(ValueError): pass

def sha256(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def _load_model():
    binding=json.loads((ROOT/'PROCESS_BINDING.json').read_text())
    p=CTX/'F_source_model.py'
    if sha256(p)!=binding['source_hashes']['F_source_model.py']:
        raise DecisionInputError('frozen model hash mismatch')
    spec=importlib.util.spec_from_file_location('f_source_model_decision',p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

SM=_load_model()
ROWS={m:np.asarray(SM.all_models()[m][0],float) for m in MODELS}

def _dec_fraction(x,name):
    if isinstance(x,bool): raise DecisionInputError(f'{name}: boolean not allowed')
    s=str(x)
    try:d=Decimal(s)
    except InvalidOperation:raise DecisionInputError(f'{name}: invalid decimal')
    if not d.is_finite():raise DecisionInputError(f'{name}: nonfinite')
    return Fraction(d)

def _positive_int(n,name='n'):
    if isinstance(n,bool):raise DecisionInputError(f'{name}: boolean not allowed')
    if isinstance(n,int):v=n
    elif isinstance(n,str) and n.isdigit():v=int(n)
    else:raise DecisionInputError(f'{name}: finite positive integer required')
    if v<=0:raise DecisionInputError(f'{name}: finite positive integer required')
    return v

def validate_counts(counts,n):
    n=_positive_int(n)
    if not isinstance(counts,(list,tuple)) or len(counts)!=K:raise DecisionInputError('counts: exactly 13 entries required')
    cc=[]
    for i,c in enumerate(counts):
        if isinstance(c,bool):raise DecisionInputError(f'counts[{i}]: boolean not allowed')
        if isinstance(c,int):v=c
        elif isinstance(c,str) and c.isdigit():v=int(c)
        else:raise DecisionInputError(f'counts[{i}]: exact nonnegative integer required')
        if v<0:raise DecisionInputError(f'counts[{i}]: negative')
        cc.append(v)
    if sum(cc)!=n:raise DecisionInputError('counts sum inconsistent with n')
    return cc,n

def _iv_rat(fr:Fraction):
    return mp.iv.mpf([fr.numerator,fr.numerator])/fr.denominator

def _binary_fraction(x:float):
    n,d=float(x).as_integer_ratio();return Fraction(n,d)

def _lo(ivx):
    return math.nextafter(float(ivx.a), -math.inf)

def _hi(ivx):
    return math.nextafter(float(ivx.b), math.inf)

def confidence_radius_interval(n,alpha='0.025'):
    n=_positive_int(n);a=_dec_fraction(alpha,'alpha')
    if not (0<a<1):raise DecisionInputError('alpha must satisfy 0<alpha<1')
    ratio=Fraction(2**K-2,1)/a
    x=mp.iv.sqrt(mp.iv.log(_iv_rat(ratio))/(2*n))
    return {'lower':_lo(x),'upper':_hi(x),'backend':'mpmath.iv','conditional_on_backend_correctness':True}

@lru_cache(maxsize=32)
def _grid_laws(model,lo_s,hi_s,intervals):
    if model not in MODELS:raise DecisionInputError('unadmitted model')
    if not isinstance(intervals,int) or isinstance(intervals,bool) or intervals<=0:raise DecisionInputError('intervals: positive integer required')
    lo=_dec_fraction(lo_s,'gain_lo');hi=_dec_fraction(hi_s,'gain_hi')
    if lo>hi:raise DecisionInputError('gain_lo>gain_hi')
    row=[_binary_fraction(x) for x in ROWS[model]]
    step=(hi-lo)/intervals
    laws=[]
    for i in range(intervals+1):
        g=lo+i*step
        ee=[mp.iv.exp(_iv_rat(g*r/2)) for r in row]
        den=sum(ee)
        laws.append(tuple(e/den for e in ee))
    # Conservative source convention: ||dp/dg||_TV <= range(row)/4, and nearest grid <= step/2.
    rng=max(row)-min(row)
    mesh_allow=rng*step/8
    return laws,mesh_allow

def lower_distance_interval(counts,n,model,detector_tv='0.005',gain_lo='2.97',gain_hi='3.03',intervals=240):
    cc,n=validate_counts(counts,n)
    det=_dec_fraction(detector_tv,'detector_tv')
    if not (0<=det<=1):raise DecisionInputError('detector_tv must be in [0,1]')
    laws,mesh_allow=_grid_laws(model,str(gain_lo),str(gain_hi),intervals)
    emp=[_iv_rat(Fraction(c,n)) for c in cc[:12]]
    emp_invalid=_iv_rat(Fraction(cc[12],n))
    tvs=[]
    for p in laws:
        s=sum(abs(emp[i]-p[i]) for i in range(12))+abs(emp_invalid)
        tvs.append(s/2)
    min_lo=min(_lo(x) for x in tvs)
    min_hi=min(_hi(x) for x in tvs)
    sub=float(mesh_allow+det)
    lower=math.nextafter(max(0.0,min_lo-sub),-math.inf) if min_lo>sub else 0.0
    upper_grid=math.nextafter(max(0.0,min_hi-float(det)),math.inf)
    return {
        'rejection_distance_lower':lower,
        'min_grid_distance_interval':{'lower':min_lo,'upper':min_hi},
        'grid_to_continuum_TV_allowance':float(mesh_allow),
        'detector_TV_radius':float(det),'gain_box':[float(_dec_fraction(gain_lo,'gain_lo')),float(_dec_fraction(gain_hi,'gain_hi'))],
        'grid_intervals':intervals,'backend':'mpmath.iv plus exact-rational binary64 source row',
        'conditional_on_backend_correctness':True,
        'upper_diagnostic_not_membership_certificate':upper_grid
    }

def classify(counts,n,alpha='0.025',hypotheses=None,numerical_guard='1e-12'):
    cc,n=validate_counts(counts,n)
    eps=confidence_radius_interval(n,alpha)
    guard=float(_dec_fraction(numerical_guard,'numerical_guard'))
    if guard<0 or not math.isfinite(guard):raise DecisionInputError('numerical_guard invalid')
    if hypotheses is None:
        hypotheses={m:{'detector_tv':'0.005','gain_lo':'2.97','gain_hi':'3.03','intervals':240} for m in MODELS}
    if set(hypotheses)!=set(MODELS):raise DecisionInputError('hypotheses must specify FIN_A7 and FLAT_P7_TRACE')
    stats={}
    excluded={}
    boundary_open={}
    for m in MODELS:
        h=hypotheses[m]
        s=lower_distance_interval(cc,n,m,h.get('detector_tv','0.005'),h.get('gain_lo','2.97'),h.get('gain_hi','3.03'),int(h.get('intervals',240)))
        stats[m]=s
        d=s['rejection_distance_lower']
        boundary_open[m]=abs(d-eps['upper'])<=guard
        excluded[m]=(d>eps['upper']+guard)
    if any(boundary_open.values()):
        out='NUMERICAL_BOUND_OPEN'; semantic='INCONCLUSIVE'
    elif all(excluded.values()):out='NEITHER_MODEL_FITS';semantic=out
    elif excluded['FIN_A7'] and not excluded['FLAT_P7_TRACE']:out='FLAT_ONLY';semantic=out
    elif excluded['FLAT_P7_TRACE'] and not excluded['FIN_A7']:out='FIN_ONLY';semantic=out
    else:out='INCONCLUSIVE';semantic=out
    return {
      'status':'CONDITIONAL_DECISION_TOOL_ONLY','outcome':out,'decision_semantic':semantic,
      'n':n,'counts':cc,'alpha_validation':float(_dec_fraction(alpha,'alpha')),
      'confidence_radius_TV_interval':eps,'lower_distance_bounds':stats,'excluded_models':excluded,
      'numerical_boundary_open':boundary_open,'numerical_guard':guard,
      'iid_certified':False,'calibration_coverage_certified':False,'physical_evidence_certified':False,
      'interpretation':'Non-rejection is not certified model membership. Tool is licensed only after external calibration/sampling premises are admitted.'
    }

def main(argv=None):
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('counts_json',help='JSON file with {n, counts[13], optional alpha,hypotheses}')
    ns=ap.parse_args(argv)
    try:
        d=json.loads(Path(ns.counts_json).read_text())
        r=classify(d['counts'],d['n'],d.get('alpha','0.025'),d.get('hypotheses'))
        print(json.dumps(r,indent=2,sort_keys=True));return 0
    except (DecisionInputError,KeyError,json.JSONDecodeError,ValueError) as e:
        print(json.dumps({'status':'FAIL_CLOSED','error':str(e),'physical_evidence_certified':False},indent=2));return 2
if __name__=='__main__':raise SystemExit(main())
