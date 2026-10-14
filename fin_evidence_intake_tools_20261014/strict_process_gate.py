#!/usr/bin/env python3
"""Model-bound process/instrument validator for READY-001.

Checks the admitted finite N2 quotient contract at software level. It is not a
physical symmetry, calibration, or provenance certificate.
"""
from __future__ import annotations
import hashlib, importlib.util, json, math, sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
CTX=ROOT/'source_context'
Q=12
KINETICS=('heat_bath','metropolis','barker')
MODELS=('FIN_A7','FULL_POTTS_TRACE','FLAT_P7_TRACE')

def sha256(path):
    h=hashlib.sha256();
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def load_module(name,path):
    spec=importlib.util.spec_from_file_location(name,path); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def _finite(a): return bool(np.isfinite(np.asarray(a,dtype=float)).all())

def _stochastic_columns(K,tol=2e-13):
    K=np.asarray(K,float)
    return bool(K.shape==(13,12) and _finite(K) and np.min(K)>=-tol and np.max(np.abs(K.sum(axis=0)-1))<=tol)

def _prob(p,tol=2e-13):
    p=np.asarray(p,float)
    return bool(p.shape==(12,) and _finite(p) and np.min(p)>=-tol and abs(float(p.sum())-1)<=tol)

def validate(model_id='FIN_A7',kind='heat_bath',g=3.0,p0=None,K=None, shortcut=False):
    binding=json.loads((ROOT/'PROCESS_BINDING.json').read_text())
    required={
        'F_source_model.py':CTX/'F_source_model.py',
        'F_positive_process.py':CTX/'F_positive_process.py',
        'F_PROCESS_CONTRACT.json':CTX/'F_PROCESS_CONTRACT.json',
        'F_INSTRUMENT_CONTRACT.json':CTX/'F_INSTRUMENT_CONTRACT.json',
    }
    hash_checks={k:(sha256(v)==binding['source_hashes'][k]) for k,v in required.items()}
    errors=[]; warnings=[]
    if not all(hash_checks.values()): errors.append('SOURCE_HASH_MISMATCH')
    if model_id not in MODELS: errors.append('UNADMITTED_MODEL_ID')
    if kind not in KINETICS: errors.append('UNADMITTED_KINETIC_ID')
    if not isinstance(g,(int,float)) or not math.isfinite(float(g)): errors.append('NONFINITE_GAIN')
    if errors:
        return {'status':'FAIL_CLOSED','errors':errors,'hash_checks':hash_checks,'physical_claim':False}
    sm=load_module('f_source_model_ready',CTX/'F_source_model.py')
    sys.modules['source_model']=sm
    pp=load_module('f_positive_process_ready',CTX/'F_positive_process.py')
    A=np.asarray(sm.all_models()[model_id],float)
    if A.shape!=(12,12) or not _finite(A):errors.append('INVALID_A_SHAPE_OR_NONFINITE')
    # Exact quotient premise: symmetric circulant, constant diagonal, zero-field model class.
    circ=max(float(np.max(np.abs(A[i]-np.roll(A[0],i)))) for i in range(Q))
    sym=float(np.max(np.abs(A-A.T)))
    diag=float(np.ptp(np.diag(A)))
    if circ>2e-12:errors.append('NONCIRCULANT_INTERACTION')
    if sym>2e-12:errors.append('NONSYMMETRIC_INTERACTION')
    if diag>2e-12:errors.append('NONCONSTANT_DIAGONAL')
    B=np.asarray(pp.difference_generator(A,float(g),kind),float)
    rowsum=float(np.max(np.abs(B.sum(axis=1))))
    minoff=float(np.min(B-np.diag(np.diag(B))))
    if not _finite(B) or rowsum>5e-12 or minoff<-5e-13:errors.append('INVALID_GENERATOR')
    pstat=np.asarray(pp.stationary_difference(A,float(g)),float)
    if p0 is None:p0=np.ones(Q)/Q
    if K is None:K=np.vstack([np.eye(Q),np.zeros((1,Q))])
    if not _prob(p0):errors.append('INVALID_INITIAL_LAW')
    if not _stochastic_columns(K):errors.append('INVALID_INSTRUMENT')
    shortcut_error=None
    if shortcut:
        if kind!='heat_bath': errors.append('RESET_SHORTCUT_ONLY_HEAT_BATH')
        reset=2*(np.ones((Q,1))*pstat[None,:]-np.eye(Q))
        shortcut_error=float(np.max(np.abs(B-reset)))
        if shortcut_error>5e-12:errors.append('RESET_SHORTCUT_INCONSISTENT_WITH_BOUND_MODEL')
    return {
        'status':'PASS_SOFTWARE_CONTRACT_ONLY' if not errors else 'FAIL_CLOSED',
        'model_id':model_id,'kinetic':kind,'gain':float(g),'source_hash_checks':hash_checks,
        'circulant_max_abs':circ,'symmetry_max_abs':sym,'diagonal_range':diag,
        'generator_row_sum_max_abs':rowsum,'generator_min_offdiag':minoff,
        'initial_law_valid':_prob(p0),'instrument_valid':_stochastic_columns(K),
        'reset_shortcut_checked':bool(shortcut),'reset_shortcut_max_abs':shortcut_error,
        'errors':errors,'warnings':warnings,'physical_symmetry_certified':False,
        'calibration_certified':False
    }

def main():
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--model',default='FIN_A7');ap.add_argument('--kinetic',default='heat_bath');ap.add_argument('--gain',type=float,default=3.0);ap.add_argument('--shortcut',action='store_true')
    ns=ap.parse_args();r=validate(ns.model,ns.kinetic,ns.gain,shortcut=ns.shortcut);print(json.dumps(r,indent=2,sort_keys=True));raise SystemExit(0 if r['status'].startswith('PASS') else 2)
if __name__=='__main__':main()
