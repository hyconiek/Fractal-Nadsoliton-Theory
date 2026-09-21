"""FR15: strengthened FR9 parity arm discovered by residual-wall search."""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
import json
from frontier_local_boxes import raw_box, _short
ROOT=Path(__file__).resolve().parents[1]

def fr15_record():
    good=raw_box(F(1,6500),F(1,2432),F(1,8192),F(1,14000))
    bad=raw_box(F(1,6500),F(1,2432),F(1,8192),F(1,13000))
    assert good['status']=='INTERVAL_CERTIFIED'
    assert bad['status']=='FAILED' and not bad['boundary_ok']
    return {
      'id':'FR15-strengthened-FR9-parity-arm','status':'INTERVAL_CERTIFIED_REPLAYED',
      'proof_type':'anisotropic second-order rational interval AD + R7P-044 characteristic/inertia disjunction',
      'domain':{
        'abs_r_minus_rstar':'<=1/6500','u=1-exp(-3J4/2)':'<=1/2432',
        'v=1-exp(-J5/2)':'<=1/8192','e=1-q_even':'<=1/14000',
        'J_fields':'J3,J4,J5,J6>=0 with the exact physical shared-field odd law'},
      'conclusion':'lambda2(Mtilde)<=sigma_* and hence lambda2(M4)<=sigma_* throughout the declared box.',
      'strict_checks':{
        'signs':good['signs'],'boundary_schur_pass':good['boundary_ok'],'P1_zero_endpoint_schur_pass':good['endpoint_ok'],
        'boundary_Nvv_interval':_short(good['boundary_schur'][3]),'endpoint_Nvv_interval':_short(good['endpoint_schur'][3]),
        'c2_interval':_short(good['c2'].v)},
      'negative_control':{
        'e=1/13000':{'status':bad['status'],'boundary_schur_pass':bad['boundary_ok'],'endpoint_pass':bad['endpoint_ok'],
                     'interpretation':'checker failure only; not a physical counterexample'}},
      'research_origin':'Fixed-seed residual search repeatedly hit the old FR9/FR11 artificial wall; the enlarged interval box removes that wall.',
      'scope':'Physical four-amplitude local equality-neighborhood theorem only; no full-seven-coordinate transfer.'}

def write_result(path=None):
    r=fr15_record(); path=Path(path or ROOT/'results/FR15_strengthened_parity_arm.json'); path.write_text(json.dumps(r,indent=2)+'\n'); return r
if __name__=='__main__': print(json.dumps(write_result(),indent=2))
