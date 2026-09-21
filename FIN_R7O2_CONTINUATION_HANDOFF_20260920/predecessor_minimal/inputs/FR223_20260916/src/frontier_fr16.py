"""FR16: wall-adapted r-v box from the compact residual search."""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
import json
from frontier_local_boxes import raw_box, _short
ROOT=Path(__file__).resolve().parents[1]

def fr16_record():
    good=raw_box(F(1,5000),F(1,8192),F(1,5600),F(1,100000))
    bad=raw_box(F(1,4800),F(1,8192),F(1,5600),F(1,100000))
    assert good['status']=='INTERVAL_CERTIFIED'
    assert bad['status']=='FAILED' and not bad['boundary_ok']
    return {
      'id':'FR16-wall-adapted-rv-box','status':'INTERVAL_CERTIFIED_REPLAYED',
      'proof_type':'anisotropic second-order rational interval AD + R7P-044 characteristic/inertia disjunction',
      'domain':{'abs_r_minus_rstar':'<=1/5000','u=1-exp(-3J4/2)':'<=1/8192',
                'v=1-exp(-J5/2)':'<=1/5600','e=1-q_even':'<=1/100000',
                'J_fields':'J3,J4,J5,J6>=0 with the exact physical shared-field odd law'},
      'conclusion':'lambda2(Mtilde)<=sigma_* and hence lambda2(M4)<=sigma_* throughout the declared wall-adapted box.',
      'strict_checks':{'signs':good['signs'],'boundary_schur_pass':good['boundary_ok'],'P1_zero_endpoint_schur_pass':good['endpoint_ok'],
                       'boundary_Nvv_interval':_short(good['boundary_schur'][3]),'endpoint_Nvv_interval':_short(good['endpoint_schur'][3]),
                       'c2_interval':_short(good['c2'].v)},
      'negative_control':{'rx=1/4800':{'status':bad['status'],'boundary_schur_pass':bad['boundary_ok'],
                                      'interpretation':'checker failure only; not a physical counterexample'}},
      'research_origin':'Residual search after FR15 converged to the FR14 r-wall near v≈1/5654. Narrowing v to 1/5600 permits a substantially longer certified r-arm.',
      'scope':'Physical four-amplitude local equality-neighborhood theorem only; no full-seven-coordinate transfer.'}

def write_result(path=None):
    r=fr16_record(); Path(path or ROOT/'results/FR16_wall_adapted_rv_box.json').write_text(json.dumps(r,indent=2)+'\n'); return r
if __name__=='__main__': print(json.dumps(write_result(),indent=2))
