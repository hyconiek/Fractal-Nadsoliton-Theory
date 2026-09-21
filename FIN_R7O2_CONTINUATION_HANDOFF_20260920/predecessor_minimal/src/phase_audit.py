from pathlib import Path
import sys,json,math,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[1];H=ROOT/'inputs/FR223_20260916';sys.path.insert(0,str(H))
from src.phase_cumulants import k4_phase_value_grad_hess
old=json.load(open(H/'results/R7P-091_quartic_complement_partial.json'))
q=json.load(open(H/'certificates/R7P-090_quartic_root_boxes.json'))
f=json.load(open(H/'certificates/R7P-092_full_phase_roots.json'))
foundation=json.load(open(H/'results/R7P-081_090_phase_foundation.json'))
audit={'task':'R7N-033','old_partition':{'processed':old['processed'],'budget':old['budget'],'initial_n':old['initial_n'],
  'safe_leaf_count':old['safe_leaf_count'],'root_neighborhood_leaf_count':old['root_neighborhood_leaf_count'],
  'unresolved_leaf_count':old['unresolved_leaf_count'],'complete':old['complete'],'root_neighborhood_radius':old['root_neighborhood_radius']},
 'local_catalogs':{'quartic_count':q['count'],'quartic_radius':q['radius'],'quartic_global_pass':q['global_pass'],
                   'full_count':f['count'],'full_global_pass':f['global_pass']},
 'fixture_source':foundation['fixture'],
 'provenance_findings':['old complement used rounded decimal fixture','old phase endpoints are binary floats derived from 2*pi','old complement removed zero root-neighborhood leaves','radius 0.05 is not licensed by radius-1e-7 quartic uniqueness boxes','396 old safe leaves require revalidation before promotion to corrected exact fixture'],
 'scientific_status':'NUMERICAL_REPRODUCED_AND_PROVENANCE_AUDIT','conclusion':'Old cover is a subdivision seed only; it does not exhaust the exact corrected fixture.'}
json.dump(audit,open(ROOT/'results/R7N-033_phase_cover_audit.json','w'),indent=2)
# R7N-034 normalized chart derivative convention diagnostic
rng=np.random.default_rng(34034);amp=tuple(float(foundation['fixture'][k]) for k in ['r3','r4','r5','z6'])
errs=[]
for z in rng.random((16,3)):
    phi=2*math.pi*z
    val,g,Hh=k4_phase_value_grad_hess(*amp,phi)
    # finite difference derivative wrt z to validate 2pi scaling
    eps=1e-6; gfd=[]
    for i in range(3):
        zp=z.copy();zm=z.copy();zp[i]+=eps;zm[i]-=eps
        vp=k4_phase_value_grad_hess(*amp,2*math.pi*zp)[0];vm=k4_phase_value_grad_hess(*amp,2*math.pi*zm)[0]
        gfd.append((vp-vm)/(2*eps))
    errs.append(float(np.max(np.abs(np.asarray(gfd)-(2*math.pi)*g))))
contract={'task':'R7N-034','chart':'z=phi/(2*pi) in [0,1]^3 with periodic seam identification',
          'derivative_rule':'grad_z=(2*pi) grad_phi; Hess_z=(2*pi)^2 Hess_phi',
          'evaluation_rule':'use rigorous interval pi in future proof cells; rational z endpoints define coverage exactly',
          'fixture':foundation['fixture'],'alternating_sign':'z6 remains negative and fixed for this fixture',
          'numeric_chain_rule_max_error':max(errs),'samples':len(errs),
          'seam_contract':'z_i=0 and z_i=1 denote the same phase; a cover must account for seam cells exactly once without a gap',
          'scientific_status':'CONDITIONAL_LEMMA_AND_NUMERICAL_DIAGNOSTIC',
          'blocked_exact_piece':'corrected exact-decimal phase interval evaluator/certificates from the 2026-09-19 consolidated audit are unavailable in supplied inputs'}
json.dump(contract,open(ROOT/'proofs/R7N-034_normalized_phase_contract.json','w'),indent=2)
print(json.dumps({'R7N-033':audit,'R7N-034':contract},indent=2))
