#!/usr/bin/env python3
import os
import json, math, pathlib, sys
import numpy as np
R7P=pathlib.Path(os.environ.get('R7P_ROOT','/mnt/data/r7p_source/unpacked/fin_rank7_followup')); sys.path.insert(0,str(R7P))
from src.model import feature_spaces
basep=pathlib.Path(os.environ.get('MP7_WORK_ROOT','/mnt/data/fin_rank7_mathphysics_next'))
q=json.load(open(basep/'results/MP7-022_023_quantitative_response.json'))
m=json.load(open(basep/'results/MP7-020_spectral_margin.json'))
cross=json.load(open(basep/'results/MP7-039_transverse_crossing.json'))
shift=json.load(open(basep/'results/MP7-035_local_finite_size_shift.json'))
# control 1: witness rescaling changes raw det but leaves normalized det/trace bound invariant
raw=2.6e-8; trG=3.0000017389; tau=67/250
bound=4*raw/(tau*tau*trG**3)
scale=2.0
raw2=raw*scale**6; trG2=trG*scale**2
bound2=4*raw2/(tau*tau*trG2**3)
# control 2: two response conventions differ. Reconstruct midpoint H at event and compare H^-1 vs microscopic response
fs=feature_spaces(); X=fs[3]; C4=fs[4]
cert=json.load(open(R7P/'certificates/R7P-026_equal_energy_event.json'))
s=np.array([(float(a)+float(b))/2 for a,b in cert['root_box'][:4]]); g=sum(map(float,cert['root_box'][4]))/2
h=C4@s; p=np.exp(h-h.max()); p/=p.sum(); mu=p@C4; Y=C4-mu; M=Y.T@(p[:,None]*Y); H=np.eye(4)/g-M
A=np.linalg.inv(H); micro=np.linalg.inv(np.eye(4)-g*M)@M
# control 3 M4 substitution breaks full determinant identity
p12=np.arange(1,13,dtype=float); p12/=p12.sum(); Sig=np.diag(p12)-np.outer(p12,p12); M7=X.T@Sig@X; M4=C4.T@Sig@C4
# full 11D identity
B0=np.column_stack([np.eye(12)[:,i]-np.eye(12)[:,11] for i in range(11)]); A7=X@X.T; gg=1.3
Hc=B0.T@(np.diag(1/p12)-gg*A7)@B0
lhs=np.linalg.det(Hc)*np.prod(p12); rhs7=np.linalg.det(np.eye(7)-gg*M7); rhs4=np.linalg.det(np.eye(4)-gg*M4)
controls=[
 {'name':'raw_PD_determinant_is_not_spectral_gap','passed':bool(abs(bound-bound2)<1e-20 and abs(raw2/raw-64)<1e-12),'detail':{'raw_det_original':raw,'raw_det_after_B_times_2':raw2,'normalized_bound_original':bound,'normalized_bound_rescaled':bound2}},
 {'name':'dual_and_microscopic_source_responses_not_interchangeable','passed':bool(float(np.linalg.norm(A-micro))>1e-2),'detail':{'frobenius_difference':float(np.linalg.norm(A-micro))}},
 {'name':'MP7_039_crossing_not_1D_pitchfork','passed':bool(cross['strict_inclusion'] and cross['critical_dimension']==2 and cross['branch_transverse_crossing_is_simple']),'detail':{'critical_dimension':cross['critical_dimension'],'D3_cubic':cross['cubic_invariant'],'dDelta_dg':cross['branch_total_derivative_delta_interval']}},
 {'name':'M4_substitution_breaks_full_simplex_determinant_identity','passed':bool(abs(lhs-rhs7)<1e-12 and abs(lhs-rhs4)>1e-4),'detail':{'lhs_full_simplex':lhs,'rhs_M7':rhs7,'rhs_wrong_M4':rhs4,'wrong_abs_error':abs(lhs-rhs4)}},
 {'name':'multiplicity_only_gives_wrong_finiteN_shift_logic','passed':bool(math.log(12)>0 and shift['A_equal_energy_log_weight_constant'][1]<0 and shift['shift_sign'].startswith('positive')),'detail':{'log12':math.log(12),'full_A':shift['A_equal_energy_log_weight_constant'],'certified_shift_sign':shift['shift_sign']}}
]
print(json.dumps({'task':'MP7-044-supplement','controls':controls,'all_passed':all(x['passed'] for x in controls)},indent=2))
