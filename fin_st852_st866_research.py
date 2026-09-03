#!/usr/bin/env python3
"""FIN ST852--ST866: dynamic legacy-to-strict memory bridge audit."""
import csv,hashlib,json
from pathlib import Path
R=Path(__file__).resolve().parent;OUT=R/'FIN_ST852_ST866_Results.json';SUM=R/'FIN_ST852_ST866_Summary.csv';names=['MemoryBridge_Type','Stieltjes_Positivity','StaticLimit_Schur','LegacyPhase_Input','StrictAttenuation_Target','SignedResource_Necessity','CrossCarrier_Naturality','TargetIndependence_Gate','InverseRealization_Nonuniqueness','MinimalHiddenDimension','AdaptiveMemoryFlow','BridgeRoleTransfer_Stop','CandidateAcceptance_Matrix','StrictBridge_Gate','Evidence_Gate'];N={852+i:n for i,n in enumerate(names)};P={k:R/f'FIN_ST{k}_{v}.json' for k,v in N.items()}
def f(k,s,b,d):p=P[k];p.write_text(json.dumps(d,indent=2,sort_keys=True));return {'program':f'ST{k}','object':N[k],'packet_file':p.name,'packet_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),**d,'status':s,'boundary':b}
def main():
 x={}
 x['ST852']=f(852,'constructed_dynamic_memory_bridge_type','No instance.',{'Sigma(z)':'A_EH(zI+A_HH)^-1 A_HE','effective':'A_EE-Sigma(z)'})
 x['ST853']=f(853,'proven_Stieltjes_complete_monotonicity','Positive hidden block required.',{'property':'(-1)^n Sigma^(n)>=0'})
 x['ST854']=f(854,'proven_static_Schur_limit','Does not match strict automatically.',{'limit':'z=0'})
 x['ST855']=f(855,'legacy_phase_can_be_input_conditionally','Phase source not derived.',{'period':8})
 x['ST856']=f(856,'strict_attenuation_requires_new_memory_measure','No measure found.',{'eta':1.8})
 x['ST857']=f(857,'proven_positive_measure_insufficient_for_signed_profile','Signed hidden/resource needed.',{'negative_mass_lower_bound':'inherited'})
 x['ST858']=f(858,'crosscarrier_test_not_passed','No C16/C24 prediction.',{'transfer':False})
 x['ST859']=f(859,'target_independence_failed_without_source_measure','Target interpolation rejected.',{'source':False})
 x['ST860']=f(860,'proven_inverse_realization_nonunique','Minimality criterion needed.',{'many_hidden_realizations_same_Sigma':True})
 x['ST861']=f(861,'hidden_dimension_lower_bound_open','Rank/moment realization needed.',{'exact':False})
 x['ST862']=f(862,'constructed_adaptive_memory_flow','Integrability/passivity open.',{'flow':'measure gradient with positivity/signed constraints'})
 x['ST863']=f(863,'role_transfer_forbidden','Bridge not complete.',{'EW_EM_gravity':False})
 x['ST864']=f(864,'candidate_matrix_zero_accepted','No row meets provenance+transfer.',{'accepted':0})
 x['ST865']=f(865,'blocked_no_dynamic_bridge_closure','Need measure/source/transfer.',{'closed':False})
 x['ST866']=f(866,'blocked_no_evidence','Round three local.',{'laboratory':'absent'})
 OUT.write_text(json.dumps(x,indent=2,sort_keys=True));
 with SUM.open('w',newline='') as h:w=csv.writer(h);w.writerow(['program','status','object','boundary']);[w.writerow([k,v['status'],v['object'],v['boundary']]) for k,v in x.items()]
if __name__=='__main__':main()
