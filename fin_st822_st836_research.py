#!/usr/bin/env python3
"""FIN ST822--ST836: recommended-frontier intake and conditional source candidates."""
import csv,hashlib,json
from pathlib import Path
R=Path(__file__).resolve().parent;OUT=R/'FIN_ST822_ST836_Results.json';SUM=R/'FIN_ST822_ST836_Summary.csv';names=['NewObject_IntakeGate','Transition_GlobalMethod_Intake','Degree8_ExactMethod_Intake','IR_CorrelatedInterval_Intake','ContinuumCone_ProofUpgrade','FormalCore_DependencyFreeze','NonlinearNadsoliton_EquationCandidate','SevenCharge_Object_Intake','DynamicLegacyStrict_MemoryCandidate','MetricResponse_Gravity_Intake','GaugeComplex_Provenance_Intake','FineTomography_TwoClock_Intake','IndependentReplay_Plan','StrictClosure_Gate','Evidence_Gate'];N={822+i:n for i,n in enumerate(names)};P={k:R/f'FIN_ST{k}_{v}.json' for k,v in N.items()}
def f(k,s,b,d):p=P[k];p.write_text(json.dumps(d,indent=2,sort_keys=True));return {'program':f'ST{k}','object':N[k],'packet_file':p.name,'packet_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),**d,'status':s,'boundary':b}
def main():
 x={}
 x['ST822']=f(822,'proven_seven_field_intake_gate','Gate does not create candidates.',{'required':['type','transformation law','source formula','nonzero witness','coupling','target independence','heldout carrier']})
 x['ST823']=f(823,'transition_lane_requires_new_interval_or_SOS_method','Multistarts gated.',{'accepted_methods':['invariant-coordinate interval cover','SOS lower bound','explicit lower competitor']})
 x['ST824']=f(824,'degree8_lane_requires_structure_aware_Q_method','More primes gated.',{'blocks':[38,199],'accepted':['representation decomposition','p-adic height-controlled lift','exact CAS checkpoint']})
 x['ST825']=f(825,'IR_lane_requires_correlated_interval_AD','Radius/log replay gated.',{'needed':['F_b interval','predictor tube','common-root links']})
 x['ST826']=f(826,'conditional_continuum_upgrade_ready','Refinement conditional.',{'theorem_chain':['consistency','stability','density','uniform tail']})
 x['ST827']=f(827,'constructed_formal_dependency_DAG','Not machine checked.',{'nodes':['Dirichlet','spectral calculus','Schur','torsors','vacuum no-go','speed fiber','current-strict no-go']})
 x['ST828']=f(828,'constructed_conditional_nonlinear_nadsoliton_candidate','No strict source.',{'equation':'i psi_dot=A psi-chi |psi|^2 psi+i(P-gamma||psi||^2)psi','missing':['chi,P,gamma source','localization theorem','orbital stability']})
 x['ST829']=f(829,'seven_charge_candidate_not_supplied','No formula.',{'charges':['R','T','D','G','C','O','S'],'candidate':None})
 x['ST830']=f(830,'constructed_dynamic_memory_bridge_schema','No target-independent kernel.',{'object':'Sigma(z)=A_EH(zI+A_HH)^-1 A_HE','requirements':['legacy phase','strict attenuation','cross-carrier prediction']})
 x['ST831']=f(831,'metric_response_gravity_intake_failed_current_A8','New functor needed.',{'required':['g->A(g)','nonproxy variation','conserved T','Bianchi-I','units','heldout geometry'],'current_rows_complete':0})
 x['ST832']=f(832,'gauge_provenance_intake_failed','Constructed cubic complex conditional.',{'strict_d0_d1_source':False})
 x['ST833']=f(833,'fine_tomography_two_clock_transfer_ready','No apparatus.',{'q_recovery':'odd shift/2','clock_exponents':[2,1]})
 x['ST834']=f(834,'constructed_independent_replay_plan','No external reviewer.',{'artifacts':['scripts','JSON','hashes','proof claims']})
 x['ST835']=f(835,'blocked_no_new_strict_object_passes_intake','No strict closure.',{'accepted_candidates':0})
 x['ST836']=f(836,'blocked_no_evidence','Round one local.',{'laboratory':'absent'})
 OUT.write_text(json.dumps(x,indent=2,sort_keys=True));
 with SUM.open('w',newline='') as h:w=csv.writer(h);w.writerow(['program','status','object','boundary']);[w.writerow([k,v['status'],v['object'],v['boundary']]) for k,v in x.items()]
if __name__=='__main__':main()
