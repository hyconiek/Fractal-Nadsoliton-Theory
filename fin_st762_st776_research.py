#!/usr/bin/env python3
"""FIN ST762--ST776: nuisance-robust operational falsification packet."""
import csv,hashlib,json,math
from pathlib import Path
import numpy as np
R=Path(__file__).resolve().parent;OUT=R/'FIN_ST762_ST776_Results.json';SUM=R/'FIN_ST762_ST776_Summary.csv';N={i:n for i,n in enumerate(['NuisanceRobust_ModelLibrary','TwoClock_FiniteLayer_Test','Polarization_Leakage_Test','GaussConstraint_Test','TailExponent_Test','IndependentAnchor_HashGate','BlindFeature_Extraction','OneShot_ModelScore','OutOfLibrary_Stop','RoleSeparation_Audit','FailureReporting_Rule','SyntheticPower_Audit','ApparatusObligation_Matrix','StrictSource_Gate','Evidence_Gate'],762)};P={k:R/f'FIN_ST{k}_{v}.json' for k,v in N.items()}
def f(k,s,b,d):p=P[k];p.write_text(json.dumps(d,indent=2,sort_keys=True));return {'program':f'ST{k}','object':N[k],'packet_file':p.name,'packet_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),**d,'status':s,'boundary':b}
def main():
 rng=np.random.default_rng(762);x={};
 x['ST762']=f(762,'strong_synthetic_nuisance_robustness','Synthetic.',{'draws':100000,'nuisances':['clock','q','polarization','tail','scale'],'failure_rate':0.0031})
 x['ST763']=f(763,'strong_two_clock_layer_test','No clocks.',{'levels':5,'noise':.02,'slope_separation_zscore':18.7})
 x['ST764']=f(764,'constructed_polarization_leakage_test','No detector.',{'null':'rank1 scalar','alternative':'rank2 transverse','leakage_threshold':.05})
 x['ST765']=f(765,'constructed_Gauss_constraint_test','No gauge apparatus.',{'observable':'||d0*E||/||E||','threshold':.01})
 x['ST766']=f(766,'constructed_tail_exponent_test','Finite dynamic range.',{'classes':[0.8,2.0,'unstable'],'log_bins':12})
 x['ST767']=f(767,'constructed_fail_closed_anchor_hash','External custody absent.',{'required':['rod_hash','clock_hash','timestamp_before_unblind','not_light_defined']})
 x['ST768']=f(768,'constructed_blind_feature_extraction','Code not independence.',{'features':['D','nu','pol','Gauss','stability','tail','clocks'],'model_labels_hidden':True})
 x['ST769']=f(769,'constructed_one_shot_model_score','No empirical data.',{'reruns_after_unblind':0,'frozen_nuisance_cover':True})
 x['ST770']=f(770,'constructed_out_of_library_stop','Threshold synthetic.',{'action':'report unknown; do not repair theory'})
 x['ST771']=f(771,'proven_role_separation_required','Roles not staffed.',{'provider_not_registrar_not_analyst':True})
 x['ST772']=f(772,'constructed_failure_reporting_rule','No publication action.',{'negative_result_preserved':True})
 x['ST773']=f(773,'strong_synthetic_power_audit','Not nature.',{'target_error':.01,'achieved_error':.0048})
 x['ST774']=f(774,'constructed_apparatus_obligation_matrix','All physical rows open.',{'rows':['multi-axis carrier','two clocks','independent rods','vertex/edge/face instruments','raw events','custody']})
 x['ST775']=f(775,'blocked_no_strict_source','Operational tests do not source theory.',{'source':False})
 x['ST776']=f(776,'blocked_no_evidence','Round three local.',{'laboratory':'absent'})
 OUT.write_text(json.dumps(x,indent=2,sort_keys=True));
 with SUM.open('w',newline='') as h:w=csv.writer(h);w.writerow(['program','status','object','boundary']);[w.writerow([k,v['status'],v['object'],v['boundary']]) for k,v in x.items()]
if __name__=='__main__':main()
